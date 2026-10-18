# 02 数据预处理与 MIL 样本构造

> **目标文献**：Tang Z, et al. *scTransMIL bridges patient-level disease states and single-cell transcriptomics for cancer screening and heterogeneity inference.* Nature Communications (2026) 17:10800. DOI: 10.1038/s41467-026-77538-5
> **模块对应原文**：Methods「Data preprocessing」；Fig. 1a（MIL 数据格式）；作者代码 `celldata_pancancer/1-process_v3.py`、`utils/tokenizer.py`
> **前置依赖**：`01_环境搭建与数据获取.md`（`raw/pancancer/*.h5ad`、Geneformer 字典、`gene_id.csv`）
> **复现难度 / 预计工时 / 计算资源**：★★★★☆ / 2–4 人日 / 内存 64–128 GB（逐文件处理）、磁盘 1–2 TB

---

## ① 科学问题与逻辑位置

- **上游**：104 个泛癌 h5ad（细胞 × 基因计数 + `Patient`/`Tissue`/`cnv_status`/`Cancer type`）。
- **下游**：
  - `processed/bags/<slide_id>/bag.h5ad`（≤256 cells）；
  - `processed/tokenized/<slide_id>.pt`（`input_ids`/`attention_mask`，N×1024）；
  - `processed/bag_index_pancancer.csv`（`case_id`/`slide_id`/`label`/`n_cells`）。
- **本模块回答**：**"如何把 4 M 个细胞切成 2 万个弱监督 bag，并让每个细胞变成一个不含绝对表达值的 token 序列？"**

这一步是整条流水线的"语义层"：bag 的边界、标签的来源、token 的顺序，任何一处与作者不同，后面的 attention 与 tumor score 都不可比。

---

## ② 背景原理

### 2.1 MIL 的数据格式

- **bag**：一次模型输入，包含 $m$ 个 instance（细胞）。作者把每个原始样本按 **256 cells** 切块；$m\le256$。
- **bag 标签**：$Y_j=1$ 若 bag 内至少一个细胞是恶性（CNV status = tumor），否则 $Y_j=0$。
- **instance 标签**：训练时不可见；模型靠 bag 标签 + 注意力/聚类损失间接学习。

> 这就是 MIL 的"OR 假设"：$Y_j=\max_i y_i$。它允许我们用样本级诊断训练，却得到细胞级分数。

### 2.2 为什么切 sub-sample

一个 sample 可能有数万细胞；一次 forward 把整样本放进 GPU 会 OOM。切块后：

- 每块 ≤256 cells × 1024 tokens，显存可控；
- 同一 donor 的多个块共享 bag 标签，增加样本量；
- 推理时再把块级预测聚合回样本级（05 号）。

### 2.3 rank value encoding（核心）

对细胞 $c$ 的基因表达向量 $\mathbf x_c\in\mathbb R^G$：

1. **文库归一化 + median scaling**：$x'_{c,g}=\dfrac{x_{c,g}/n_c\times10^4}{\text{median}_g}$（先按细胞总 counts $n_c$ 归一化到 1e4，再除以 `gene_median_dictionary`），避免文库大小与高表达基因垄断排序。
2. **去零**：只保留 $x'_{c,g}>0$ 的基因，长度自适应（$G_c\le G$）。
3. **降序排序**：按 $x'$ 从大到小，得到基因序列 $(g_{(1)},\dots,g_{(G_c)})$。
4. **token 映射**：$t_k=\text{vocab}[g_{(k)}]$（Ensembl → token ID）。
5. **截断/补齐**：取前 `max_len=1024` 个 token；不足补 `<pad>`，`attention_mask` 标记有效位。

**为什么 rank 能缓解批次效应**：不同平台/批次的绝对表达值差异大，但"哪些基因排在前 1000"相对稳定。rank 编码丢弃了数值尺度，只保留**相对次序**，因此对文库大小、归一化方法、平台差异更鲁棒。

> ⚠️ **导师批注**：rank 编码也丢掉了"表达量级"信息——一个基因从 100 降到 10，如果仍排在同一名次，模型看不见。这与 09 号"BLK 表达无差异但仍被 attention 找回"是同一枚硬币的两面：优势是鲁棒，代价是对"幅度变化"不敏感。

### 2.4 尾部回绕补齐（代码细节，正文未写）

作者 `split_ids()`：

```python
for i in range(0, len(ids), 256):
    subset = ids[i:i+256]
    if len(subset) < 256:
        subset += ids[:256-len(subset)]   # ★ 回绕复制开头细胞
```

最后一块不足 256 时，从样本开头复制细胞补齐。**副作用**：这些细胞在一个 bag 里出现两次（权重被放大）。复现时必须照做（D5），并在敏感性分析里试"丢弃最后一块"。

### 2.5 bag 过滤与多分类标签

- step2：`input_data.size(0) < 17` 的 bag 被删除（细胞数 <17 的 bag 信息量太低）。
- step3（多分类）：只保留 tumor bag；标签取 bag 内 tumor 细胞的 `Cancer type` **众数**。
- 二分类：`tumor_tissue` / `normal_tissue`。

---

## ③ 原文逐句拆解

| 原文句 | 中文转述 | 对应操作 | 原文是否给参数 | 缺失项决策 |
|---|---|---|---|---|
| "We utilized gene expression count matrices acquired from public repositories rather than processing raw sequencing data." | 用已定量矩阵 | 不跑 CellRanger | — | — |
| "genes within each individual cell were converted into a rank-based representation to mitigate technical variation" | rank 编码 | `TranscriptomeTokenizer` | 未给截断长度 | D4 |
| "we computationally partitioned each original sample into multiple sub-samples, with a maximum size of 256 cells per sub-sample." | 256 cells/bag | `split_ids` | 给了 256 | D5 |
| "a sub-sample was labeled as a 'tumor sub-sample' if it contained at least one malignant cell" | OR 标签 | `if 'tumor' in cnv_status` | 给了规则 | — |
| "For the consistent description, we use the bag to represent samples or sub-samples." | 术语 | — | — | — |
| 作者代码 step2 | 检查 <17 cells 的 bag | `input_data.size(0) < 17` | 代码给了 17 | ★ 记录：`split_ids` 已把每块补到 256，故该检查**实际不触发**（D5） |
| 作者代码 step3 | 多分类标签 | tumor 细胞的 `Cancer type` 众数 | 代码 | D16 |

---

## ④ 逐步可运行代码

### 4.1 基因 ID 映射表

```python
# code/02_00_gene_map.py
# 输入：作者仓库 celldata_pancancer/gene_id.csv（列：Gene name, Gene stable ID）
# 输出：processed/gene_id_map.tsv（gene_name, ensembl_id, token_id, mapped）
import pickle, pandas as pd

TOKEN_PKL = "raw/pretrained/geneformer/token_dictionary_gc95M.pkl"
with open(TOKEN_PKL, "rb") as f:
    token_dict = pickle.load(f)          # Ensembl ID -> token id

g = pd.read_csv("raw/pretrained/reference_code/celldata_pancancer/gene_id.csv")
# ★ 按列名读取（真实表头是 Gene stable ID,Gene name；不要按位置赋值）
g = g.rename(columns={"Gene stable ID": "ensembl_id", "Gene name": "gene_name"})
g["token_id"] = g["ensembl_id"].map(token_dict)
g["mapped"] = g["token_id"].notna()
assert g["token_id"].notna().sum() > 0, "token 映射为 0：检查列名/token 字典"
g.to_csv("processed/gene_id_map.tsv", sep="\t", index=False)
print("genes:", len(g), "mapped:", int(g.mapped.sum()),
      "lost:", int((~g.mapped).sum()))
```

**判据**：mapped 比例应 >90%；丢失的基因在 `logs/02_gene_map_lost.tsv` 记录。

### 4.2 切 bag 与标签

```python
# code/02_01_make_bags.py
# 用法：python code/02_01_make_bags.py --h5 raw/pancancer/xxx.h5ad --out processed/bags/pancancer
import argparse, os, numpy as np, pandas as pd, scanpy as sc

SUBSET = 256          # D5
MIN_CELLS = 17        # 作者 step2

def safe(x):
    # 文件名安全：空格/斜杠/冒号等替换为下划线
    import re
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(x))

def split_ids(ids, subset_size=SUBSET):
    subsets = []
    for i in range(0, len(ids), subset_size):
        s = ids[i:i+subset_size]
        if len(s) < subset_size:
            s = s + ids[:subset_size-len(s)]   # ★ 回绕复制（D5）
        subsets.append(s)
    return subsets

ap = argparse.ArgumentParser()
ap.add_argument("--h5", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()

ad = sc.read_h5ad(a.h5)
# 1) 基因映射到词表（用 4.1 的 map）
gmap = pd.read_csv("processed/gene_id_map.tsv", sep="\t")
name2ens = dict(zip(gmap.gene_name, gmap.ensembl_id))
name2ens = {k: v for k, v in name2ens.items() if isinstance(v, str)}

# 2) 只保留能映射的基因（兼容 var_names 是 SYMBOL 或 Ensembl 两种情况）
if str(ad.var_names[0]).startswith("ENSG"):
    ad.var["ensembl_id"] = list(ad.var_names)
else:
    ad.var["ensembl_id"] = [name2ens.get(g) for g in ad.var_names]
ad = ad[:, [x is not None for x in ad.var["ensembl_id"]]].copy()
ad.var_names = ad.var["ensembl_id"]
ad.var_names_make_unique()

rows = []
for patient in ad.obs["Patient"].unique():
    ap_ = ad[ad.obs["Patient"] == patient]
    for tissue in ap_.obs["Tissue"].unique():
        at = ap_[ap_.obs["Tissue"] == tissue]
        # 作者不在这里预过滤：split_ids 会把不足 256 的块回绕补到 256；
        # step2 的 `<17` 检查因此实际不触发（D5）。
        for idx, ids in enumerate(split_ids(list(at.obs_names))):
            bag = at[ids].copy()
            slide_id = f"{safe(patient)}_{safe(tissue)}_{idx}"
            outdir = os.path.join(a.out, slide_id)
            os.makedirs(outdir, exist_ok=True)
            bag.write(os.path.join(outdir, "bag.h5ad"))
            label = ("tumor_tissue" if "tumor" in set(bag.obs["cnv_status"])
                     else "normal_tissue")
            rows.append(dict(case_id=safe(patient), slide_id=slide_id,
                             label=label, n_cells=bag.n_obs,
                             dataset=os.path.basename(a.out)))
pd.DataFrame(rows).to_csv(f"{a.out}/bag_index_part.csv", index=False)
```

> 真实运行时对 104 个文件循环，最后用 `02_05_merge_index.py` 把所有 `bag_index_part.csv` 合并为 `processed/bag_index_pancancer.csv`（加 `dataset` 列、去重、校验 slide_id 唯一）。**注意**：不同文件的 `Patient`/`Tissue` 列名可能不同，需在 01 号的 inventory 里核对。

### 4.3 多分类标签

```python
# code/02_03_multiclass_labels.py
import pandas as pd, scanpy as sc
idx = pd.read_csv("processed/bag_index_pancancer.csv")
tumor = idx[idx.label == "tumor_tissue"].copy()
labels = []
for sid in tumor.slide_id:
    ad = sc.read_h5ad(f"processed/bags/pancancer/{sid}/bag.h5ad")
    t = ad[ad.obs["cnv_status"] == "tumor"]
    labels.append(t.obs["Cancer type"].value_counts().index[0])
tumor["label"] = labels
tumor.to_csv("processed/bag_index_pancancer_multicls.csv", index=False)
print(tumor.label.value_counts())
```

**判据**：类别数应为 29（代码 `label_dict`）或 30（正文）；缺类要登记（D16）。

### 4.3b 合并 104 个 part

```python
# code/02_05_merge_index.py
import glob, pandas as pd
parts = sorted(glob.glob("processed/bags/pancancer/*/bag_index_part.csv"))
df = pd.concat([pd.read_csv(f) for f in parts], ignore_index=True)
assert df.slide_id.is_unique, "slide_id 重复：检查 patient/tissue 命名"
df["dataset"] = "pancancer"
df.to_csv("processed/bag_index_pancancer.csv", index=False)
print("bags:", len(df), "donors:", df.case_id.nunique())
print(df.label.value_counts())
```

### 4.4 Tokenize（复用 Geneformer）

```python
# code/02_02_tokenize.py
# 用法：python code/02_02_tokenize.py --bag processed/bags/pancancer/P001_T_0 \
#        --out processed/tokenized/pancancer/P001_T_0.pt
import argparse, os, sys, pickle, numpy as np, torch
from datasets import load_from_disk

# ★ 作者 tokenizer 在 raw/pretrained/reference_code/utils/
sys.path.insert(0, "raw/pretrained/reference_code")
from utils.tokenizer import TranscriptomeTokenizer

TOKEN_PKL   = "raw/pretrained/geneformer/token_dictionary_gc95M.pkl"
MEDIAN_PKL  = "raw/pretrained/geneformer/gene_median_dictionary_gc95M.pkl"
MAPPING_PKL = "raw/pretrained/geneformer/ensembl_mapping_dict_gc95M.pkl"

with open(TOKEN_PKL, "rb") as f:
    token_dict = pickle.load(f)

def tokenize_one(bag_dir, out_dir, prefix):
    # ★ 必须显式传 4 个路径 + model_input_size/special_token；
    #   作者 utils/tokenizer.py 的默认值是 /aaa/... 与 4096/True
    tk = TranscriptomeTokenizer(
        custom_attr_name_dict={"Celltype": "Celltype",
                               "Cancer type": "Cancer type",
                               "cnv_status": "cnv_status"},
        nproc=4,
        model_input_size=4096,      # 作者实际值（与其注释的 30M 建议相反，见 D10）
        special_token=True,         # 作者实际值；会插入 <cls>/<eos>
        gene_median_file=MEDIAN_PKL,
        token_dictionary_file=TOKEN_PKL,
        gene_mapping_file=MAPPING_PKL,
    )
    tk.tokenize_data(data_directory=bag_dir, output_directory=out_dir,
                     output_prefix=prefix, file_format="h5ad")

def pad_and_truncate(token_dir, prefix, max_len=1024):
    ds = load_from_disk(os.path.join(token_dir, prefix))
    pad_id = int(token_dict.get("<pad>", 0))
    n = len(ds)
    ids = np.full((n, max_len), pad_id, dtype=np.int64)
    mask = np.zeros((n, max_len), dtype=np.int64)
    for i in range(n):
        seq = list(ds[i]["input_ids"])[:max_len]
        ids[i, :len(seq)] = seq
        mask[i, :len(seq)] = 1
    return torch.from_numpy(ids), torch.from_numpy(mask)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--bag", required=True)      # .../<slide_id>/bag.h5ad 的目录
    ap.add_argument("--out", required=True)      # .../<slide_id>.pt
    ap.add_argument("--tmp", default="processed/tokenized/tmp")
    a = ap.parse_args()
    prefix = os.path.basename(a.bag.rstrip("/"))
    os.makedirs(a.tmp, exist_ok=True)
    tokenize_one(a.bag, a.tmp, prefix)
    ids, mask = pad_and_truncate(a.tmp, prefix, max_len=1024)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    # ★ 唯一契约格式：dict(input_ids, attention_mask, slide_id)
    torch.save({"input_ids": ids, "attention_mask": mask, "slide_id": prefix}, a.out)
    print(prefix, ids.shape)
```

> **契约**：`processed/tokenized/<dataset>/<slide_id>.pt` 是 **dict**，含 `input_ids`（int64 N×1024）、`attention_mask`（int64 N×1024）、`slide_id`（str）。03/04/05 全部按此解包。作者原脚本保存两个独立 `.pt`；本项目改为单 dict 是为了避免两个文件错配，**内容与作者一致**（在 D5 记录格式差异）。

### 4.5 质量断言

```python
# code/02_04_assert_bags.py
import glob, torch, pandas as pd
rows = []
for f in glob.glob("processed/tokenized/pancancer/*.pt"):
    d = torch.load(f)
    ids, mask = d["input_ids"], d["attention_mask"]
    assert ids.shape == mask.shape
    assert ids.max() < 50000
    assert mask.sum(1).min() > 0
    rows.append(dict(slide_id=f.split("/")[-1][:-3], n_cells=ids.shape[0],
                     mean_len=float(mask.sum(1).float().mean())))
df = pd.DataFrame(rows)
print(df.describe())
assert df.n_cells.max() <= 256
df.to_csv("results/02_token_stats.csv", index=False)
```

---

## ⑤ 参数解释与调参空间

| 参数 | 本项目取值 | 可选 | 改变后果 |
|---|---|---|---|
| `subset_size` | 256 | 128–512 | 越大 bag 信息越多、显存越高；影响 tumor score 粒度 |
| `MIN_CELLS` | 17 | 1–50 | 过滤极小 bag；太小会引入噪声 |
| `max_len` | 1024 | 512–2048 | 序列越长信息越多、显存/时间越高（D4） |
| median scaling | Geneformer 字典 | 不用 | 排序被高表达基因主导 |
| 零值过滤 | 是 | 否 | 序列含大量 pad，注意力被稀释 |
| 尾部补齐 | 回绕复制 | 丢弃最后一块 | 影响最后一块的细胞权重（D5） |
| 基因映射 | `gene_id.csv` | 新版 Ensembl | 丢失率变化（D2） |

---

## ⑥ 预期结果与 QC 检查点

| 检查点 | 预期（原文/代码） | 允许波动 | Go/No-Go |
|---|---|---|---|
| 总 bag 数 | ≈20,000（原文报告值；4M/256≈15.6k 是下界，差异需在 gap 表解释） | ±10% | ❌ |
| 每 bag 细胞数 | ≤256 | 必须 | ❌ |
| bag 细胞数中位 | 200–256 | — | ⚠️ |
| 映射基因数 | ≈30,000 | ±5% | ⚠️ |
| 基因映射丢失率 | <10% | — | ⚠️ |
| 二分类标签比例 | tumor 占多数 | — | ⚠️ 记录 |
| 多分类类别数 | 29/30 | 见 D16 | ⚠️ |
| token 平均长度 | 数百–1024 | — | ⚠️ |
| `input_ids` 最大值 | < 词表大小 | 必须 | ❌ |

---

## ⑦ 常见陷阱与排错表

| 现象 | 根因 | 解决 |
|---|---|---|
| 映射基因为 0 | `gene_id.csv` 的列名/大小写不同 | 打印 `g.head()`；按实际列名改 |
| token 全是 pad | 序列被截断到 0 / token 字典 key 类型不符 | 检查 `token_dict` 的 key 是 str 还是 bytes |
| bag 数远少于预期 | `Patient`/`Tissue` 列名不同 | 用 01 号 inventory 核对 |
| 最后一块细胞重复 | 回绕补齐（设计如此） | 记录；敏感性分析可改丢弃 |
| 多分类缺类 | 某些癌种在过滤后无 tumor bag | 记录缺哪类（D16） |
| `load_from_disk` 字段名不同 | Geneformer 版本差异 | `print(ds.features)` 后调整 |
| 内存爆掉 | 一次读 4 M cells | 逐文件处理，及时 `del ad; gc.collect()` |

---

## ⑧ 与原文图表的对应关系

| 原文图 | 本模块产物 | 说明 |
|---|---|---|
| Fig. 1a（MIL 数据格式） | `processed/bag_index_pancancer.csv` | bag/instance 结构 |
| Fig. 1b（MIL 范式） | `processed/tokenized/*.pt` | 模型输入 |
| Fig. 2a（训练量 50/75/100%） | bag 索引 + splits | 05 号按比例抽样 |
| Fig. S3a（Donor setting） | 标签列 | 02 号可生成 Donor setting 版本 |
| Methods「Data preprocessing」 | 全部 | — |

> **交付物**：`processed/gene_id_map.tsv`、`processed/bag_index_pancancer.csv`、`processed/bags/`、`processed/tokenized/`、`results/02_token_stats.csv`。完成后进入 `03_模型架构_rank编码_Transformer_门控与聚合.md`。
