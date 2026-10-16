#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
spatial_cn_kmeans.py —— 空间邻域组成（Cellular Neighborhood, CN）的无监督划分
================================================================================
【用途】
  读入一张细胞表（每行一个细胞：x, y, celltype，可选分组列），对每个细胞取其
  空间窗口内的**细胞类型组成向量**（比例），按 K-means 聚类，得到每个细胞的
  CN 标签；并输出 CN × 细胞类型组成矩阵与每样本 CN 频率。

  本脚本**与具体文献/数据无关**：不硬编码任何项目路径、队列名、marker 面板、
  细胞类型名称，全部走参数（见 REPRO_SPEC.md §1.4 脚本分层）。

【方法学出处】
  after Wang et al. 2026 bioRxiv, Methods「CN definition」：
    70 μm 半径窗口 → 窗口内细胞类型组成比例 → K-means（k=15）→ 每细胞 CN 标签。
  本脚本把 70 μm / k=15 作为**参数默认值**保留（`--radius-um 70`、`--k 15`），
  但允许按数据集调整；窗口模式另提供 k-NN 窗口（`--window-mode knn`）作为
  半径窗口在**没有像素物理尺寸**时的等价替代。

【关键约定 / 常见坑（见 REPRO_SPEC.md 附录 C「坐标错配」）】
  * `--window-mode radius` **必须**同时给出 `--pixel-size-um`（像素的物理边长，
    μm/px），脚本据此把 μm 换算为像素。缺参数直接报错退出（不猜、不默认 1.0）。
  * `--window-mode knn` 不需要像素尺寸：k-NN 窗口按"邻居个数"定义，与单位无关。
  * 输出表行序稳定（按分组名、CN 标签号排序后写出），随机种子默认 20260101。

【用法示例】
  # 1) 半径窗口（推荐，复现原文 70 μm）
  python3 tools/spatial_cn_kmeans.py \\
      --cells processed/cells.csv --group-col sample \\
      --window-mode radius --radius-um 70 --pixel-size-um 0.5 \\
      --k 15 --out-prefix results/cn

  # 2) k-NN 窗口（无像素尺寸时）
  python3 tools/spatial_cn_kmeans.py --cells processed/cells.tsv --sep '\\t' \\
      --window-mode knn --k 70 --kmeans-k 15 --out-prefix results/cn

  # 3) 自检（合成组织 + ARI 断言）
  python3 tools/spatial_cn_kmeans.py --selftest

【输出】（前缀由 `--out-prefix` 指定）
  <prefix>_cell_cn.csv           细胞级表：原表全部列 + `cn`（CN 标签）+ `cn_raw`
                                 + `n_cells_in_window`（窗口内邻居数，不含自身）
                                 + `window_valid`（是否通过 --min-cells-in-window）
  <prefix>_cn_composition.csv    CN × 细胞类型组成矩阵（均值比例）+ n_cells + frac_cells
  <prefix>_cn_freq_by_sample.csv 每样本 CN 频率（宽表）+ n_cells

【依赖】
  必需：numpy, pandas, scipy, scikit-learn
  无（本脚本不使用任何可选依赖）

【退出码】
  0 = 成功
  1 = 未预期异常
  2 = 参数/输入校验失败（缺列、坐标非数值、缺少单位换算参数等）
  3 = 质量断言失败（产物为空、全部细胞被过滤、聚类退化）
  4 = --selftest 断言未通过
  5 = 缺少必需依赖
================================================================================
"""
from __future__ import annotations

import argparse
import os
import sys
import time

EXIT_OK = 0
EXIT_UNEXPECTED = 1
EXIT_USAGE = 2
EXIT_QUALITY = 3
EXIT_SELFTEST = 4
EXIT_DEP = 5


class UsageError(Exception):
    """参数/输入校验失败（退出码 2）。"""


class QualityError(Exception):
    """质量断言失败（退出码 3）。"""


def _require_deps():
    try:
        import numpy  # noqa: F401
        import pandas  # noqa: F401
        import scipy  # noqa: F401
        import sklearn  # noqa: F401
    except ImportError as exc:  # pragma: no cover
        sys.stderr.write(
            "[FATAL] 缺少必需依赖：%s\n"
            "        安装：python3 -m pip install --user numpy pandas scipy scikit-learn\n" % exc
        )
        raise SystemExit(EXIT_DEP)


# ----------------------------------------------------------------------------
# 通用 IO / 校验
# ----------------------------------------------------------------------------
def read_table(path, sep=None):
    import pandas as pd

    if not os.path.exists(path):
        raise UsageError("输入文件不存在：%s" % path)
    if os.path.isdir(path):
        raise UsageError("输入是目录而非文件：%s" % path)
    if sep is None:
        sep = "\t" if path.lower().endswith((".tsv", ".tab", ".txt")) else ","
    try:
        df = pd.read_csv(path, sep=sep)
    except Exception as exc:
        raise UsageError("读取失败（--sep 是否为 '%s'？）：%s" % (sep, exc))
    if df.shape[0] == 0:
        raise UsageError("输入表 0 行：%s" % path)
    df.columns = [str(c).strip() for c in df.columns]
    return df


def check_columns(df, cols, path):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise UsageError(
            "输入表缺少列 %s：%s\n        实际列：%s"
            % (missing, path, list(df.columns)[:30])
        )


def numeric_coords(df, x_col, y_col):
    import numpy as np

    out = []
    for col in (x_col, y_col):
        v = df[col]
        if v.dtype.kind not in "iuf":
            v = __import__("pandas").to_numeric(v, errors="coerce")
        bad = v.isna().to_numpy()
        if bad.any():
            raise UsageError(
                "列 '%s' 含 %d 个非数值/NA（示例行号：%s）—— 坐标必须为数值"
                % (col, int(bad.sum()), np.flatnonzero(bad)[:5].tolist())
            )
        out.append(v.to_numpy(dtype=float))
    return out[0], out[1]


# ----------------------------------------------------------------------------
# 核心：窗口组成向量
# ----------------------------------------------------------------------------
def window_composition(coords, codes, n_types, mode, radius_px=None, knn=None,
                       min_cells_in_window=1):
    """返回 (composition[n, n_types] 比例矩阵, n_eff[n] 窗口内邻居数, valid[n] bool)。

    `codes` 为 0..n_types-1 的整数细胞类型编码。窗口**不含自身**。
    """
    import numpy as np
    from scipy.spatial import cKDTree

    n = coords.shape[0]
    tree = cKDTree(coords)
    comp = np.zeros((n, n_types), dtype=float)
    n_eff = np.zeros(n, dtype=int)

    if mode == "radius":
        assert radius_px is not None
        neighbors = tree.query_ball_point(coords, r=float(radius_px))
    else:
        kq = min(int(knn) + 1, n)
        _, idx = tree.query(coords, k=kq)
        if idx.ndim == 1:  # n == 1
            idx = idx.reshape(1, -1)
        neighbors = [row for row in idx]

    for i in range(n):
        idxs = np.asarray(neighbors[i], dtype=int)
        if idxs.size:
            idxs = idxs[idxs != i]
        n_eff[i] = idxs.size
        if idxs.size:
            counts = np.bincount(codes[idxs], minlength=n_types)
            tot = counts.sum()
            if tot > 0:
                comp[i] = counts / float(tot)

    valid = n_eff >= int(min_cells_in_window)
    return comp, n_eff, valid


def relabel_by_size(labels, n_clusters):
    """把 K-means 原始簇号按簇内细胞数降序重映射为 1..k（输出稳定、可读）。"""
    import numpy as np

    counts = np.bincount(labels, minlength=n_clusters)
    order = np.argsort(-counts, kind="stable")
    remap = np.zeros(n_clusters, dtype=int)
    for new_id, old_id in enumerate(order, start=1):
        remap[old_id] = new_id
    return remap[labels], counts


def fit_cn(df, x_col, y_col, ct_col, group_col, mode, radius_um, pixel_size_um,
           knn, k, seed, scale, min_cells_in_window, relabel):
    """完整流程：窗口组成 → K-means → 细胞级 CN 标签。返回 dict。"""
    import numpy as np
    import pandas as pd
    from sklearn.cluster import KMeans
    from sklearn.preprocessing import StandardScaler

    x, y = numeric_coords(df, x_col, y_col)
    coords = np.column_stack([x, y])

    ct = df[ct_col].astype(str)
    if ct.isna().any() or (ct.str.strip() == "").any():
        raise UsageError("列 '%s' 含空/NA 细胞类型标签" % ct_col)
    types = sorted(ct.unique().tolist())
    type_to_code = {t: i for i, t in enumerate(types)}
    codes_all = ct.map(type_to_code).to_numpy(dtype=int)

    if group_col is not None:
        if group_col not in df.columns:
            raise UsageError(
                "输入表缺少分组列 '%s'；实际列：%s" % (group_col, list(df.columns)[:30])
            )
        groups_all = df[group_col].astype(str).to_numpy()
    else:
        groups_all = np.full(len(df), "__all__", dtype=object)

    radius_px = None
    if mode == "radius":
        if pixel_size_um is None:
            raise UsageError(
                "--window-mode radius 需要 --pixel-size-um（像素物理边长，μm/px）"
                "才能把 --radius-um 换算为像素。\n"
                "        若确实没有像素物理尺寸，请改用 --window-mode knn --k <邻居数>。"
                "（REPRO_SPEC 附录 C「坐标错配」高危项，本脚本拒绝猜测默认值）"
            )
        radius_px = radius_um / float(pixel_size_um)

    comp = np.zeros((len(df), len(types)), dtype=float)
    n_eff = np.zeros(len(df), dtype=int)
    valid = np.zeros(len(df), dtype=bool)

    group_order = sorted(pd.unique(groups_all).tolist())
    n_windows_by_group = {}
    for g in group_order:
        m = groups_all == g
        if m.sum() == 0:
            continue
        c, ne, v = window_composition(
            coords[m], codes_all[m], len(types), mode,
            radius_px=radius_px, knn=knn, min_cells_in_window=min_cells_in_window,
        )
        comp[m] = c
        n_eff[m] = ne
        valid[m] = v
        n_windows_by_group[g] = int(m.sum())

    if valid.sum() < k:
        raise QualityError(
            "通过 --min-cells-in-window=%d 的细胞只有 %d 个，少于 --k=%d —— 无法聚类。"
            "请放宽 --min-cells-in-window 或检查坐标单位/像素尺寸。"
            % (min_cells_in_window, int(valid.sum()), k)
        )

    X = comp[valid]
    if scale:
        X = StandardScaler().fit_transform(X)

    km = KMeans(n_clusters=k, random_state=seed, n_init=10)
    raw = km.fit_predict(X)

    if relabel:
        labels, cluster_sizes = relabel_by_size(raw, k)
    else:
        labels = raw + 1
        cluster_sizes = np.bincount(raw, minlength=k)

    cn_full = np.full(len(df), None, dtype=object)
    cn_raw_full = np.full(len(df), np.nan)
    cn_full[valid] = ["cn_%d" % v for v in labels]
    cn_raw_full[valid] = labels

    out = df.copy()
    out["cn"] = cn_full
    out["cn_raw"] = cn_raw_full
    out["n_cells_in_window"] = n_eff
    out["window_valid"] = valid

    return {
        "df": out,
        "types": types,
        "codes": codes_all,
        "groups": groups_all,
        "valid": valid,
        "labels": cn_full,
        "k": k,
        "comp": comp,
        "cluster_sizes": cluster_sizes,
        "n_windows_by_group": n_windows_by_group,
        "radius_px": radius_px,
        "inertia": float(km.inertia_),
    }


def build_outputs(res):
    import numpy as np
    import pandas as pd

    types = res["types"]
    df = res["df"]
    valid = res["valid"]
    labels = res["labels"]

    lab_valid = np.asarray(labels[valid], dtype=object)
    cn_order = sorted(pd.unique(lab_valid).tolist(),
                      key=lambda s: int(str(s).split("_")[1]))

    # --- CN × 细胞类型组成矩阵（各 CN 内细胞窗口组成比例的均值）----------------
    comp = res["comp"][valid]
    rows = []
    for cn in cn_order:
        m = lab_valid == cn
        rows.append(
            dict([("cn", cn)] + [(t, float(comp[m, i].mean())) for i, t in enumerate(types)]
                 + [("n_cells", int(m.sum())),
                    ("frac_cells", float(m.sum()) / float(valid.sum()))])
        )
    cn_comp = pd.DataFrame(rows, columns=["cn"] + types + ["n_cells", "frac_cells"])

    # --- 每样本 CN 频率（宽表）------------------------------------------------
    groups = np.asarray(res["groups"], dtype=object)
    freq_rows = []
    for g in sorted(pd.unique(groups[valid]).tolist()):
        m = valid & (groups == g)
        labs = np.asarray(labels[m], dtype=object)
        tot = labs.size
        row = {"group": g, "n_cells": int(tot),
               "n_unassigned": int((~valid & (groups == g)).sum())}
        for cn in cn_order:
            row[cn] = float((labs == cn).sum()) / float(tot) if tot else float("nan")
        freq_rows.append(row)
    cn_freq = pd.DataFrame(freq_rows, columns=["group", "n_cells", "n_unassigned"] + cn_order)

    return cn_comp, cn_freq


def write_outputs(res, out_prefix):
    import pandas as pd

    cn_comp, cn_freq = build_outputs(res)
    outdir = os.path.dirname(os.path.abspath(out_prefix))
    os.makedirs(outdir, exist_ok=True)

    cell_path = "%s_cell_cn.csv" % out_prefix
    comp_path = "%s_cn_composition.csv" % out_prefix
    freq_path = "%s_cn_freq_by_sample.csv" % out_prefix

    cell = res["df"].sort_values(["cn"], na_position="last", kind="stable")
    cell.to_csv(cell_path, index=False)
    cn_comp.to_csv(comp_path, index=False)
    cn_freq.to_csv(freq_path, index=False)

    for p in (cell_path, comp_path, freq_path):
        if (not os.path.exists(p)) or os.path.getsize(p) == 0:
            raise QualityError("产物为空：%s（REPRO_SPEC 附录 C 静默失败）" % p)
    if cn_comp.shape[0] != res["k"]:
        raise QualityError(
            "CN 组成表只有 %d 行，应为 --k=%d 行（有 CN 未获得任何细胞）"
            % (cn_comp.shape[0], res["k"])
        )
    return cell_path, comp_path, freq_path, cn_comp, cn_freq


# ----------------------------------------------------------------------------
# 自检：合成组织（已知 CN 结构）+ ARI
# ----------------------------------------------------------------------------
def make_synthetic_tissue(seed=20260101, cells_per_zone=300, n_samples=2,
                          zone_size=400.0, gap=400.0, sparse_cells=40):
    """构造 4 个已知组成的细胞邻域区 + 1 个稀疏区。

    区块之间留 gap（> 2×窗口半径），使窗口不跨区污染；稀疏区用于检验
    `--min-cells-in-window` 的稀疏窗口过滤确实生效。

    返回 (DataFrame[celltype,x,y,sample,true_cn], 组成字典)
    """
    import numpy as np
    import pandas as pd

    mixtures = {
        "core":     {"tumor": 0.80, "stroma": 0.06, "tcell": 0.10, "macro": 0.04},
        "boundary": {"tumor": 0.45, "stroma": 0.35, "tcell": 0.15, "macro": 0.05},
        "stroma":   {"tumor": 0.10, "stroma": 0.80, "tcell": 0.06, "macro": 0.04},
        "immune":   {"tumor": 0.08, "stroma": 0.12, "tcell": 0.60, "macro": 0.20},
    }
    types_all = ["tumor", "stroma", "tcell", "macro"]
    names = sorted(mixtures.keys())
    frames = []
    for s in range(n_samples):
        for zi, zn in enumerate(names):
            gx, gy = divmod(zi, 2)
            x0, y0 = gx * (zone_size + gap), gy * (zone_size + gap)
            rng_z = np.random.default_rng(seed + 1000 * s + zi)
            probs = np.array([mixtures[zn][t] for t in types_all])
            pick = rng_z.choice(len(types_all), size=cells_per_zone, p=probs)
            frames.append(pd.DataFrame({
                "celltype": [types_all[i] for i in pick],
                "x": x0 + rng_z.uniform(0, zone_size, cells_per_zone),
                "y": y0 + rng_z.uniform(0, zone_size, cells_per_zone),
                "sample": "sample_%d" % (s + 1),
                "true_cn": zn,
            }))
        # 稀疏区（第 3 行）：密度低到窗口内邻居数常 < 5，应被 --min-cells-in-window 过滤
        rng_s = np.random.default_rng(seed + 9999 + s)
        frames.append(pd.DataFrame({
            "celltype": [types_all[i] for i in rng_s.choice(len(types_all), sparse_cells)],
            "x": rng_s.uniform(0, zone_size, sparse_cells),
            "y": 2 * (zone_size + gap) + rng_s.uniform(0, zone_size, sparse_cells),
            "sample": "sample_%d" % (s + 1),
            "true_cn": "sparse",
        }))
    return pd.concat(frames, ignore_index=True), mixtures


def selftest(seed=20260101, pixel_size_um=1.0):
    import numpy as np
    import pandas as pd
    from sklearn.metrics import adjusted_rand_score

    print("=" * 78)
    print("spatial_cn_kmeans.py --selftest")
    print("=" * 78)
    df, mixtures = make_synthetic_tissue(seed=seed)
    print("[合成组织] 细胞数=%d，样本=%s，真实 CN=%s"
          % (len(df), sorted(df["sample"].unique()), sorted(df["true_cn"].unique())))
    print("[合成组织] 真实组成：")
    for zn in sorted(mixtures.keys()):
        print("           %-9s %s" % (zn, {k: round(v, 2) for k, v in mixtures[zn].items()}))

    results = {}
    for mode, kwargs in (
        ("radius", dict(mode="radius", radius_um=70.0, pixel_size_um=pixel_size_um, knn=None)),
        ("knn",    dict(mode="knn", radius_um=70.0, pixel_size_um=pixel_size_um, knn=70)),
    ):
        t0 = time.time()
        res = fit_cn(df, "x", "y", "celltype", "sample", k=4, seed=seed, scale=False,
                     min_cells_in_window=5, relabel=True, **kwargs)
        valid = res["valid"]
        # ARI 只在真实 4 区上计算（稀疏区属于"应被过滤"的细胞，不参与 ARI）
        in_truth4 = df["true_cn"].to_numpy() != "sparse"
        keep = valid & in_truth4
        ari = adjusted_rand_score(
            df["true_cn"].to_numpy()[keep], np.asarray(res["labels"], dtype=object)[keep]
        )
        dt = time.time() - t0
        n_drop = int((~valid).sum())
        n_sparse = int((~valid & ~in_truth4).sum())
        print("\n[模式 %s] 用时 %.2f s；有效细胞 %d/%d（%d 个窗口邻居数 < 5 被过滤，"
              "其中稀疏区 %d/%d）"
              % (mode, dt, int(valid.sum()), len(df), n_drop, n_sparse, int((~in_truth4).sum())))
        print("[模式 %s] K-means inertia = %.1f" % (mode, res["inertia"]))
        # CN × 真实区的交叉表（仅真实 4 区）
        ct = pd.crosstab(np.asarray(res["labels"], dtype=object)[keep],
                         df["true_cn"].to_numpy()[keep])
        print("[模式 %s] CN × 真实 CN 交叉表：" % mode)
        for line in ct.to_string().splitlines():
            print("           " + line)
        print("[模式 %s] >>> 调整兰德指数 ARI = %.4f <<<" % (mode, ari))
        results[mode] = dict(ari=ari, seconds=dt, n_valid=int(valid.sum()),
                             n_drop=n_drop, n_sparse=n_sparse)

    # --- 汇总 + 断言 ---------------------------------------------------------
    print("\n" + "-" * 78)
    print("汇总：")
    for mode, r in results.items():
        print("  %-7s ARI=%.4f  用时=%.2fs  有效细胞=%d  过滤=%d（稀疏区 %d）"
              % (mode, r["ari"], r["seconds"], r["n_valid"], r["n_drop"], r["n_sparse"]))

    failed = []
    for mode, r in results.items():
        if not (r["ari"] > 0.7):
            failed.append("%s: ARI=%.4f <= 0.7" % (mode, r["ari"]))
        if r["n_valid"] < 0.5 * len(df):
            failed.append("%s: 有效细胞过少（%d/%d）" % (mode, r["n_valid"], len(df)))
        if mode == "radius" and r["n_sparse"] < 1:
            # 稀疏区过滤只在"半径窗口"下有定义：k-NN 窗口总能凑够 k 个邻居
            failed.append("%s: --min-cells-in-window 未过滤掉任何稀疏区细胞（过滤逻辑未生效）"
                          % mode)

    if failed:
        print("\n[FAIL] selftest 未通过：")
        for f in failed:
            print("  - " + f)
        return EXIT_SELFTEST
    print("\n[PASS] selftest 通过：半径窗口与 k-NN 窗口的 ARI 均 > 0.7")
    return EXIT_OK


# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------
def build_parser():
    p = argparse.ArgumentParser(
        prog="spatial_cn_kmeans.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "空间邻域组成（CN）K-means 划分：窗口组成比例 → K-means → 每细胞 CN 标签。\n"
            "after Wang et al. 2026 bioRxiv, Methods「CN definition」"
            "（70 μm 半径窗口 → 组成比例 → K-means k=15）。"
        ),
        epilog="""\
单位与符号约定
  --window-mode radius : 窗口半径由 --radius-um（μm）给出，必须配合 --pixel-size-um
                         （μm/px）换算为像素；缺 --pixel-size-um 直接报错退出（不猜默认值）。
  --window-mode knn    : 窗口 = 空间上最近的 --k 个细胞（不含自身），与单位无关，
                         不需要 --pixel-size-um。
  窗口内计数不含细胞自身；组成向量按比例归一（和为 1）。

输出表
  <prefix>_cell_cn.csv            细胞级：原列 + cn / cn_raw / n_cells_in_window / window_valid
  <prefix>_cn_composition.csv     CN × 细胞类型：各 CN 内窗口组成比例的均值 + n_cells + frac_cells
  <prefix>_cn_freq_by_sample.csv  每样本 CN 频率（宽表）+ n_cells + n_unassigned

示例
  python3 tools/spatial_cn_kmeans.py --cells processed/cells.csv --group-col sample \\
      --window-mode radius --radius-um 70 --pixel-size-um 0.5 --k 15 --out-prefix results/cn
  python3 tools/spatial_cn_kmeans.py --cells processed/cells.tsv --sep '\\t' \\
      --window-mode knn --k 70 --out-prefix results/cn
  python3 tools/spatial_cn_kmeans.py --selftest

退出码：0=成功 1=未预期异常 2=参数/输入校验失败 3=质量断言失败 4=selftest 失败 5=缺必需依赖
""",
    )
    io = p.add_argument_group("输入输出")
    io.add_argument("--cells", metavar="FILE", help="细胞表（CSV/TSV，每行一个细胞）")
    io.add_argument("--sep", default=None, help="列分隔符（默认按扩展名：.tsv/.tab/.txt 用制表符，否则逗号）")
    io.add_argument("--x-col", default="x", help="x 坐标列名（默认 x）")
    io.add_argument("--y-col", default="y", help="y 坐标列名（默认 y）")
    io.add_argument("--celltype-col", default="celltype", help="细胞类型列名（默认 celltype）")
    io.add_argument("--group-col", default="sample",
                    help="分组列名（默认 sample）；置为 '' 表示不分组（全部细胞视为一组）")
    io.add_argument("--out-prefix", default="cn_out", metavar="PREFIX",
                    help="输出前缀，产物为 <PREFIX>_cell_cn.csv / _cn_composition.csv / _cn_freq_by_sample.csv")

    w = p.add_argument_group("窗口定义")
    w.add_argument("--window-mode", choices=["radius", "knn"], default="radius",
                   help="窗口模式：radius（--radius-um + --pixel-size-um）或 knn（--k 个邻居）")
    w.add_argument("--radius-um", type=float, default=70.0,
                   help="半径窗口半径，单位 μm（默认 70，复现原文）")
    w.add_argument("--pixel-size-um", type=float, default=None,
                   help="像素物理边长，μm/px（radius 模式**必需**；knn 模式忽略）")
    w.add_argument("--knn", "--knn-k", "--window-knn-k", dest="knn", type=int, default=70,
                   help="knn 模式下窗口邻居数（默认 70）。注意：本脚本 --k 指 K-means 簇数，"
                        "k-NN 窗口大小请用 --knn，避免二者混淆")
    w.add_argument("--min-cells-in-window", type=int, default=1, metavar="N",
                   help="窗口内邻居数 < N 的细胞标记为 window_valid=False 且不参与聚类（稀疏区过滤，默认 1）")

    c = p.add_argument_group("K-means")
    c.add_argument("--k", "--kmeans-k", "--kmeans_k", dest="kmeans_k", type=int, default=15,
                   help="K-means 簇数 k（默认 15，复现原文）。--k / --kmeans-k 等价；"
                        "k-NN 窗口大小用 --knn，二者不共用")
    c.add_argument("--seed", type=int, default=20260101, help="随机种子（默认 20260101）")
    c.add_argument("--scale", dest="scale", action="store_true",
                   help="聚类前对组成向量做 z-score 标准化")
    c.add_argument("--no-scale", dest="scale", action="store_false",
                   help="聚类前不标准化，直接用原始比例向量（默认）")
    p.set_defaults(scale=False)
    c.add_argument("--relabel-by-size", dest="relabel", action="store_true",
                   help="把簇号按簇内细胞数降序重编号为 cn_1..cn_k（默认开启，输出稳定）")
    c.add_argument("--no-relabel-by-size", dest="relabel", action="store_false",
                   help="保留 K-means 原始簇号 +1")
    p.set_defaults(relabel=True)

    p.add_argument("--selftest", action="store_true",
                   help="合成组织自检：跑完整流程并报告调整兰德指数 ARI（断言 ARI > 0.7）")
    return p


def validate(args):
    if args.selftest:
        return
    if not args.cells:
        raise UsageError("缺少 --cells（或使用 --selftest）")
    if args.x_col == args.y_col:
        raise UsageError("--x-col 与 --y-col 不能相同（%s）" % args.x_col)
    if args.window_mode == "radius":
        if args.radius_um is None or not (args.radius_um > 0):
            raise UsageError("--radius-um 必须为正数，当前 %r" % args.radius_um)
        if args.pixel_size_um is None:
            raise UsageError(
                "--window-mode radius 需要 --pixel-size-um（μm/px）。\n"
                "        缺此参数会把 μm 当像素用，导致尺度错配（REPRO_SPEC 附录 C）。\n"
                "        若没有像素物理尺寸，请用 --window-mode knn --k <邻居数>。"
            )
    if args.pixel_size_um is not None and not (args.pixel_size_um > 0):
        raise UsageError("--pixel-size-um 必须为正数，当前 %r" % args.pixel_size_um)
    if args.window_mode == "knn" and args.knn < 1:
        raise UsageError("--knn 必须 >= 1，当前 %d" % args.knn)
    if args.kmeans_k < 2:
        raise UsageError("--k/--kmeans-k 必须 >= 2，当前 %d" % args.kmeans_k)
    if args.min_cells_in_window < 0:
        raise UsageError("--min-cells-in-window 必须 >= 0，当前 %d" % args.min_cells_in_window)


def main(argv=None):
    args = build_parser().parse_args(argv)
    _require_deps()

    if args.selftest:
        return selftest()
    try:
        validate(args)
        group_col = args.group_col if args.group_col not in (None, "") else None
        df = read_table(args.cells, args.sep)
        check_columns(df, [args.x_col, args.y_col, args.celltype_col], args.cells)
        if group_col is not None:
            check_columns(df, [group_col], args.cells)

        t0 = time.time()
        res = fit_cn(df, args.x_col, args.y_col, args.celltype_col, group_col,
                     args.window_mode, args.radius_um, args.pixel_size_um, args.knn,
                     args.kmeans_k, args.seed, args.scale, args.min_cells_in_window,
                     args.relabel)
        cell_p, comp_p, freq_p, cn_comp, cn_freq = write_outputs(res, args.out_prefix)
        dt = time.time() - t0

        valid = res["valid"]
        print("[INFO] 细胞数=%d，细胞类型数=%d，分组数=%d"
              % (len(df), len(res["types"]), len(res["n_windows_by_group"])))
        print("[INFO] 窗口模式=%s%s，有效细胞 %d/%d（--min-cells-in-window=%d）"
              % (args.window_mode,
                 ("，radius=%.1f μm = %.2f px（pixel_size=%.4g μm/px）"
                  % (args.radius_um, res["radius_px"], args.pixel_size_um))
                 if args.window_mode == "radius" else "，knn_k=%d" % args.knn,
                 int(valid.sum()), len(df), args.min_cells_in_window))
        print("[INFO] K-means k=%d，inertia=%.2f，标准化=%s，用时 %.2fs"
              % (args.kmeans_k, res["inertia"], args.scale, dt))
        print("[INFO] 各 CN 细胞数：%s"
              % {r["cn"]: int(r["n_cells"]) for _, r in cn_comp.iterrows()})

        tiny = [r["cn"] for _, r in cn_comp.iterrows() if r["frac_cells"] < 0.01]
        if tiny:
            print("[WARN] 以下 CN 占比 < 1%%，结果解释需谨慎（可能为噪声簇）：%s" % tiny)
        print("[OK] 写出：%s\n[OK] 写出：%s\n[OK] 写出：%s" % (cell_p, comp_p, freq_p))
        return EXIT_OK
    except UsageError as exc:
        sys.stderr.write("[USAGE] %s\n" % exc)
        return EXIT_USAGE
    except QualityError as exc:
        sys.stderr.write("[QUALITY] %s\n" % exc)
        return EXIT_QUALITY
    except SystemExit:
        raise
    except Exception as exc:  # pragma: no cover
        import traceback
        traceback.print_exc()
        sys.stderr.write("[FATAL] %s: %s\n" % (type(exc).__name__, exc))
        return EXIT_UNEXPECTED


if __name__ == "__main__":
    sys.exit(main())
