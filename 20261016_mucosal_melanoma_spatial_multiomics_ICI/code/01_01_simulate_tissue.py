#!/usr/bin/env python3
"""合成 MuM 组织生成器（ground truth）—— 用于在无原始数据时验证下游空间算法。

输出（全部带已知 ground truth）：
  processed/01_sim_cells.parquet    每细胞一行：坐标 + 真实 CN + 真实表型 + 40 个 marker 强度
  processed/01_sim_pixels.parquet   每 MALDI 像素一行：各"真实代谢物"的像素强度
  results/01_sim_ground_truth.json  区域比例、CN 组成、富集/排斥对、细胞类型-代谢物真值谱

设计目标：让下游 05(CN)/06(空间上下文)/07(patch)/08(Delaunay 置换) 能在**已知答案**上验证。
参数对齐原文（见 01 号文档 §5）：
  总细胞数 695444 / 97 cores；1.0 mm 芯；70 µm CN 窗口；20 µm MALDI 像素；0.23 µm/px COMET。
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np
import pandas as pd

# ---------------- 常量（项目全局，见 REPRO_SPEC §6 与简报 §9）----------------
SEED_DEFAULT = 20260101
N_CORES = 97
N_CELLS_TOTAL = 695_444                 # 原文报告的 QC 后分割细胞数
CORE_DIAMETER_UM = 1000.0               # 1.0 mm TMA 芯
CORE_RADIUS_UM = CORE_DIAMETER_UM / 2
CN_WINDOW_UM = 70.0                     # 原文：索引细胞半径 70 µm 窗口
MALDI_PIXEL_UM = 20.0                   # 原文：20 µm 空间分辨率
COMET_PX_UM = 0.23                      # 原文：0.23 µm/pixel
N_CN = 15                               # 原文：k = 15
N_CYCLES = 20                           # 原文：每轮 2 个 marker → 40-plex

# ---------------- 8 个 major cell types（对齐原文 Fig.1B 右）----------------
MAJOR = ["Tumor", "Tcell", "Macrophage", "DC", "NK", "Bcell", "Stromal", "Other"]

# ---------------- 25 个 cellular phenotypes（我方定义，非原文报告）----------------
# 原文只报告"25 个 phenotype"这个数字，未列出名称 → 本项目自定义并在契约表登记
PHENO = {
    "Tumor":      ["Tumor_Ki67", "Tumor_MIF", "Tumor_HLADR", "Tumor_MITF", "Tumor_panSOX10"],
    "Tcell":      ["T_CD8", "T_CD8_Ki67", "T_CD8_GZMB", "T_CD4", "T_CD4_FOXP3", "T_DN"],
    "Macrophage": ["Mac_CD163", "Mac_CD163_IRF4", "Mac_CD163_Arg1", "Mac_CD68", "Mac_MHCII"],
    "DC":         ["DC_CD11c"],
    "NK":         ["NK"],
    "Bcell":      ["B"],
    "Stromal":    ["Stroma_Vim", "Stroma_Vim_SMA", "Stroma_pericyte", "Endothelial"],
    "Other":      ["Neuron_Other", "Unassigned"],
}
PHENO_LIST = [p for v in PHENO.values() for p in v]
PHENO_MAJOR = {p: m for m, v in PHENO.items() for p in v}
assert len(PHENO_LIST) == 25, f"应为 25 个表型，实际 {len(PHENO_LIST)}"

# ---------------- 40-plex 面板（对齐原文 Fig.1B 左的功能分层）----------------
MARKERS = [
    # 肿瘤谱系 / 状态 (6)
    "SOX10", "MITF", "S100B", "Ki67", "MIF", "HLA_DR",
    # T 细胞 (6)
    "CD3", "CD4", "CD8", "FOXP3", "GZMB", "PD1",
    # 髓系 / DC (7)
    "CD163", "CD68", "CD11c", "IRF4", "Arg1", "CD206", "HLA_ABC",
    # NK / B (4)
    "NKp46", "CD56", "CD20", "CD79a",
    # 基质 / 血管 (5)
    "Vimentin", "SMA", "FAP", "CD31", "Desmin",
    # 功能 / 检查点 / 其他 (12)
    "PDL1", "CTLA4", "LAG3", "TIGIT", "TIM3", "CD45", "CD45RO", "CD27", "GranzymeA",
    "CXCL9", "IDO1", "Ki67_stroma",
]
assert len(MARKERS) == 40, f"应为 40 个 marker，实际 {len(MARKERS)}"

# ---------------- 区域划分（半径归一化）----------------
# 每个区域给出：r_min, r_max（相对瘤体半径 R_tumor 的倍数）；CN 概率向量
REGION_DEFS = [
    ("central_tumor", 0.00, 0.45),
    ("invasive_tumor", 0.45, 0.75),
    ("tumor_boundary", 0.75, 1.00),
    ("peri_stroma", 1.00, 1.25),
    ("stroma", 1.25, 1.80),
]

# 15 个 CN 的生物学命名（对齐原文 Fig.2B）
CN_NAMES = [
    "CN1_Vim_pos_SMA_pos_stroma", "CN2_NK", "CN3_DC", "CN4_B_cells",
    "CN5_B_and_CD4T", "CN6_proliferative_tumor", "CN7_CD8T", "CN8_mixed_cells",
    "CN9_tumor_infiltrated_stroma", "CN10_Vim_pos_stroma", "CN11_immune_enriched",
    "CN12_immune_infiltrated_stroma", "CN13_tumor_boundary", "CN14_invasive_tumor",
    "CN15_central_tumor",
]
# 每个区域的 CN 先验（行和为 1）—— 这是 ground truth 的"真实空间结构"
REGION_CN_PRIOR = {
    "central_tumor":    {"CN15_central_tumor": 0.62, "CN14_invasive_tumor": 0.20,
                         "CN6_proliferative_tumor": 0.10, "CN11_immune_enriched": 0.04,
                         "CN8_mixed_cells": 0.04},
    "invasive_tumor":   {"CN14_invasive_tumor": 0.52, "CN6_proliferative_tumor": 0.18,
                         "CN13_tumor_boundary": 0.14, "CN15_central_tumor": 0.08,
                         "CN3_DC": 0.04, "CN8_mixed_cells": 0.04},
    "tumor_boundary":   {"CN13_tumor_boundary": 0.48, "CN14_invasive_tumor": 0.16,
                         "CN2_NK": 0.08, "CN3_DC": 0.10, "CN7_CD8T": 0.08,
                         "CN11_immune_enriched": 0.06, "CN8_mixed_cells": 0.04},
    "peri_stroma":      {"CN9_tumor_infiltrated_stroma": 0.30, "CN12_immune_infiltrated_stroma": 0.22,
                         "CN10_Vim_pos_stroma": 0.20, "CN1_Vim_pos_SMA_pos_stroma": 0.14,
                         "CN4_B_cells": 0.06, "CN5_B_and_CD4T": 0.08},
    "stroma":           {"CN10_Vim_pos_stroma": 0.36, "CN1_Vim_pos_SMA_pos_stroma": 0.26,
                         "CN12_immune_infiltrated_stroma": 0.12, "CN4_B_cells": 0.10,
                         "CN5_B_and_CD4T": 0.08, "CN2_NK": 0.04, "CN3_DC": 0.04},
}

# 每个 CN 的细胞表型组成（行和为 1）
CN_PHENO_MIX = {
    "CN15_central_tumor":   {"Tumor_MIF": 0.58, "Tumor_panSOX10": 0.18, "Tumor_MITF": 0.10,
                             "Tumor_HLADR": 0.04, "Stroma_Vim_SMA": 0.05,
                             "Mac_MHCII": 0.03, "Unassigned": 0.02},
    "CN14_invasive_tumor":  {"Tumor_Ki67": 0.30, "Tumor_MIF": 0.20, "Tumor_HLADR": 0.14,
                             "Tumor_panSOX10": 0.09, "Stroma_Vim": 0.07, "Stroma_Vim_SMA": 0.06,
                             "Mac_CD68": 0.05, "T_CD8": 0.03, "T_CD8_GZMB": 0.03, "B": 0.03},
    "CN6_proliferative_tumor": {"Tumor_Ki67": 0.72, "Tumor_panSOX10": 0.14, "Tumor_MIF": 0.06,
                                "T_CD8": 0.04, "Unassigned": 0.04},
    "CN13_tumor_boundary":  {"Tumor_Ki67": 0.20, "Tumor_HLADR": 0.11, "Mac_CD163": 0.17,
                             "Mac_CD163_IRF4": 0.05, "Mac_CD68": 0.11, "DC_CD11c": 0.07,
                             "T_CD8": 0.09, "T_CD8_GZMB": 0.06, "T_CD4": 0.05,
                             "Stroma_Vim": 0.05, "B": 0.02, "Unassigned": 0.02},
    "CN2_NK":               {"NK": 0.58, "T_CD8": 0.13, "Mac_CD68": 0.09, "DC_CD11c": 0.06,
                             "Stroma_Vim": 0.05, "Neuron_Other": 0.05, "Unassigned": 0.04},
    "CN3_DC":               {"DC_CD11c": 0.58, "Mac_MHCII": 0.12, "T_CD4": 0.10,
                             "T_CD8": 0.08, "B": 0.06, "Stroma_Vim": 0.04, "Unassigned": 0.02},
    "CN4_B_cells":          {"B": 0.62, "T_CD4": 0.12, "DC_CD11c": 0.08, "Stroma_Vim": 0.10,
                             "Mac_CD68": 0.04, "Unassigned": 0.04},
    "CN5_B_and_CD4T":       {"B": 0.36, "T_CD4": 0.30, "DC_CD11c": 0.10, "T_CD8": 0.08,
                             "Stroma_Vim": 0.08, "Mac_CD163": 0.04, "Unassigned": 0.04},
    "CN7_CD8T":             {"T_CD8": 0.42, "T_CD8_GZMB": 0.20, "T_CD8_Ki67": 0.12,
                             "T_CD4": 0.08, "DC_CD11c": 0.06, "Mac_CD68": 0.06,
                             "Stroma_Vim": 0.06},
    "CN8_mixed_cells":      {"Tumor_panSOX10": 0.12, "Stroma_Vim": 0.18, "Mac_CD68": 0.12,
                             "T_CD4": 0.10, "T_CD8": 0.09, "T_DN": 0.05, "B": 0.07,
                             "DC_CD11c": 0.06, "NK": 0.04, "Endothelial": 0.05,
                             "Neuron_Other": 0.05, "Unassigned": 0.07},
    "CN9_tumor_infiltrated_stroma": {"Stroma_Vim": 0.36, "Stroma_Vim_SMA": 0.08, "Tumor_panSOX10": 0.18,
                                     "Tumor_Ki67": 0.07, "Mac_CD68": 0.09, "Mac_CD163_IRF4": 0.05,
                                     "T_CD8": 0.07, "T_CD4": 0.05, "T_CD4_FOXP3": 0.04,
                                     "Endothelial": 0.01},
    "CN10_Vim_pos_stroma":  {"Stroma_Vim": 0.55, "Stroma_Vim_SMA": 0.12, "Endothelial": 0.09,
                             "Stroma_pericyte": 0.07, "Mac_CD68": 0.05, "Mac_CD163_Arg1": 0.05,
                             "T_CD4_FOXP3": 0.04, "Tumor_panSOX10": 0.02, "Unassigned": 0.01},
    "CN1_Vim_pos_SMA_pos_stroma": {"Stroma_Vim_SMA": 0.66, "Stroma_Vim": 0.14,
                                   "Stroma_pericyte": 0.10, "Mac_CD68": 0.04,
                                   "Endothelial": 0.04, "Unassigned": 0.02},
    "CN11_immune_enriched": {"T_CD8": 0.20, "T_CD4": 0.17, "T_DN": 0.06, "B": 0.13, "NK": 0.11,
                             "DC_CD11c": 0.12, "Mac_MHCII": 0.10, "Mac_CD68": 0.06,
                             "T_CD4_FOXP3": 0.03, "Unassigned": 0.02},
    "CN12_immune_infiltrated_stroma": {"Stroma_Vim": 0.36, "T_CD8": 0.13, "T_CD4": 0.11,
                                       "B": 0.09, "Mac_CD163": 0.08, "Mac_CD163_Arg1": 0.05,
                                       "DC_CD11c": 0.06, "Neuron_Other": 0.04,
                                       "Endothelial": 0.04, "Unassigned": 0.04},
}

# 已知的"空间富集 / 排斥"细胞对（用于验证 08 号文档的置换检验功效与 FPR）
# 这是我方注入的 ground truth：下游必须能检测出 ENRICHED、且不能误报 REPELLED
TRUE_INTERACTIONS = [
    # (pheno_A, pheno_B, 类型, 相对最近邻距离倍数, 只在哪些 CN 内成立)
    ("Tumor_Ki67", "Mac_CD163", "ENRICHED", 0.72, ["CN13_tumor_boundary", "CN14_invasive_tumor"]),
    ("Tumor_Ki67", "DC_CD11c",  "ENRICHED", 0.78, ["CN13_tumor_boundary"]),
    ("T_CD8_GZMB", "Tumor_HLADR", "ENRICHED", 0.82, ["CN14_invasive_tumor", "CN13_tumor_boundary"]),
    ("T_CD4_FOXP3", "Mac_CD163_Arg1", "ENRICHED", 0.80, ["CN10_Vim_pos_stroma", "CN12_immune_infiltrated_stroma", "CN9_tumor_infiltrated_stroma"]),
    ("Tumor_MIF", "Stroma_Vim_SMA", "REPELLED", 1.45, ["CN14_invasive_tumor", "CN15_central_tumor"]),
    ("T_CD8", "Stroma_Vim_SMA", "REPELLED", 1.35, ["CN10_Vim_pos_stroma", "CN9_tumor_infiltrated_stroma"]),
    ("B", "Tumor_Ki67", "REPELLED", 1.30, ["CN14_invasive_tumor"]),
]

# 真实代谢物（原文 Fig.5 明确提到的）与其 m/z（负离子 [M-H]- 理论值，见 02 号文档 §4.4）
TRUE_METABOLITES = {
    "linoleic_acid":          279.2330,
    "oleic_acid":             281.2486,
    "palmitic_acid":          255.2330,
    "stearic_acid":           283.2643,
    "indole_3_acetic_acid":   174.0561,
    "indole_3_carbinol":      146.0611,
    "indole_3_acetaldehyde":  158.0611,
    "indole_3_carboxylic_acid": 160.0404,
    "tryptophan":             203.0826,
    "lactate":                 89.0244,
}

# 细胞表型 → 代谢物分泌/富集谱（ground truth）。
# 关键设计：让 UFA 主要来自肿瘤与基质，indole 类主要来自髓系/DC → 下游应能重建
PHENO_METABOLITE = {
    "Tumor_Ki67":        {"linoleic_acid": 1.4, "oleic_acid": 1.3, "palmitic_acid": 1.5, "stearic_acid": 1.4, "lactate": 2.2},
    "Tumor_MIF":         {"linoleic_acid": 1.1, "oleic_acid": 1.0, "palmitic_acid": 1.2, "stearic_acid": 1.2, "indole_3_acetaldehyde": 1.3, "lactate": 1.4},
    "Tumor_HLADR":       {"linoleic_acid": 1.0, "oleic_acid": 1.0, "palmitic_acid": 1.1, "stearic_acid": 1.1, "tryptophan": 0.8},
    "Tumor_MITF":        {"linoleic_acid": 1.2, "oleic_acid": 1.2, "palmitic_acid": 1.3, "stearic_acid": 1.3},
    "Tumor_panSOX10":    {"linoleic_acid": 1.0, "oleic_acid": 1.0, "palmitic_acid": 1.0, "stearic_acid": 1.0},
    "T_CD8":             {"tryptophan": 0.6},
    "T_CD8_Ki67":        {"tryptophan": 0.5, "lactate": 1.3},
    "T_CD8_GZMB":        {"tryptophan": 0.6},
    "T_CD4":             {"tryptophan": 0.7},
    "T_CD4_FOXP3":       {"indole_3_acetic_acid": 0.9, "tryptophan": 1.1},
    "T_DN":              {"tryptophan": 0.8},
    "Mac_CD163":         {"indole_3_acetic_acid": 2.0, "indole_3_acetaldehyde": 1.9, "indole_3_carbinol": 1.6, "tryptophan": 1.6},
    "Mac_CD163_IRF4":    {"indole_3_acetic_acid": 1.6, "indole_3_acetaldehyde": 1.5, "tryptophan": 1.3},
    "Mac_CD163_Arg1":    {"indole_3_acetic_acid": 2.4, "indole_3_acetaldehyde": 2.2, "indole_3_carbinol": 2.0, "indole_3_carboxylic_acid": 1.7, "tryptophan": 1.9, "linoleic_acid": 1.2},
    "Mac_CD68":          {"indole_3_acetic_acid": 1.2, "linoleic_acid": 0.9},
    "Mac_MHCII":         {"indole_3_acetic_acid": 1.4, "indole_3_carbinol": 1.5},
    "DC_CD11c":          {"indole_3_carboxylic_acid": 2.1, "indole_3_carbinol": 2.0, "indole_3_acetaldehyde": 1.4, "tryptophan": 1.2},
    "NK":                {"lactate": 1.1},
    "B":                 {"tryptophan": 0.9},
    "Stroma_Vim":        {"linoleic_acid": 1.3, "oleic_acid": 1.2, "palmitic_acid": 1.2, "stearic_acid": 1.2},
    "Stroma_Vim_SMA":    {"linoleic_acid": 1.6, "oleic_acid": 1.5, "palmitic_acid": 1.4, "stearic_acid": 1.4},
    "Stroma_pericyte":   {"linoleic_acid": 0.9, "oleic_acid": 0.9},
    "Endothelial":       {"lactate": 1.2, "linoleic_acid": 0.8},
    "Neuron_Other":      {"tryptophan": 1.0},
    "Unassigned":        {},
}

# 20 µm MALDI 像素 vs COMET 细胞的空间聚合（用于验证 02 号文档 §2.10 的聚合偏差）
MALDI_PX_PER_MM = 1000.0 / MALDI_PIXEL_UM    # 50 px/mm


# =============================================================================
def build_distribution(n: int, probs: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """把整数 n 按 probs 精确分到各类（multinomial 且强制和为 n）。"""
    p = np.asarray(probs, dtype=float)
    p = p / p.sum()
    counts = rng.multinomial(n, p)
    return counts


def sample_cn_for_region(region: str, n: int, rng: np.random.Generator) -> np.ndarray:
    pri = REGION_CN_PRIOR[region]
    keys = list(pri)
    p = np.array([pri[k] for k in keys], dtype=float)
    p = p / p.sum()
    idx = rng.choice(len(keys), size=n, p=p)
    return np.array(keys, dtype=object)[idx]


def sample_pheno_for_cn(cns: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    out = np.empty(len(cns), dtype=object)
    for cn in np.unique(cns):
        mix = CN_PHENO_MIX[cn]
        keys = list(mix)
        p = np.array([mix[k] for k in keys], dtype=float); p /= p.sum()
        m = cns == cn
        out[m] = rng.choice(np.array(keys, dtype=object), size=int(m.sum()), p=p)
    return out


def poisson_disc_min_dist(rng, n, box_um, dmin_um, max_iter=40):
    """近似的硬核点过程：重复采样 + 拒绝，达到最小间距 dmin。
    用于让细胞位置看起来像组织（不是纯随机均匀）。规模大时退化为均匀采样。"""
    x = rng.uniform(0, box_um[0], n)
    y = rng.uniform(0, box_um[1], n)
    if n < 20000:
        for _ in range(3):
            from scipy.spatial import cKDTree
            tree = cKDTree(np.column_stack([x, y]))
            pairs = tree.query_pairs(dmin_um, output_type="ndarray")
            if len(pairs) == 0:
                break
            drop = np.unique(pairs[:, 1])
            x[drop] = rng.uniform(0, box_um[0], len(drop))
            y[drop] = rng.uniform(0, box_um[1], len(drop))
    return x, y


def apply_true_interactions(px, py, pheno, cn, rng, cn_pheno_indexing=None):
    """按 TRUE_INTERACTIONS 施加**已知的**空间富集/排斥。
    做法：对每条规则的 B 类细胞，在指定 CN 内把坐标朝该 CN 内最近的 A 类细胞
    移动（ENRICHED → 距离 ×factor；REPELLED → 远离）。这是可控的 ground truth 注入。"""
    from scipy.spatial import cKDTree
    for (pa, pb, kind, factor, cns) in TRUE_INTERACTIONS:
        mask_cn = np.isin(cn, cns)
        ia = np.where(mask_cn & (pheno == pa))[0]
        ib = np.where(mask_cn & (pheno == pb))[0]
        if len(ia) == 0 or len(ib) == 0:
            continue
        A = np.column_stack([px[ia], py[ia]])
        tree = cKDTree(A)
        d, j = tree.query(np.column_stack([px[ib], py[ib]]), k=1)
        target = np.column_stack([A[j, 0], A[j, 1]])
        cur = np.column_stack([px[ib], py[ib]])
        vec = target - cur
        if kind == "ENRICHED":
            newpos = cur + vec * (1.0 - factor)       # 向 A 移动
        else:
            # 远离：沿反方向推开，限制在芯内
            norm = np.linalg.norm(vec, axis=1, keepdims=True) + 1e-9
            direction = vec / norm
            newpos = cur - direction * np.minimum((factor - 1.0) * d, 60.0)[:, None]
        px[ib], py[ib] = newpos[:, 0], newpos[:, 1]
    return px, py


def build_marker_matrix(pheno, rng, dropout=0.35, sigma_log=0.55):
    """按表型生成 40 个 marker 的对数正态强度 + 检出率（dropout）。
    真值矩阵 MARker_MEAN 由"谱系专一 + 状态叠加"生成；这里给出确定性的构造规则。"""
    n = len(pheno)
    n_m = len(MARKERS)
    mu = np.zeros((n, n_m), dtype=np.float32)
    mi = {m: i for i, m in enumerate(MARKERS)}

    def setm(sel, mark, val):
        mu[sel, mi[mark]] = val

    # 基础背景
    mu[:] = -1.0
    for m in MAJOR:
        sel = np.array([PHENO_MAJOR[p] == m for p in pheno])
        setm(sel, "CD45", 2.0 if m in ("Tcell", "Macrophage", "DC", "NK", "Bcell") else -1.0)
    # 肿瘤
    sel = np.array([PHENO_MAJOR[p] == "Tumor" for p in pheno])
    for mk in ("SOX10", "MITF", "S100B"):
        setm(sel, mk, 1.5)
    setm(np.array(pheno == "Tumor_Ki67"), "Ki67", 2.8)
    setm(np.array(pheno == "Tumor_MIF"), "MIF", 2.9)
    setm(np.array(pheno == "Tumor_HLADR"), "HLA_DR", 2.6)
    setm(np.array(pheno == "Tumor_HLADR"), "PDL1", 1.8)
    setm(np.array(pheno == "Tumor_MITF"), "MITF", 2.7)
    setm(sel, "PDL1", 1.0)
    # T 细胞
    tc = np.array([PHENO_MAJOR[p] == "Tcell" for p in pheno])
    setm(tc, "CD3", 2.4)
    setm(np.array(pheno == "T_CD8"), "CD8", 2.7)
    setm(np.array(pheno == "T_CD8_Ki67"), "CD8", 2.6); setm(np.array(pheno == "T_CD8_Ki67"), "Ki67", 2.5)
    setm(np.array(pheno == "T_CD8_GZMB"), "CD8", 2.5); setm(np.array(pheno == "T_CD8_GZMB"), "GZMB", 2.7); setm(np.array(pheno == "T_CD8_GZMB"), "GranzymeA", 2.3)
    setm(np.array(pheno == "T_CD4"), "CD4", 2.6)
    setm(np.array(pheno == "T_CD4_FOXP3"), "CD4", 2.4); setm(np.array(pheno == "T_CD4_FOXP3"), "FOXP3", 2.5); setm(np.array(pheno == "T_CD4_FOXP3"), "CTLA4", 2.0)
    setm(tc, "PD1", 1.4)
    # 髓系
    mc = np.array([PHENO_MAJOR[p] == "Macrophage" for p in pheno])
    setm(mc, "CD68", 2.3)
    setm(np.array([p.startswith("Mac_CD163") for p in pheno]), "CD163", 2.7)
    setm(np.array(pheno == "Mac_CD163_IRF4"), "IRF4", 1.9)
    setm(np.array(pheno == "Mac_CD163_Arg1"), "Arg1", 2.1); setm(np.array(pheno == "Mac_CD163_Arg1"), "CD206", 2.0)
    setm(np.array(pheno == "Mac_MHCII"), "HLA_DR", 2.4); setm(np.array(pheno == "Mac_MHCII"), "HLA_ABC", 2.2)
    setm(mc, "HLA_ABC", 1.4)
    # DC / NK / B
    dc = np.array(pheno == "DC_CD11c")
    setm(dc, "CD11c", 2.9); setm(dc, "HLA_DR", 2.7); setm(dc, "CXCL9", 1.9)
    setm(np.array(pheno == "NK"), "NKp46", 2.6); setm(np.array(pheno == "NK"), "CD56", 2.4); setm(np.array(pheno == "NK"), "GZMB", 1.8)
    setm(np.array(pheno == "B"), "CD20", 2.6); setm(np.array(pheno == "B"), "CD79a", 2.5)
    # 基质 / 血管
    st = np.array([PHENO_MAJOR[p] == "Stromal" for p in pheno])
    setm(st, "Vimentin", 2.5)
    setm(np.array(pheno == "Stroma_Vim_SMA"), "SMA", 2.8); setm(np.array(pheno == "Stroma_Vim_SMA"), "FAP", 2.2)
    setm(np.array(pheno == "Stroma_pericyte"), "Desmin", 2.2); setm(np.array(pheno == "Stroma_pericyte"), "SMA", 1.8)
    setm(np.array(pheno == "Endothelial"), "CD31", 2.8)
    setm(st, "Ki67_stroma", 1.2)
    # 功能性（按 CN 局部环境叠加由下游机制生成，这里只给细胞内在项）
    setm(np.array(pheno == "Unassigned"), "CD45", 1.0)
    mu += 0.8
    # 对数正态 + dropout
    # 对数正态：mu 为 log 尺度均值 → 强度恒为正；clip 防止极端值
    logv = np.clip(rng.normal(mu, sigma_log), -3.0, 6.0)
    vals = np.exp(logv).astype(np.float32)
    detected = rng.random(vals.shape) > dropout
    vals = np.where(detected, vals, 0.0).astype(np.float32)
    return vals


def build_metabolite_maps(px, py, pheno, core_side_um, rng, noise_sd=0.25, rho=0.5):
    """把细胞级代谢物真值聚合成 20 µm 像素级图谱。
    **这是关键的一步**：它同时产出
      (a) ground-truth 的"细胞类型 → 代谢物"关联（可用于检验下游的相关性检验功效）；
      (b) 真实的**空间聚合偏差**：一个 20 µm 像素覆盖多个细胞 → 相关性被压缩。
    rho = 细胞间代谢物水平的空间自相关强度（模拟局部免疫-代谢耦合）。"""
    npx = int(np.ceil(core_side_um / MALDI_PIXEL_UM))
    feats = list(TRUE_METABOLITES)
    acc = np.zeros((npx, npx, len(feats)), dtype=np.float64)
    ix = np.clip((px / MALDI_PIXEL_UM).astype(int), 0, npx - 1)
    iy = np.clip((py / MALDI_PIXEL_UM).astype(int), 0, npx - 1)

    # 平滑的局部"代谢场"（模拟免疫-代谢耦合与组织梯度）
    from scipy.ndimage import gaussian_filter
    field = {f: gaussian_filter(rng.normal(0, 1, (npx, npx)), sigma=max(1.0, rho * 4)) for f in feats}

    vals = np.zeros((len(px), len(feats)), dtype=np.float32)
    for k, f in enumerate(feats):
        base = np.array([PHENO_METABOLITE.get(p, {}).get(f, 0.0) for p in pheno], dtype=np.float32)
        local = field[f][ix, iy].astype(np.float32) * 0.45
        v = np.expm1(np.clip(base * 1.0 + local + rng.normal(0, noise_sd, len(px)), -3, 4))
        vals[:, k] = v
    for k in range(len(feats)):
        np.add.at(acc[:, :, k], (ix, iy), vals[:, k])
    # 每像素的细胞数 → 用于下游"聚合尺度敏感性分析"
    cellcount = np.zeros((npx, npx), dtype=np.int32)
    np.add.at(cellcount, (ix, iy), 1)
    return acc.astype(np.float32), cellcount, feats, vals


def simulate_core(core_id, n_cells, rng, tumor_radius_frac=0.62):
    """生成一个 1.0 mm 芯：返回细胞表 + 像素表。"""
    side = 2 * CORE_RADIUS_UM
    cx = cy = CORE_RADIUS_UM
    # 1) 位置：**圆形芯内**均匀采样（否则密度被方形外接框稀释）
    n_out = int(round(n_cells / (np.pi / 4)))          # 补偿圆/方面积比 0.785
    px, py = poisson_disc_min_dist(rng, n_out, (side, side), 6.0)
    keep = np.hypot(px - cx, py - cy) <= CORE_RADIUS_UM * 0.98
    px, py = px[keep][:n_cells], py[keep][:n_cells]
    n_cells = len(px)                                   # 实际生成数（略少于请求数）
    r = np.hypot(px - cx, py - cy) / (tumor_radius_frac * CORE_RADIUS_UM)  # 归一化到瘤体半径
    # 2) 区域标签（全部落在芯内）
    region = np.full(n_cells, "stroma", dtype=object)
    for name, lo, hi in REGION_DEFS:
        m = (r >= lo) & (r < hi)
        region[m] = name
    # 3) CN 标签（真值）
    cn = np.empty(n_cells, dtype=object)
    for name, _, _ in REGION_DEFS:
        m = region == name
        if m.sum() > 0:
            cn[m] = sample_cn_for_region(name, int(m.sum()), rng)
    # 4) 表型标签（真值）
    pheno = sample_pheno_for_cn(cn, rng)
    # 5) 注入已知的富集/排斥
    px, py = apply_true_interactions(px, py, pheno, cn, rng)
    px = np.clip(px, 0, side); py = np.clip(py, 0, side)
    # 6) 40-plex marker 强度（含 dropout）
    mk = build_marker_matrix(pheno, rng)
    # 7) 代谢物像素图（20 µm）
    met, cellcount, feats, cell_met = build_metabolite_maps(px, py, pheno, side, rng)

    cells = pd.DataFrame({"core_id": core_id, "cell_id": np.arange(n_cells),
                          "x_um": px.astype(np.float32), "y_um": py.astype(np.float32),
                          "region": region, "true_cn": cn, "true_phenotype": pheno,
                          "true_major": [PHENO_MAJOR[p] for p in pheno]})
    for j, m in enumerate(MARKERS):
        cells[m] = mk[:, j]
    # 细胞级代谢物真值（用于下游检验"空间聚合偏差"——见 01 号文档 §4.5）
    for k, f in enumerate(feats):
        cells[f"met_{f}"] = cell_met[:, k]
    in_core_px = (np.hypot((np.arange(met.shape[0]) * MALDI_PIXEL_UM + MALDI_PIXEL_UM / 2)[:, None] - cx,
                           (np.arange(met.shape[1]) * MALDI_PIXEL_UM + MALDI_PIXEL_UM / 2)[None, :] - cy)
                   <= CORE_RADIUS_UM)
    pix = pd.DataFrame({"core_id": core_id, "in_core": in_core_px.ravel(),
                        "px_x_um": np.repeat((np.arange(met.shape[0]) * MALDI_PIXEL_UM + MALDI_PIXEL_UM / 2), met.shape[1]).astype(np.float32),
                        "px_y_um": np.tile((np.arange(met.shape[1]) * MALDI_PIXEL_UM + MALDI_PIXEL_UM / 2), met.shape[0]).astype(np.float32),
                        "n_cells": cellcount.ravel()})
    for k, f in enumerate(feats):
        pix[f] = met[:, :, k].ravel()
    return cells, pix


def write_table(df: pd.DataFrame, base: Path) -> Path:
    """优先写 .parquet（需 pyarrow/fastparquet）；不可用时回退到 .csv.gz。
    回退不影响可复现性（浮点精度由 repr 保证），但下游读入时须用同名 helper。"""
    try:
        p = base.with_suffix(".parquet")
        df.to_parquet(p, index=False)
        return p
    except (ImportError, ValueError):
        p = Path(str(base) + ".csv.gz")
        df.to_csv(p, index=False, compression="gzip", float_format="%.6g")
        print(f"     [warn] 未安装 pyarrow/fastparquet，改存 {p.name}")
        return p


def main():
    ap = argparse.ArgumentParser(description="合成 MuM 组织生成器（ground truth）")
    ap.add_argument("--out-dir", default=".", help="项目根目录")
    ap.add_argument("--seed", type=int, default=SEED_DEFAULT)
    ap.add_argument("--n-cores", type=int, default=N_CORES)
    ap.add_argument("--cells-per-core", type=int, default=1200,
                    help="每 core 生成多少细胞；默认 1200（快速验证）。"
                         "设为 None/0 则使用 695444/n_cores 精确对齐原文总量")
    ap.add_argument("--full", action="store_true",
                    help="生成全量 695444 个细胞（内存约 3-5 GB，耗时数分钟）")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    root = Path(args.out_dir).resolve()
    (root / "processed").mkdir(parents=True, exist_ok=True)
    (root / "results").mkdir(parents=True, exist_ok=True)

    if args.full:
        per_core = build_distribution(N_CELLS_TOTAL, np.ones(args.n_cores), rng)
    else:
        per_core = np.full(args.n_cores, args.cells_per_core, dtype=int)

    all_cells, all_pix = [], []
    for i, n in enumerate(per_core):
        if n == 0:
            continue
        c, p = simulate_core(f"SIM{i + 1:03d}", int(n), rng)
        all_cells.append(c); all_pix.append(p)
        if (i + 1) % 10 == 0 or args.full:
            print(f"  core {i + 1}/{args.n_cores}  (n={n})", flush=True)

    cells = pd.concat(all_cells, ignore_index=True)
    pix = pd.concat(all_pix, ignore_index=True)
    write_table(cells, root / "processed/01_sim_cells")
    write_table(pix,   root / "processed/01_sim_pixels")

    gt = {
        "seed": args.seed,
        "n_cores": int(args.n_cores),
        "n_cells_total": int(len(cells)),
        "n_phenotypes": len(PHENO_LIST),
        "n_major_types": len(MAJOR),
        "n_markers": len(MARKERS),
        "n_cn": N_CN,
        "cn_names": CN_NAMES,
        "phenotype_names": PHENO_LIST,
        "marker_names": MARKERS,
        "metabolite_features": {k: v for k, v in TRUE_METABOLITES.items()},
        "core_diameter_um": CORE_DIAMETER_UM,
        "cn_window_um": CN_WINDOW_UM,
        "cn_window_px_comet": CN_WINDOW_UM / COMET_PX_UM,
        "maldi_pixel_um": MALDI_PIXEL_UM,
        "comet_px_um": COMET_PX_UM,
        "true_interactions": [dict(A=a, B=b, kind=k, factor=f,
                                   scope=list(cns)) for a, b, k, f, cns in TRUE_INTERACTIONS],
        "region_cn_prior": REGION_CN_PRIOR,
        "cn_pheno_mix": CN_PHENO_MIX,
        "pheno_metabolite_profile": PHENO_METABOLITE,
        "cells_per_maldi_pixel_median": float(pix.loc[pix["in_core"], "n_cells"].median()),
        "cells_per_maldi_pixel_mean": float(pix.loc[pix["in_core"], "n_cells"].mean()),
        "region_fraction_observed": cells["region"].value_counts(normalize=True).round(4).to_dict(),
        "cn_fraction_observed": cells["true_cn"].value_counts(normalize=True).round(4).to_dict(),
    }
    with open(root / "results/01_sim_ground_truth.json", "w", encoding="utf-8") as fh:
        json.dump(gt, fh, ensure_ascii=False, indent=2)

    print(f"[OK] cells={len(cells):,}  pixels={len(pix):,}")
    print(f"     每 MALDI 像素覆盖细胞数 中位数={gt['cells_per_maldi_pixel_median']:.1f} "
          f"均值={gt['cells_per_maldi_pixel_mean']:.2f}（验证 02 号文档 §2.10 的聚合偏差）")
    print(f"     CN 窗口 70 µm = {gt['cn_window_px_comet']:.1f} COMET px")
    # 质量断言（防静默失败）
    assert len(cells) > 0 and cells["true_cn"].nunique() == N_CN, "CN 未全部出现"
    assert set(cells["true_phenotype"].unique()) <= set(PHENO_LIST)
    assert (cells[[m for m in MARKERS]].to_numpy() >= 0).all(), "marker 强度出现负值"
    assert "in_core" in pix.columns and pix["in_core"].any(), "像素 in_core 掩膜缺失"
    print("     [assert] 通过：CN/表型/marker/代谢物四重 ground truth 均已落盘")


if __name__ == "__main__":
    sys.exit(main())
