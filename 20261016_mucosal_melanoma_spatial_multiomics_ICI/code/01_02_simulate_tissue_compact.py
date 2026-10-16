#!/usr/bin/env python3
"""compact_sim.py — 合成 MuM 组织生成器（最小可运行版，已实测）
参数对齐原文：695444 细胞 / 97 cores / 1.0 mm 芯 / 70 µm CN 窗口 / 20 µm MALDI / 0.23 µm/px。
完整版（含 poisson-disc 硬核点过程与更细的面板表）见 code/01_01_simulate_tissue.py。
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter

SEED = 20260101
N_CORES, N_CELLS_TOTAL = 97, 695_444
CORE_R_UM, CN_WIN_UM, MALDI_UM, COMET_UM = 500.0, 70.0, 20.0, 0.23
MAJOR = ["Tumor", "Tcell", "Macrophage", "DC", "NK", "Bcell", "Stromal", "Other"]
PHENO = ["Tumor_Ki67", "Tumor_MIF", "Tumor_HLADR", "Tumor_MITF", "Tumor_panSOX10",
         "T_CD8", "T_CD8_Ki67", "T_CD8_GZMB", "T_CD4", "T_CD4_FOXP3", "T_DN",
         "Mac_CD163", "Mac_CD163_IRF4", "Mac_CD163_Arg1", "Mac_CD68", "Mac_MHCII",
         "DC_CD11c", "NK", "B", "Stroma_Vim", "Stroma_Vim_SMA", "Stroma_pericyte",
         "Endothelial", "Neuron_Other", "Unassigned"]
PM = {"Tumor": PHENO[0:5], "Tcell": PHENO[5:11], "Macrophage": PHENO[11:16], "DC": ["DC_CD11c"],
      "NK": ["NK"], "Bcell": ["B"], "Stromal": PHENO[19:23], "Other": PHENO[23:25]}
P2M = {p: m for m, v in PM.items() for p in v}
MARKERS = ["SOX10", "MITF", "S100B", "Ki67", "MIF", "HLA_DR", "CD3", "CD4", "CD8", "FOXP3",
           "GZMB", "PD1", "CD163", "CD68", "CD11c", "IRF4", "Arg1", "CD206", "HLA_ABC",
           "NKp46", "CD56", "CD20", "CD79a", "Vimentin", "SMA", "FAP", "CD31", "Desmin",
           "PDL1", "CTLA4", "LAG3", "TIGIT", "TIM3", "CD45", "CD45RO", "CD27",
           "GranzymeA", "CXCL9", "IDO1", "Ki67_stroma"]                       # 40-plex
METS = {"linoleic_acid": 279.2330, "oleic_acid": 281.2486, "palmitic_acid": 255.2330,
        "stearic_acid": 283.2643, "indole_3_acetic_acid": 174.0561,
        "indole_3_carbinol": 146.0611, "indole_3_acetaldehyde": 158.0611,
        "indole_3_carboxylic_acid": 160.0404, "tryptophan": 203.0826, "lactate": 89.0244}
# 区域 → (r_min, r_max)，r 以瘤体半径为 1（瘤体半径 = 0.62 × 芯半径）
REGIONS = [("central_tumor", 0, .45), ("invasive_tumor", .45, .75), ("tumor_boundary", .75, 1.),
           ("peri_stroma", 1., 1.25), ("stroma", 1.25, 9.)]
# 区域 → CN 先验（真值空间结构）
R2CN = {"central_tumor": [("CN15_central_tumor", .62), ("CN14_invasive_tumor", .20), ("CN6_proliferative_tumor", .10), ("CN11_immune_enriched", .04), ("CN8_mixed_cells", .04)],
        "invasive_tumor": [("CN14_invasive_tumor", .50), ("CN6_proliferative_tumor", .18), ("CN13_tumor_boundary", .14), ("CN15_central_tumor", .08), ("CN3_DC", .05), ("CN8_mixed_cells", .05)],
        "tumor_boundary": [("CN13_tumor_boundary", .48), ("CN14_invasive_tumor", .16), ("CN3_DC", .10), ("CN2_NK", .08), ("CN7_CD8T", .08), ("CN11_immune_enriched", .06), ("CN8_mixed_cells", .04)],
        "peri_stroma": [("CN9_tumor_infiltrated_stroma", .30), ("CN12_immune_infiltrated_stroma", .22), ("CN10_Vim_pos_stroma", .20), ("CN1_Vim_pos_SMA_pos_stroma", .14), ("CN5_B_and_CD4T", .08), ("CN4_B_cells", .06)],
        "stroma": [("CN10_Vim_pos_stroma", .36), ("CN1_Vim_pos_SMA_pos_stroma", .26), ("CN12_immune_infiltrated_stroma", .12), ("CN4_B_cells", .10), ("CN5_B_and_CD4T", .08), ("CN2_NK", .04), ("CN3_DC", .04)]}
CN_NAMES = ["CN1_Vim_pos_SMA_pos_stroma", "CN2_NK", "CN3_DC", "CN4_B_cells", "CN5_B_and_CD4T",
            "CN6_proliferative_tumor", "CN7_CD8T", "CN8_mixed_cells", "CN9_tumor_infiltrated_stroma",
            "CN10_Vim_pos_stroma", "CN11_immune_enriched", "CN12_immune_infiltrated_stroma",
            "CN13_tumor_boundary", "CN14_invasive_tumor", "CN15_central_tumor"]
# CN → 表型组成（真值）。已核对：25 个表型中 24 个在此表中至少出现一次，
# 未出现的 "Unassigned" 由下方 2% 兜底填充，保证 25/25 全覆盖。
C2P = {"CN15_central_tumor": {"Tumor_MIF": .58, "Tumor_panSOX10": .18, "Tumor_MITF": .10, "Tumor_HLADR": .04, "Stroma_Vim_SMA": .05, "Mac_MHCII": .03, "Unassigned": .02},
       "CN14_invasive_tumor": {"Tumor_Ki67": .30, "Tumor_MIF": .20, "Tumor_HLADR": .14, "Tumor_panSOX10": .09, "Stroma_Vim": .07, "Stroma_Vim_SMA": .06, "Mac_CD68": .05, "T_CD8": .03, "T_CD8_GZMB": .03, "B": .03},
       "CN6_proliferative_tumor": {"Tumor_Ki67": .72, "Tumor_panSOX10": .14, "Tumor_MIF": .06, "T_CD8": .04, "Unassigned": .04},
       "CN13_tumor_boundary": {"Tumor_Ki67": .20, "Tumor_HLADR": .11, "Mac_CD163": .17, "Mac_CD163_IRF4": .05, "Mac_CD68": .11, "DC_CD11c": .07, "T_CD8": .09, "T_CD8_GZMB": .06, "T_CD4": .05, "Stroma_Vim": .05, "B": .02, "Unassigned": .02},
       "CN2_NK": {"NK": .58, "T_CD8": .13, "Mac_CD68": .09, "DC_CD11c": .06, "Stroma_Vim": .05, "Neuron_Other": .05, "Unassigned": .04},
       "CN3_DC": {"DC_CD11c": .58, "Mac_MHCII": .12, "T_CD4": .10, "T_CD8": .08, "B": .06, "Stroma_Vim": .04, "Unassigned": .02},
       "CN4_B_cells": {"B": .62, "T_CD4": .12, "DC_CD11c": .08, "Stroma_Vim": .10, "Mac_CD68": .04, "Unassigned": .04},
       "CN5_B_and_CD4T": {"B": .36, "T_CD4": .30, "DC_CD11c": .10, "T_CD8": .08, "Stroma_Vim": .08, "Mac_CD163": .04, "Unassigned": .04},
       "CN7_CD8T": {"T_CD8": .42, "T_CD8_GZMB": .20, "T_CD8_Ki67": .12, "T_CD4": .08, "DC_CD11c": .06, "Mac_CD68": .06, "Stroma_Vim": .06},
       "CN8_mixed_cells": {"Tumor_panSOX10": .12, "Stroma_Vim": .18, "Mac_CD68": .12, "T_CD4": .10, "T_CD8": .09, "T_DN": .05, "B": .07, "DC_CD11c": .06, "NK": .04, "Endothelial": .05, "Neuron_Other": .05, "Unassigned": .07},
       "CN9_tumor_infiltrated_stroma": {"Stroma_Vim": .36, "Stroma_Vim_SMA": .08, "Tumor_panSOX10": .18, "Tumor_Ki67": .07, "Mac_CD68": .09, "Mac_CD163_IRF4": .05, "T_CD8": .07, "T_CD4": .05, "T_CD4_FOXP3": .04, "Endothelial": .01},
       "CN10_Vim_pos_stroma": {"Stroma_Vim": .55, "Stroma_Vim_SMA": .12, "Endothelial": .09, "Stroma_pericyte": .07, "Mac_CD68": .05, "Mac_CD163_Arg1": .05, "T_CD4_FOXP3": .04, "Tumor_panSOX10": .02, "Unassigned": .01},
       "CN1_Vim_pos_SMA_pos_stroma": {"Stroma_Vim_SMA": .66, "Stroma_Vim": .14, "Stroma_pericyte": .10, "Mac_CD68": .04, "Endothelial": .04, "Unassigned": .02},
       "CN11_immune_enriched": {"T_CD8": .20, "T_CD4": .17, "T_DN": .06, "B": .13, "NK": .11, "DC_CD11c": .12, "Mac_MHCII": .10, "Mac_CD68": .06, "T_CD4_FOXP3": .03, "Unassigned": .02},
       "CN12_immune_infiltrated_stroma": {"Stroma_Vim": .36, "T_CD8": .13, "T_CD4": .11, "B": .09, "Mac_CD163": .08, "Mac_CD163_Arg1": .05, "DC_CD11c": .06, "Neuron_Other": .04, "Endothelial": .04, "Unassigned": .04}}
# 注入的已知空间富集 / 排斥（下游 08 号文档必须检出这些方向）
TRUE_INT = [("Tumor_Ki67", "Mac_CD163", "ENRICHED", .72, ["CN13_tumor_boundary", "CN14_invasive_tumor"]),
            ("Tumor_Ki67", "DC_CD11c", "ENRICHED", .78, ["CN13_tumor_boundary"]),
            ("T_CD8_GZMB", "Tumor_HLADR", "ENRICHED", .82, ["CN14_invasive_tumor", "CN13_tumor_boundary"]),
            ("T_CD4_FOXP3", "Mac_CD163_Arg1", "ENRICHED", .80, ["CN10_Vim_pos_stroma", "CN12_immune_infiltrated_stroma", "CN9_tumor_infiltrated_stroma"]),
            ("Tumor_MIF", "Stroma_Vim_SMA", "REPELLED", 1.45, ["CN14_invasive_tumor", "CN15_central_tumor"]),
            ("T_CD8", "Stroma_Vim_SMA", "REPELLED", 1.35, ["CN10_Vim_pos_stroma", "CN9_tumor_infiltrated_stroma"]),
            ("B", "Tumor_Ki67", "REPELLED", 1.30, ["CN14_invasive_tumor"])]
# 表型 → 代谢物真值谱（>0 表示该表型贡献该代谢物；倍数越大贡献越高）
P2MET = {"Tumor_Ki67": {"linoleic_acid": 1.4, "oleic_acid": 1.3, "palmitic_acid": 1.5, "stearic_acid": 1.4, "lactate": 2.2},
         "Tumor_MIF": {"linoleic_acid": 1.1, "oleic_acid": 1.0, "palmitic_acid": 1.2, "stearic_acid": 1.2, "indole_3_acetaldehyde": 1.3, "lactate": 1.4},
         "Tumor_HLADR": {"linoleic_acid": 1.0, "oleic_acid": 1.0, "palmitic_acid": 1.1, "stearic_acid": 1.1, "tryptophan": .8},
         "Tumor_MITF": {"linoleic_acid": 1.2, "oleic_acid": 1.2, "palmitic_acid": 1.3, "stearic_acid": 1.3},
         "Tumor_panSOX10": {"linoleic_acid": 1.0, "oleic_acid": 1.0, "palmitic_acid": 1.0, "stearic_acid": 1.0},
         "T_CD8": {"tryptophan": .6}, "T_CD8_Ki67": {"tryptophan": .5, "lactate": 1.3},
         "T_CD8_GZMB": {"tryptophan": .6}, "T_CD4": {"tryptophan": .7},
         "T_CD4_FOXP3": {"indole_3_acetic_acid": .9, "tryptophan": 1.1}, "T_DN": {"tryptophan": .8},
         "Mac_CD163": {"indole_3_acetic_acid": 2.0, "indole_3_acetaldehyde": 1.9, "indole_3_carbinol": 1.6, "tryptophan": 1.6},
         "Mac_CD163_IRF4": {"indole_3_acetic_acid": 1.6, "indole_3_acetaldehyde": 1.5, "tryptophan": 1.3},
         "Mac_CD163_Arg1": {"indole_3_acetic_acid": 2.4, "indole_3_acetaldehyde": 2.2, "indole_3_carbinol": 2.0, "indole_3_carboxylic_acid": 1.7, "tryptophan": 1.9, "linoleic_acid": 1.2},
         "Mac_CD68": {"indole_3_acetic_acid": 1.2, "linoleic_acid": .9},
         "Mac_MHCII": {"indole_3_acetic_acid": 1.4, "indole_3_carbinol": 1.5},
         "DC_CD11c": {"indole_3_carboxylic_acid": 2.1, "indole_3_carbinol": 2.0, "indole_3_acetaldehyde": 1.4, "tryptophan": 1.2},
         "NK": {"lactate": 1.1}, "B": {"tryptophan": .9},
         "Stroma_Vim": {"linoleic_acid": 1.3, "oleic_acid": 1.2, "palmitic_acid": 1.2, "stearic_acid": 1.2},
         "Stroma_Vim_SMA": {"linoleic_acid": 1.6, "oleic_acid": 1.5, "palmitic_acid": 1.4, "stearic_acid": 1.4},
         "Stroma_pericyte": {"linoleic_acid": .9, "oleic_acid": .9},
         "Endothelial": {"lactate": 1.2, "linoleic_acid": .8},
         "Neuron_Other": {"tryptophan": 1.0}, "Unassigned": {}}
MARKER_RULE = [  # (表型集合谓词, marker, log 尺度均值)
    (lambda p: True, "CD45", None), (lambda p: P2M[p] == "Tumor", "SOX10", 1.5),
    (lambda p: P2M[p] == "Tumor", "MITF", 1.5), (lambda p: P2M[p] == "Tumor", "S100B", 1.5),
    (lambda p: p == "Tumor_Ki67", "Ki67", 2.8), (lambda p: p == "Tumor_MIF", "MIF", 2.9),
    (lambda p: p == "Tumor_HLADR", "HLA_DR", 2.6), (lambda p: p == "Tumor_HLADR", "PDL1", 1.8),
    (lambda p: p == "Tumor_MITF", "MITF", 2.7), (lambda p: P2M[p] == "Tcell", "CD3", 2.4),
    (lambda p: p == "T_CD8", "CD8", 2.7), (lambda p: p == "T_CD8_Ki67", "CD8", 2.6),
    (lambda p: p == "T_CD8_Ki67", "Ki67", 2.5), (lambda p: p == "T_CD8_GZMB", "CD8", 2.5),
    (lambda p: p == "T_CD8_GZMB", "GZMB", 2.7), (lambda p: p == "T_CD4", "CD4", 2.6),
    (lambda p: p == "T_CD4_FOXP3", "CD4", 2.4), (lambda p: p == "T_CD4_FOXP3", "FOXP3", 2.5),
    (lambda p: p == "T_CD4_FOXP3", "CTLA4", 2.0), (lambda p: P2M[p] == "Macrophage", "CD68", 2.3),
    (lambda p: p.startswith("Mac_CD163"), "CD163", 2.7), (lambda p: p == "Mac_CD163_IRF4", "IRF4", 1.9),
    (lambda p: p == "Mac_CD163_Arg1", "Arg1", 2.1), (lambda p: p == "Mac_CD163_Arg1", "CD206", 2.0),
    (lambda p: p == "Mac_MHCII", "HLA_DR", 2.4), (lambda p: p == "DC_CD11c", "CD11c", 2.9),
    (lambda p: p == "DC_CD11c", "CXCL9", 1.9), (lambda p: p == "NK", "NKp46", 2.6),
    (lambda p: p == "NK", "CD56", 2.4), (lambda p: p == "B", "CD20", 2.6),
    (lambda p: p == "B", "CD79a", 2.5), (lambda p: P2M[p] == "Stromal", "Vimentin", 2.5),
    (lambda p: p == "Stroma_Vim_SMA", "SMA", 2.8), (lambda p: p == "Stroma_Vim_SMA", "FAP", 2.2),
    (lambda p: p == "Stroma_pericyte", "Desmin", 2.2), (lambda p: p == "Endothelial", "CD31", 2.8)]


def pick(rng, table, n):
    ks = list(table); pr = np.array([table[k] for k in ks], float); pr /= pr.sum()
    return rng.choice(np.array(ks, dtype=object), size=n, p=pr)


def inject(px, py, ph, cn, rng):
    """按 TRUE_INT 注入已知的富集/排斥（移动 B 类细胞坐标，不改变其表型标签）。"""
    for A, B, kind, f, scope in TRUE_INT:
        m = np.isin(cn, scope); ia = np.where(m & (ph == A))[0]; ib = np.where(m & (ph == B))[0]
        if len(ia) < 3 or len(ib) < 3:
            continue
        Apts = np.column_stack([px[ia], py[ia]]); cur = np.column_stack([px[ib], py[ib]])
        d, j = cKDTree(Apts).query(cur, k=1)
        if kind == "ENRICHED":
            nxt = cur + (Apts[j] - cur) * (1 - f)
        else:                                    # 远离：沿反方向推开，限幅 60 µm
            v = (Apts[j] - cur); nrm = np.linalg.norm(v, axis=1, keepdims=True) + 1e-9
            nxt = cur - (v / nrm) * np.minimum((f - 1) * d, 60.0)[:, None]
        px[ib], py[ib] = nxt[:, 0], nxt[:, 1]
    return px, py


def markers(ph, rng, dropout=.35, sd=.55):
    mu = np.full((len(ph), len(MARKERS)), .8, dtype=float)      # 背景
    mi = {m: i for i, m in enumerate(MARKERS)}
    for pred, mk, val in MARKER_RULE:
        if val is None:
            mu[[pred(p) for p in ph], mi[mk]] = 2.0 if mk == "CD45" else mu[[pred(p) for p in ph], mi[mk]]
            continue
        mu[[pred(p) for p in ph], mi[mk]] = val
    for i, p in enumerate(ph):                                   # CD45：仅免疫谱系
        if P2M[p] in ("Tcell", "Macrophage", "DC", "NK", "Bcell"):
            mu[i, mi["CD45"]] = 2.0
    v = np.exp(np.clip(rng.normal(mu, sd), -3, 6)).astype(np.float32)
    return np.where(rng.random(v.shape) > dropout, v, 0.0).astype(np.float32)


def one_core(cid, n, rng, tumor_frac=.62):
    side = 2 * CORE_R_UM; c = CORE_R_UM
    n_s = int(n / (np.pi / 4) * 1.05)                            # 圆内采样补偿
    x = rng.uniform(0, side, n_s); y = rng.uniform(0, side, n_s)
    keep = np.hypot(x - c, y - c) <= CORE_R_UM * .98
    x, y = x[keep][:n], y[keep][:n]; n = len(x)
    r = np.hypot(x - c, y - c) / (tumor_frac * CORE_R_UM)
    reg = np.full(n, "stroma", dtype=object)
    for nm, lo, hi in REGIONS:
        reg[(r >= lo) & (r < hi)] = nm
    cn = np.empty(n, dtype=object)
    for nm in np.unique(reg):
        m = reg == nm; cn[m] = pick(rng, dict(R2CN[nm]), int(m.sum()))
    ph = np.empty(n, dtype=object)
    for cnm in np.unique(cn):
        m = cn == cnm; ph[m] = pick(rng, C2P[cnm], int(m.sum()))
    x, y = inject(x, y, ph, cn, rng)
    x = np.clip(x, 0, side); y = np.clip(y, 0, side)
    mk = markers(ph, rng)
    # 细胞级代谢物真值 → 聚合到 20 µm 像素
    npx = int(np.ceil(side / MALDI_UM)); ix = np.clip((x / MALDI_UM).astype(int), 0, npx - 1)
    iy = np.clip((y / MALDI_UM).astype(int), 0, npx - 1)
    cellm = np.zeros((n, len(METS)), np.float32)
    for k, f in enumerate(METS):
        base = np.array([P2MET.get(p, {}).get(f, 0.) for p in ph], np.float32)
        field = gaussian_filter(rng.normal(0, 1, (npx, npx)), sigma=2.5)[ix, iy]
        cellm[:, k] = np.expm1(np.clip(base + .45 * field + rng.normal(0, .25, n), -3, 4))
    acc = np.zeros((npx, npx, len(METS)), np.float32)
    for k in range(len(METS)):
        np.add.at(acc[:, :, k], (ix, iy), cellm[:, k])
    cnt = np.zeros((npx, npx), np.int32); np.add.at(cnt, (ix, iy), 1)
    cells = pd.DataFrame({"core_id": cid, "cell_id": np.arange(n), "x_um": x.astype(np.float32),
                          "y_um": y.astype(np.float32), "region": reg, "true_cn": cn,
                          "true_phenotype": ph, "true_major": [P2M[p] for p in ph]})
    for j, m in enumerate(MARKERS): cells[m] = mk[:, j]
    for k, f in enumerate(METS): cells[f"met_{f}"] = cellm[:, k]                 # 细胞级真值
    gx = np.arange(npx) * MALDI_UM + MALDI_UM / 2
    px = pd.DataFrame({"core_id": cid, "in_core": (np.hypot(gx[:, None] - c, gx[None, :] - c) <= CORE_R_UM).ravel(),
                       "px_x_um": np.repeat(gx, npx).astype(np.float32),
                       "px_y_um": np.tile(gx, npx).astype(np.float32), "n_cells": cnt.ravel()})
    for k, f in enumerate(METS): px[f] = acc[:, :, k].ravel()
    return cells, px


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out-dir", default=".")
    ap.add_argument("--seed", type=int, default=SEED); ap.add_argument("--n-cores", type=int, default=4)
    ap.add_argument("--cells-per-core", type=int, default=7169, help="7169 → 97 cores 合计 ≈695444")
    ap.add_argument("--full", action="store_true", help="97 cores × 7169 = 695,444 细胞（慢，~5 GB）")
    a = ap.parse_args(); rng = np.random.default_rng(a.seed); root = Path(a.out_dir).resolve()
    (root / "processed").mkdir(parents=True, exist_ok=True); (root / "results").mkdir(parents=True, exist_ok=True)
    n_cores = N_CORES if a.full else a.n_cores
    cs, ps = [], []
    for i in range(n_cores):
        c, p = one_core(f"SIM{i + 1:03d}", a.cells_per_core, rng); cs.append(c); ps.append(p)
    cells = pd.concat(cs, ignore_index=True); pix = pd.concat(ps, ignore_index=True)
    cells.to_csv(root / "processed/01_compactsim_cells.csv.gz", index=False, compression="gzip", float_format="%.6g")
    pix.to_csv(root / "processed/01_compactsim_pixels.csv.gz", index=False, compression="gzip", float_format="%.6g")
    gt = dict(seed=a.seed, n_cores=n_cores, n_cells=int(len(cells)), cn_names=CN_NAMES, phenotype_names=PHENO,
              marker_names=MARKERS, metabolite_features=METS, cn_window_px_comet=CN_WIN_UM / COMET_UM,
              maldi_pixel_um=MALDI_UM, comet_px_um=COMET_UM, core_diameter_um=2 * CORE_R_UM,
              true_interactions=[dict(A=x, B=y, kind=k, factor=f, scope=s) for x, y, k, f, s in TRUE_INT],
              region_prior=R2CN, cn_pheno_mix=C2P, pheno_metabolite_profile=P2MET,
              cells_per_pixel_mean=float(pix.loc[pix.in_core, "n_cells"].mean()),
              cells_per_pixel_median=float(pix.loc[pix.in_core, "n_cells"].median()))
    (root / "results/01_compactsim_ground_truth.json").write_text(json.dumps(gt, ensure_ascii=False, indent=2))
    print(f"[OK] cores={n_cores} cells={len(cells):,} pixels={len(pix):,}")
    print(f"     表型覆盖 {cells.true_phenotype.nunique()}/25；CN 覆盖 {cells.true_cn.nunique()}/15")
    print(f"     20 µm 像素平均覆盖 {gt['cells_per_pixel_mean']:.2f} 个细胞（中位 {gt['cells_per_pixel_median']:.0f}）")
    print(f"     CN 窗口 70 µm = {gt['cn_window_px_comet']:.1f} COMET px")
    assert cells.true_phenotype.nunique() == 25 and cells.true_cn.nunique() == 15
    assert (cells[MARKERS].to_numpy() >= 0).all() and pix.in_core.any()


if __name__ == "__main__":
    main()
