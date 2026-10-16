#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
spatial_interaction_test.py —— 细胞-细胞空间相互作用的邻接-置换检验
================================================================================
【用途】
  读入一张细胞表（每行一个细胞：x, y, celltype，可含分组列如 fov/sample），
  用 **Delaunay 三角剖分** 建立空间邻接 → 距离阈值过滤 → 对每一对细胞类型统计
  邻接边数与边距离 → **在每个分组内做位置/标签置换**构造零分布 → 用
  **Mann–Whitney U** 比较观测距离分布与零分布 → 输出 logFC 与 BH-FDR。

  本脚本**与具体文献/数据无关**：不硬编码任何项目路径、队列名、细胞类型名，
  全部走参数（见 REPRO_SPEC.md §1.4 脚本分层）。

【方法学出处】
  after Wang et al. 2026 bioRxiv, Methods「Cell-cell interaction analysis」：
    Delaunay 三角剖分 → 距离阈值 → 100 次位置置换 → Mann–Whitney U → logFC。

【符号约定（务必注意，已写入 --help）】
  logFC = log2( mean(观测边距离) / mean(置换边距离) )
      logFC < 0  →  观测距离比置换更近  →  **该细胞对空间富集（enriched / co-localized）**
      logFC > 0  →  观测距离比置换更远  →  **该细胞对空间排斥（depleted / segregated）**
      logFC ≈ 0  →  与零模型无异

【零模型】（置换**只在分组内进行**，不同 fov/sample 之间不混）
  position_shuffle（默认）：分组内随机重排细胞坐标后**重建三角剖分**。
      注意：Delaunay 只依赖点集，而点集在坐标重排下不变，因此该零模型在统计上
      与 label_shuffle **等价**（两者都随机化"细胞↔几何位置"的对应）；保留为默认
      以复现原文描述的"位置置换"。
  label_shuffle：分组内随机重排细胞类型标签（几何与三角剖分不变，更快）。
  toroidal：分组内随机挑选一个细胞类型，把它的坐标做环形平移（wrap-around，
      越界按包围盒回卷）；**保持该类型内部的空间结构**，只随机化类型间相对位置 ——
      与前两者**不同**，用于检验"跨类型共定位"是否超出各自内部聚集所能解释。

【p 值与置换分辨率（方法学要点）】
  p 来自 Mann–Whitney U（大样本正态近似），其下界远低于 1/n_perm；因此当
  p < 1/n_perm 时脚本会打印明确警告：**该量级的 p 反映的是距离分布的大样本近似，
  而不是置换分辨率**。增加 --n-perm 不能降低这个下限。q 值为全体细胞对上的 BH-FDR。

【用法示例】
  python3 tools/spatial_interaction_test.py \\
      --cells processed/cells.csv --group-col fov \\
      --max-dist 30 --dist-unit um --pixel-size-um 0.5 \\
      --null-model position_shuffle --n-perm 100 --out-prefix results/cci

  python3 tools/spatial_interaction_test.py --cells processed/cells.tsv --sep '\\t' \\
      --group-col sample --max-dist 100 --dist-unit pixels --n-perm 200 --out-prefix results/cci

  python3 tools/spatial_interaction_test.py --selftest

【输出】
  <prefix>_interaction_stats.csv
      pair, type_a, type_b, n_obs, mean_dist_obs, mean_dist_null, logFC,
      p_mwu, q_bh, n_perm, mean_n_null, n_groups_obs
      （按 pair 字符串升序 → 输出行序稳定）
  <prefix>_run_meta.csv
      运行元信息（null_model、n_perm、距离单位与换算、分组数、警告项、logFC 符号约定、
      `seconds` 运行耗时）。**`seconds` 每次不同，故该文件不是逐字节可复现的**；
      `_interaction_stats.csv` 的全部统计量在固定 --seed 下完全确定（行序亦固定）。

【依赖】
  必需：numpy, pandas, scipy, scikit-learn

【退出码】
  0 = 成功
  1 = 未预期异常
  2 = 参数/输入校验失败
  3 = 质量断言失败（产物为空、无有效边等）
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
        raise UsageError("输入表缺少列 %s：%s\n        实际列：%s"
                         % (missing, path, list(df.columns)[:30]))


def numeric_coords(df, x_col, y_col):
    import numpy as np
    import pandas as pd

    out = []
    for col in (x_col, y_col):
        v = df[col]
        if v.dtype.kind not in "iuf":
            v = pd.to_numeric(v, errors="coerce")
        bad = v.isna().to_numpy()
        if bad.any():
            raise UsageError("列 '%s' 含 %d 个非数值/NA（示例行号：%s）—— 坐标必须为数值"
                             % (col, int(bad.sum()), np.flatnonzero(bad)[:5].tolist()))
        out.append(v.to_numpy(dtype=float))
    return out[0], out[1]


# ----------------------------------------------------------------------------
# 邻接边（Delaunay + 距离阈值）
# ----------------------------------------------------------------------------
def delaunay_edges(coords, max_dist=None):
    """返回 (edges[m,2] 去重无向边（已排序），dists[m])，仅保留长度 <= max_dist 的边。"""
    import numpy as np
    from scipy.spatial import Delaunay

    n = coords.shape[0]
    if n < 3:
        return np.zeros((0, 2), dtype=int), np.zeros(0, dtype=float)
    try:
        tri = Delaunay(coords)
    except Exception as exc:
        raise QualityError(
            "Delaunay 三角剖分失败（%s）；请检查是否有大量重复坐标或全部点共线" % exc
        )
    s = tri.simplices
    e = np.vstack([s[:, [0, 1]], s[:, [1, 2]], s[:, [2, 0]]])
    e = np.unique(np.sort(e, axis=1), axis=0)
    d = np.linalg.norm(coords[e[:, 0]] - coords[e[:, 1]], axis=1)
    if max_dist is not None:
        keep = d <= float(max_dist)
        e, d = e[keep], d[keep]
    return e, d


def edge_pair_stats(codes_a, codes_b, dists, all_pairs):
    """把边按无序类型对归类。

    返回 { (i, j): {"n": int, "sum": float, "dists": np.ndarray} }，覆盖 all_pairs 全部项。
    """
    import numpy as np

    out = {}
    if len(dists) == 0:
        for pk in all_pairs:
            out[pk] = {"n": 0, "sum": 0.0, "dists": np.zeros(0)}
        return out
    lo = np.minimum(codes_a, codes_b).astype(np.int64)
    hi = np.maximum(codes_a, codes_b).astype(np.int64)
    flat = lo * 100000 + hi
    for (i, j) in all_pairs:
        m = flat == (i * 100000 + j)
        c = int(m.sum())
        out[(i, j)] = {
            "n": c,
            "sum": float(dists[m].sum()) if c else 0.0,
            "dists": dists[m] if c else np.zeros(0),
        }
    return out


# ----------------------------------------------------------------------------
# 置换检验主流程
# ----------------------------------------------------------------------------
def permutation_edges(P, C, null_model, rng, max_dist_px, n_types, e_obs, d_obs):
    """生成一次置换下的 (边端类型编码, 边距离)。

    e_obs/d_obs 为观测三角剖分的边与距离（供 label_shuffle 复用，无需重建）。
    """
    import numpy as np

    if null_model == "label_shuffle":
        lab = rng.permutation(C)
        return lab[e_obs[:, 0]], lab[e_obs[:, 1]], d_obs

    if null_model == "position_shuffle":
        # 分组内随机重排坐标 → 重建三角剖分（Delaunay 只依赖点集）
        perm = rng.permutation(P.shape[0])
        e, d = delaunay_edges(P[perm], max_dist=max_dist_px)
        return C[e[:, 0]], C[e[:, 1]], d

    if null_model == "toroidal":
        shift_type = int(rng.integers(n_types))
        sel = C == shift_type
        if not sel.any():
            return C[e_obs[:, 0]], C[e_obs[:, 1]], d_obs
        lo = P.min(axis=0)
        span = P.max(axis=0) - lo
        span = np.where(span > 0, span, 1.0)
        delta = rng.uniform(0.0, 1.0, size=2) * span
        Q = P.copy()
        Q[sel] = lo + np.mod(P[sel] - lo + delta, span)
        e, d = delaunay_edges(Q, max_dist=max_dist_px)
        return C[e[:, 0]], C[e[:, 1]], d

    raise UsageError("未知 --null-model：%s" % null_model)


def run_interaction(df, x_col, y_col, ct_col, group_col, max_dist_px, null_model,
                    n_perm, seed):
    """主流程。返回 (stats DataFrame, types list, group_order list, meta dict)。"""
    import numpy as np
    import pandas as pd
    from scipy.stats import mannwhitneyu

    x, y = numeric_coords(df, x_col, y_col)
    coords_all = np.column_stack([x, y])
    ct = df[ct_col].astype(str)
    types = sorted(ct.unique().tolist())
    n_types = len(types)
    tcode = {t: i for i, t in enumerate(types)}
    codes_all = ct.map(tcode).to_numpy(dtype=int)

    if group_col is None:
        groups_all = np.full(len(df), "__all__", dtype=object)
    else:
        groups_all = df[group_col].astype(str).to_numpy()
    group_order = sorted(pd.unique(groups_all).tolist())

    all_pairs = [(i, j) for i in range(n_types) for j in range(i, n_types)]

    obs = {pk: {"n": 0, "sum": 0.0, "dists": np.zeros(0)} for pk in all_pairs}
    obs_groups = {pk: set() for pk in all_pairs}
    null_n = {pk: [] for pk in all_pairs}
    null_dists = {pk: [] for pk in all_pairs}
    skipped_groups = []
    n_delaunay_builds = 0

    rng = np.random.default_rng(seed)

    for g in group_order:
        m = groups_all == g
        n_g = int(m.sum())
        if n_g < 3:
            skipped_groups.append((g, n_g))
            continue
        P = coords_all[m]
        C = codes_all[m]

        e0, d0 = delaunay_edges(P, max_dist=max_dist_px)
        n_delaunay_builds += 1

        st = edge_pair_stats(C[e0[:, 0]], C[e0[:, 1]], d0, all_pairs)
        for pk, s in st.items():
            if s["n"]:
                obs[pk]["n"] += s["n"]
                obs[pk]["sum"] += s["sum"]
                obs[pk]["dists"] = np.concatenate([obs[pk]["dists"], s["dists"]])
                obs_groups[pk].add(g)

        for _ in range(n_perm):
            ca, cb, dp = permutation_edges(P, C, null_model, rng, max_dist_px,
                                           n_types, e0, d0)
            n_delaunay_builds += 0 if null_model == "label_shuffle" else 1
            stn = edge_pair_stats(ca, cb, dp, all_pairs)
            for pk, s in stn.items():
                null_n[pk].append(s["n"])
                if s["n"]:
                    null_dists[pk].append(s["dists"])

    rows = []
    n_p_lt = 0
    for (i, j) in all_pairs:
        pk = (i, j)
        n_obs = obs[pk]["n"]
        mean_obs = obs[pk]["sum"] / n_obs if n_obs else float("nan")
        nn = np.asarray(null_n[pk], dtype=float)
        mean_n_null = float(nn.mean()) if nn.size else float("nan")
        pool = np.concatenate(null_dists[pk]) if null_dists[pk] else np.zeros(0)
        mean_null = float(pool.mean()) if pool.size else float("nan")

        logfc, p = float("nan"), float("nan")
        if n_obs >= 2 and pool.size >= 2:
            logfc = float(np.log2(mean_obs / mean_null)) if mean_null > 0 else float("nan")
            try:
                p = float(mannwhitneyu(obs[pk]["dists"], pool,
                                       alternative="two-sided").pvalue)
            except ValueError:
                p = float("nan")
        if p == p and p < 1.0 / n_perm:
            n_p_lt += 1

        rows.append({
            "pair": "%s|%s" % (types[i], types[j]),
            "type_a": types[i], "type_b": types[j],
            "n_obs": int(n_obs),
            "mean_dist_obs": mean_obs,
            "mean_dist_null": mean_null,
            "logFC": logfc,
            "p_mwu": p,
            "q_bh": float("nan"),
            "n_perm": int(n_perm),
            "mean_n_null": mean_n_null,
            "n_groups_obs": int(len(obs_groups[pk])),
        })

    stats = pd.DataFrame(rows)
    stats = bh_fdr(stats, "p_mwu", "q_bh")
    stats = stats.sort_values(["pair"], kind="stable").reset_index(drop=True)

    meta = {
        "n_types": n_types,
        "group_order": group_order,
        "skipped_groups": skipped_groups,
        "n_pairs_p_lt_1_over_nperm": n_p_lt,
        "n_delaunay_builds": n_delaunay_builds,
        "n_pairs_with_obs_edges": int((stats["n_obs"] > 0).sum()),
    }
    return stats, types, group_order, meta


def bh_fdr(df, p_col, q_col):
    """Benjamini–Hochberg FDR；p 为 NaN 的行 q 保持 NaN（不参与校正）。"""
    import numpy as np

    p = df[p_col].to_numpy(dtype=float)
    q = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    m = int(ok.sum())
    if m:
        idx = np.flatnonzero(ok)
        order = idx[np.argsort(p[ok], kind="stable")]
        ranked = p[order]
        adj = ranked * m / np.arange(1, m + 1)
        adj = np.minimum.accumulate(adj[::-1])[::-1]
        q[order] = np.clip(adj, 0.0, 1.0)
    df[q_col] = q
    return df


# ----------------------------------------------------------------------------
# 自检：功效 / 假阳性率
# ----------------------------------------------------------------------------
def make_synthetic_interaction(seed, enrich=True, n_blobs=3, n_pairs=50, sigma=30.0,
                               pair_dist=(1.5, 5.0), ab_frac=0.5):
    """合成组织（每个 fov）。两种模式：

    几何（两种模式相同）：在 n_blobs 个团块内生成"紧邻点对"（间距 1.5~5 px），
    即点过程本身是聚集 + 成对的，**所有细胞共享同一套几何**。

    enrich=True（**功效用**）：随机一半紧邻点对被赋予 A/B 标签 → A 与 B 始终互为
        紧邻配对（真实空间富集），其边距离显著短于置换零分布。

    enrich=False（**假阳性率用，真正的零假设数据集**）：四种标签在**全部细胞上完全
        随机**分配，与几何无关 → 不存在任何真实的空间富集，任何显著结果都是假阳性。

    为什么 FPR 必须在**单独的、无植入结构**的数据集上估计：本脚本默认的
    标签/位置置换零模型，其零分布是"随机标签下的**通用细胞对**"的距离分布，
    因此每个细胞对都是与"整体点过程"比较。在 enrich=True 的组织里，A-B 被刻意
    拉近后，其余细胞对相对整体点过程不再"可交换"（exchangeable），它们被判显著
    反映的是"与整体点过程的差异"，**不是 H0 意义上的假阳性**。详见 --selftest 输出。

    返回 DataFrame[x, y, celltype, fov, true_pair_id]
    """
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(seed)
    frames = []
    for fov in (1, 2):
        centers = [(120.0, 120.0), (330.0, 200.0), (200.0, 380.0)][:n_blobs]
        xs, ys, pair_id = [], [], []
        pid = 0
        for (bx, by) in centers:
            cx = bx + rng.normal(0.0, sigma, n_pairs)
            cy = by + rng.normal(0.0, sigma, n_pairs)
            ang = rng.uniform(0.0, 2 * np.pi, n_pairs)
            dist = rng.uniform(pair_dist[0], pair_dist[1], n_pairs)
            for k in range(n_pairs):
                xs.append(cx[k]); ys.append(cy[k]); pair_id.append(pid)
                xs.append(cx[k] + dist[k] * np.cos(ang[k]))
                ys.append(cy[k] + dist[k] * np.sin(ang[k]))
                pair_id.append(pid)
                pid += 1
        n = len(xs)
        labels = np.empty(n, dtype=object)
        if enrich:
            is_ab = rng.random(pid) < ab_frac
            rest = []
            for k in range(pid):
                i0, i1 = 2 * k, 2 * k + 1
                if is_ab[k]:
                    labels[i0], labels[i1] = "A", "B"
                else:
                    rest.extend([i0, i1])
            rest = np.asarray(rest, dtype=int)
            if rest.size >= 2:
                shuffled = rng.permutation(rest)
                half = shuffled.size // 2
                labels[shuffled[:half]] = "C"
                labels[shuffled[half:]] = "D"
        else:
            labels = rng.choice(["A", "B", "C", "D"], size=n)
        frames.append(pd.DataFrame({"x": xs, "y": ys, "celltype": labels,
                                    "fov": "fov%d" % fov,
                                    "true_pair_id": pair_id}))
    return pd.concat(frames, ignore_index=True)


def selftest(seed=20260101, n_rep=12, n_perm=100):
    import numpy as np

    print("=" * 78)
    print("spatial_interaction_test.py --selftest")
    print("=" * 78)
    print("[几何] 每个 fov：3 个团块 × 50 个紧邻点对（间距 1.5~5 px）——所有细胞共享同一套几何")
    print("[数据集 1：enrich=True] 随机一半点对赋予 A/B → A-B 真实空间富集 → 用于估计**功效**")
    print("[数据集 2：enrich=False] 四种标签在全部细胞上完全随机分配（真 H0）→ 用于估计**假阳性率**")
    print("[判阳] q_bh < 0.05 且 logFC < 0（logFC<0 = 观测距离更近 = 空间富集）")
    print("[参数] 复本数=%d，每复本 --n-perm=%d，--group-col fov，--max-dist 60 pixels"
          % (n_rep, n_perm))

    # ---- 数据集 1：功效 -----------------------------------------------------
    print("\n" + "-" * 78)
    print("[数据集 1] A/B 紧邻配对（真实富集）")
    det_ab, det_cd_enriched, rows = 0, 0, []
    t0 = time.time()
    for rep in range(n_rep):
        df = make_synthetic_interaction(seed + rep, enrich=True)
        stats, _t, _g, _m = run_interaction(
            df, "x", "y", "celltype", "fov", max_dist_px=60.0,
            null_model="position_shuffle", n_perm=n_perm, seed=seed + rep)
        s_ = stats.set_index("pair")
        ab, cd = s_.loc["A|B"], s_.loc["C|D"]
        hit_ab = bool((ab["q_bh"] < 0.05) and (ab["logFC"] < 0))
        hit_cd = bool((cd["q_bh"] < 0.05) and (cd["logFC"] < 0))
        det_ab += int(hit_ab)
        det_cd_enriched += int(hit_cd)
        rows.append(dict(rep=rep, ab_logFC=ab["logFC"], ab_p=ab["p_mwu"], ab_q=ab["q_bh"],
                         ab_n=int(ab["n_obs"]), ab_detected=hit_ab,
                         cd_logFC=cd["logFC"], cd_q=cd["q_bh"], cd_n=int(cd["n_obs"]),
                         cd_detected=hit_cd,
                         n_sig_pairs=int((stats["q_bh"] < 0.05).sum()),
                         n_pairs=int(stats.shape[0])))
        print("  复本 %2d | A|B logFC=%+.3f p=%.3g q=%.3g n_obs=%3d %s | "
              "显著对 %d/%d；最具富集的 3 对：%s"
              % (rep, ab["logFC"], ab["p_mwu"], ab["q_bh"], int(ab["n_obs"]),
                 "检出" if hit_ab else "未检出",
                 int((stats["q_bh"] < 0.05).sum()), int(stats.shape[0]),
                 ", ".join(stats.sort_values("logFC")["pair"].head(3).tolist())))
    dt1 = time.time() - t0
    power = det_ab / float(n_rep)

    # ---- 数据集 2：假阳性率 -------------------------------------------------
    print("\n" + "-" * 78)
    print("[数据集 2] 标签完全随机（真 H0，无任何真实富集）→ 任何显著结果都是假阳性")
    fp_cd, fp_all, n_tested = 0, 0, 0
    t0 = time.time()
    for rep in range(n_rep):
        df = make_synthetic_interaction(seed + 5000 + rep, enrich=False)
        stats, _t, _g, _m = run_interaction(
            df, "x", "y", "celltype", "fov", max_dist_px=60.0,
            null_model="position_shuffle", n_perm=n_perm, seed=seed + 5000 + rep)
        cd = stats.set_index("pair").loc["C|D"]
        hit = bool((cd["q_bh"] < 0.05) and (cd["logFC"] < 0))
        sig = int((stats["q_bh"] < 0.05).sum())
        fp_cd += int(hit)
        fp_all += sig
        n_tested += int(stats.shape[0])
        rows.append(dict(rep=rep, cd_logFC=cd["logFC"], cd_q=cd["q_bh"], cd_detected=hit))
        print("  复本 %2d | C|D logFC=%+.3f p=%.3g q=%.3g n_obs=%3d %s | 显著对 %d/%d"
              % (rep, cd["logFC"], cd["p_mwu"], cd["q_bh"], int(cd["n_obs"]),
                 "检出错判" if hit else "未检出", sig, int(stats.shape[0])))
    dt2 = time.time() - t0
    fpr_cd = fp_cd / float(n_rep)
    fpr_pair = fp_all / float(n_tested) if n_tested else float("nan")

    ab_fc = float(np.mean([r["ab_logFC"] for r in rows[:n_rep]]))
    print("\n[结果] 数据集 1（真实富集）用时 %.2f s" % dt1)
    print("[结果] 功效 power（A|B 检出率）= %d/%d = %.3f" % (det_ab, n_rep, power))
    print("[结果] A|B logFC 均值 = %+.3f（应显著为负 = 富集）" % ab_fc)
    print("[结果] 数据集 2（真 H0）用时 %.2f s" % dt2)
    print("[结果] 假阳性率 FPR（C|D 误检率）= %d/%d = %.3f" % (fp_cd, n_rep, fpr_cd))
    print("[结果] 逐细胞对假阳性率 = %d/%d = %.3f" % (fp_all, n_tested, fpr_pair))
    print("[说明] 数据集 1 中 C|D 的检出率 = %d/%d：该组织里 A-B 被刻意拉近，其余细胞对"
          "相对整体点过程不再" % (det_cd_enriched, n_rep))
    print("       可交换，故其在数据集 1 中被检出反映的是与整体点过程的差异，"
          "**不是 H0 假阳性**；")
    print("       这正是本零模型（标签/位置置换 = 随机标签下的通用细胞对）的性质："
          "解释 logFC 必须结合整体点过程，")
    print("       不能把它当作 CSR（完全随机）检验。因此 FPR 只在数据集 2（无植入结构）估计。")
    min_p = float(np.nanmin([r["ab_p"] for r in rows[:n_rep]]))
    if min_p < 1.0 / n_perm:
        print("[方法学要点] 数据集 1 最小 p = %.3g < 1/n_perm = %.3g → 该量级 p 由 "
              "Mann–Whitney U 的大样本近似给出，**不是置换分辨率**；" % (min_p, 1.0 / n_perm))
        print("             正式分析须在报告中写明。")

    failed = []
    if not (power > 0.8):
        failed.append("功效 %.3f <= 0.8" % power)
    if not (fpr_cd < 0.1):
        failed.append("假阳性率 %.3f >= 0.1" % fpr_cd)
    if not (fpr_pair < 0.1):
        failed.append("逐细胞对假阳性率 %.3f >= 0.1" % fpr_pair)
    if ab_fc >= 0:
        failed.append("A|B logFC 均值 %.3f 不为负（符号约定异常）" % ab_fc)
    if failed:
        print("\n[FAIL] selftest 未通过：" + "；".join(failed))
        return EXIT_SELFTEST
    print("\n[PASS] selftest 通过：功效 %.3f > 0.8；真 H0 数据集上 FPR %.3f < 0.1"
          "（逐细胞对 %.3f）" % (power, fpr_cd, fpr_pair))
    return EXIT_OK


# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------
def build_parser():
    p = argparse.ArgumentParser(
        prog="spatial_interaction_test.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "细胞-细胞空间相互作用：Delaunay 邻接 → 距离阈值 → 分组内置换 → "
            "Mann–Whitney U → logFC。\n"
            "after Wang et al. 2026 bioRxiv, Methods「Cell-cell interaction analysis」。"
        ),
        epilog="""\
符号约定（重要）
  logFC = log2( mean(观测边距离) / mean(置换边距离) )
    logFC < 0 → 观测更近 → 该细胞对**空间富集 / 共定位（enriched）**
    logFC > 0 → 观测更远 → 该细胞对**空间排斥（depleted / segregated）**
    logFC ≈ 0 → 与零模型无异
  p_mwu：Mann–Whitney U（two-sided），观测边距离分布 vs 置换边距离池。
  q_bh ：全体细胞对上的 Benjamini–Hochberg FDR。
  n_perm：置换次数（写入输出表）。当 p_mwu < 1/n_perm 时打印警告：该量级的 p 来自
          U 检验的大样本近似，**不代表置换分辨率**（增加 n_perm 不会降低该下限）。

单位
  --max-dist 的数值单位由 --dist-unit 决定（pixels 或 um）。--dist-unit um 时**必须**
  提供 --pixel-size-um（μm/px）做换算（REPRO_SPEC 附录 C「坐标错配」高危项）。

零模型
  position_shuffle（默认）分组内重排坐标并重建三角剖分；label_shuffle 分组内重排
  细胞类型标签（更快，与前者统计等价）；toroidal 分组内对随机一个类型做环形平移
  （保持类型内部结构，与前两者不同）。置换只在分组内进行。

性能
  每次置换都要重建三角剖分（label_shuffle 除外），故耗时约为
  `n_perm × 分组数` 次 Delaunay；运行时会把"实际构建次数"打印出来。经验值：单分组
  2 万细胞 × 100 置换约 13 s（本机），可按此估算规模。

示例
  python3 tools/spatial_interaction_test.py --cells processed/cells.csv --group-col fov \\
      --max-dist 30 --dist-unit um --pixel-size-um 0.5 --n-perm 100 --out-prefix results/cci
  python3 tools/spatial_interaction_test.py --selftest

退出码：0=成功 1=未预期异常 2=参数/输入校验失败 3=质量断言失败 4=selftest 失败 5=缺必需依赖
""",
    )
    io = p.add_argument_group("输入输出")
    io.add_argument("--cells", metavar="FILE", help="细胞表（CSV/TSV，每行一个细胞）")
    io.add_argument("--sep", default=None, help="列分隔符（默认按扩展名）")
    io.add_argument("--x-col", default="x", help="x 坐标列名（默认 x）")
    io.add_argument("--y-col", default="y", help="y 坐标列名（默认 y）")
    io.add_argument("--celltype-col", default="celltype", help="细胞类型列名（默认 celltype）")
    io.add_argument("--group-col", default="fov",
                    help="分组列名（默认 fov；常用 sample/fov/roi）。置为 '' 表示不分组"
                         "（全部细胞视为一组，置换在全体上进行）")
    io.add_argument("--out-prefix", default="interaction_out", metavar="PREFIX",
                    help="输出前缀（<PREFIX>_interaction_stats.csv、<PREFIX>_run_meta.csv）")

    d = p.add_argument_group("邻接与距离")
    d.add_argument("--max-dist", type=float, default=None, metavar="D",
                   help="邻接边距离阈值（单位见 --dist-unit）；默认 None = 不设阈值，"
                        "使用全部 Delaunay 边")
    d.add_argument("--dist-unit", choices=["pixels", "um"], default="pixels",
                   help="--max-dist 的单位（默认 pixels）")
    d.add_argument("--pixel-size-um", type=float, default=None,
                   help="像素物理边长 μm/px（--dist-unit um 时必需）")

    w = p.add_argument_group("零模型与置换")
    w.add_argument("--null-model", choices=["position_shuffle", "label_shuffle", "toroidal"],
                   default="position_shuffle",
                   help="零模型（默认 position_shuffle，复现原文的位置置换）；"
                        "详见 --help 顶部说明")
    w.add_argument("--n-perm", type=int, default=100, help="置换次数（默认 100，复现原文）")
    w.add_argument("--seed", type=int, default=20260101, help="随机种子（默认 20260101）")

    p.add_argument("--selftest", action="store_true",
                   help="合成组织自检：报告功效（A-B 刻意富集）与假阳性率（C-D 零假设）"
                        "（断言 power > 0.8 且 FPR < 0.1）")
    return p


def validate(args):
    if args.selftest:
        return
    if not args.cells:
        raise UsageError("缺少 --cells（或使用 --selftest）")
    if args.x_col == args.y_col:
        raise UsageError("--x-col 与 --y-col 不能相同（%s）" % args.x_col)
    if args.max_dist is not None and not (args.max_dist > 0):
        raise UsageError("--max-dist 必须为正数，当前 %r" % args.max_dist)
    if args.dist_unit == "um" and args.max_dist is not None:
        if args.pixel_size_um is None:
            raise UsageError("--dist-unit um 需要 --pixel-size-um（μm/px）换算；缺参数会把 μm "
                             "当像素用（REPRO_SPEC 附录 C 尺度错配）")
        if not (args.pixel_size_um > 0):
            raise UsageError("--pixel-size-um 必须为正数，当前 %r" % args.pixel_size_um)
    if args.n_perm < 1:
        raise UsageError("--n-perm 必须 >= 1，当前 %d" % args.n_perm)


def main(argv=None):
    args = build_parser().parse_args(argv)
    _require_deps()

    if args.selftest:
        return selftest()
    try:
        validate(args)
        import pandas as pd

        group_col = args.group_col if args.group_col not in (None, "") else None
        df = read_table(args.cells, args.sep)
        check_columns(df, [args.x_col, args.y_col, args.celltype_col], args.cells)
        if group_col is not None:
            check_columns(df, [group_col], args.cells)

        max_dist_px = args.max_dist
        if max_dist_px is not None and args.dist_unit == "um":
            max_dist_px = args.max_dist / float(args.pixel_size_um)

        t0 = time.time()
        stats, types, group_order, meta = run_interaction(
            df, args.x_col, args.y_col, args.celltype_col, group_col, max_dist_px,
            args.null_model, args.n_perm, args.seed)
        dt = time.time() - t0

        if stats.shape[0] == 0:
            raise QualityError("没有任何细胞对统计量（细胞类型数=%d）" % len(types))
        if int(stats["n_obs"].sum()) == 0:
            raise QualityError(
                "所有细胞对的观测边数都是 0 —— 距离阈值过小、分组内细胞过少或坐标尺度可疑"
                "（REPRO_SPEC 附录 C 静默失败）")
        for g, n_g in meta["skipped_groups"]:
            print("[WARN] 分组 '%s' 只有 %d 个细胞（<3），已跳过" % (g, n_g))

        outdir = os.path.dirname(os.path.abspath(args.out_prefix))
        os.makedirs(outdir, exist_ok=True)
        stats_path = "%s_interaction_stats.csv" % args.out_prefix
        meta_path = "%s_run_meta.csv" % args.out_prefix
        stats.to_csv(stats_path, index=False)

        n_lt = meta["n_pairs_p_lt_1_over_nperm"]
        nsig = int((stats["q_bh"] < 0.05).sum())
        warns = []
        if n_lt:
            warns.append("%d 个细胞对 p_mwu < 1/n_perm（大样本近似，非置换分辨率）" % n_lt)
        if nsig == 0:
            warns.append("BH-FDR < 0.05 的细胞对为 0")
        pd.DataFrame([{
            "n_cells": len(df), "n_celltypes": len(types), "n_pairs": stats.shape[0],
            "celltypes": ";".join(types),
            "n_groups": len(group_order), "group_col": group_col or "__all__",
            "groups": ";".join(group_order[:50]),
            "null_model": args.null_model, "n_perm": args.n_perm, "seed": args.seed,
            "max_dist_input": args.max_dist, "dist_unit": args.dist_unit,
            "pixel_size_um": args.pixel_size_um,
            "max_dist_px": max_dist_px if max_dist_px is not None else "Inf",
            "n_pairs_with_obs_edges": meta["n_pairs_with_obs_edges"],
            "n_pairs_p_lt_1_over_nperm": n_lt, "n_pairs_q_lt_0.05": nsig,
            "logFC_convention": "log2(mean_dist_obs/mean_dist_null); <0 = enriched",
            "seconds": round(dt, 3), "warnings": " ; ".join(warns),
        }]).to_csv(meta_path, index=False)

        for pp in (stats_path, meta_path):
            if (not os.path.exists(pp)) or os.path.getsize(pp) == 0:
                raise QualityError("产物为空：%s（REPRO_SPEC 附录 C 静默失败）" % pp)

        print("[INFO] 细胞数=%d，细胞类型=%s，分组数=%d（%s），细胞对=%d"
              % (len(df), types, len(group_order), ",".join(group_order[:8]),
                 stats.shape[0]))
        print("[INFO] 零模型=%s，n_perm=%d，seed=%d，max_dist=%s %s（=%.6g px），Delaunay 构建 %d 次"
              % (args.null_model, args.n_perm, args.seed, args.max_dist, args.dist_unit,
                 max_dist_px if max_dist_px is not None else float("inf"),
                 meta["n_delaunay_builds"]))
        print("[INFO] 用时 %.2f s；有观测边的细胞对 %d/%d"
              % (dt, meta["n_pairs_with_obs_edges"], stats.shape[0]))
        print("[INFO] 统计表（按 pair 排序）：")
        for line in stats[["pair", "n_obs", "mean_dist_obs", "mean_dist_null",
                           "logFC", "p_mwu", "q_bh"]].to_string(index=False).splitlines():
            print("       " + line)
        if n_lt:
            print("[WARN] %d 个细胞对的 p_mwu < 1/n_perm = %.3g。" % (n_lt, 1.0 / args.n_perm))
            print("[WARN] 该量级的 p 由 Mann–Whitney U 大样本近似给出，**不是置换分辨率**；")
            print("[WARN] 报告中必须说明（增加 --n-perm 不能降低该下限）。")
        if nsig == 0:
            print("[WARN] BH-FDR < 0.05 的细胞对为 0，请核对距离阈值/分组是否符合预期。")
        print("[OK] 写出：%s\n[OK] 写出：%s" % (stats_path, meta_path))
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
