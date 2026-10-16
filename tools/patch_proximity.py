#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
patch_proximity.py —— Patch（细胞团块）识别 → 边缘细胞 → 环带内异类细胞计数/密度
================================================================================
【用途】
  从一张细胞表中按**指定通道**取出 patch（例如某个 CN、某个细胞类型）→ 用 HDBSCAN
  把该通道细胞聚成若干 patch → 提取每个 patch 的**边缘细胞**（凹包边界附近）→ 把
  边缘细胞与其 k 个最近邻合并为 **anchor 集合** → 统计距 anchor 不同**半径环带**内
  **不属于该 patch** 的细胞（按细胞类型分别计数，并给出**密度 = 计数 / 环带面积**）。

  本脚本**与具体文献/数据无关**：不硬编码任何项目路径、队列名、细胞类型名，
  全部走参数（见 REPRO_SPEC.md §1.4 脚本分层）。

【方法学出处】
  after Wang et al. 2026 bioRxiv, Methods「Patch proximity analysis」：
    HDBSCAN 聚类 → 凹包（concave hull）边缘 → 3-NN 扩展 → 5/10/15/20/25 μm 环带计数。
  本脚本把 `--radii-um 5,10,15,20,25`、`--edge-knn 3`、`--min-cluster-size` 等作为参数
  默认值保留，可按数据集调整。

【关键实现细节（会写进输出，便于审阅）】
  1. 聚类：`--cluster-method hdbscan` 优先用 `hdbscan` 包，缺失时回退
     `sklearn.cluster.HDBSCAN`（sklearn ≥ 1.3），两者都缺失则明确报错；
     HDBSCAN 的噪声点（label = -1）不计入任何 patch。
     `--cluster-method none` 表示不做聚类，把该通道全部细胞当作 1 个 patch。
  2. 凹包：优先 `shapely.concave_hull(MultiPoint, ratio=...)`（shapely ≥ 2.1）；
     不可用时回退**自实现 alpha-shape**（Delaunay 三角形按外接圆半径阈值筛选，
     取其只被一个保留三角形共享的边作为边界）；再不可用则回退凸包
     （`--edge-method convex_fallback` 可强制）。**`hull_method` 会写入输出表。**
  3. 边缘细胞判据（确定性）：patch 内细胞到凹包**边界折线**的距离 < `--edge-band-um`
     （默认 None → 自适应取 2 × patch 内最近邻距离的中位数）；到边界距离用
     "点到线段距离" 逐段取最小值计算，不依赖 shapely 的几何运算。
  4. anchor 集合 = 边缘细胞 ∪ 每个边缘细胞的 `--edge-knn` 个最近邻（在 patch 内）。
  5. 环带：非 patch 细胞按 `--ring-origin` 定义的"距离"划入 (r_{i-1}, r_i] 环带；
     面积随半径平方增长，因此**同时输出 density = count / 环带面积**，不要只看 count。
     `--ring-origin anchor`（默认）：距离 = 到**最近 anchor 细胞**的距离。
        注意这是**有系统偏大**的度量：anchor 位于凹包边界**内侧**（边缘带宽 ≤
        `--edge-band-um`），因此外部细胞到 anchor 的距离比到边界的真实距离大
        0 ~ edge_band；最内环带会少计、外环带会多计。输出表给出
        `anchor_mean_depth_um`（anchor 到边界的平均深度）便于判断偏倚量级。
     `--ring-origin boundary`：距离 = 到**凹包边界折线**的距离（无偏，便于与
        几何理论值比较；注意落在凹包内部的非 patch 细胞距离≈0，会进入最内环带）。

【单位】
  `--radii-um` / `--edge-band-um` / `--alpha-um`（若给）/ `--min-cluster-size` 中的
  长度量均以 μm 计，**必须**提供 `--pixel-size-um`（μm/px）做换算；缺参数直接报错退出。
  `--min-cluster-size` 是**细胞数**（无量纲）。

【用法示例】
  python3 tools/patch_proximity.py \\
      --cells processed/cells_with_cn.csv --patch-by cn --patch-value cn_1 \\
      --cluster-method hdbscan --min-cluster-size 50 --edge-knn 3 \\
      --radii-um 5,10,15,20,25 --pixel-size-um 0.5 \\
      --group-col sample --out-prefix results/patch

  python3 tools/patch_proximity.py --cells processed/cells.tsv --sep '\\t' \\
      --patch-by celltype --patch-value Tumor --pixel-size-um 0.5 --out-prefix results/patch

  python3 tools/patch_proximity.py --selftest

【输出】
  <prefix>_patch_summary.csv
      group, patch_id, n_patch_cells, n_noise_cells_group, n_edge_cells, n_anchor_cells,
      hull_method, cluster_method, edge_method, edge_band_um, edge_band_source,
      anchor_mean_depth_um, hull_area_um2, perimeter_um, pixel_size_um, radii_um,
      ring_origin
  <prefix>_ring_counts.csv
      group, patch_id, radius_um, ring_inner_um, ring_outer_um, celltype,
      count, ring_area_um2, density_per_um2, n_anchor_cells, n_patch_cells
      （排序：group, patch_id, radius_um, celltype → 行序稳定）

【依赖】
  必需：numpy, pandas, scipy, scikit-learn
  可选：hdbscan（缺失→sklearn.cluster.HDBSCAN），shapely ≥ 2.1（缺失→自实现 alpha-shape）

【退出码】
  0 = 成功
  1 = 未预期异常
  2 = 参数/输入校验失败
  3 = 质量断言失败（产物为空、无 patch、无 anchor 等）
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
# 聚类
# ----------------------------------------------------------------------------
def spatial_cluster(coords, method, min_cluster_size, min_samples=None):
    """返回 (labels[-1 表示噪声], cluster_method 说明字符串)。"""
    import numpy as np

    if method == "none":
        return np.zeros(coords.shape[0], dtype=int), "none(single_patch)"
    n = coords.shape[0]
    if n < max(2, int(min_cluster_size)):
        raise QualityError(
            "patch 通道细胞只有 %d 个，少于 --min-cluster-size=%d，无法聚类"
            % (n, min_cluster_size))

    try:
        import hdbscan as _hdbscan
        cl = _hdbscan.HDBSCAN(
            min_cluster_size=int(min_cluster_size),
            min_samples=(None if min_samples is None else int(min_samples)),
            cluster_selection_method="eom", allow_single_cluster=True)
        lab = cl.fit_predict(coords)
        return lab, "hdbscan.HDBSCAN(min_cluster_size=%d,min_samples=%s,eom)" % (
            min_cluster_size, min_samples)
    except ImportError:
        pass

    try:
        from sklearn.cluster import HDBSCAN as SkHDBSCAN
        kw = dict(min_cluster_size=int(min_cluster_size),
                  min_samples=(None if min_samples is None else int(min_samples)),
                  cluster_selection_method="eom", allow_single_cluster=True)
        import inspect
        if "copy" in inspect.signature(SkHDBSCAN).parameters:
            kw["copy"] = True          # sklearn >= 1.9 显式传值，避免 FutureWarning
        cl = SkHDBSCAN(**kw)
        lab = cl.fit_predict(coords)
        return lab, "sklearn.cluster.HDBSCAN(min_cluster_size=%d,min_samples=%s,eom)" % (
            min_cluster_size, min_samples)
    except ImportError:
        raise QualityError(
            "缺少 HDBSCAN 实现：`hdbscan` 与 `sklearn.cluster.HDBSCAN` 都不可用。\n"
            "        安装：python3 -m pip install --user hdbscan\n"
            "        或升级 scikit-learn ≥ 1.3；也可用 --cluster-method none 跳过聚类。")


# ----------------------------------------------------------------------------
# 边界：凹包 / alpha-shape / 凸包 → 边界折线（线段数组）
# ----------------------------------------------------------------------------
def segments_from_ring(ring_xy):
    """闭合环 -> 线段数组 (m,2,2)。"""
    import numpy as np

    p = np.asarray(ring_xy, dtype=float)
    if p.shape[0] < 3:
        return np.zeros((0, 2, 2))
    if np.allclose(p[0], p[-1]):
        p = p[:-1]
    nxt = np.roll(np.arange(len(p)), -1)
    return np.stack([p, p[nxt]], axis=1)


def boundary_length(segs):
    """边界折线总长（用于量化 alpha-shape 的"碎裂"程度）。"""
    import numpy as np

    if segs.shape[0] == 0:
        return 0.0
    return float(np.linalg.norm(segs[:, 1, :] - segs[:, 0, :], axis=1).sum())


def alpha_shape_segments(coords, alpha, return_area=False):
    """自实现 alpha-shape：Delaunay 三角形外接圆半径 <= alpha 保留；
    边界边 = 只被一个保留三角形共享的边。

    alpha 单位与 coords 相同（px）。返回 (segments, area_px2) 或 segments。
    """
    import numpy as np
    from scipy.spatial import Delaunay

    m = coords.shape[0]
    if m < 3:
        return (np.zeros((0, 2, 2)), 0.0) if return_area else np.zeros((0, 2, 2))
    tri = Delaunay(coords)
    s = tri.simplices
    a, b, c = s[:, 0], s[:, 1], s[:, 2]
    p0, p1, p2 = coords[a], coords[b], coords[c]
    # 边长
    la = np.linalg.norm(p1 - p2, axis=1)
    lb = np.linalg.norm(p0 - p2, axis=1)
    lc = np.linalg.norm(p0 - p1, axis=1)
    # 面积（有向）
    area2 = (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1]) - \
            (p2[:, 0] - p0[:, 0]) * (p1[:, 1] - p0[:, 1])
    area = np.abs(area2) / 2.0
    with np.errstate(divide="ignore", invalid="ignore"):
        circum = (la * lb * lc) / (4.0 * area)
    keep = circum <= float(alpha)
    keep &= np.isfinite(circum)
    if not keep.any():
        return (np.zeros((0, 2, 2)), 0.0) if return_area else np.zeros((0, 2, 2))

    # 边界边：在保留三角形中出现 1 次的边
    edges = np.vstack([s[keep][:, [0, 1]], s[keep][:, [1, 2]], s[keep][:, [2, 0]]])
    edges = np.sort(edges, axis=1)
    uniq, counts = np.unique(edges, axis=0, return_counts=True)
    bnd = uniq[counts == 1]
    segs = np.stack([coords[bnd[:, 0]], coords[bnd[:, 1]]], axis=1) if bnd.size \
        else np.zeros((0, 2, 2))
    if return_area:
        return segs, float(area[keep].sum())
    return segs


def concave_hull_segments(coords, ratio=0.3, method="auto"):
    """返回 (segments, hull_method, area_px2, perimeter_px)。

    method: auto / concave / convex_fallback / alpha
    """
    import numpy as np

    if method in ("auto", "concave"):
        try:
            import shapely
            fn = getattr(shapely, "concave_hull", None)
            if fn is None:
                raise AttributeError("shapely.concave_hull 不可用（需 shapely ≥ 2.1）")
            hull = fn(shapely.MultiPoint(coords), ratio=float(ratio))
            if hull.is_empty or hull.geom_type != "Polygon":
                raise ValueError("concave_hull 未返回有效 Polygon（ratio 过小？）")
            segs = segments_from_ring(np.asarray(hull.exterior.coords))
            return segs, "shapely.concave_hull(ratio=%.2f)" % ratio, float(hull.area), \
                float(hull.length)
        except Exception as exc:
            if method == "concave":
                raise QualityError("shapely 凹包失败：%s；可改用 --edge-method alpha/convex_fallback"
                                   % exc)
            print("[WARN] shapely.concave_hull 不可用（%s）→ 回退自实现 alpha-shape" % exc)

    if method in ("auto", "concave", "alpha"):
        # 自适应 alpha：从 2×最近邻距离中位数起逐步放大，直到边界总长接近凸包周长。
        # 理由：alpha 过小时 Delaunay 三角形被大量剔除，边界边在图元之间"碎裂"，
        # 会几乎把 patch 内每一个细胞都判成边缘细胞（假边缘）。边界总长是"碎裂程度"
        # 的可计算代理量，凸包周长是任何边界的下界，故以 1.5×凸包周长为收敛判据。
        from scipy.spatial import cKDTree, ConvexHull
        tree = cKDTree(coords)
        d, _ = tree.query(coords, k=2)
        nn = float(np.median(d[:, 1])) if d.shape[1] > 1 else 1.0
        alpha = 2.0 * max(nn, 1e-6)
        hull_len = float(ConvexHull(coords).area)
        segs, area, used = np.zeros((0, 2, 2)), 0.0, None
        for _ in range(14):
            cand, cand_area = alpha_shape_segments(coords, alpha, return_area=True)
            used = alpha
            if cand.shape[0] >= 3:
                segs, area = cand, cand_area
                if boundary_length(cand) <= 1.5 * hull_len:
                    break
            alpha *= 1.3
        if segs.shape[0] >= 3:
            return segs, ("alpha_shape(alpha=%.6g px=%.2f×median NN)" % (used, used / nn)), \
                float(area), float("nan")
        print("[WARN] alpha-shape 边界边不足 3 条 → 回退凸包")

    from scipy.spatial import ConvexHull
    hull = ConvexHull(coords)
    pts = coords[hull.vertices]
    segs = segments_from_ring(np.vstack([pts, pts[:1]]))
    return segs, "convex_fallback(ConvexHull)", float(hull.volume), float(hull.area)


def point_to_segments_min_dist(points, segs, chunk=2000):
    """点到线段集合的最近距离（逐点取所有线段最小值）。

    points: (n,2)；segs: (m,2,2) → 返回 (n,)。
    """
    import numpy as np

    if segs.shape[0] == 0:
        return np.full(points.shape[0], np.inf)
    out = np.full(points.shape[0], np.inf)
    A = segs[:, 0, :]
    B = segs[:, 1, :]
    AB = B - A
    denom = np.einsum("ij,ij->i", AB, AB)
    denom = np.where(denom > 0, denom, 1.0)
    for s in range(0, points.shape[0], chunk):
        P = points[s:s + chunk]
        AP = P[:, None, :] - A[None, :, :]
        t = np.einsum("nmi,mi->nm", AP, AB) / denom[None, :]
        t = np.clip(t, 0.0, 1.0)
        proj = A[None, :, :] + t[:, :, None] * AB[None, :, :]
        d = np.linalg.norm(P[:, None, :] - proj, axis=2)
        out[s:s + chunk] = d.min(axis=1)
    return out


# ----------------------------------------------------------------------------
# 主流程
# ----------------------------------------------------------------------------
def run_patch(df, x_col, y_col, ct_col, group_col, patch_col, patch_value,
              cluster_method, min_cluster_size, edge_method, edge_band_um,
              edge_knn, radii_um, pixel_size_um, ring_origin="anchor",
              concave_ratio=0.3, hdbscan_min_samples=None):
    import numpy as np
    import pandas as pd
    from scipy.spatial import cKDTree

    x, y = numeric_coords(df, x_col, y_col)
    coords = np.column_stack([x, y])
    ct = df[ct_col].astype(str).to_numpy()

    if group_col is None:
        groups_all = np.full(len(df), "__all__", dtype=object)
    else:
        groups_all = df[group_col].astype(str).to_numpy()
    is_patch = (df[patch_col].astype(str) == str(patch_value)).to_numpy()
    if not is_patch.any():
        raise UsageError(
            "列 '%s' 中没有值为 '%s' 的细胞；该列取值示例：%s"
            % (patch_col, patch_value, sorted(pd.unique(df[patch_col].astype(str)))[:20]))
    if not (ct[~is_patch].size):
        raise QualityError("除 patch 之外没有任何其他细胞，无法统计环带")

    seg_of_um = 1.0 / float(pixel_size_um)
    radii_um = sorted(float(r) for r in radii_um)
    summary_rows, ring_rows = [], []
    group_order = sorted(pd.unique(groups_all).tolist())

    for g in group_order:
        gm = groups_all == g
        pm = gm & is_patch
        if pm.sum() < 3:
            print("[WARN] 分组 '%s' 的 patch 通道细胞只有 %d 个（<3），跳过" % (g, int(pm.sum())))
            continue
        P = coords[pm]
        labs, cluster_desc = spatial_cluster(P, cluster_method, min_cluster_size,
                                             hdbscan_min_samples)
        noise = labs == -1
        if labs[~noise].size == 0:
            print("[WARN] 分组 '%s'：HDBSCAN 未形成任何簇（全部为噪声），跳过" % g)
            continue

        for cid in sorted(set(labs[~noise].tolist())):
            Q = P[labs == cid]                           # 该 patch 的细胞坐标
            if Q.shape[0] < 3:
                continue
            segs, hull_desc, h_area_px2, h_perim_px = concave_hull_segments(
                Q, ratio=concave_ratio, method=edge_method)

            # 边缘带宽（μm → px）；未指定时自适应 = 2 × patch 内最近邻距离中位数
            tree_q = cKDTree(Q)
            if edge_band_um is None:
                dq, _ = tree_q.query(Q, k=min(2, Q.shape[0]))
                band_px = 2.0 * float(np.median(dq[:, -1]))
                band_um = band_px * pixel_size_um
                band_desc = "auto(2*median_NN)"
            else:
                band_px = float(edge_band_um) * seg_of_um
                band_um = float(edge_band_um)
                band_desc = "explicit"

            dist_to_bnd = point_to_segments_min_dist(Q, segs)
            edge_mask = dist_to_bnd < band_px
            if edge_mask.mean() > 0.8 and Q.shape[0] > 20:
                print("[WARN] 分组 '%s' patch %s：%.0f%% 的 patch 细胞落在边缘带宽内，"
                      "边界可能已碎裂（凹包/alpha-shape 内切或参数过小），"
                      % (g, cid, 100 * edge_mask.mean()))
                print("       边缘细胞集合会接近整个 patch；建议减小 --edge-band-um，"
                      "或调整 --concave-ratio / --edge-method。")
            if not edge_mask.any():                      # 兜底：取距边界最近的 5%
                kmin = max(3, int(0.05 * Q.shape[0]))
                idx = np.argsort(dist_to_bnd, kind="stable")[:kmin]
                edge_mask = np.zeros(Q.shape[0], dtype=bool)
                edge_mask[idx] = True
                print("[WARN] 分组 '%s' patch %s：无细胞落在边缘带宽内，"
                      "改取距边界最近的 %d 个细胞" % (g, cid, kmin))

            # anchor = 边缘细胞 ∪ 每个边缘细胞的 edge_knn 个最近邻（patch 内）
            kk = min(int(edge_knn) + 1, Q.shape[0])
            _, nb = tree_q.query(Q[edge_mask], k=kk)
            nb = np.atleast_2d(nb)
            anchor = np.zeros(Q.shape[0], dtype=bool)
            anchor[np.unique(nb.ravel())] = True
            A = Q[anchor]
            anchor_depth_um = float(dist_to_bnd[anchor].mean()) * pixel_size_um

            perim_um = h_perim_px * pixel_size_um
            area_um2 = h_area_px2 * (pixel_size_um ** 2)

            # 环带坐标：非 patch 细胞 → 到 anchor 的最近距离 或 到凹包边界的距离
            om = gm & (~is_patch)
            if om.sum() == 0:
                continue
            O = coords[om]
            oct_ = ct[om]
            if ring_origin == "boundary":
                d_px = point_to_segments_min_dist(O, segs)
            else:
                d_px, _ = cKDTree(A).query(O, k=1)
            d_um = d_px * pixel_size_um

            other_types_g = sorted(set(oct_.tolist()))
            for i, r_um in enumerate(radii_um):
                lo_um = radii_um[i - 1] if i > 0 else 0.0
                m = (d_um > lo_um) & (d_um <= r_um)
                ring_area = np.pi * (r_um ** 2 - lo_um ** 2)
                vals, cnts = np.unique(oct_[m], return_counts=True)
                counts = dict(zip(vals.tolist(), cnts.tolist()))
                # 该分组的全部非 patch 类型都出零行，保证表格完整（避免静默缺格）
                for t in other_types_g:
                    c = int(counts.get(t, 0))
                    ring_rows.append({
                        "group": g, "patch_id": "patch_%d" % cid,
                        "radius_um": r_um, "ring_inner_um": lo_um, "ring_outer_um": r_um,
                        "celltype": t, "count": c,
                        "ring_area_um2": float(ring_area),
                        "density_per_um2": c / ring_area,
                        "n_anchor_cells": int(anchor.sum()),
                        "n_patch_cells": int(Q.shape[0]),
                        "ring_origin": ring_origin,
                    })

            n_channel_group = int(pm.sum())
            noise_frac = float(noise.sum()) / float(n_channel_group) if n_channel_group else 0.0
            if noise_frac > 0.2:
                print("[WARN] 分组 '%s'：patch 通道 %d/%d 个细胞（%.0f%%）被 HDBSCAN 判为噪声，"
                      "patch 只覆盖该通道的一部分。"
                      % (g, int(noise.sum()), n_channel_group, 100 * noise_frac))
                print("       HDBSCAN 基于密度对比找簇，均匀密度的连续区域可能只被取出一"
                      "小核心 + 大量噪声。")
                print("       若该通道本身就是一个连续团块，请用 --cluster-method none"
                      "（整体当作 1 个 patch），或调整 --min-cluster-size / "
                      "--hdbscan-min-samples。")
            summary_rows.append({
                "group": g, "patch_id": "patch_%d" % cid,
                "n_patch_cells": int(Q.shape[0]),
                "n_noise_cells_group": int(noise.sum()),
                "noise_frac_group": noise_frac,
                "n_edge_cells": int(edge_mask.sum()),
                "n_anchor_cells": int(anchor.sum()),
                "hull_method": hull_desc, "cluster_method": cluster_desc,
                "edge_method": edge_desc(band_px),
                "edge_band_um": float(band_um), "edge_band_source": band_desc,
                "anchor_mean_depth_um": anchor_depth_um,
                "hull_area_um2": float(area_um2), "perimeter_um": float(perim_um),
                "pixel_size_um": float(pixel_size_um),
                "radii_um": ",".join("%g" % r for r in radii_um),
                "ring_origin": ring_origin,
            })

    if not summary_rows:
        raise QualityError(
            "没有产生任何 patch 产物。可能原因：patch 通道细胞过少、--min-cluster-size 过大、"
            "或每个分组内 HDBSCAN 都把全部细胞判为噪声（可尝试减小 --min-cluster-size，"
            "或用 --cluster-method none 把该通道当作单一 patch）。\n"
            "        （REPRO_SPEC 附录 C：禁止静默产出空文件）")
    summary = pd.DataFrame(summary_rows).sort_values(
        ["group", "patch_id"], kind="stable").reset_index(drop=True)
    rings = pd.DataFrame(ring_rows)
    if rings.shape[0]:
        rings = rings.sort_values(["group", "patch_id", "radius_um", "celltype"],
                                  kind="stable").reset_index(drop=True)
    return summary, rings


def edge_desc(band_px):
    return "hull_boundary(dist_to_boundary<%.6g px)" % band_px


# ----------------------------------------------------------------------------
# 自检：环带计数 vs 理论值
# ----------------------------------------------------------------------------
def make_synthetic_patch(seed=20260101, a_um=100.0, b_um=60.0, n_patch=6000,
                         density=0.5, pixel_size_um=0.5):
    """构造椭圆 patch + 按"到边界距离"均匀（恒定面密度）分布的外部细胞。

    外部细胞用椭圆外平行体精确采样：δ 的密度 ∝ (L + 2πδ)（Steiner 公式：
    外平行体面积 = L·δ + πδ²），沿外法线偏移 δ，故其**到边界的真实距离恰为 δ**。
    返回 (DataFrame[x_px,y_px,celltype], 真值字典)
    """
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(seed)
    radii = [5.0, 10.0, 15.0, 20.0, 25.0]
    dmax = max(radii)

    # --- patch 内细胞（椭圆盘内均匀）-----------------------------------------
    r = np.sqrt(rng.uniform(0, 1, n_patch))
    th = rng.uniform(0, 2 * np.pi, n_patch)
    px = a_um * r * np.cos(th)
    py = b_um * r * np.sin(th)

    # --- 外部细胞：delta 采样 + 沿外法线偏移 ---------------------------------
    L = float(np.pi * (3 * (a_um + b_um) -
                       np.sqrt((3 * a_um + b_um) * (a_um + 3 * b_um))))   # Ramanujan 周长
    area_band = L * dmax + np.pi * dmax ** 2
    n_out = int(round(density * area_band))
    u = rng.uniform(0, 1, n_out)
    # 解 πδ² + Lδ - u·area_band = 0
    delta = (-L + np.sqrt(L ** 2 + 4 * np.pi * u * area_band)) / (2 * np.pi)
    t = rng.uniform(0, 2 * np.pi, n_out)
    bx, by = a_um * np.cos(t), b_um * np.sin(t)
    nx, ny = b_um * np.cos(t), a_um * np.sin(t)        # 椭圆外法线方向
    nn = np.hypot(nx, ny)
    nx, ny = nx / nn, ny / nn
    ox = bx + delta * nx
    oy = by + delta * ny

    other_types = ["immune", "stroma", "tcell"]
    otype = rng.choice(other_types, size=n_out, p=[0.5, 0.3, 0.2])

    df = pd.concat([
        pd.DataFrame({"x_px": px / pixel_size_um, "y_px": py / pixel_size_um,
                      "celltype": "patchcell", "true_delta_um": np.nan}),
        pd.DataFrame({"x_px": ox / pixel_size_um, "y_px": oy / pixel_size_um,
                      "celltype": otype, "true_delta_um": delta}),
    ], ignore_index=True)
    truth = {"perimeter_um": L, "radii": radii, "density": density,
             "pixel_size_um": pixel_size_um}
    return df, truth


def selftest_hdbscan_path(seed=20260101):
    """结构性检查：走 HDBSCAN 聚类 + 多 patch + 分组 + 输出表完整性。

    几何：把 selftest 的椭圆 patch 细胞复制并平移，得到 2 个分离的 patch；
    另加一批均匀背景细胞（非 patch 类型）。不做理论值断言，只检查：
      * HDBSCAN 路径被记录（cluster_method 字符串）
      * 找到的 patch 数 >= 2
      * 环带表完整（patch 数 × 环带数 × 该分组非 patch 类型数，含 count=0 的补零行）
    """
    import numpy as np
    import pandas as pd

    pixel = 0.5
    df, truth = make_synthetic_patch(seed=seed, pixel_size_um=pixel, n_patch=3000)
    patch_cells = df[df["celltype"] == "patchcell"].copy()
    shifted = patch_cells.copy()
    shifted["x_px"] = shifted["x_px"] + 400.0 / pixel      # 平移 +400 μm，形成第二个 patch
    rng = np.random.default_rng(seed + 1)
    n_bg = 4000
    bg = pd.DataFrame({
        "x_px": rng.uniform(-200.0 / pixel, 900.0 / pixel, n_bg),
        "y_px": rng.uniform(-200.0 / pixel, 300.0 / pixel, n_bg),
        "celltype": rng.choice(["immune", "stroma", "tcell"], n_bg, p=[0.5, 0.3, 0.2]),
        "true_delta_um": np.nan,
    })
    d2 = pd.concat([patch_cells, shifted, bg], ignore_index=True)
    d2["channel"] = np.where(d2["celltype"] == "patchcell", "yes", "no")

    radii = truth["radii"]
    t0 = time.time()
    summary, rings = run_patch(
        d2, "x_px", "y_px", "celltype", None, "channel", "yes",
        cluster_method="hdbscan", min_cluster_size=200, edge_method="auto",
        edge_band_um=8.0, edge_knn=3, radii_um=radii, pixel_size_um=pixel,
        ring_origin="boundary")
    dt = time.time() - t0

    n_other = d2.loc[d2["channel"] == "no", "celltype"].nunique()
    expect_rows = summary.shape[0] * len(radii) * n_other
    print("\n" + "=" * 78)
    print("[HDBSCAN 路径检查] 细胞=%d（2 个真实 patch + 背景）；cluster_method=%s"
          % (len(d2), summary["cluster_method"].iloc[0]))
    print("[HDBSCAN 路径检查] 找到 patch 数=%d；各 patch 细胞数=%s；噪声占通道比例=%.2f"
          % (summary.shape[0], summary["n_patch_cells"].tolist(),
             float(summary["noise_frac_group"].iloc[0])))
    print("[HDBSCAN 路径检查] hull_method=%s" % summary["hull_method"].iloc[0])
    print("[HDBSCAN 路径检查] 环带表行数=%d（期望 patch%d × 环带%d × 类型%d = %d，含补零行）；用时 %.2f s"
          % (rings.shape[0], summary.shape[0], len(radii), n_other, expect_rows, dt))

    failed = []
    if summary.shape[0] < 2:
        failed.append("HDBSCAN 只找到 %d 个 patch（应为 2）" % summary.shape[0])
    if "HDBSCAN" not in summary["cluster_method"].iloc[0]:
        failed.append("cluster_method 未记录 HDBSCAN 路径")
    if rings.shape[0] != expect_rows:
        failed.append("环带表行数 %d != 期望 %d（补零缺失？）" % (rings.shape[0], expect_rows))
    if set(rings["radius_um"].unique()) != set(radii):
        failed.append("环带表缺少部分半径：%s" % sorted(rings["radius_um"].unique()))
    return failed, dict(n_patch=int(summary.shape[0]),
                        n_ring_rows=int(rings.shape[0]),
                        noise_frac=float(summary["noise_frac_group"].iloc[0]))


def selftest(seed=20260101):
    import numpy as np

    print("=" * 78)
    print("patch_proximity.py --selftest")
    print("=" * 78)
    pixel = 0.5
    df, truth = make_synthetic_patch(seed=seed, pixel_size_um=pixel)
    radii = truth["radii"]
    L, dens = truth["perimeter_um"], truth["density"]
    print("[合成组织] 椭圆 patch a=100 μm,b=60 μm；周长 L≈%.1f μm；"
          "外部细胞面密度=%.3f /μm²" % (L, dens))
    print("[合成组织] 外部细胞沿外法线偏移 δ 放置（其到边界真实距离恰为 δ），"
          "δ 密度 ∝ (L+2πδ)（Steiner 公式）")
    print("[合成组织] 像素尺寸=%.2f μm/px，坐标写入 x_px/y_px；细胞总数=%d；"
          "环带=%s μm" % (pixel, len(df), ",".join("%g" % r for r in radii)))
    print("[理论值] 环带 (r_lo, r_hi] 期望细胞数 = density × (L·Δr + π(r_hi²-r_lo²))")

    # 三组配置：两种 ring 坐标原点 × 两种凹包 ratio
    #   tolerance = 该配置允许的最大相对误差（None = 仅报告不设断言）
    configs = [
        ("boundary", 0.3, 0.15,
         "默认凹包 ratio=0.3 + 到凹包边界距离（ratio 越小边界越向内切）"),
        ("boundary", 1.0, 0.15,
         "凹包 ratio=1.0（近凸，贴合椭圆）+ 到凹包边界距离"),
        ("anchor", 1.0, 0.25,
         "凹包 ratio=1.0 + 到最近 anchor 距离（--ring-origin anchor，默认）"),
    ]

    results = []
    for origin, cratio, tol, label in configs:
        tag = "origin=%s,ratio=%.1f" % (origin, cratio)
        t0 = time.time()
        summary, rings = run_patch(
            df, "x_px", "y_px", "celltype", None, "celltype", "patchcell",
            cluster_method="none", min_cluster_size=10, edge_method="auto",
            edge_band_um=8.0, edge_knn=3, radii_um=radii, pixel_size_um=pixel,
            ring_origin=origin, concave_ratio=cratio)
        dt = time.time() - t0
        if rings.shape[0] == 0:
            print("[FAIL] %s 未产生任何环带计数" % tag)
            return EXIT_SELFTEST
        s = summary.iloc[0]
        print("\n" + "-" * 78)
        print("[%s] %s" % (tag, label))
        print("[%s] patch 细胞=%d，边缘细胞=%d，anchor=%d（anchor 平均深度 %.2f μm）；用时 %.2f s"
              % (tag, s["n_patch_cells"], s["n_edge_cells"], s["n_anchor_cells"],
                 s["anchor_mean_depth_um"], dt))
        print("[%s] hull_method=%s | cluster_method=%s | edge_band=%.3f μm(%s)"
              % (tag, s["hull_method"], s["cluster_method"], s["edge_band_um"],
                 s["edge_band_source"]))
        agg = rings.groupby("radius_um", as_index=False)["count"].sum()
        print("[%s] %-10s %10s %12s %12s %10s"
              % (tag, "环带(μm)", "观测count", "理论count", "绝对误差", "相对误差"))
        worst, per_band = 0.0, {}
        for i, r_um in enumerate(radii):
            lo = radii[i - 1] if i > 0 else 0.0
            expected = dens * (np.pi * (r_um ** 2 - lo ** 2) + L * (r_um - lo))
            obs = int(agg.loc[agg["radius_um"] == r_um, "count"].sum())
            rel = abs(obs - expected) / expected
            per_band[r_um] = rel
            worst = max(worst, rel)
            print("[%s] %-10s %10d %12.1f %12.1f %9.2f%%"
                  % (tag, "%.0f-%.0f" % (lo, r_um), obs, expected,
                     obs - expected, 100 * rel))
        print("[%s] >>> 最大相对误差 = %.2f%% <<<" % (tag, 100 * worst))
        results.append(dict(tag=tag, origin=origin, ratio=cratio, tol=tol, worst=worst,
                            seconds=dt, n_edge=int(s["n_edge_cells"]),
                            n_anchor=int(s["n_anchor_cells"]),
                            depth=float(s["anchor_mean_depth_um"]), per_band=per_band))

    print("\n" + "=" * 78)
    print("[结论] 汇总（环带计数 vs 几何理论值，最大相对误差）：")
    for r in results:
        lim = "<%.0f%%" % (100 * r["tol"]) if r["tol"] else "不设断言"
        print("  %-24s 最大相对误差 %6.2f%%（要求 %s）" % (r["tag"], 100 * r["worst"], lim))
    ab = [r for r in results if r["origin"] == "boundary"]
    an = [r for r in results if r["origin"] == "anchor"][0]
    best_b = min(ab, key=lambda r: r["worst"])
    print("[结论] 到凹包边界的距离（--ring-origin boundary）在凹包贴合真实形状时无偏："
          "最大相对误差 %.2f%%" % (100 * best_b["worst"]))
    print("[结论] 到最近 anchor 的距离（--ring-origin anchor，默认）系统性偏大：anchor 平均"
          "位于边界内侧 %.2f μm，" % an["depth"])
    print("       最内环带相对误差 %.1f%%、最外环带 %.1f%%（可预期的几何偏倚，非实现错误）；"
          % (100 * an["per_band"][radii[0]], 100 * an["per_band"][radii[-1]]))
    print("       需与几何理论值/文献环带半径严格对齐时请用 --ring-origin boundary。")
    print("[结论] 凹包 ratio 影响边界内切程度：ratio=0.3 时 boundary 最大相对误差 %.2f%%，"
          "ratio=1.0 时为 %.2f%%"
          % (100 * [r for r in ab if r["ratio"] == 0.3][0]["worst"],
             100 * [r for r in ab if r["ratio"] == 1.0][0]["worst"]))

    failed = []
    for r in results:
        if r["tol"] is not None and not (r["worst"] < r["tol"]):
            failed.append("%s 最大相对误差 %.2f%% >= %.0f%%"
                          % (r["tag"], 100 * r["worst"], 100 * r["tol"]))
    if results[0]["n_anchor"] <= results[0]["n_edge"]:
        failed.append("anchor 数未超过边缘细胞数（--edge-knn 扩展未生效）")

    hdb_failed, hdb_info = selftest_hdbscan_path(seed=seed)
    failed.extend(hdb_failed)
    if failed:
        print("\n[FAIL] selftest 未通过：" + "；".join(failed))
        return EXIT_SELFTEST
    print("\n[PASS] selftest 通过：ring_origin=boundary 全部环带相对误差 < 15%%"
          "（ratio=0.3 时 %.2f%%，ratio=1.0 时 %.2f%%）；"
          "ring_origin=anchor 偏倚在预期内（最大 %.2f%%）；"
          "HDBSCAN 路径找到 %d 个 patch，环带表 %d 行完整"
          % (100 * [r for r in results if r["origin"] == "boundary" and r["ratio"] == 0.3][0]["worst"],
             100 * [r for r in results if r["origin"] == "boundary" and r["ratio"] == 1.0][0]["worst"],
             100 * [r for r in results if r["origin"] == "anchor"][0]["worst"],
             hdb_info["n_patch"], hdb_info["n_ring_rows"]))
    return EXIT_OK


# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------
def build_parser():
    p = argparse.ArgumentParser(
        prog="patch_proximity.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "Patch 邻近分析：指定通道细胞 → HDBSCAN 聚类 → 凹包边缘 → k-NN 扩展 → "
            "半径环带内异类细胞计数与密度。\n"
            "after Wang et al. 2026 bioRxiv, Methods「Patch proximity analysis」。"
        ),
        epilog="""\
流程与判据
  1) 取 --patch-by/--patch-value 指定通道的细胞 → --cluster-method hdbscan 聚成 patch
     （噪声点 label=-1 不计入任何 patch；--cluster-method none = 全部当作 1 个 patch）。
  2) 每个 patch 求凹包（shapely.concave_hull，回退自实现 alpha-shape，再回退凸包）；
     hull_method 写入输出表，可核对用了哪条路径。
  3) 边缘细胞 = patch 内"到凹包边界折线距离 < --edge-band-um"的细胞（默认自适应
     2×patch 内最近邻距离中位数）。
  4) anchor = 边缘细胞 ∪ 每个边缘细胞在 patch 内的 --edge-knn 个最近邻。
  5) 环带 = (r_{i-1}, r_i]，按非 patch 细胞的 --ring-origin 距离划分：
     anchor（默认）= 到最近 anchor 细胞的距离（anchor 位于边界内侧 0~edge_band，
     故最内环带会少计、外环带多计）；boundary = 到凹包边界折线的距离（无偏）；
     同时输出 count 与 density = count/环带面积（**面积随半径平方增长，勿只看 count**）。

单位
  --radii-um / --edge-band-um 为 μm，需 --pixel-size-um（μm/px）换算；缺参数报错退出。
  --min-cluster-size / --edge-knn 为单位无关的**细胞数**。

示例
  python3 tools/patch_proximity.py --cells processed/cells_with_cn.csv \\
      --patch-by cn --patch-value cn_1 --cluster-method hdbscan --min-cluster-size 50 \\
      --edge-knn 3 --radii-um 5,10,15,20,25 --pixel-size-um 0.5 \\
      --group-col sample --out-prefix results/patch
  python3 tools/patch_proximity.py --selftest

退出码：0=成功 1=未预期异常 2=参数/输入校验失败 3=质量断言失败 4=selftest 失败 5=缺必需依赖
""",
    )
    io = p.add_argument_group("输入输出")
    io.add_argument("--cells", metavar="FILE", help="细胞表（CSV/TSV，每行一个细胞）")
    io.add_argument("--sep", default=None, help="列分隔符（默认按扩展名）")
    io.add_argument("--x-col", default="x", help="x 坐标列名（默认 x）")
    io.add_argument("--y-col", default="y", help="y 坐标列名（默认 y）")
    io.add_argument("--celltype-col", default="celltype", help="细胞类型列名（默认 celltype）")
    io.add_argument("--group-col", default=None,
                    help="分组列名（如 sample/fov；默认 None = 全部细胞视为一组）")
    io.add_argument("--out-prefix", default="patch_proximity_out", metavar="PREFIX",
                    help="输出前缀（<PREFIX>_patch_summary.csv、<PREFIX>_ring_counts.csv）")

    t = p.add_argument_group("patch 通道")
    t.add_argument("--patch-by", choices=["cn", "celltype"], default="cn",
                   help="用哪一列定义 patch 通道（cn 用 --cn-col；celltype 用 --celltype-col）")
    t.add_argument("--cn-col", default="cn", help="--patch-by cn 时的列名（默认 cn）")
    t.add_argument("--patch-value", metavar="VAL",
                   help="patch 通道的取值（如 cn_1 / Tumor）；**必需**")

    c = p.add_argument_group("聚类与边缘")
    c.add_argument("--cluster-method", choices=["hdbscan", "none"], default="hdbscan",
                   help="patch 聚类方法（默认 hdbscan：优先 hdbscan 包，回退 sklearn HDBSCAN）")
    c.add_argument("--min-cluster-size", type=int, default=50,
                   help="HDBSCAN 最小簇大小，单位=细胞数（默认 50）")
    c.add_argument("--hdbscan-min-samples", type=int, default=None, metavar="N",
                   help="HDBSCAN 的 min_samples（细胞数；默认 None = 由 HDBSCAN 自行取 "
                        "min_cluster_size）。调大更保守、噪声更多；调小更易把连续区域连成一簇")
    c.add_argument("--edge-method", choices=["auto", "concave", "alpha", "convex_fallback"],
                   default="auto",
                   help="边界提取路径：auto（凹包→alpha-shape→凸包，默认）/ concave（仅凹包）"
                        "/ alpha（自实现 alpha-shape）/ convex_fallback（凸包）")
    c.add_argument("--edge-band-um", type=float, default=None,
                   help="边缘带宽（μm）：到凹包边界距离小于该值的 patch 细胞算边缘细胞；"
                        "默认 None = 自适应 2×patch 内最近邻距离中位数")
    c.add_argument("--edge-knn", type=int, default=3,
                   help="每个边缘细胞额外纳入的最近邻个数（默认 3，复现原文）")
    c.add_argument("--concave-ratio", type=float, default=0.3, metavar="R",
                   help="shapely.concave_hull 的 ratio 参数，0~1（默认 0.3；越大越接近凸包）。"
                        "ratio 越小边界越向内切，环带半径的几何含义会随之偏移")

    r = p.add_argument_group("环带")
    r.add_argument("--radii-um", default="5,10,15,20,25",
                   help="环带外半径，逗号分隔，单位 μm（默认 5,10,15,20,25，复现原文）")
    r.add_argument("--pixel-size-um", type=float, default=None,
                   help="像素物理边长 μm/px（长度类参数换算**必需**）")
    r.add_argument("--ring-origin", choices=["anchor", "boundary"], default="anchor",
                   help="环带坐标原点：anchor（默认）= 到最近 anchor 细胞的距离"
                        "（系统性偏大 0~edge_band）；boundary = 到凹包边界的距离（无偏，"
                        "便于与几何理论值比较）。详见 --help 顶部说明")
    r.add_argument("--seed", type=int, default=20260101, help="随机种子（默认 20260101）")

    p.add_argument("--selftest", action="store_true",
                   help="合成椭圆 patch 自检：报告环带计数与理论值的相对误差"
                        "（断言 < 15%%）")
    return p


def parse_radii(s):
    try:
        vals = [float(v) for v in str(s).split(",") if str(v).strip() != ""]
    except ValueError:
        raise UsageError("--radii-um 解析失败：%r（应为逗号分隔的正数，如 5,10,15,20,25）" % s)
    if not vals or any(v <= 0 for v in vals):
        raise UsageError("--radii-um 必须为正数列表，当前 %r" % s)
    return sorted(set(vals))


def validate(args):
    if args.selftest:
        return
    if not args.cells:
        raise UsageError("缺少 --cells（或使用 --selftest）")
    if args.patch_value in (None, ""):
        raise UsageError("缺少 --patch-value（patch 通道取值，如 cn_1 或某细胞类型名）")
    if args.x_col == args.y_col:
        raise UsageError("--x-col 与 --y-col 不能相同（%s）" % args.x_col)
    if args.pixel_size_um is None:
        raise UsageError(
            "--radii-um/--edge-band-um 均为 μm，需要 --pixel-size-um（μm/px）换算；"
            "缺参数会把 μm 当像素用（REPRO_SPEC 附录 C 尺度错配）")
    if not (args.pixel_size_um > 0):
        raise UsageError("--pixel-size-um 必须为正数，当前 %r" % args.pixel_size_um)
    if args.edge_band_um is not None and not (args.edge_band_um > 0):
        raise UsageError("--edge-band-um 必须为正数，当前 %r" % args.edge_band_um)
    if args.edge_knn < 0:
        raise UsageError("--edge-knn 必须 >= 0，当前 %d" % args.edge_knn)
    if args.concave_ratio is not None and not (0.0 <= args.concave_ratio <= 1.0):
        raise UsageError("--concave-ratio 必须在 [0, 1]，当前 %r" % args.concave_ratio)
    if args.min_cluster_size < 2:
        raise UsageError("--min-cluster-size 必须 >= 2，当前 %d" % args.min_cluster_size)


def main(argv=None):
    args = build_parser().parse_args(argv)
    _require_deps()

    if args.selftest:
        return selftest()
    try:
        validate(args)
        import pandas as pd

        group_col = args.group_col if args.group_col not in (None, "") else None
        patch_col = args.cn_col if args.patch_by == "cn" else args.celltype_col
        radii = parse_radii(args.radii_um)

        df = read_table(args.cells, args.sep)
        need = [args.x_col, args.y_col, args.celltype_col, patch_col]
        if group_col is not None:
            need.append(group_col)
        check_columns(df, need, args.cells)

        t0 = time.time()
        summary, rings = run_patch(
            df, args.x_col, args.y_col, args.celltype_col, group_col, patch_col,
            args.patch_value, args.cluster_method, args.min_cluster_size,
            args.edge_method, args.edge_band_um, args.edge_knn, radii,
            args.pixel_size_um, args.ring_origin, args.concave_ratio,
            args.hdbscan_min_samples)
        dt = time.time() - t0

        if summary.shape[0] == 0:
            raise QualityError(
                "没有产生任何 patch（可能：patch 通道细胞过少、--min-cluster-size 过大、"
                "或分组内全为 HDBSCAN 噪声）—— REPRO_SPEC 附录 C 静默失败")
        if rings.shape[0] == 0:
            raise QualityError(
                "patch 已找到但没有任何环带计数（环带内无其他细胞？请检查 --radii-um 与坐标单位）")

        outdir = os.path.dirname(os.path.abspath(args.out_prefix))
        os.makedirs(outdir, exist_ok=True)
        s_path = "%s_patch_summary.csv" % args.out_prefix
        r_path = "%s_ring_counts.csv" % args.out_prefix
        summary.to_csv(s_path, index=False)
        rings.to_csv(r_path, index=False)
        for pp in (s_path, r_path):
            if (not os.path.exists(pp)) or os.path.getsize(pp) == 0:
                raise QualityError("产物为空：%s（REPRO_SPEC 附录 C 静默失败）" % pp)

        hull_methods = sorted(summary["hull_method"].unique().tolist())
        cl_methods = sorted(summary["cluster_method"].unique().tolist())
        print("[INFO] 细胞数=%d；patch 通道=%s=='%s'；patch 数=%d；分组数=%d"
              % (len(df), patch_col, args.patch_value, summary.shape[0],
                 summary["group"].nunique()))
        print("[INFO] 聚类路径：%s" % " | ".join(cl_methods))
        print("[INFO] 边界路径（hull_method）：%s" % " | ".join(hull_methods))
        print("[INFO] patch 细胞数范围 %d–%d；边缘细胞数范围 %d–%d；anchor 数范围 %d–%d"
              % (summary["n_patch_cells"].min(), summary["n_patch_cells"].max(),
                 summary["n_edge_cells"].min(), summary["n_edge_cells"].max(),
                 summary["n_anchor_cells"].min(), summary["n_anchor_cells"].max()))
        print("[INFO] 环带=%s μm；density 单位 = 计数/μm²；用时 %.2f s"
              % (",".join("%g" % r for r in radii), dt))
        print("[INFO] patch 汇总：")
        for line in summary[["group", "patch_id", "n_patch_cells", "n_edge_cells",
                             "n_anchor_cells", "hull_method"]].to_string(
                index=False).splitlines():
            print("       " + line)
        print("[OK] 写出：%s\n[OK] 写出：%s" % (s_path, r_path))
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
