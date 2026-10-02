#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
step1_atlas.py — 单细胞 PBMC 图谱 + 轴基因细胞类型分布（GSE285335）
================================================================================
标准单细胞流程的**可扩展**实现（不依赖 scanpy/leidenalg 的原生扩展）：

  (1) 全 282,463 细胞：HVG（去技术基因）z 标准化 → 截断 SVD(PCA, 50) →
      MiniBatchKMeans 微观簇(600) → 微观簇质心 kNN 图 → Louvain 社区(res=2)
  (2) 对每个社区用规范标记表达谱注释主要 PBMC 谱系 → **全部细胞**获得细胞类型标签
  (3) B 区室（B + Plasma）内再细分 naive / memory / plasma / other
  (4) 量化 5 个轴基因在 各细胞类型 与 B 亚群 的检出率与表达量
  (5) 分层抽样 30,000 细胞做 2-D 嵌入（UMAP，失败退 t-SNE）用于图谱图

红线：只回答"轴基因在何种细胞状态表达"，**不用于疾病方向宣称**。
"""

# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os
import re
import json
import argparse

import numpy as np
import pandas as pd
import scipy.sparse as sp

BASE = _paths.LEGACY
DATA_DIR = BASE("阶段4.5_复现包", "reproducibility", "data", "sc_regulatory")
RES_DIR = BASE("阶段4.5_复现包", "reproducibility", "results", "sc_regulatory")
FIG_DIR = _paths.results("sc_regulatory")  # 原目的地为投稿材料 additional_file_1
for d in (RES_DIR, FIG_DIR):
    os.makedirs(d, exist_ok=True)

AXIS = ["C1GALT1", "C1GALT1C1", "GALNT2", "GALNT12", "ST6GALNAC2"]

# 从 HVG 中剔除技术性/广谱管家基因（核糖体、线粒体、血红蛋白、应激即刻基因等），
# 否则 PCA 会被其主导，产生与谱系无关的人为结构。
TECH = re.compile(
    r"^(RP[LS]\d|MRP[LS]\d|MT-|MALAT1$|NEAT1$|HB[ABDEGQZ]\d|TMSB\d|EEF\d|EIF\d|"
    r"FTL$|FTH1$|B2M$|ACTB$|GAPDH$|TUBB$|TUBA1B$|FOS$|FOSB$|JUN$|JUNB$|EGR1$|"
    r"HSP\d|DNAJ|SERPINA1$|LGALS1$|S100A|HLA-|IGHV|IGKV|IGLV|TRBV|TRGV|"
    r"XIST$|MTRNR|RPSA$|PTMA$|NPM1$|UBA52$|UBB$|UBC$)")

# 更特异的谱系标记（避免 GZMB/IRF7 这类广谱基因造成 pDC 假富集）
LINEAGE = {
    "B":        ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1", "CD37", "TNFRSF13C"],
    "Plasma":   ["MZB1", "JCHAIN", "XBP1", "DERL3", "TNFRSF17", "IGKC", "CD38"],
    "T":        ["CD3D", "CD3E", "CD3G", "TRAC"],
    "NK":       ["NKG7", "GNLY", "KLRF1", "FGFBP2", "KLRD1"],
    "Mono":     ["CD14", "LYZ", "FCN1", "S100A8", "S100A9", "VCAN", "MNDA", "CSF1R"],
    "DC":       ["FCER1A", "CD1C", "CLEC10A", "CLEC9A", "CD1E"],
    "pDC":      ["LILRA4", "CLEC4C", "IL3RA", "SERPINF1", "PLD4"],
    "Platelet": ["PPBP", "PF4", "ITGA2B", "GP9", "TUBB1"],
    "HSPC":     ["CD34", "SPINK2", "AVP", "CRHBP"],
}

B_SUBSET = {
    "B_naive":  ["IGHD", "IGHM", "TCL1A", "FCER2"],
    "B_memory": ["CD27", "TNFRSF13B", "IGHG1", "IGHA1"],
    "B_plasma": ["MZB1", "JCHAIN", "XBP1", "DERL3", "TNFRSF17", "IGKC", "CD38"],
}


def load_panel():
    M = sp.load_npz(os.path.join(DATA_DIR, "panel_matrix.npz")).tocsr()
    genes = [g for g in open(os.path.join(DATA_DIR, "panel_genes.txt"), encoding="utf-8").read().split("\n") if g]
    meta = pd.read_csv(os.path.join(DATA_DIR, "panel_meta.tsv"), sep="\t")
    hvg = [g for g in open(os.path.join(DATA_DIR, "hvg.txt"), encoding="utf-8").read().split("\n") if g]
    return M, genes, meta, hvg


def mean_of_sets(X, genes, sets):
    gi = {g: i for i, g in enumerate(genes)}
    cols, names = [], []
    for name, gl in sets.items():
        idx = [gi[g] for g in gl if g in gi]
        names.append(name)
        cols.append(np.asarray(X[idx, :].mean(axis=0)).ravel() if idx
                    else np.zeros(X.shape[1], dtype=np.float32))
    return np.vstack(cols).T, names          # cells x sets


def louvain_small(adj, seed=0, resolution=2.0):
    """仅用于 ~数百个微观簇的小图。resolution 越高，社区数越多。"""
    try:
        import igraph as ig, leidenalg as la
        g = ig.Graph(n=adj.shape[0], edges=list(zip(*adj.nonzero())))
        return np.array(la.find_partition(
            g, la.RBConfigurationVertexPartition, seed=seed,
            resolution_parameter=resolution).membership)
    except Exception:
        pass
    import networkx as nx
    G = nx.from_scipy_sparse_array(adj)
    comms = nx.community.louvain_communities(G, seed=seed, resolution=resolution)
    lab = np.zeros(adj.shape[0], dtype=int)
    for c, nodes in enumerate(comms):
        lab[list(nodes)] = c
    return lab


def embed(coords, seed=0):
    try:
        import umap
        print("  [embed] UMAP")
        return umap.UMAP(n_neighbors=15, min_dist=0.3, random_state=seed).fit_transform(coords), "UMAP"
    except Exception as e:
        print(f"  [embed] t-SNE（UMAP 不可用：{type(e).__name__}）")
        from sklearn.manifold import TSNE
        return TSNE(n_components=2, perplexity=30, learning_rate="auto", max_iter=750,
                    init="pca", random_state=seed).fit_transform(coords), "t-SNE"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nsub", type=int, default=30000)
    ap.add_argument("--n-micro", type=int, default=600)
    ap.add_argument("--resolution", type=float, default=2.0)
    ap.add_argument("--no-fig", action="store_true")
    a = ap.parse_args()

    print("=" * 74)
    print("step1 单细胞 PBMC 图谱（全细胞社区 + 标记注释）")
    print("=" * 74)
    M, genes, meta, hvg = load_panel()
    gi = {g: i for i, g in enumerate(genes)}
    ncell = M.shape[1]
    print(f"面板：{M.shape[0]} 基因 x {ncell} 细胞 | HVG={len(hvg)}")

    # ---------- 全细胞 PCA ----------
    hvg_use = [g for g in hvg if not TECH.match(g)]
    print(f"HVG 去技术基因后：{len(hvg_use)} / {len(hvg)}")
    hvg_i = [gi[g] for g in hvg_use if g in gi]
    Xh = M[hvg_i, :].T.tocsr().astype(np.float32)          # cells x HVG
    mu = np.asarray(Xh.mean(axis=0)).ravel()
    sq = np.asarray(Xh.multiply(Xh).mean(axis=0)).ravel()
    sd = np.sqrt(np.maximum(sq - mu ** 2, 1e-8))
    coo = Xh.tocoo()
    Xz = sp.csr_matrix(((coo.data - mu[coo.col]) / sd[coo.col], (coo.row, coo.col)),
                       shape=Xh.shape)
    np.clip(Xz.data, -10, 10, out=Xz.data)
    del Xh
    from sklearn.decomposition import TruncatedSVD
    svd = TruncatedSVD(n_components=50, random_state=0)
    pcs = svd.fit_transform(Xz)
    print(f"  PCA：{pcs.shape}，前 50 PC 解释方差比 = {svd.explained_variance_ratio_.sum():.3f}")
    del svd

    # ---------- 微观簇 → 社区 ----------
    from sklearn.cluster import MiniBatchKMeans
    from sklearn.neighbors import NearestNeighbors
    km = MiniBatchKMeans(n_clusters=a.n_micro, random_state=0, n_init=5, batch_size=8192)
    micro = km.fit_predict(pcs)
    cent = km.cluster_centers_
    nn = NearestNeighbors(n_neighbors=13).fit(cent)
    dist, ind = nn.kneighbors(cent)
    rows, cols, vals = [], [], []
    for i in range(a.n_micro):
        for j, d in zip(ind[i, 1:], dist[i, 1:]):
            rows.append(i); cols.append(j); vals.append(1.0 / (1.0 + d))
    Am = sp.csr_matrix((vals, (rows, cols)), shape=(a.n_micro, a.n_micro))
    Am = Am.maximum(Am.T)
    lab_micro = louvain_small(Am, resolution=a.resolution)
    cell_cluster = lab_micro[micro]
    ncomm = len(set(cell_cluster))
    print(f"  微观簇={a.n_micro} → 社区={ncomm}")

    # ---------- 细胞类型：透明门控规则（互斥；稀有/特异谱系优先） ----------
    def col(g):
        return (np.asarray(M[gi[g], :].todense()).ravel() if g in gi
                else np.zeros(ncell, dtype=np.float32))

    CD3D, MS4A1, CD79A = col("CD3D"), col("MS4A1"), col("CD79A")
    MZB1, JCHAIN, DERL3 = col("MZB1"), col("JCHAIN"), col("DERL3")
    CD14, LYZ, FCN1 = col("CD14"), col("LYZ"), col("FCN1")
    FCER1A, CD1C, LILRA4, CLEC4C = col("FCER1A"), col("CD1C"), col("LILRA4"), col("CLEC4C")
    PPBP, PF4, NKG7, GNLY, CD34 = col("PPBP"), col("PF4"), col("NKG7"), col("GNLY"), col("CD34")

    ct = np.array(["Other"] * ncell, dtype=object)
    ct[(NKG7 > 0) | (GNLY > 0)] = "NK"                                  # 细胞毒效应
    ct[CD3D > 0] = "T"                                                  # CD3 复合体
    ct[(MS4A1 > 0) & (CD79A > 0)] = "B"
    ct[(CD14 > 0) | ((LYZ > 0) & (FCN1 > 0))] = "Mono"
    ct[(FCER1A > 0) | (CD1C > 1)] = "DC"
    ct[(LILRA4 > 0) | (CLEC4C > 0)] = "pDC"
    ct[(MZB1 > 0) | (JCHAIN > 1) | (DERL3 > 0)] = "Plasma"
    ct[(PPBP > 2) | ((PPBP > 0) & (PF4 > 0))] = "Platelet"
    ct[CD34 > 0] = "HSPC"
    cell_type = ct
    print("\n细胞类型（门控规则，互斥）：")
    print(pd.Series(cell_type).value_counts().to_string())

    # ---------- B 亚群 ----------
    sb, bnames = mean_of_sets(M, genes, B_SUBSET)
    b_subset = np.array(["None"] * ncell, dtype=object)
    bmask = np.isin(cell_type, ["B", "Plasma"])
    if bmask.sum():
        sub_scores = sb[bmask]
        bidx = sub_scores.argmax(axis=1)
        b_subset[bmask] = [bnames[i] if sub_scores[k, i] > 0 else "B_other"
                           for k, i in enumerate(bidx)]

    meta = meta.copy()
    meta["cell_type"] = cell_type
    meta["cluster"] = cell_cluster
    meta["b_subset"] = b_subset
    print("\n细胞类型计数：")
    print(meta["cell_type"].value_counts().to_string())
    print("\nB 区室内亚群：")
    print(meta.loc[bmask, "b_subset"].value_counts().to_string())

    # ---------- 轴基因量化 ----------
    rows = []
    for g in AXIS:
        if g not in gi:
            continue
        col = np.asarray(M[gi[g], :].todense()).ravel()
        for kind, key in [("cell_type", "cell_type"), ("b_subset", "b_subset")]:
            for lvl in pd.unique(meta[key]):
                if lvl in ("None",):
                    continue
                m = (meta[key] == lvl).to_numpy()
                pos = col[m] > 0
                rows.append({"gene": g, "group_by": kind, "level": lvl, "n_cells": int(m.sum()),
                             "detection_rate": float(pos.mean()),
                             "mean_log1pCPM": float(col[m].mean()),
                             "mean_nonzero_log1pCPM": float(col[m][pos].mean()) if pos.any() else 0.0})
    ax_df = pd.DataFrame(rows)
    ax_df.to_csv(os.path.join(RES_DIR, "axis_expression_by_celltype.tsv"), sep="\t", index=False)
    print("\n轴基因 × 细胞类型（检出率）：")
    print(ax_df[ax_df.group_by == "cell_type"].pivot(index="level", columns="gene",
                                                     values="detection_rate").round(4).to_string())

    comp = pd.crosstab(meta["group"], meta["cell_type"], normalize="index")
    comp.to_csv(os.path.join(RES_DIR, "celltype_composition_by_group.tsv"), sep="\t")
    meta.to_csv(os.path.join(RES_DIR, "cell_type_assignment.tsv"), sep="\t", index=False)
    pd.crosstab(pd.Series(cell_cluster, name="community"),
                pd.Series(cell_type, name="cell_type")).to_csv(
        os.path.join(RES_DIR, "community_marker_profile.tsv"), sep="\t")

    # ---------- 图（抽样嵌入） ----------
    rng = np.random.default_rng(0)
    nsub = min(a.nsub, ncell)
    sub_idx = np.sort(rng.choice(ncell, size=nsub, replace=False))
    emb_name, coords = "none", None
    if not a.no_fig:
        coords, emb_name = embed(pcs[sub_idx])
    sub = meta.iloc[sub_idx].reset_index(drop=True).copy()
    if coords is not None:
        sub["e1"] = coords[:, 0]; sub["e2"] = coords[:, 1]
        sub.to_csv(os.path.join(RES_DIR, "atlas_cells.tsv"), sep="\t", index=False)

    if not a.no_fig:
        make_fig(sub, ax_df, emb_name)

    summary = {"n_cells": int(ncell), "n_subsample": int(nsub), "n_communities": int(ncomm),
               "embedding": emb_name, "n_hvg_used": len(hvg_i),
               "pca_var_explained_50": float(np.var(pcs, axis=0).sum()),
               "celltype_counts": {k: int(v) for k, v in meta["cell_type"].value_counts().items()},
               "b_subset_counts": {k: int(v) for k, v in meta.loc[bmask, "b_subset"].value_counts().items()}}
    with open(os.path.join(RES_DIR, "step1_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print("\n[done] step1")


def make_fig(sub, ax_df, emb_name):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    order = ["T", "NK", "B", "Mono", "DC", "pDC", "Plasma", "Platelet", "HSPC", "Other"]
    palette = {"T": "#2ca02c", "NK": "#9467bd", "B": "#1f77b4", "Mono": "#ff7f0e",
               "DC": "#8c564b", "pDC": "#e377c2", "Plasma": "#d62728",
               "Platelet": "#7f7f7f", "HSPC": "#bcbd22", "Other": "#d9d9d9"}
    fig = plt.figure(figsize=(15, 5.2))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.05, 1.05], wspace=0.34)

    ax = fig.add_subplot(gs[0, 0])
    for ct in order:
        m = (sub["cell_type"] == ct).to_numpy()
        if m.sum() == 0:
            continue
        ax.scatter(sub.loc[m, "e1"], sub.loc[m, "e2"], s=1.3, c=palette[ct],
                   label=f"{ct} ({int(m.sum()):,})", linewidths=0, rasterized=True)
    ax.set_title(f"PBMC atlas across 26 donors\n(282,463 cells; {emb_name} of 30,000-cell subset)",
                 fontsize=10.5, fontweight="bold")
    ax.set_xlabel(f"{emb_name}-1"); ax.set_ylabel(f"{emb_name}-2")
    ax.legend(fontsize=6.8, markerscale=3, loc="best", frameon=False)
    ax.set_xticks([]); ax.set_yticks([])

    ax = fig.add_subplot(gs[0, 1])
    d = ax_df[ax_df.group_by == "cell_type"]
    cts = [c for c in order if (d["level"] == c).any()]
    cc = [c for c in AXIS if (d["gene"] == c).any()]
    for xi, c in enumerate(cts):
        for yi, g in enumerate(cc):
            r = d[(d["level"] == c) & (d["gene"] == g)]
            if r.empty:
                continue
            dr = float(r["detection_rate"].iloc[0]); ex = float(r["mean_nonzero_log1pCPM"].iloc[0])
            ax.scatter(xi, yi, s=20 + dr * 420, c=[ex], cmap="viridis",
                       vmin=0.7, vmax=1.25, edgecolors="k", linewidths=0.4)
    ax.set_xticks(range(len(cts))); ax.set_xticklabels(cts, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(cc))); ax.set_yticklabels(cc, fontsize=8)
    ax.set_title("Axis-gene expression across cell types\n(dot size = detection rate)",
                 fontsize=10.5, fontweight="bold")
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(0.7, 1.25))
    cb = fig.colorbar(sm, ax=ax, fraction=0.045, pad=0.02)
    cb.set_label("log1p CPM in +cells", fontsize=7.5); cb.ax.tick_params(labelsize=7)

    ax = fig.add_subplot(gs[0, 2])
    d2 = ax_df[ax_df.group_by == "b_subset"]
    bo = [o for o in ["B_naive", "B_memory", "B_plasma", "B_other"] if (d2["level"] == o).any()]
    xw = 0.8 / max(len(cc), 1)
    for k, g in enumerate(cc):
        vals = []
        for o in bo:
            r = d2[(d2["level"] == o) & (d2["gene"] == g)]
            vals.append(float(r["detection_rate"].iloc[0]) if not r.empty else np.nan)
        ax.bar(np.arange(len(bo)) + k * xw - 0.4 + xw / 2, vals, width=xw, label=g)
    ax.set_xticks(range(len(bo)))
    ax.set_xticklabels([o.replace("B_", "") for o in bo], fontsize=8.5)
    ax.set_ylabel("detection rate"); ax.legend(fontsize=7, frameon=False)
    ax.set_title("B-lineage trajectory\naxis-gene detection by B state", fontsize=10.5, fontweight="bold")
    ax.grid(axis="y", alpha=0.3)

    fig.tight_layout()
    for ext in ("png", "pdf"):
        out = os.path.join(FIG_DIR, f"Figure_S5_scRNA_atlas_axis.{ext}")
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(f"  图已写：{out}")
    try:
        from PIL import Image
        im = Image.open(os.path.join(FIG_DIR, "Figure_S5_scRNA_atlas_axis.png"))
        im.save(os.path.join(FIG_DIR, "Figure_S5_scRNA_atlas_axis.tiff"),
                format="TIFF", compression="tiff_lzw", dpi=(600, 600))
    except Exception as e:
        print(f"  [警告] TIFF 失败：{e}")


if __name__ == "__main__":
    main()
