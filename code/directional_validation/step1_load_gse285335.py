#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
方向性验证：GSE285335 IgAN PBMC 单细胞 B 细胞亚群中 C1GALT1/C1GALT1C1 表达轨迹
================================================================================
目的：
  在公开可复现的 IgAN 患者外周血 PBMC 单细胞数据中，验证半乳糖化关键酶
  C1GALT1（及伴侣 C1GALT1C1）是否在 B 细胞亚群中存在特异表达 / 差异表达，
  把"机制轴驱动 Gd-IgA1 定量"落到细胞类型层面（B 细胞是 Gd-IgA1 产生者）。

数据：GSE285335（Kim G, Sci Rep 2025），26 个 IgAN PBMC 10x 样本
  分组：L1-L6 晚期 IgAN；E1-E11 早期 IgAN；H1-H9 健康对照

实现说明：
  纯 scipy + pandas + matplotlib，不依赖 scanpy/numba（规避 numba DLL 问题）。

输出：
  - 各 B 细胞亚群 / 全 PBMC 中 C1GALT1、C1GALT1C1 表达统计（CSV）
  - 表达轨迹图（PNG/PDF/TIFF）
"""

# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os
import gzip
import tarfile
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

BASE = _paths.LEGACY
RAW_TAR = BASE("00_rawdata", "GSE285335", "GSE285335_RAW.tar")
EXTRACT_DIR = BASE("00_rawdata", "GSE285335", "extracted")
OUT_DIR = _paths.results("directional_validation")
FIG_DIR = _paths.results("directional_validation")  # 原目的地为投稿材料 figures/
for d in (EXTRACT_DIR, OUT_DIR, FIG_DIR):
    os.makedirs(d, exist_ok=True)


def read_10x_triplet(dirpath, stem):
    """读取某样本的 10x 三件套，返回 gene_names(list), barcodes(list), csr(genes x cells)."""
    import scipy.io as sio
    import scipy.sparse as sp

    bf = os.path.join(dirpath, stem + "_barcodes.tsv.gz")
    ff = None
    mf = os.path.join(dirpath, stem + "_matrix.mtx.gz")
    # features 可能是 features.tsv.gz 或 genes.tsv.gz
    for cand in (stem + "_features.tsv.gz", stem + "_genes.tsv.gz"):
        p = os.path.join(dirpath, cand)
        if os.path.exists(p):
            ff = p
            break
    if ff is None:
        raise FileNotFoundError(f"{stem} 缺 features/genes 文件")

    with gzip.open(bf, "rt") as f:
        barcodes = [ln.strip() for ln in f]

    feats = pd.read_csv(ff, sep="\t", header=None)
    gene_names = feats.iloc[:, 1].astype(str).tolist()

    with gzip.open(mf, "rt") as f:
        mat = sio.mmread(f)
    return gene_names, barcodes, sp.csc_matrix(mat)


def load_all():
    """解压并加载所有样本，返回 (gene_names, cell_meta_df, X csr genes x cells)."""
    if not os.path.isdir(EXTRACT_DIR) or not os.listdir(EXTRACT_DIR):
        print(f"解压 {os.path.basename(RAW_TAR)} ...")
        with tarfile.open(RAW_TAR, "r") as tar:
            tar.extractall(EXTRACT_DIR)

    # 枚举样本 stem
    stems = set()
    for fn in os.listdir(EXTRACT_DIR):
        if fn.endswith("_barcodes.tsv.gz"):
            stems.add(fn.replace("_barcodes.tsv.gz", ""))

    import scipy.sparse as sp

    gene_global = None
    mats = []
    meta_rows = []
    for stem in sorted(stems):
        gene_names, barcodes, mat = read_10x_triplet(EXTRACT_DIR, stem)
        if gene_global is None:
            gene_global = gene_names
        elif gene_names != gene_global:
            # 理论上 10x 多样本共用 features；若不一致则中止以保正确
            raise ValueError(f"{stem} 基因列表与其他样本不一致，需人工对齐")

        grp = stem.split("_")[-1]  # L1/E1/H1
        if grp.startswith("L"):
            group = "Late"
        elif grp.startswith("E"):
            group = "Early"
        elif grp.startswith("H"):
            group = "HC"
        else:
            group = "Unknown"

        mats.append(mat)
        for bc in barcodes:
            meta_rows.append({"sample": stem, "group": group, "barcode": bc})

    X = sp.hstack(mats, format="csr")  # genes x cells
    meta = pd.DataFrame(meta_rows)
    return gene_global, meta, X


def main():
    print("=" * 72)
    print("方向性验证：GSE285335 B 细胞 C1GALT1/C1GALT1C1 表达轨迹")
    print("=" * 72)

    if not os.path.exists(RAW_TAR):
        print(f"\n[错误] 找不到 {RAW_TAR}")
        print("请先按 00_rawdata/GSE285335/下载说明.md 下载数据并放到该路径。")
        return

    gene_names, meta, X = load_all()
    n_genes, n_cells = X.shape
    print(f"\n数据加载完成：{n_genes} 基因 x {n_cells} 细胞")
    print(f"分组分布：\n{meta['group'].value_counts().to_string()}")

    gene_idx = {g: i for i, g in enumerate(gene_names)}

    # ---- 目标基因 + B 细胞 marker ----
    target = ["C1GALT1", "C1GALT1C1"]
    b_markers = {
        "pan-B": ["MS4A1", "CD19", "CD79A", "CD79B"],
        "naive": ["IGHD", "IGHM", "TCL1A", "FCER2"],
        "memory": ["CD27", "TNFRSF13B", "CD40"],
        "plasma": ["SDC1", "MZB1", "JCHAIN", "XBP1"],
        "germinal": ["BCL6", "AICDA", "TNFRSF13C"],
    }
    print("\n目标/标记基因是否存在：")
    for g in target + [m for v in b_markers.values() for m in v]:
        print(f"  {g:12s}: {'在' if g in gene_idx else '缺失'}")

    # ---- 粗分 B 细胞（任一 pan-B marker 检出即判为 B 细胞）----
    def gene_col(g):
        return np.asarray(X[gene_idx[g], :].todense()).ravel() if g in gene_idx else None

    panB_genes = [g for g in b_markers["pan-B"] if g in gene_idx]
    panB_mat = np.zeros(n_cells, dtype=bool)
    for g in panB_genes:
        panB_mat |= (gene_col(g) > 0)
    meta["is_B"] = panB_mat

    n_b = int(panB_mat.sum())
    print(f"\nB 细胞（任一 pan-B marker 检出）：{n_b} / {n_cells} "
          f"({100*n_b/n_cells:.1f}%)")

    # ---- B 细胞亚群划分（基于 marker 检出，允许重叠）----
    for sub, gs in b_markers.items():
        if sub == "pan-B":
            continue
        gs_ok = [g for g in gs if g in gene_idx]
        mask = np.zeros(n_cells, dtype=bool)
        for g in gs_ok:
            mask |= (gene_col(g) > 0)
        meta["B_" + sub] = mask & panB_mat

    # ---- 目标基因在各分组/亚群的表达统计 ----
    rows = []
    for g in target:
        if g not in gene_idx:
            continue
        col = gene_col(g)
        # 全 PBMC 按分组
        for grp in ["HC", "Early", "Late"]:
            m = (meta["group"] == grp).to_numpy()
            if m.sum() == 0:
                continue
            rows.append({
                "gene": g, "subset": f"PBMC-{grp}", "group": grp,
                "n_cells": int(m.sum()),
                "frac_detected": float((col[m] > 0).mean()),
                "mean_expr": float(col[m].mean()),
            })
        # B 细胞整体 + 亚群
        b_mask = meta["is_B"].to_numpy()
        rows.append({
            "gene": g, "subset": "B_all", "group": "all",
            "n_cells": int(b_mask.sum()),
            "frac_detected": float((col[b_mask] > 0).mean()),
            "mean_expr": float(col[b_mask].mean()),
        })
        for sub in ["B_naive", "B_memory", "B_plasma", "B_germinal"]:
            m = meta[sub].to_numpy()
            if m.sum() == 0:
                continue
            rows.append({
                "gene": g, "subset": sub, "group": "all",
                "n_cells": int(m.sum()),
                "frac_detected": float((col[m] > 0).mean()),
                "mean_expr": float(col[m].mean()),
            })

    df = pd.DataFrame(rows)
    csv_path = os.path.join(OUT_DIR, "C1GALT1_Bcell_expression_summary.csv")
    df.to_csv(csv_path, index=False)
    print(f"\n表达统计已写：{csv_path}")
    print(df.to_string(index=False))

    # ---- 绘图 ----
    plot_trajectory(df, meta, gene_col, gene_idx, target)

    print("\n完成。")


def plot_trajectory(df, meta, gene_col, gene_idx, target):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # 子集顺序：PBMC（各分组）+ B 细胞亚群
    order = ["PBMC-HC", "PBMC-Early", "PBMC-Late",
             "B_all", "B_naive", "B_memory", "B_plasma", "B_germinal"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for gi, g in enumerate(target):
        ax = axes[gi]
        sub = df[df["gene"] == g].set_index("subset")
        labels = []
        fracs = []
        means = []
        for s in order:
            if s in sub.index:
                labels.append(s)
                fracs.append(sub.loc[s, "frac_detected"] * 100)
                means.append(sub.loc[s, "mean_expr"])

        x = np.arange(len(labels))
        bars = ax.bar(x, fracs, color="#4C72B0", alpha=0.85)
        ax.set_ylabel("% cells expressing", fontsize=10)
        ax.set_title(f"{g} detection across PBMC and B-cell subsets",
                     fontsize=11, fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
        # 标注 mean expr
        for i, m in enumerate(means):
            ax.text(x[i], fracs[i] + 0.5, f"{m:.2f}", ha="center",
                    va="bottom", fontsize=7, color="#333333")
        ax.set_ylim(0, max(fracs) * 1.2 + 2)
        ax.grid(axis="y", alpha=0.3)

    fig.suptitle("Directional validation: C1GALT1 / C1GALT1C1 expression "
                 "in PBMCs and B-cell subsets (GSE285335)",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    for ext in ["png", "pdf"]:
        out = os.path.join(FIG_DIR, f"FigS4_C1GALT1_Bcell_trajectory.{ext}")
        fig.savefig(out, dpi=300 if ext == "png" else None,
                    bbox_inches="tight")
        print(f"  图已写：{out}")

    # TIFF（PIL 600dpi LZW）
    try:
        from PIL import Image
        png = os.path.join(FIG_DIR, "FigS4_C1GALT1_Bcell_trajectory.png")
        im = Image.open(png)
        tiff = os.path.join(FIG_DIR, "FigS4_C1GALT1_Bcell_trajectory.tiff")
        im.save(tiff, format="TIFF", compression="tiff_lzw",
                dpi=(600, 600))
        print(f"  图已写：{tiff}")
    except Exception as e:
        print(f"  [警告] TIFF 生成失败：{e}")


if __name__ == "__main__":
    main()
