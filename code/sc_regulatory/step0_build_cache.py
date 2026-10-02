#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
step0_build_cache.py — 单细胞调控分析：缓存构建（GSE285335）
================================================================================
目的
----
把 26 个 10x 样本（gene x cell 稀疏阵）一次性处理为**可复用的精简缓存**，供后续
  ① 无监督图谱（step1）、② 上游调控子推断（step2）、③ motif×CpG 交叉（step3）
使用。全程 scipy/pandas，不依赖 scanpy/numba。

两遍流式读取（每遍 26 文件），控制内存
----
  Pass 1：逐样本 mmread → 每细胞 QC（总 counts / 检出基因数 / 线粒体比例）
          → 用**该样本内**的每细胞 library size 做 CPM+log1p → 累积每基因
          sum / sumsq（用于 HVG）
  Pass 2：用全局选定的基因面板（HVG ∪ 标记 ∪ 轴基因 ∪ CollecTRI TF）重建
          genes x cells 稀疏阵

输出（reproducibility/data/sc_regulatory/）
----
  cell_qc.tsv           每细胞：sample, group, barcode, total_counts, n_genes, mito_frac, keep
  gene_stats.tsv        每基因：gene, mean, var, dispersion, mean_bin, disp_z, is_hvg
  hvg.txt               选定的 HVG 基因名
  panel_genes.txt       最终面板基因名
  panel_matrix.npz      面板 genes x cells 稀疏阵（float32，log1p CPM）
  panel_meta.tsv        面板矩阵的细胞顺序元数据
  build_log.json        运行摘要（细胞数、基因数、参数）

依赖
----
  第三方资源 **CollecTRI**（OmniPath 的 TF–靶基因网络）属"不随仓库重分发"项：
  缺失时由 `code/lib/collectri.py` 自动从 OmniPath 获取（该模块会打印与论文
  快照的版本差异告警）。
"""

# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
_sys.path.insert(0, _paths.LIB)                  # 随包分发的辅助模块
from collectri import load as load_collectri  # noqa: E402  （CollecTRI 缺失时自动获取）
# --------------------------------------------------------------------------
import os
import gzip
import json
import time
import argparse

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.io as sio

BASE = _paths.LEGACY
EXTRACT = BASE("00_rawdata", "GSE285335", "extracted")
DATA_DIR = BASE("阶段4.5_复现包", "reproducibility", "data", "sc_regulatory")
os.makedirs(DATA_DIR, exist_ok=True)

MIN_COUNTS = 500          # 与稿件 Methods 一致：library size <= 500 移除
MIN_GENES = 200           # 追加最低检出基因数（近乎不删，记录用）
N_HVG = 2000

# ---------------------------------------------------------------------------
# 基因面板（标记 + 轴基因）
# ---------------------------------------------------------------------------
AXIS_GENES = ["C1GALT1", "C1GALT1C1", "GALNT2", "GALNT12", "ST6GALNAC2"]

MARKERS = {
    "B_cell":     ["MS4A1", "CD19", "CD79A", "CD79B", "CD37", "BANK1", "CD74",
                   "TNFRSF13C", "HLA-DRA"],
    "naive_B":    ["IGHD", "IGHM", "TCL1A", "FCER2", "IL4R"],
    "memory_B":   ["CD27", "TNFRSF13B", "IGHG1", "IGHG2", "IGHA1"],
    "plasma":     ["MZB1", "JCHAIN", "XBP1", "PRDM1", "DERL3", "TNFRSF17", "CD38", "IGKC"],
    "plasmablast":["SDC1", "AICDA", "CD38"],
    "T_cell":     ["CD3D", "CD3E", "CD3G", "TRAC", "IL7R", "CD2", "CD28", "LTB"],
    "CD4_T":      ["CD4", "IL7R", "CD40LG", "CTLA4", "FOXP3"],
    "CD8_T":      ["CD8A", "CD8B", "GZMK", "CCL5", "GZMH"],
    "NK":         ["NKG7", "GNLY", "KLRD1", "KLRF1", "KLRB1", "FGFBP2", "PRF1", "NCAM1"],
    "monocyte":   ["CD14", "LYZ", "FCN1", "S100A8", "S100A9", "VCAN", "MNDA",
                   "CSF1R", "CD68", "FCGR3A", "MS4A7"],
    "DC":         ["FCER1A", "CST3", "CD1C", "CLEC10A", "CLEC9A", "CD1E", "IRF8", "BATF3"],
    "pDC":        ["LILRA4", "CLEC4C", "IL3RA", "SERPINF1", "PLD4", "TCF4"],
    "platelet":   ["PPBP", "PF4", "ITGA2B", "GP9", "TUBB1", "NRGN"],
    "HSPC":       ["CD34", "SPINK2", "AVP", "CRHBP", "PRSS57"],
    "prolifer":   ["MKI67", "TOP2A", "STMN1", "TYMS"],
}


def sample_stems():
    return sorted(f[: -len("_barcodes.tsv.gz")]
                  for f in os.listdir(EXTRACT) if f.endswith("_barcodes.tsv.gz"))


def group_of(stem):
    tag = stem.split("_")[-1]
    if tag.startswith("L"):
        return "Late"
    if tag.startswith("E"):
        return "Early"
    if tag.startswith("H"):
        return "HC"
    return "Unknown"


def read_genes_and_barcodes(stem):
    gsm = stem.split("_")[0]
    ff = os.path.join(EXTRACT, gsm + "_features.tsv.gz")
    bf = os.path.join(EXTRACT, stem + "_barcodes.tsv.gz")
    feats = pd.read_csv(ff, sep="\t", header=None, dtype=str)
    genes = feats.iloc[:, 1].astype(str).tolist()
    first = feats.iloc[:, 0].astype(str).tolist()   # ensembl id（诊断用）
    with gzip.open(bf, "rt") as f:
        barcodes = [ln.strip() for ln in f]
    return genes, first, barcodes


def read_matrix(stem):
    mf = os.path.join(EXTRACT, stem + "_matrix.mtx.gz")
    with gzip.open(mf, "rt") as f:
        mat = sio.mmread(f)
    return sp.csc_matrix(mat)      # genes x cells, csc 便于列操作


def col_sums_csc(mat):
    return np.asarray(mat.sum(axis=0)).ravel()


def build(limit=None):
    stems = sample_stems()
    if limit:
        stems = stems[:limit]
    t0 = time.time()
    print(f"[step0] 样本数 = {len(stems)}")

    genes_ref = None
    # ---- Pass 1 ----
    qc_rows = []
    g_sum = g_sq = None
    n_total = 0
    for k, stem in enumerate(stems, 1):
        ts = time.time()
        genes, ens, barcodes = read_genes_and_barcodes(stem)
        if genes_ref is None:
            genes_ref = genes
        elif genes != genes_ref:
            raise ValueError(f"{stem} 基因列表与其他样本不一致")

        mat = read_matrix(stem)                       # genes x cells
        tot = col_sums_csc(mat)                       # 每细胞 library size
        ncell = mat.shape[1]
        ngenes = np.diff(mat.indptr)                  # csc: 每列 nnz = 检出基因数

        is_mito = np.array([g.startswith("MT-") for g in genes])
        mito_cnt = np.asarray(mat[is_mito, :].sum(axis=0)).ravel()
        mito_frac = np.divide(mito_cnt, tot, out=np.zeros_like(mito_cnt, dtype=float),
                              where=tot > 0)
        keep = tot >= MIN_COUNTS
        # 该样本内 CPM + log1p（显式按列缩放，避免稀疏广播歧义）
        inv = np.divide(1e6, tot, out=np.zeros_like(tot, dtype=float), where=tot > 0)
        coo = mat.tocoo()
        norm = sp.csr_matrix(
            (np.log1p(coo.data * inv[coo.col]), (coo.row, coo.col)),
            shape=mat.shape)

        if g_sum is None:
            g_sum = np.zeros(len(genes), dtype=np.float64)
            g_sq = np.zeros(len(genes), dtype=np.float64)
        # 只累计 keep 的细胞
        idx = np.where(keep)[0]
        sub = norm[:, idx]
        g_sum += np.asarray(sub.sum(axis=1)).ravel()
        g_sq += np.asarray(sub.multiply(sub).sum(axis=1)).ravel()
        n_total += int(keep.sum())

        for i in range(ncell):
            qc_rows.append((stem, group_of(stem), barcodes[i], float(tot[i]),
                            int(ngenes[i]), float(mito_frac[i]), bool(keep[i])))
        print(f"  [{k:2d}/{len(stems)}] {stem:18s} {ncell:6d} cells  "
              f"keep={int(keep.sum()):6d}  {time.time()-ts:5.1f}s")

    qc = pd.DataFrame(qc_rows, columns=["sample", "group", "barcode",
                                        "total_counts", "n_genes", "mito_frac", "keep"])
    qc.to_csv(os.path.join(DATA_DIR, "cell_qc.tsv"), sep="\t", index=False)
    print(f"[pass1] 总细胞 {len(qc):,}  保留 {n_total:,}  ({100*n_total/len(qc):.1f}%)  "
          f"{time.time()-t0:.0f}s")

    # ---- HVG ----
    mean = g_sum / n_total
    var = np.maximum(g_sq / n_total - mean ** 2, 0.0)
    disp = np.divide(var, mean, out=np.zeros_like(var), where=mean > 0)
    gs = pd.DataFrame({"gene": genes_ref, "mean": mean, "var": var, "dispersion": disp})
    # 分箱 z-score（seurat 风格）
    gs["mean_bin"] = pd.cut(gs["mean"], bins=20, labels=False)
    z = np.zeros(len(gs))
    for b, sub in gs.groupby("mean_bin"):
        d = sub["dispersion"].to_numpy()
        sd = d.std()
        z[sub.index] = (d - d.mean()) / sd if sd > 0 else 0.0
    gs["disp_z"] = z
    gs["is_hvg"] = False
    top = gs.sort_values("disp_z", ascending=False).head(N_HVG).index
    gs.loc[top, "is_hvg"] = True
    gs.to_csv(os.path.join(DATA_DIR, "gene_stats.tsv"), sep="\t", index=False)

    hvg = gs.loc[gs["is_hvg"], "gene"].tolist()
    with open(os.path.join(DATA_DIR, "hvg.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(hvg))

    # ---- 面板 ----
    # CollecTRI（OmniPath；第三方资源、不随仓库分发）：缺失时经 code/lib/collectri.py
    # 自动获取；文件已存在时读取路径与结果与原先完全一致。
    ctri = load_collectri(data_dir=DATA_DIR)
    tf_set = set(ctri["source_genesymbol"].dropna().astype(str)) | \
             set(ctri["target_genesymbol"].dropna().astype(str))
    marker_set = {m for v in MARKERS.values() for m in v}
    want = set(hvg) | marker_set | set(AXIS_GENES) | tf_set
    panel = [g for g in genes_ref if g in want]
    with open(os.path.join(DATA_DIR, "panel_genes.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(panel))
    print(f"[panel] HVG={len(hvg)} 标记={len(marker_set)} 轴={len(AXIS_GENES)} "
          f"CollecTRI基因={len(tf_set)}  ->  面板={len(panel)} 基因")

    # ---- Pass 2：重建面板矩阵 ----
    gidx = {g: i for i, g in enumerate(genes_ref)}
    rows_all, cols_all, data_all = [], [], []
    meta_rows = []
    col_off = 0
    for k, stem in enumerate(stems, 1):
        genes, ens, barcodes = read_genes_and_barcodes(stem)
        mat = read_matrix(stem)
        tot = col_sums_csc(mat)
        keep = tot >= MIN_COUNTS
        idx = np.where(keep)[0]
        sub = mat[:, idx].tocsc()
        inv = np.divide(1e6, tot[idx], out=np.zeros(int(keep.sum())), where=tot[idx] > 0)

        # 选面板行
        row_sel = np.array([gidx[g] for g in panel])
        subp = sub[row_sel, :].tocoo()
        vals = subp.data * inv[subp.col]
        vals = np.log1p(vals)
        rows_all.append(subp.row.astype(np.int32))
        cols_all.append((subp.col + col_off).astype(np.int32))
        data_all.append(vals.astype(np.float32))
        for j in idx:
            meta_rows.append((stem, group_of(stem), barcodes[j], col_off + j))
        col_off += int(keep.sum())
        print(f"  [p2 {k:2d}/{len(stems)}] {stem:18s} panel-nnz={subp.nnz:,}")

    R = np.concatenate(rows_all); C = np.concatenate(cols_all); D = np.concatenate(data_all)
    M = sp.csr_matrix((D, (R, C)), shape=(len(panel), int(col_off)), dtype=np.float32)
    sp.save_npz(os.path.join(DATA_DIR, "panel_matrix.npz"), M)
    pd.DataFrame(meta_rows, columns=["sample", "group", "barcode", "cell_idx"]).to_csv(
        os.path.join(DATA_DIR, "panel_meta.tsv"), sep="\t", index=False)

    log = {
        "n_samples": len(stems), "n_cells_total": int(len(qc)),
        "n_cells_kept": int(n_total), "n_genes_total": len(genes_ref),
        "n_hvg": len(hvg), "n_panel": len(panel), "panel_nnz": int(M.nnz),
        "min_counts": MIN_COUNTS, "min_genes": MIN_GENES,
        "axis_present": [g for g in AXIS_GENES if g in gidx],
        "axis_missing": [g for g in AXIS_GENES if g not in gidx],
        "elapsed_sec": round(time.time() - t0, 1),
    }
    with open(os.path.join(DATA_DIR, "build_log.json"), "w", encoding="utf-8") as f:
        json.dump(log, f, indent=2, ensure_ascii=False)
    print("[done]", json.dumps(log, ensure_ascii=False))
    return log


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()
    build(limit=a.limit)
