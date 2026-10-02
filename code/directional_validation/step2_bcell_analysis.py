#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Directional validation — Step 2: B-cell subset analysis of GSE285335 (IgAN PBMC scRNA-seq).

Purpose:
  Test whether C1GALT1 / C1GALT1C1 (the core of the Gd-IgA1 O-glycosylation axis)
  show a graded expression trajectory along the B-cell maturation axis
  (naive -> memory -> plasma/plasmablast), and whether their expression differs
  between late-stage IgAN, early-stage IgAN and healthy controls within B cells.

  This is the cell-biology landing point for the paper's mechanism axis
  ("mechanism axis drives Gd-IgA1 quantity in B cells").

Design (fully reproducible, no scanpy/numba needed):
  1. Load 26 samples (L1-L6 late, E1-E11 early, H1-H9 healthy) 10x matrices.
  2. Merge, then CPM-normalise + log1p (sparse-safe, in-memory).
  3. Identify B cells by marker expression (MS4A1/CD19/CD79A/CD79B).
  4. Subset B cells: naive (IGHD+IGHM+, CD27-), memory (CD27+), plasma/blast
     (SDC1+/MZB1+/XBP1+/JCHAIN+).
  5. Plot C1GALT1 / C1GALT1C1 expression across subsets and groups.
  6. Statistical tests (Mann-Whitney U) for group comparisons within B cells.

Data source: GSE285335 (Kim G, Sci Rep 2025), 26 IgAN PBMC samples, 10x 5' v1.
Accession matrix: GSE285335_RAW.tar (barcodes/features/matrix triplets per sample).
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
import glob
import numpy as np
import scipy.sparse as sp
import scipy.io as sio
import scipy.stats as stats

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
RAW_DIR = _paths.at("00_rawdata/GSE285335/extracted")
OUT_DIR = _paths.at("阶段4.5_复现包/reproducibility/results/directional_validation")
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Sample metadata
# ---------------------------------------------------------------------------
GROUP_OF = {}
for i in range(1, 7):
    GROUP_OF[f"L{i}"] = "Late IgAN"
for i in range(1, 12):
    GROUP_OF[f"E{i}"] = "Early IgAN"
for i in range(1, 10):
    GROUP_OF[f"H{i}"] = "Healthy"

# ---------------------------------------------------------------------------
# Marker genes
# ---------------------------------------------------------------------------
GENE_OF_INTEREST = ["C1GALT1", "C1GALT1C1"]
B_MARKERS = ["MS4A1", "CD19", "CD79A", "CD79B"]
NAIVE_MARKERS = ["IGHD", "IGHM"]
MEMORY_MARKER = "CD27"
PLASMA_MARKERS = ["SDC1", "MZB1", "XBP1", "JCHAIN"]


def load_triplet(stem):
    """Load one sample's 10x triplet -> (gene_names, barcodes, sparse counts).

    GSE285335 naming: barcodes/matrix carry the sample suffix (e.g.
    GSM8700986_L1_barcodes.tsv.gz), but the features file is shared per-GSM
    (GSM8700986_features.tsv.gz, identical across samples). The GSM id is the
    numeric part of the stem before the "_L1"/"_E1"/"_H1" suffix.
    """
    barc = [l.rstrip("\n") for l in gzip.open(RAW_DIR(f"{stem}_barcodes.tsv.gz"), "rt")]
    gsm = stem.split("_")[0]  # GSM8700986
    feat = [l.rstrip("\n").split("\t")[0] for l in gzip.open(RAW_DIR(f"{gsm}_features.tsv.gz"), "rt")]
    M = sio.mmread(RAW_DIR(f"{stem}_matrix.mtx.gz")).tocsr()
    return feat, barc, M


def main():
    # ------------------------------------------------------------------
    # 1. Load all 26 samples
    # ------------------------------------------------------------------
    stems = sorted(set(
        os.path.basename(f).split("_barcodes.tsv.gz")[0]
        for f in glob.glob(RAW_DIR("*_barcodes.tsv.gz"))
    ))
    print(f"Loading {len(stems)} samples...")

    all_genes = None
    blocks = []          # list of (sample_stem, sparse_cols)
    sample_label = []    # per-sample stem (repeated per cell)
    ncell_per_sample = []

    for stem in stems:
        feat, barc, M = load_triplet(stem)
        if all_genes is None:
            all_genes = feat
        else:
            assert all_genes == feat, f"Gene order mismatch at {stem}"
        blocks.append(M)
        sample_label.append(stem.split("_")[-1])  # L1 / E3 / H7 ...
        ncell_per_sample.append(M.shape[1])
        print(f"  {stem}: {M.shape[1]} cells")

    # Merge into one genes x cells sparse matrix (stack columns)
    X = sp.hstack(blocks, format="csr")  # genes x cells
    n_genes, n_cells = X.shape
    sample_per_cell = np.repeat(sample_label, ncell_per_sample)
    print(f"\nMerged: {n_genes} genes x {n_cells} cells")

    # ------------------------------------------------------------------
    # 2. Normalise: CPM + log1p  (only for the genes we need — memory-safe)
    # ------------------------------------------------------------------
    # library size per cell (sum of counts), from the integer sparse matrix
    libsize = np.asarray(X.sum(axis=0)).ravel()
    keep_cell = libsize > 500  # basic QC: drop empty/broken droplets
    libsize = libsize[keep_cell]
    sample_per_cell = sample_per_cell[keep_cell]
    n_cells = int(keep_cell.sum())
    print(f"After QC (libsize>500): {n_cells} cells")

    # We only need ~25 genes (C1GALT1/C1GALT1C1 + B/plasma markers).
    # Extract those rows from the integer sparse matrix, drop QC'd cells,
    # then CPM-normalise + log1p ONLY on that small sub-matrix.
    all_markers = sorted(set(
        GENE_OF_INTEREST + B_MARKERS + NAIVE_MARKERS + [MEMORY_MARKER] + PLASMA_MARKERS
    ))
    gene_idx = {g: i for i, g in enumerate(all_genes)}
    row_idx = np.array([gene_idx[g] for g in all_markers if g in gene_idx])
    row_names = [all_genes[i] for i in row_idx]

    # Subset rows and QC'd columns (keep integer type — counts are int32/64)
    Xsub = X[row_idx, :][:, keep_cell].tocsr().astype(np.float64)

    # CPM + log1p per column (cell). For a small matrix, multiply by a
    # dense row-vector broadcast is cheap and memory-safe.
    Xsub = Xsub.multiply((1e4 / libsize)[None, :]).tocsr()
    Xsub.data = np.log1p(Xsub.data)

    # Build a dense genes x cells array (small: ~25 x 282k = ~7M floats = 56 MB)
    G = np.asarray(Xsub.todense())          # shape (n_markers, n_cells)
    expr_map = {g: G[i] for i, g in enumerate(row_names)}

    def expr(gene):
        """Return per-cell log1p-CPM expression vector for a gene."""
        v = expr_map.get(gene)
        if v is None:
            return np.zeros(n_cells)
        return v

    # ------------------------------------------------------------------
    # 3. Identify B cells
    # ------------------------------------------------------------------
    bscore = np.zeros(n_cells)
    for g in B_MARKERS:
        bscore += expr(g)
    # B cell: expresses >=2 of the 4 B markers at log1p CPM >= 0.5
    b_positive_count = sum((expr(g) >= 0.5).astype(int) for g in B_MARKERS)
    is_b = b_positive_count >= 2

    print(f"\nB cells identified: {is_b.sum()} / {n_cells} ({100*is_b.mean():.2f}%)")

    # ------------------------------------------------------------------
    # 4. Subset B cells into naive / memory / plasma
    # ------------------------------------------------------------------
    igd = expr("IGHD"); ighm = expr("IGHM"); cd27 = expr("CD27")
    plasma_score = sum(expr(g) for g in PLASMA_MARKERS)
    naive_mask = (igd >= 0.5) | (ighm >= 0.5)
    plasma_mask = plasma_score >= 2.0

    def assign_subset():
        sub = np.full(n_cells, "other_B", dtype=object)
        sub[is_b & plasma_mask] = "plasma"
        sub[is_b & (~plasma_mask) & naive_mask] = "naive"
        sub[is_b & (~plasma_mask) & (~naive_mask) & (cd27 >= 0.5)] = "memory"
        sub[is_b & (~plasma_mask) & (~naive_mask) & (cd27 < 0.5)] = "naive"  # CD27- naive
        return sub

    subset = assign_subset()
    b_cells = is_b

    # Count summary
    from collections import Counter
    cnt = Counter(subset[b_cells])
    print("B-cell subset counts:", dict(cnt))

    # ------------------------------------------------------------------
    # 5. Expression summary of C1GALT1/C1GALT1C1 by subset
    # ------------------------------------------------------------------
    print("\n=== C1GALT1 / C1GALT1C1 mean expression by B-cell subset ===")
    subset_order = ["naive", "memory", "plasma", "other_B"]
    for g in GENE_OF_INTEREST:
        e = expr(g)
        row = []
        for s in subset_order:
            m = e[b_cells & (subset == s)]
            row.append(f"{s}={m.mean():.3f}" if len(m) else f"{s}=NA")
        print(f"  {g}: " + " | ".join(row))

    # ------------------------------------------------------------------
    # 6. Statistical tests
    # ------------------------------------------------------------------
    results = []
    for g in GENE_OF_INTEREST:
        e = expr(g)
        # within B cells: naive vs plasma (trajectory gradient test)
        naive_e = e[b_cells & (subset == "naive")]
        plasma_e = e[b_cells & (subset == "plasma")]
        memory_e = e[b_cells & (subset == "memory")]
        if len(naive_e) > 5 and len(plasma_e) > 5:
            u, p = stats.mannwhitneyu(naive_e, plasma_e, alternative="two-sided")
            results.append((g, "naive_vs_plasma", u, p, naive_e.mean(), plasma_e.mean()))
        if len(naive_e) > 5 and len(memory_e) > 5:
            u, p = stats.mannwhitneyu(naive_e, memory_e, alternative="two-sided")
            results.append((g, "naive_vs_memory", u, p, naive_e.mean(), memory_e.mean()))

        # group comparison within B cells (Late vs Healthy)
        late_e = e[b_cells & np.isin(sample_per_cell, ["L1","L2","L3","L4","L5","L6"])]
        hc_e = e[b_cells & np.isin(sample_per_cell, ["H1","H2","H3","H4","H5","H6","H7","H8","H9"])]
        early_e = e[b_cells & np.isin(sample_per_cell, ["E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11"])]
        for lbl, a, b in [("late_vs_hc", late_e, hc_e), ("late_vs_early", late_e, early_e), ("early_vs_hc", early_e, hc_e)]:
            if len(a) > 5 and len(b) > 5:
                u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
                results.append((g, lbl, u, p, a.mean(), b.mean()))

    print("\n=== Statistical tests ===")
    for r in results:
        g, lbl, u, p, ma, mb = r
        print(f"  {g} {lbl}: mean {ma:.3f} vs {mb:.3f}, U={u:.1f}, P={p:.3e}")

    # ------------------------------------------------------------------
    # 7. Save results
    # ------------------------------------------------------------------
    import csv
    with open(OUT_DIR("bcell_c1galt1_summary.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["gene", "comparison", "U", "P", "mean_A", "mean_B"])
        for r in results:
            w.writerow(r)

    # subset counts
    with open(OUT_DIR("bcell_subset_counts.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["subset", "n_cells"])
        for s in subset_order:
            w.writerow([s, int((b_cells & (subset == s)).sum())])

    print(f"\nResults saved to {OUT_DIR}")


if __name__ == "__main__":
    main()
