#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Directional validation — Step 5: doublet-contamination sensitivity analysis of the
B-lineage maturation result in GSE285335.

Motivation (reviewer point): in the marker-gated B-lineage partition a large share
of cells fell into the plasma/plasmablast-like class, and plasma cells are rare in
PBMC. High-RNA cells (many detected genes / many counts) are the classic doublet
proxy, and a doublet enriched in the plasma-like gate could inflate the apparent
naive -> plasma detection-rate shift of C1GALT1/C1GALT1C1 (a doublet carries two
transcriptomes, so more genes are detected per droplet).

This script does NOT add a new result to the paper; it is a robustness check:
  1. per-cell library size (counts) and detected-gene number (nFeature) for every cell;
  2. the naive vs plasma-like nFeature / library-size distributions within the
     marker-gated B-lineage pool (is the plasma-like gate enriched in high-RNA cells?);
  3. the naive -> plasma-like C1GALT1 / C1GALT1C1 detection-rate shift recomputed
     (a) on all B-lineage cells, and (b) after removing a transparent high-RNA
     doublet proxy (cells above the per-sample 99th percentile of nFeature), both
     pooled (Fisher exact) and donor-level (paired Wilcoxon; donor x state pseudobulk).

Gating is identical to step4_pseudobulk.py (B: >=2 of 4 canonical markers at
log1p CPM > 0; plasma-like: sum of plasma-marker log1p CPM >= log1p(3)).

Outputs (results/directional_validation/):
  bcell_doublet_sensitivity.json
  bcell_doublet_sensitivity.tsv
"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402
# --------------------------------------------------------------------------
import os
import gzip
import glob
import json
import csv
import numpy as np
import scipy.io as sio
import scipy.stats as stats

RAW_DIR = _paths.at("00_rawdata/GSE285335/extracted")
OUT_DIR = _paths.at("阶段4.5_复现包/reproducibility/results/directional_validation")
os.makedirs(OUT_DIR, exist_ok=True)

GENE_OF_INTEREST = ["C1GALT1", "C1GALT1C1"]
B_MARKERS = ["MS4A1", "CD19", "CD79A", "CD79B"]
NAIVE_MARKERS = ["IGHD", "IGHM"]
MEMORY_MARKER = "CD27"
PLASMA_MARKERS = ["XBP1", "PRDM1", "JCHAIN", "MZB1"]
NEEDED = sorted(set(GENE_OF_INTEREST + B_MARKERS + NAIVE_MARKERS + [MEMORY_MARKER] + PLASMA_MARKERS))

NFeature_PCT = 99.0   # per-sample high-nFeature doublet proxy (top 1%)


def load_triplet(stem):
    barc = [l.rstrip("\n") for l in gzip.open(RAW_DIR(f"{stem}_barcodes.tsv.gz"), "rt")]
    gsm = stem.split("_")[0]
    feat = [l.rstrip("\n").split("\t")[0] for l in gzip.open(RAW_DIR(f"{gsm}_features.tsv.gz"), "rt")]
    M = sio.mmread(RAW_DIR(f"{stem}_matrix.mtx.gz")).tocsr()
    return feat, barc, M


def main():
    stems = sorted(set(
        os.path.basename(f).split("_barcodes.tsv.gz")[0]
        for f in glob.glob(RAW_DIR("*_barcodes.tsv.gz"))
    ))

    gene_vals = {g: [] for g in NEEDED}
    sample_of_cell, nfeat_of_cell, libsize_of_cell, doublet_proxy = [], [], [], []

    for stem in stems:
        feat, barc, M = load_triplet(stem)
        gi = {g: i for i, g in enumerate(feat)}
        libsize = np.asarray(M.sum(axis=0)).ravel()
        keep = libsize > 500
        # detected genes per cell = nonzeros per column
        nfeat = np.diff(M.tocsc().indptr)  # len == n_cells
        libsize_k = libsize[keep]
        nfeat_k = nfeat[keep]
        thr = np.percentile(nfeat_k, NFeature_PCT)
        doublet_proxy.append((nfeat_k > thr))
        sample_of_cell.append(np.array([stem.split("_")[-1]] * int(keep.sum()), dtype=object))
        nfeat_of_cell.append(nfeat_k)
        libsize_of_cell.append(libsize_k)
        scale = 1e4 / libsize_k
        for g in NEEDED:
            if g in gi:
                row = np.asarray(M[gi[g], :].todense()).ravel()[keep]
                vals = np.log1p(row * scale)
            else:
                vals = np.zeros(int(keep.sum()))
            gene_vals[g].append(vals)
        print(f"  loaded {stem}: {int(keep.sum())} cells (nFeature 99th pct = {thr:.0f})")

    sample_of_cell = np.concatenate(sample_of_cell)
    nfeat_of_cell = np.concatenate(nfeat_of_cell)
    libsize_of_cell = np.concatenate(libsize_of_cell)
    doublet_proxy = np.concatenate(doublet_proxy)
    expr_map = {g: np.concatenate(gene_vals[g]) for g in NEEDED}
    n_cells = sample_of_cell.shape[0]

    def expr(gene):
        return expr_map.get(gene, np.zeros(n_cells))

    b_pos = sum((expr(g) > 0).astype(int) for g in B_MARKERS)
    is_b = b_pos >= 2

    igd = expr("IGHD"); ighm = expr("IGHM"); cd27 = expr("CD27")
    plasma_score = sum(expr(g) for g in PLASMA_MARKERS)
    naive_mask = (igd > 0) | (ighm > 0)
    plasma_mask = plasma_score >= np.log1p(3)

    subset = np.full(n_cells, "other", dtype=object)
    subset[is_b & plasma_mask] = "plasma"
    subset[is_b & (~plasma_mask) & naive_mask] = "naive"
    subset[is_b & (~plasma_mask) & (~naive_mask) & (cd27 > 0)] = "memory"
    subset[is_b & (~plasma_mask) & (~naive_mask) & (cd27 == 0)] = "naive"

    summary = {"n_cells_total": int(n_cells), "n_blineage": int(is_b.sum()),
               "doublet_proxy_rule": f"nFeature > per-sample {NFeature_PCT}th percentile",
               "n_doublet_proxy_all": int(doublet_proxy.sum()),
               "n_doublet_proxy_in_blineage": int((doublet_proxy & is_b).sum())}

    # --- high-RNA enrichment in the plasma-like gate ---
    for state in ("naive", "memory", "plasma"):
        m = is_b & (subset == state)
        summary[f"{state}_n"] = int(m.sum())
        summary[f"{state}_nfeat_median"] = float(np.median(nfeat_of_cell[m])) if m.sum() else float("nan")
        summary[f"{state}_libsize_median"] = float(np.median(libsize_of_cell[m])) if m.sum() else float("nan")
        summary[f"{state}_doublet_proxy_frac"] = float(doublet_proxy[m].mean()) if m.sum() else float("nan")

    # Fisher: plasma-like vs naive doublet-proxy enrichment
    tab = [[int((doublet_proxy & is_b & (subset == "plasma")).sum()),
            int((~doublet_proxy & is_b & (subset == "plasma")).sum())],
           [int((doublet_proxy & is_b & (subset == "naive")).sum()),
            int((~doublet_proxy & is_b & (subset == "naive")).sum())]]
    odds, p_enr = stats.fisher_exact(tab)
    summary["plasma_vs_naive_doublet_oddsratio"] = float(odds)
    summary["plasma_vs_naive_doublet_fisher_P"] = float(p_enr)

    rows = []

    def detection_report(tag, mask_cells, note):
        """Pooled (Fisher) naive vs plasma detection shift on the selected cells."""
        rec = {"analysis": tag, "note": note, "n_cells": int(mask_cells.sum())}
        for gene in GENE_OF_INTEREST:
            e = expr(gene)
            a = e[mask_cells & is_b & (subset == "naive")] > 0
            b = e[mask_cells & is_b & (subset == "plasma")] > 0
            rec[f"{gene}_naive_det"] = float(a.mean()) if a.size else float("nan")
            rec[f"{gene}_plasma_det"] = float(b.mean()) if b.size else float("nan")
            table = [[int(b.sum()), int((~b).sum())], [int(a.sum()), int((~a).sum())]]
            _, p = stats.fisher_exact(table)
            rec[f"{gene}_fisher_P"] = float(p)
            # donor-level paired Wilcoxon on donor x state units
            by = {}
            for d in sorted(set(sample_of_cell)):
                for s in ("naive", "plasma"):
                    m = mask_cells & is_b & (subset == s) & (sample_of_cell == d)
                    if m.sum() > 0:
                        by[(d, s)] = float((e[m] > 0).mean())
            donors = [d for d in sorted({k[0] for k in by}) if (d, "naive") in by and (d, "plasma") in by]
            av = np.array([by[(d, "naive")] for d in donors])
            bv = np.array([by[(d, "plasma")] for d in donors])
            try:
                _, pw = stats.wilcoxon(av, bv, alternative="two-sided")
            except Exception:
                pw = float("nan")
            rec[f"{gene}_donor_naive_mean"] = float(av.mean()) if av.size else float("nan")
            rec[f"{gene}_donor_plasma_mean"] = float(bv.mean()) if bv.size else float("nan")
            rec[f"{gene}_donor_wilcoxon_P"] = float(pw)
            rec[f"{gene}_n_donors"] = int(len(donors))
        rows.append(rec)

    all_cells = np.ones(n_cells, dtype=bool)
    detection_report("all_blineage", all_cells, "marker-gated B-lineage, no doublet removal")
    detection_report("doublet_removed", ~doublet_proxy, "high-nFeature (top 1% per sample) cells removed")

    with open(OUT_DIR("bcell_doublet_sensitivity.json"), "w", encoding="utf-8") as f:
        json.dump({"summary": summary, "reports": rows}, f, indent=2)
    with open(OUT_DIR("bcell_doublet_sensitivity.tsv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)

    print(json.dumps(summary, indent=2))
    for r in rows:
        print(r)


if __name__ == "__main__":
    main()
