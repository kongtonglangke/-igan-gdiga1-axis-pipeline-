#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Directional validation — Step 4: donor-level pseudobulk analysis of the B-cell
maturation trajectory in GSE285335.

Motivation (reviewer point): C1GALT1C1 is detected in only ~3-5% of B cells by the
10x 5' assay, so per-cell tests (Fisher exact on detection; Mann-Whitney on
per-cell level) treat individual cells as independent replicates and are sensitive
to drop-out. A donor-level pseudobulk analysis aggregates cells within each
(donor x B-cell state) unit and tests the naive -> plasma shift with the donor as
the unit of replication (paired Wilcoxon signed-rank), which is the standard,
conservative read-out for sparse single-cell data.

Outputs (results/directional_validation/):
  bcell_pseudobulk_by_donor.tsv   per donor x state: n_cells, detection rate, mean log1p CPM
  bcell_pseudobulk_summary.json   paired test results + group means

Implementation note: per-sample processing (no full genes-by-cells hstack) keeps
peak memory low; only the ~25 marker genes are densified.

Gating is identical to step3_robust_analysis.py (B: >=2 of 4 canonical markers,
log1p CPM > 0; plasma: sum of plasma-marker log1p CPM >= log1p(3)).
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
import scipy.sparse as sp
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

    gene_vals = {g: [] for g in NEEDED}   # per-gene list of per-cell arrays
    sample_of_cell = []

    for stem in stems:
        feat, barc, M = load_triplet(stem)
        gi = {g: i for i, g in enumerate(feat)}
        libsize = np.asarray(M.sum(axis=0)).ravel()
        keep = libsize > 500
        libsize = libsize[keep]
        sample_of_cell.append(np.array([stem.split("_")[-1]] * int(keep.sum()), dtype=object))
        scale = 1e4 / libsize
        for g in NEEDED:
            if g in gi:
                # rows are 1 x ncells, convert the single row to a dense 1-D vector
                row = np.asarray(M[gi[g], :].todense()).ravel()[keep]
                vals = np.log1p(row * scale)
            else:
                vals = np.zeros(int(keep.sum()))
            gene_vals[g].append(vals)
        print(f"  loaded {stem}: {int(keep.sum())} cells")

    sample_of_cell = np.concatenate(sample_of_cell)
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

    print(f"B cells: {int(is_b.sum())} / {n_cells}")

    rows = []
    for donor in sorted(set(sample_of_cell)):
        for state in ("naive", "memory", "plasma"):
            m = is_b & (subset == state) & (sample_of_cell == donor)
            n = int(m.sum())
            if n == 0:
                continue
            rec = {"donor": donor, "state": state, "n_cells": n}
            for gene in GENE_OF_INTEREST:
                e = expr(gene)[m]
                rec[f"{gene}_det"] = float((e > 0).mean())
                rec[f"{gene}_mean"] = float(e.mean())
            rows.append(rec)
    rows.sort(key=lambda r: (r["donor"], r["state"]))

    with open(OUT_DIR("bcell_pseudobulk_by_donor.tsv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)

    by = {(r["donor"], r["state"]): r for r in rows}
    donors_both = [d for d in sorted({r["donor"] for r in rows})
                   if (d, "naive") in by and (d, "plasma") in by]

    summary = {"n_donors_total": len(stems), "n_cells": int(n_cells),
               "n_bcells": int(is_b.sum()), "n_donors_naive_and_plasma": len(donors_both)}
    for gene in GENE_OF_INTEREST:
        for metric in ("det", "mean"):
            a = np.array([by[(d, "naive")][f"{gene}_{metric}"] for d in donors_both])
            b = np.array([by[(d, "plasma")][f"{gene}_{metric}"] for d in donors_both])
            try:
                wstat, p = stats.wilcoxon(a, b, alternative="two-sided")
            except Exception:
                wstat, p = float("nan"), float("nan")
            summary[f"{gene}_{metric}_naive_mean"] = float(np.mean(a))
            summary[f"{gene}_{metric}_plasma_mean"] = float(np.mean(b))
            summary[f"{gene}_{metric}_naive_median"] = float(np.median(a))
            summary[f"{gene}_{metric}_plasma_median"] = float(np.median(b))
            summary[f"{gene}_{metric}_n_up"] = int((b > a).sum())
            summary[f"{gene}_{metric}_n_down"] = int((b < a).sum())
            summary[f"{gene}_{metric}_wilcoxon_W"] = float(wstat)
            summary[f"{gene}_{metric}_wilcoxon_P"] = float(p)

    for gene in GENE_OF_INTEREST:
        for state in ("naive", "memory", "plasma"):
            m = is_b & (subset == state)
            summary[f"pooled_{gene}_{state}_det"] = float((expr(gene)[m] > 0).mean())
            summary[f"pooled_{gene}_{state}_n"] = int(m.sum())

    with open(OUT_DIR("bcell_pseudobulk_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
