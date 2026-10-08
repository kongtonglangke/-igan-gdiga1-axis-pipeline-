#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Directional validation — Step 6: audit of the marker-gated plasma-like class and
definition-robustness of the naive -> plasma detection-rate shift in GSE285335.

Motivation (reviewer point 3): in the marker-gated B-lineage partition the
plasma/plasmablast-like class holds ~39% of the B-lineage pool, which is far above
the expected plasma-cell fraction of PBMC. Two questions follow:
  (Q-A) is the broad "plasma-like" gate capturing real plasma cells, or is it a
        low-threshold marker rule that also admits activated B cells?
  (Q-B) does the naive -> plasma detection-rate shift of C1GALT1 / C1GALT1C1 hold
        when the B-lineage pool is split by a *different* definition?

For (Q-A) this script reports, within the marker-gated plasma-like class, the
positivity rate for each canonical B-cell marker (MS4A1, CD19, CD79A, CD79B) and for
each plasma marker individually, plus the distribution of the combined plasma score.

For (Q-B) it re-uses the unsupervised, cluster-based partition of the same dataset
(reproducibility/results/sc_regulatory/cell_type_assignment.tsv; produced by
code/sc_regulatory/step1_atlas.py on the identical 282,463 cells) and recomputes the
axis detection rates for two of its resolutions:
  * compartment level  : cell_type "B" versus cell_type "Plasma";
  * sub-partition level: b_subset "B_naive" versus "B_plasma".
Both are reported with a pooled Fisher exact test and a donor-level paired
Wilcoxon signed-rank test (unit = donor x state pseudobulk).

This script adds no new claim to the paper by itself; it quantifies a robustness
property of the single-cell cross-check.

Outputs (results/directional_validation/):
  bcell_gate_audit.json
  bcell_gate_audit.tsv
  bcell_definition_crosscheck.tsv
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
ATLAS = _paths.at("阶段4.5_复现包/reproducibility/results/sc_regulatory/cell_type_assignment.tsv")
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


def wilcoxon_paired(a, b):
    try:
        _, p = stats.wilcoxon(a, b, alternative="two-sided")
        return float(p)
    except Exception:
        return float("nan")


def donor_detection(sample_of_cell, state_of_cell, expr_pos, want_a, want_b):
    """donor-level paired means of the detection rate for states want_a / want_b."""
    by = {}
    for d in sorted(set(sample_of_cell)):
        for s in (want_a, want_b):
            m = (sample_of_cell == d) & (state_of_cell == s)
            if m.sum() > 0:
                by[(d, s)] = float(expr_pos[m].mean())
    donors = [d for d in sorted({k[0] for k in by}) if (d, want_a) in by and (d, want_b) in by]
    av = np.array([by[(d, want_a)] for d in donors])
    bv = np.array([by[(d, want_b)] for d in donors])
    return av, bv, wilcoxon_paired(av, bv), donors


def main():
    stems = sorted(set(
        os.path.basename(f).split("_barcodes.tsv.gz")[0]
        for f in glob.glob(RAW_DIR("*_barcodes.tsv.gz"))
    ))

    gene_vals = {g: [] for g in NEEDED}
    sample_of_cell, barcode_of_cell = [], []

    for stem in stems:
        feat, barc, M = load_triplet(stem)
        gi = {g: i for i, g in enumerate(feat)}
        libsize = np.asarray(M.sum(axis=0)).ravel()
        keep = libsize > 500
        scale = 1e4 / libsize[keep]
        lbl = stem.split("_")[-1]
        sample_of_cell.append(np.array([lbl] * int(keep.sum()), dtype=object))
        barcode_of_cell.append(np.array(barc, dtype=object)[keep])
        for g in NEEDED:
            if g in gi:
                row = np.asarray(M[gi[g], :].todense()).ravel()[keep]
                vals = np.log1p(row * scale)
            else:
                vals = np.zeros(int(keep.sum()))
            gene_vals[g].append(vals)
        print(f"  loaded {stem}: {int(keep.sum())} cells")

    sample_of_cell = np.concatenate(sample_of_cell)
    barcode_of_cell = np.concatenate(barcode_of_cell)
    expr_map = {g: np.concatenate(gene_vals[g]) for g in NEEDED}
    n_cells = sample_of_cell.shape[0]

    def expr(gene):
        return expr_map.get(gene, np.zeros(n_cells))

    # ---------------- marker gate (identical to step4/step5) ----------------
    b_pos = sum((expr(g) > 0).astype(int) for g in B_MARKERS)
    is_b = b_pos >= 2
    igd = expr("IGHD"); ighm = expr("IGHM"); cd27 = expr("CD27")
    plasma_score = sum(expr(g) for g in PLASMA_MARKERS)
    naive_mask = (igd > 0) | (ighm > 0)
    plasma_mask = plasma_score >= np.log1p(3)

    state = np.full(n_cells, "other", dtype=object)
    state[is_b & plasma_mask] = "plasma"
    state[is_b & (~plasma_mask) & naive_mask] = "naive"
    state[is_b & (~plasma_mask) & (~naive_mask) & (cd27 > 0)] = "memory"
    state[is_b & (~plasma_mask) & (~naive_mask) & (cd27 == 0)] = "naive"

    pg = is_b & (state == "plasma")
    ng = is_b & (state == "naive")

    summary = {
        "n_cells": int(n_cells),
        "n_blineage": int(is_b.sum()),
        "n_naive": int(ng.sum()),
        "n_memory": int((is_b & (state == "memory")).sum()),
        "n_plasma_like": int(pg.sum()),
        "plasma_like_frac_of_blineage": float(pg.sum() / is_b.sum()),
        "gate_rule_plasma": "sum of log1p CPM of XBP1/PRDM1/JCHAIN/MZB1 >= log1p(3)",
        "gate_rule_plasma_threshold_natural": float(np.log1p(3)),
        "gate_priority": "B-lineage gate first; plasma-like overrides naive/memory",
    }

    # ---------------- Q-A: purity of the plasma-like gate ----------------
    for name, mask in (("plasma_like", pg), ("naive", ng)):
        for mk in B_MARKERS:
            summary[f"{name}_frac_{mk}_pos"] = float((expr(mk)[mask] > 0).mean())
        for mk in PLASMA_MARKERS:
            summary[f"{name}_frac_{mk}_pos"] = float((expr(mk)[mask] > 0).mean())
        summary[f"{name}_frac_all4_plasma_pos"] = float(
            np.logical_and.reduce([expr(mk)[mask] > 0 for mk in PLASMA_MARKERS]).mean())
        summary[f"{name}_frac_any_plasma_pos"] = float(
            np.logical_or.reduce([expr(mk)[mask] > 0 for mk in PLASMA_MARKERS]).mean())
        summary[f"{name}_plasma_score_median"] = float(np.median(plasma_score[mask]))
        summary[f"{name}_plasma_score_p90"] = float(np.percentile(plasma_score[mask], 90))
        # MS4A1+ fraction: CD79A/B are near-ubiquitous on B-lineage, MS4A1 is lost on plasma cells
        summary[f"{name}_frac_MS4A1_or_CD19_pos"] = float(
            np.logical_or(expr("MS4A1")[mask] > 0, expr("CD19")[mask] > 0).mean())

    # ---------------- Q-B: unsupervised atlas as a second definition ----------------
    atlas_ct, atlas_bs = {}, {}
    with open(ATLAS, encoding="utf-8") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for row in rd:
            key = row["barcode"]          # already carries the sample-label prefix, e.g. "L1_AAAC..."
            atlas_ct[key] = row["cell_type"]
            atlas_bs[key] = row["b_subset"]
    ct = np.array([atlas_ct.get(b, "UNMAPPED") for b in barcode_of_cell], dtype=object)
    bs = np.array([atlas_bs.get(b, "UNMAPPED") for b in barcode_of_cell], dtype=object)
    summary["n_mapped_to_atlas"] = int((ct != "UNMAPPED").sum())

    xs_rows = []
    def crosscheck(tag, state_arr, a_level, b_level, note):
        rec = {"definition": tag, "level_a": a_level, "level_b": b_level, "note": note}
        ma = state_arr == a_level
        mb = state_arr == b_level
        rec["n_a"] = int(ma.sum()); rec["n_b"] = int(mb.sum())
        for gene in GENE_OF_INTEREST:
            pa = expr(gene)[ma] > 0
            pb = expr(gene)[mb] > 0
            rec[f"{gene}_det_{a_level}"] = float(pa.mean())
            rec[f"{gene}_det_{b_level}"] = float(pb.mean())
            table = [[int(pb.sum()), int((~pb).sum())], [int(pa.sum()), int((~pa).sum())]]
            _, p = stats.fisher_exact(table)
            rec[f"{gene}_fisher_P"] = float(p)
            av, bv, pw, _ = donor_detection(sample_of_cell, state_arr, expr(gene) > 0, a_level, b_level)
            rec[f"{gene}_donor_{a_level}_mean"] = float(av.mean()) if av.size else float("nan")
            rec[f"{gene}_donor_{b_level}_mean"] = float(bv.mean()) if bv.size else float("nan")
            rec[f"{gene}_donor_wilcoxon_P"] = pw
            rec[f"{gene}_n_donors"] = int(av.size)
        xs_rows.append(rec)

    crosscheck("unsupervised atlas — compartment", ct, "B", "Plasma",
               "cell_type labels of the 15-community atlas partition")
    crosscheck("unsupervised atlas — sub-partition", bs, "B_naive", "B_plasma",
               "argmax marker score within the B/Plasma compartment")

    # marker gate, same two genes, for reference in the same table shape
    ref = {"definition": "marker gate (primary)", "level_a": "naive", "level_b": "plasma",
           "note": "sum-of-markers gate used in the paper",
           "n_a": int(ng.sum()), "n_b": int(pg.sum())}
    for gene in GENE_OF_INTEREST:
        pa = expr(gene)[ng] > 0
        pb = expr(gene)[pg] > 0
        ref[f"{gene}_det_naive"] = float(pa.mean())
        ref[f"{gene}_det_plasma"] = float(pb.mean())
        table = [[int(pb.sum()), int((~pb).sum())], [int(pa.sum()), int((~pa).sum())]]
        ref[f"{gene}_fisher_P"] = float(stats.fisher_exact(table)[1])
        av, bv, pw, _ = donor_detection(sample_of_cell, state, expr(gene) > 0, "naive", "plasma")
        ref[f"{gene}_donor_naive_mean"] = float(av.mean()) if av.size else float("nan")
        ref[f"{gene}_donor_plasma_mean"] = float(bv.mean()) if bv.size else float("nan")
        ref[f"{gene}_donor_wilcoxon_P"] = pw
        ref[f"{gene}_n_donors"] = int(av.size)
    xs_rows.append(ref)

    with open(OUT_DIR("bcell_gate_audit.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    with open(OUT_DIR("bcell_gate_audit.tsv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["metric", "value"])
        for k, v in summary.items():
            w.writerow([k, v])
    cols = sorted({k for r in xs_rows for k in r})
    with open(OUT_DIR("bcell_definition_crosscheck.tsv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader(); w.writerows(xs_rows)

    print(json.dumps(summary, indent=2, ensure_ascii=False))
    for r in xs_rows:
        print({k: r[k] for k in sorted(r)})


if __name__ == "__main__":
    main()
