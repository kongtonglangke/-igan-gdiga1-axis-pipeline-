#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Directional validation — Step 3: robust B-cell trajectory analysis + figures.

Improvements over step2:
  - Plasma markers restricted to those with real signal (XBP1/PRDM1/JCHAIN/MZB1),
    dropping near-zero-detection SDC1/AICDA.
  - Two complementary read-outs, both reported honestly:
      (a) detection rate (fraction of cells expressing the gene) — Fisher exact test
      (b) expression level among expressing cells — Mann-Whitney U
    This avoids the drop-out-driven artefact where a low-detection gene's "mean"
    is dominated by the zero fraction.
  - Outputs a publication-ready figure (trajectory violin/bar) + a summary TSV.
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
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RAW_DIR = _paths.at("00_rawdata/GSE285335/extracted")
OUT_DIR = _paths.at("阶段4.5_复现包/reproducibility/results/directional_validation")
os.makedirs(OUT_DIR, exist_ok=True)

# Sample -> group
LATE = [f"L{i}" for i in range(1, 7)]
EARLY = [f"E{i}" for i in range(1, 12)]
HEALTHY = [f"H{i}" for i in range(1, 10)]

GENE_OF_INTEREST = ["C1GALT1", "C1GALT1C1"]
B_MARKERS = ["MS4A1", "CD19", "CD79A", "CD79B"]
NAIVE_MARKERS = ["IGHD", "IGHM"]
MEMORY_MARKER = "CD27"
PLASMA_MARKERS = ["XBP1", "PRDM1", "JCHAIN", "MZB1"]  # SDC1/AICDA dropped (near-zero detection)


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

    all_genes = None
    blocks = []
    sample_label = []
    ncell_per_sample = []
    for stem in stems:
        feat, barc, M = load_triplet(stem)
        if all_genes is None:
            all_genes = feat
        else:
            assert all_genes == feat
        blocks.append(M)
        sample_label.append(stem.split("_")[-1])
        ncell_per_sample.append(M.shape[1])

    X = sp.hstack(blocks, format="csr")
    n_cells = X.shape[1]
    sample_per_cell = np.repeat(sample_label, ncell_per_sample)

    # libsize + QC
    libsize = np.asarray(X.sum(axis=0)).ravel()
    keep_cell = libsize > 500
    libsize = libsize[keep_cell]
    sample_per_cell = sample_per_cell[keep_cell]
    n_cells = int(keep_cell.sum())

    # extract only marker genes -> dense
    all_markers = sorted(set(GENE_OF_INTEREST + B_MARKERS + NAIVE_MARKERS + [MEMORY_MARKER] + PLASMA_MARKERS))
    gene_idx = {g: i for i, g in enumerate(all_genes)}
    row_idx = np.array([gene_idx[g] for g in all_markers if g in gene_idx])
    row_names = [all_genes[i] for i in row_idx]

    Xsub = X[row_idx, :][:, keep_cell].tocsr().astype(np.float64)
    Xsub = Xsub.multiply((1e4 / libsize)[None, :]).tocsr()
    Xsub.data = np.log1p(Xsub.data)
    G = np.asarray(Xsub.todense())
    expr_map = {g: G[i] for i, g in enumerate(row_names)}

    def expr(gene):
        v = expr_map.get(gene)
        return v if v is not None else np.zeros(n_cells)

    # ---------- B cell identification ----------
    b_pos = sum((expr(g) > 0).astype(int) for g in B_MARKERS)
    is_b = b_pos >= 2
    n_b = int(is_b.sum())

    # ---------- subset assignment ----------
    igd = expr("IGHD"); ighm = expr("IGHM"); cd27 = expr("CD27")
    plasma_score = sum(expr(g) for g in PLASMA_MARKERS)
    naive_mask = (igd > 0) | (ighm > 0)
    plasma_mask = plasma_score >= np.log1p(3)   # >=3 counts-equivalent among plasma markers

    subset = np.full(n_cells, "other", dtype=object)
    subset[is_b & plasma_mask] = "plasma"
    subset[is_b & (~plasma_mask) & naive_mask] = "naive"
    subset[is_b & (~plasma_mask) & (~naive_mask) & (cd27 > 0)] = "memory"
    subset[is_b & (~plasma_mask) & (~naive_mask) & (cd27 == 0)] = "naive"

    # ---------- group labels ----------
    grp = np.full(n_cells, "", dtype=object)
    for g in LATE: grp[sample_per_cell == g] = "Late"
    for g in EARLY: grp[sample_per_cell == g] = "Early"
    for g in HEALTHY: grp[sample_per_cell == g] = "HC"

    # ---------- output summary ----------
    subset_order = ["naive", "memory", "plasma"]
    print(f"B cells: {n_b} / {n_cells} ({100*n_b/n_cells:.2f}%)")
    for s in subset_order:
        n = int((is_b & (subset == s)).sum())
        print(f"  {s}: {n}")

    # ---------- two read-outs per gene ----------
    rows = []
    for gene in GENE_OF_INTEREST:
        e = expr(gene)
        for s in subset_order:
            mask = is_b & (subset == s)
            det = float((e[mask] > 0).mean())          # detection rate
            mean_expr = float(e[mask].mean())           # log1p CPM mean (incl zeros)
            mean_nonzero = float(e[mask][e[mask] > 0].mean()) if (e[mask] > 0).any() else 0.0
            rows.append((gene, s, int(mask.sum()), det, mean_expr, mean_nonzero))

    # trajectory tests (naive vs plasma) — detection (Fisher) + level (MWU)
    print("\n=== Trajectory: naive -> plasma ===")
    tests = []
    for gene in GENE_OF_INTEREST:
        e = expr(gene)
        naive_m = is_b & (subset == "naive")
        plas_m = is_b & (subset == "plasma")
        mem_m = is_b & (subset == "memory")
        # detection rate Fisher
        a = int((e[naive_m] > 0).sum()); b = int(naive_m.sum()) - a
        c = int((e[plas_m] > 0).sum()); d = int(plas_m.sum()) - c
        or_, p_fisher = stats.fisher_exact([[a, b], [c, d]])
        # level MWU among expressing cells
        nz_n = e[naive_m][e[naive_m] > 0]; nz_p = e[plas_m][e[plas_m] > 0]
        if len(nz_n) > 5 and len(nz_p) > 5:
            u, p_mwu = stats.mannwhitneyu(nz_n, nz_p, alternative="two-sided")
        else:
            u, p_mwu = np.nan, np.nan
        det_n = (e[naive_m] > 0).mean(); det_p = (e[plas_m] > 0).mean()
        tests.append((gene, "naive", "plasma", det_n, det_p, or_, p_fisher,
                      nz_n.mean() if len(nz_n) else np.nan, nz_p.mean() if len(nz_p) else np.nan, p_mwu))
        print(f"  {gene}: det naive={det_n:.3f} vs plasma={det_p:.3f} (OR={or_:.2f}, P_fisher={p_fisher:.2e}); "
              f"level_nonzero {nz_n.mean():.3f} vs {nz_p.mean():.3f} (P_mwu={p_mwu:.2e})")

    # group comparison (Late / Early / HC) on the same two read-outs — reported in
    # the reproducibility package only (the dataset is not powered for clinical
    # staging; see Methods "Single-cell B-cell reanalysis")
    print("\n=== Group comparison (B cells only) ===")
    g_rows = []
    for gene in GENE_OF_INTEREST:
        e = expr(gene)
        for ga, gb in [("Late", "Early"), ("Late", "HC"), ("Early", "HC")]:
            ma = is_b & (grp == ga); mb = is_b & (grp == gb)
            a = int((e[ma] > 0).sum()); b = int(ma.sum()) - a
            c = int((e[mb] > 0).sum()); d = int(mb.sum()) - c
            or_, p_fisher = stats.fisher_exact([[a, b], [c, d]])
            nz_a = e[ma][e[ma] > 0]; nz_b = e[mb][e[mb] > 0]
            if len(nz_a) > 5 and len(nz_b) > 5:
                _, p_mwu = stats.mannwhitneyu(nz_a, nz_b, alternative="two-sided")
            else:
                p_mwu = np.nan
            g_rows.append((gene, ga, gb, int(ma.sum()), int(mb.sum()),
                           float((e[ma] > 0).mean()), float((e[mb] > 0).mean()),
                           p_fisher, p_mwu))
            print(f"  {gene} {ga} vs {gb}: det {float((e[ma]>0).mean()):.3f} vs "
                  f"{float((e[mb]>0).mean()):.3f} (P_fisher={p_fisher:.3g}); P_mwu={p_mwu:.3g}")

    # save TSVs
    import csv
    with open(OUT_DIR("bcell_expression_by_subset.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["gene", "subset", "n_cells", "detection_rate", "mean_log1pCPM", "mean_nonzero_log1pCPM"])
        w.writerows(rows)
    with open(OUT_DIR("trajectory_tests.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["gene", "subset_A", "subset_B", "det_A", "det_B", "OR", "P_fisher",
                    "level_A_nonzero", "level_B_nonzero", "P_mwu"])
        w.writerows(tests)
    with open(OUT_DIR("bcell_group_comparison.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["gene", "group_A", "group_B", "n_A", "n_B", "det_A", "det_B",
                    "P_fisher", "P_mwu"])
        w.writerows(g_rows)
    # definitive per-subset counts under the final gating rules (supersedes the
    # step2 snapshot; written here so the shipped counts can never drift from the
    # counts used in the figure and in the manuscript text)
    with open(OUT_DIR("bcell_subset_counts.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["subset", "n_cells"])
        for s in subset_order:
            w.writerow([s, int((is_b & (subset == s)).sum())])
        w.writerow(["total_B", n_b])
        w.writerow(["total_cells_qc", n_cells])

    # ---------- figure ----------
    plot_figure(G, row_names, is_b, subset, subset_order, grp, OUT_DIR, n_b=n_b)
    print("\nDone. Results + figure saved to", OUT_DIR)


def plot_figure(G, row_names, is_b, subset, subset_order, grp, outdir, n_b=None):
    """2x2 panel: detection rate + expression level, for C1GALT1 & C1GALT1C1."""
    expr_map = {g: G[i] for i, g in enumerate(row_names)}
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8})
    fig, axes = plt.subplots(2, 2, figsize=(8.5, 6.5))
    colors = {"naive": "#3B82F6", "memory": "#8B5CF6", "plasma": "#EF4444"}

    for gi, gene in enumerate(GENE_OF_INTEREST):
        e = expr_map[gene]
        # left panel: detection rate (bar)
        ax1 = axes[gi, 0]
        dets, errs, labels = [], [], []
        for s in subset_order:
            m = is_b & (subset == s)
            det = (e[m] > 0).mean()
            n = int(m.sum())
            se = np.sqrt(det * (1 - det) / n)
            dets.append(det * 100); errs.append(se * 100); labels.append(s)
        ax1.bar(range(3), dets, yerr=errs, capsize=3, color=[colors[s] for s in subset_order],
                edgecolor="black", lw=0.5)
        ax1.set_xticks(range(3)); ax1.set_xticklabels(["naive", "memory", "plasma"])
        ax1.set_ylabel(f"{gene} detection rate (%)")
        ax1.set_title(f"{gene} — fraction of B cells expressing", fontsize=9, weight="bold")
        ax1.set_ylim(0, max(dets) * 1.30 + 1)
        for i, d in enumerate(dets):
            ax1.text(i, d + errs[i] + max(dets) * 0.03, f"{d:.1f}%", ha="center", fontsize=7)
        ax1.spines[["top", "right"]].set_visible(False)

        # right panel: expression level (among expressing cells) — box + strip
        ax2 = axes[gi, 1]
        data = []; pos = []
        for si, s in enumerate(subset_order):
            m = is_b & (subset == s)
            vals = e[m][e[m] > 0]
            data.append(vals); pos.append(si)
        bp = ax2.boxplot(data, positions=pos, widths=0.5, patch_artist=True,
                         showfliers=False, medianprops=dict(color="black"))
        for patch, s in zip(bp["boxes"], subset_order):
            patch.set_facecolor(colors[s]); patch.set_alpha(0.7)
        ax2.set_xticks(pos); ax2.set_xticklabels(["naive", "memory", "plasma"])
        ax2.set_ylabel(f"{gene} log1p-CPM (expr. cells)")
        ax2.set_title(f"{gene} — level among expressing B cells", fontsize=9, weight="bold")
        ax2.spines[["top", "right"]].set_visible(False)

    n_b_txt = f"{n_b:,}" if n_b is not None else "n/a"
    fig.suptitle(f"C1GALT1 / C1GALT1C1 along the B-cell maturation axis "
                 f"(GSE285335, IgAN PBMC, {n_b_txt} B cells)",
                 fontsize=10, weight="bold", y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(os.path.join(outdir, "Figure_S4_Bcell_C1GALT1_trajectory.png"), dpi=300)
    fig.savefig(os.path.join(outdir, "Figure_S4_Bcell_C1GALT1_trajectory.pdf"))
    plt.close(fig)
    print("Figure saved: Figure_S4_Bcell_C1GALT1_trajectory.png/pdf")


if __name__ == "__main__":
    main()
