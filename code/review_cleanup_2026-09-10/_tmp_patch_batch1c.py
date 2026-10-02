#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""_tmp_patch_batch1c.py — 定点改写：ch3 ch4 ch5（去内部编号 / Q1-Q4 对接 / P0-P1 各项）"""
import os

MS = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                  "submission", "manuscript")


def patch(fn, pairs):
    p = os.path.join(MS, fn)
    with open(p, encoding="utf-8", newline="") as fh:
        t = fh.read()
    for old, new, exp in pairs:
        n = t.count(old)
        if n != exp:
            raise SystemExit(f"[FAIL] {fn}: found {n}, expected {exp}\n  >> {old[:140]}")
        t = t.replace(old, new)
    with open(p, "w", encoding="utf-8", newline="") as fh:
        fh.write(t)
    print(f"OK  {fn}  ({len(pairs)} pairs)")


# ============================== ch3 Methods ==============================
patch("03_Methods.md", [
    # P0-1 内部层编号
    ("We built a layered causal-inference framework (evidence layers L0\u2013L7; Fig. 1c) to dissect",
     "We built a layered causal-inference framework (Fig. 1c) to dissect", 1),
    ("### Differential expression analysis (evidence layer L0)", "### Differential expression analysis", 1),
    ("### eQTL context mapping (evidence layer L1a)", "### eQTL context mapping", 1),
    ("### Replication of published Gd-IgA1 lead variants (evidence layer L1b)",
     "### Replication of published Gd-IgA1 lead variants", 1),
    ("### Colocalisation and fine-mapping (evidence layer L2)", "### Colocalisation and fine-mapping", 1),
    ("### Allelic scoring and direction consistency (evidence layer L2)",
     "### Allelic scoring and direction consistency", 1),
    ("### Two-step methylation mediation framework (evidence layer L3)",
     "### Two-step methylation mediation framework", 1),
    ("### Direct lookup of axis variants in the IgAN GWAS (evidence layer L4)",
     "### Direct lookup of axis variants in the IgAN GWAS", 1),
    ("### East Asian heterogeneity and allele frequency (evidence layer L5)",
     "### East Asian heterogeneity and allele frequency", 1),
    ("### Druggability and intervention scoring (evidence layer L6)",
     "### Druggability and intervention scoring", 1),
    ("### Single-cell atlas and upstream regulator inference (evidence layer L7)",
     "### Single-cell atlas and upstream regulator inference", 1),
    ("### Methylation-sensitive motif analysis at axis CpGs (evidence layer L7)",
     "### Methylation-sensitive motif analysis at axis CpGs", 1),
    ("strongest in East Asians (R6); this ancestry mismatch",
     "strongest in East Asians (Results); this ancestry mismatch", 1),
    ("This cross-reference informs the East Asian layer (L5/R6) and the ancestry caveats of the Limitations; it introduces no new evidence layer.",
     "This cross-reference informs the East Asian analysis and the ancestry caveats of the Limitations.", 1),
    ("To corroborate the B-cell-state specificity of the cis-eQTL signal (R2) at single-cell resolution",
     "To corroborate the B-cell-state specificity of the cis-eQTL signal at single-cell resolution", 1),
    ("nominal P values for the five-gene hypothesis-driven panel (L0), permutation-based p_perm from the uniformly recomputed eQTL Catalogue (L1), and region-level Bayesian posteriors (L2/L3).",
     "nominal P values for the five-gene hypothesis-driven expression panel, permutation-based p_perm from the uniformly recomputed eQTL Catalogue for the eQTL layers, and region-level Bayesian posteriors for the colocalisation layers.", 1),
    # P0-3 单细胞过滤前后细胞数
    ("Per-sample 10x count matrices (barcodes/features/matrix triplets) were downloaded from the GEO supplementary archive, loaded with scipy.io.mmread and merged by column with scipy.sparse.hstack, yielding a single genes-by-cells sparse matrix. Per-cell library sizes were computed as the column sums; cells with library size \u2264 500 counts were removed, leaving 282,463 cells for analysis.",
     "Per-sample 10x count matrices (CellRanger-filtered barcodes, features and matrix triplets) were downloaded from the GEO supplementary archive, loaded with scipy.io.mmread and merged by column with scipy.sparse.hstack, yielding a single genes-by-cells sparse matrix. Per-cell library sizes were computed as the column sums and a minimum library size of 500 counts was applied; because the deposited matrices are already CellRanger-filtered, this threshold excluded no further cells, and all 282,463 cells (minimum library size 501 counts) were retained for analysis.", 1),
    # P1-4 decoupleR 方法口径
    ("the per-cell activity of each regulator was computed as the mean z-score of its target genes (univariate linear model, as implemented in decoupleR)",
     "the per-cell activity of each regulator was computed with decoupleR's univariate linear model (ULM), whose per-cell test statistic is the mean z-score of the regulon's target genes", 1),
    # P1-9 英文
    ("the kidney series served as directional separation analyses.",
     "the kidney series served as directional-separation datasets.", 1),
    # P1-2 IVW 首用展开
    ("combined in an inverse-variance fixed-effect meta-analysis",
     "combined in an inverse-variance-weighted (IVW) fixed-effect meta-analysis", 1),
    # P2-3 共定位先验依据
    ("Priors were set to p1 = p2 = 1 \u00d7 10\u207b\u2074 and p12 = 1 \u00d7 10\u207b\u2075 with a shared prior effect variance W = 0.04; effect alleles were aligned between eQTL and GWAS beta signs.",
     "Priors were set to p1 = p2 = 1 \u00d7 10\u207b\u2074 and p12 = 1 \u00d7 10\u207b\u2075, with a shared prior effect variance W = 0.04 (a prior standard deviation of 0.2, within the range conventional for cis-regulatory effect sizes) and aligned effect alleles between eQTL and GWAS beta signs; all priors were fixed before analysis and held identical across the 25 planned comparisons.", 1),
    # P2-1 L4 FDR 敏感性
    ("the complete lookup is tabulated in Additional file 1: Tables S4a and S4b.",
     "the complete lookup is tabulated in Additional file 1: Tables S4a and S4b. Reported without correction by design, the look-up was additionally subjected to a Benjamini\u2013Hochberg sensitivity analysis across the 13 testable look-ups (three evaluable lead variants plus the ten gene\u2013context window look-ups): the lead-variant look-ups remained null (smallest adjusted q = 0.28), and the only look-up reaching q < 0.05 was the GALNT12 regional window minimum (rs187211065, effect-allele frequency 0.0024, P = 9.4 \u00d7 10\u207b\u2076, q = 6.1 \u00d7 10\u207b\u2075) \u2014 a rare-variant window scan that remains approximately three orders of magnitude from genome-wide significance and does not involve any axis lead variant.", 1),
    # P2-8 Table 1 标题体例
    ("## Table 1. Data resources used in this study", "**Table 1.** Data resources used in this study", 1),
    # P0-4 引文年份同步（表 1 数据集行）
    ("| GSE285335 (Kim 2025) |", "| GSE285335 (Kim 2026) |", 1),
    # P1-8 LLM/AI 首用展开
    ("### LLM/AI assistance disclosure",
     "### Large language model (LLM) and artificial intelligence (AI) assistance disclosure", 1),
    # 版本日志
    ("> Reference numbers are final: [12]\u2013[24] follow the whole-manuscript order of first appearance (renumbered 2026-09-10); [25]\u2013[27] extend the sequence in Results/Discussion.",
     "> v1.8 (2026-09-10, reviewer-perspective clean-up): all internal analysis codes (L0\u2013L7, R1\u2013R8) removed from the sub-section titles and body (including the L5 cross-reference and the statistical-software paragraph); the single-cell section now states explicitly that the deposited matrices are CellRanger-filtered so that the 500-count threshold removed no further cells (resolving an apparent before/after cell-count contradiction); the decoupleR description now names the univariate linear model and its mean-z-score statistic consistently; the IVW abbreviation is expanded at first use; the colocalisation priors are justified; and an L4 Benjamini\u2013Hochberg sensitivity analysis is reported.\n> Reference numbers are final: [12]\u2013[24] follow the whole-manuscript order of first appearance (renumbered 2026-09-10); [25]\u2013[27] extend the sequence in Results/Discussion.", 1),
])

# ============================== ch4 Results ==============================
patch("04_Results.md", [
    ("(evidence layer L0; Fig. 1b)", "(Fig. 1b)", 1),
    ("supports a peripheral origin of the aberrantly glycosylated IgA1.",
     "supports a peripheral origin of the aberrantly glycosylated IgA1. Quantifying the compartment contrast: all five axis genes were down-regulated in the PBMC series (5 of 5), whereas only one of the ten gene\u2013series measurements across the two kidney datasets was significantly down-regulated (C1GALT1 in microdissected glomeruli) and four were significantly up-regulated (GALNT2 and ST6GALNAC2 in whole kidney; GALNT2 and C1GALT1C1 in glomeruli), so the two compartments differ in the sign of the axis response and not merely in its magnitude.", 1),
    ("Bulk PBMC data cannot resolve the cell type that carries the regulatory signal.",
     "Bulk PBMC data cannot resolve the cell type that carries the regulatory signal (Q1).", 1),
    ("(evidence layer L1a; Fig. 2a)", "(Fig. 2a)", 1),
    ("(evidence layer L1b; Fig. 2b\u2013c; Additional file 1: Figures S2 and S3; Table 3)",
     "(Fig. 2b\u2013c; Additional file 1: Figures S2 and S3; Table 3)", 1),
    ("(GCST90018866; evidence layer L2; Fig. 3a\u2013c; Table 4)", "(GCST90018866; Fig. 3a\u2013c; Table 4; Q3)", 1),
    ("No comparison supported a shared causal variant: posterior probability of hypothesis 4 (PPH4) reached 0.5 in none of the 14 tests (maximum 0.031), and posterior mass concentrated on PPH1\u2014an eQTL signal without a disease signal\u2014in the tissue eQTLs.",
     "No comparison supported a shared causal variant: the posterior probability of a shared causal variant (PPH4; PPH0\u2013PPH4 denote the five coloc hypotheses, with PPH1 an eQTL signal without a disease signal and PPH3 distinct causal variants) did not reach 0.5 in any of the 14 tests (maximum 0.031), and posterior mass concentrated on PPH1 in the tissue eQTLs.", 1),
    ("(evidence layer L3; Fig. 4a\u2013d)", "(Fig. 4a\u2013d; Q2)", 1),
    ("clumped inverse-variance weighting where at least three independent instruments were available",
     "clumped inverse-variance weighting (IVW) where at least three independent instruments were available", 1),
    ("(evidence layer L4), we examined", "(Q3), we examined", 1),
    ("(evidence layer L5; Fig. 5a\u2013b; Table 6)", "(Fig. 5a\u2013b; Table 6; Q4)", 1),
    ("(evidence layer L6; Fig. 6a\u2013c; Table 7)", "(Fig. 6a\u2013c; Table 7; Q4)", 1),
    ("which the L3 results show to be transcription-coupled but not disease-mediating",
     "which the methylation results reported above show to be transcription-coupled but not disease-mediating", 1),
    ("To identify the regulatory programme that sets axis-gene expression\u2014and to locate the epigenetic window of R7 within a concrete transcriptional circuit\u2014we built an unsupervised single-cell atlas of the IgAN PBMC dataset analysed above (GSE285335; 26 donors, 282,463 cells; evidence layer L7).",
     "To identify the regulatory programme that sets axis-gene expression\u2014and to locate the epigenetic window identified above within a concrete transcriptional circuit (Q2)\u2014we built an unsupervised single-cell atlas of the IgAN PBMC dataset analysed above (GSE285335; 26 donors, 282,463 cells).", 1),
    ("had the lowest detection rate of any lineage in B cells (3.2%).",
     "had the lowest detection rate of any lineage in B cells in this unsupervised atlas partition (3.2%); the marker-gated B-cell partition described above gives a higher estimate (4.2%), because the two gates resolve different B-cell pools (Fig. S4).", 1),
    ("reinforces the cell-type-precise intervention argument of R7.",
     "reinforces the cell-type-precise intervention argument developed above.", 1),
    ("(evidence layer L7; Additional file 1: Fig. S6 and Table S10)",
     "(Additional file 1: Fig. S6 and Table S10)", 1),
    ("Because the significant axis CpGs of R4 are gene-body marks",
     "Because the significant axis CpGs identified above are gene-body marks", 1),
    ("the anchor of the actionable window identified in R7",
     "the anchor of the actionable window identified above", 1),
    ("> Reference numbers are final (renumbered 2026-09-10): Results introduces [25]",
     "> v1.7 (2026-09-10, reviewer-perspective clean-up): internal analysis codes (L0\u2013L7, R1\u2013R8) removed throughout R1\u2013R8; the four Background questions are now signposted in the Results (Q1: naive-B-cell cis-regulation; Q2: epigenetic control and its upstream circuit; Q3: production phenotype versus disease susceptibility; Q4: East Asian specificity and the intervention window); PPH1 and PPH3 are defined at first use; IVW is defined at first use; the compartment contrast in R1 is quantified; and the two single-cell B-cell partitions (3.2% versus 4.2% C1GALT1C1 detection) are now explicitly distinguished at the point of use.\n> Reference numbers are final (renumbered 2026-09-10): Results introduces [25]", 1),
])

# ============================== ch5 Discussion + Limitations ==============================
patch("05_Discussion_Limitations.md", [
    ("(R2)", "(Fig. 2)", 3),
    ("(R3, R5)", "(Fig. 3)", 1),
    ("(R4)", "(Fig. 4)", 2),
    ("(R6)", "(Fig. 5a\u2013b)", 3),
    ("(R8)", "(Additional file 1: Figures S5\u2013S6)", 2),
    ("The population dimension of the study is more qualified.",
     "The population dimension of the study is more heavily qualified.", 1),
    ("Seven limitations qualify these conclusions. First, the Han-Chinese Gd-IgA1 GWAS is released with P values only and the full-effect-size dataset (Kiryluk et al. [7], dbGaP phs000431) was not accessed;",
     "Seven limitations qualify these conclusions. First, and most fundamentally, this is a causal-inference framework rather than a causal-validation framework: it integrates and bounds causal evidence derived from observational summary data and rests on no gene-perturbation, chromatin or reporter assay, so the disease-level nulls bound the shared-causal-variant evidence available from the current GWAS rather than excluding any relationship. More specifically, the Han-Chinese Gd-IgA1 GWAS is released with P values only and the full-effect-size dataset (Kiryluk et al. [7], Database of Genotypes and Phenotypes [dbGaP] phs000431) was not accessed;", 1),
    ("C1GALT1C1 is detected in only \u2248 4% of B cells",
     "C1GALT1C1 is detected in only \u2248 4% of B cells in the marker-gated partition (3.2% in the unsupervised atlas partition)", 1),
    ("The logical next steps are therefore functional: to test whether candidate interventions restore C1GALT1C1 and galactosylation capacity along the naive-to-plasma B-cell trajectory, and to determine whether the acquired promoter hypermethylation we identify as the core of the window is itself reversible in primary IgAN B cells. Such validation\u2014including the effect of candidate interventions on the naive-to-plasma B-cell trajectory\u2014is being pursued in our laboratory with orthogonal experimental systems and will be reported separately.",
     "The logical next steps are therefore functional, and we state them as explicit tests for the orthogonal systems under way in our laboratory: (i) whether demethylating agents of the 5-azacytidine class, and DNMT1-directed compounds, restore C1GALT1C1 expression and galactosylation capacity in primary B cells from patients with IgAN; (ii) whether IL-4- and IL-17-driven hypermethylation of the C1GALT1C1 promoter is reversible in those cells, and how the methylation set point moves along the naive-to-plasma B-cell trajectory; and (iii) whether the predicted methylation-sensitive KLF2 and SP1 motifs at the C1GALT1C1 promoter CpG island are indeed occupancy-blocked by methylation, tested by chromatin immunoprecipitation and reporter assays. These experiments address the one substantive gap that a design built entirely on public summary data cannot close\u2014experimental functional validation of the circuit\u2014and will be reported separately.", 1),
    ("> Reference numbers are final (renumbered 2026-09-10): Discussion introduces [27]",
     "> v1.6 (2026-09-10, reviewer-perspective clean-up): internal analysis codes (R1\u2013R8) replaced by figure and additional-file cross-references throughout the Discussion and Limitations; the Limitations now open with an explicit statement that this is a causal-inference rather than a causal-validation framework; dbGaP is expanded at first use; the two single-cell B-cell partitions are distinguished where the 4% detection figure is quoted; and the closing paragraph now states the three functional validation experiments as an explicit list.\n> Reference numbers are final (renumbered 2026-09-10): Discussion introduces [27]", 1),
])

print("batch-1c done")
