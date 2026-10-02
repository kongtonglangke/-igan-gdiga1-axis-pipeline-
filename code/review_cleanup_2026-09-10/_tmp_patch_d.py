#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""_tmp_patch_d.py — Additional file 1 图注：去内部编号 + 修图号互指错误"""
import os

PKG = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P = os.path.join(PKG, "submission", "additional_file_1", "Additional_File_1_Figure_Legends.md")

PAIRS = [
    # 图注正文
    ("This figure supports Fig. 1b (evidence layer L0; Results R1).", "This figure supports Fig. 1b.", 1),
    ("This figure supports Fig. 2 (evidence layers L1a/L1b; Results R2).", "This figure supports Fig. 2a.", 1),
    ("This figure supports Fig. 2b (evidence layer L1b; Results R2); the complete replication matrix is in Table 3.",
     "This figure supports Fig. 2b; the complete replication matrix is in Table 3.", 1),
    ("Independent cell-level cross-check of the B-cell-state specificity of the cis-eQTL signal at R2.",
     "Independent cell-level cross-check of the B-cell-state specificity of the cis-eQTL signal (Fig. 2).", 1),
    ("in line with the cell-state-specific cis-eQTL signal reported in R2.",
     "in line with the cell-state-specific cis-eQTL signal (Fig. 2).", 1),
    ("(C1GALT1C1 is detected in \u2248 4% of B cells)",
     "(C1GALT1C1 is detected in \u2248 4% of B cells in this marker-gated partition)", 1),
    ("282,463 cells) supporting evidence layer L7 (Results R8). Cells were embedded",
     "282,463 cells). Cells were embedded", 1),
    ("This figure supports Fig. 1c (evidence layer L7; Results R8). UMAP,",
     "This figure supports Fig. 1c. UMAP,", 1),
    ("Regulator inference within the B-cell compartment of GSE285335, supporting evidence layer L7 (Results R8).",
     "Upstream regulator inference within the B-cell compartment of GSE285335.", 1),
    ("Regulator activity was estimated per cell as the mean z-score of the target genes of each CollecTRI regulon restricted to the expressed panel (317 regulators with \u2265 5 panel targets; univariate linear model)",
     "Regulator activity was estimated per cell with decoupleR's univariate linear model (ULM) applied to the target genes of each CollecTRI regulon restricted to the expressed panel (317 regulators with \u2265 5 panel targets)", 1),
    ("Because target-matched random gene sets attain \u03c1 of 0.73\u20130.78, \u03c1 alone does not establish significance",
     "Because target-matched random gene sets attain a median \u03c1 of 0.73\u20130.78 (95th percentile 0.84\u20130.88), \u03c1 alone does not establish significance", 1),
    ("This figure supports Fig. 1c (evidence layer L7; Results R8).", "This figure supports Fig. 1c.", 1),
    ("B-cell-state dependence of the *C1GALT1* cis-eQTL, supporting evidence layer L1b (Results R2).",
     "B-cell-state dependence of the *C1GALT1* cis-eQTL (Fig. 2b).", 1),
    ("This figure supports Fig. 2b (evidence layer L1b; Results R2). eQTL,",
     "This figure supports Fig. 2b. eQTL,", 1),
    # 索引表：去 R 码 + 修错指图号
    ("Corresponds to R2, Table 3 and Additional file 1: Fig. S2.",
     "Corresponds to Fig. 2, Table 3 and Additional file 1: Fig. S2.", 1),
    ("Corresponds to R3, Table 4 and Additional file 1: Fig. S3.", "Corresponds to Fig. 3 and Table 4.", 1),
    ("SuSiE-versus-ABF posterior difference and the robust-no-colocalization flag. Corresponds to R3 and Table 4.",
     "SuSiE-versus-ABF posterior difference and the robust-no-colocalization flag. Corresponds to Fig. 3 and Table 4.", 1),
    ("signal-by-signal, official credible sets). Corresponds to R3 and Table 4.",
     "signal-by-signal, official credible sets). Corresponds to Fig. 3 and Table 4.", 1),
    ("effect-allele frequency and P by ancestry. Corresponds to R5.",
     "effect-allele frequency and P by ancestry. Corresponds to the direct disease-level lookup (Results).", 1),
    ("with the exact-hit flag. Corresponds to R5.",
     "with the exact-hit flag. Corresponds to the direct disease-level lookup (Results).", 1),
    ("lead-variant hits and annotation. Corresponds to R4 and Table 5.",
     "lead-variant hits and annotation. Corresponds to Fig. 4 and Table 5.", 1),
    ("MR estimate and P (Wald ratio or clumped inverse-variance weighting). Corresponds to R4 and Table 5.",
     "MR estimate and P (Wald ratio or clumped inverse-variance weighting). Corresponds to Fig. 4 and Table 5.", 1),
    ("Corresponds to R4, Table 5 and Additional file 1: Fig. S4c.", "Corresponds to Fig. 4c and Table 5.", 2),
    ("Corresponds to R7 and Table 7.", "Corresponds to Fig. 6 and Table 7.", 2),
    ("Corresponds to R2 and Additional file 1: Fig. S4.", "Corresponds to Fig. 2 and Additional file 1: Fig. S4.", 2),
    ("242 associations, 178 factors). Corresponds to R8.",
     "242 associations, 178 factors). Corresponds to Fig. 1c and Additional file 1: Table S9.", 1),
    ("950 combinations). Corresponds to R8 and Additional file 1: Fig. S6.",
     "950 combinations). Corresponds to Fig. 1c and Additional file 1: Fig. S6.", 1),
    ("heterogeneity and meta-regression tests. Corresponds to R2 and Additional file 1: Fig. S7.",
     "heterogeneity and meta-regression tests. Corresponds to Fig. 2 and Additional file 1: Fig. S7.", 1),
    ("(Japanese; E-GEAD-398). Corresponds to R6 (evidence layer L5).",
     "(Japanese; E-GEAD-398). Corresponds to Fig. 5.", 1),
    ("nominal P value and effect slope. Corresponds to R6.",
     "nominal P value and effect slope. Corresponds to Fig. 5.", 1),
    # 版本日志
    ("> v1.4 (2026-09-10): the table index was extended",
     "> v1.5 (2026-09-10, reviewer-perspective clean-up): internal analysis codes (L0\u2013L7, R1\u2013R8) removed from the figure legends and from the table index; three wrong supplementary-figure cross-references in the index corrected (Tables S2a/S2b/S3 pointed to \"Additional file 1: Fig. S3\", which is the lead-variant eQTL panel \u2014 now Fig. 3; Tables S5c/S5d pointed to \"Additional file 1: Fig. S4c\", which is the B-cell trajectory panel \u2014 now Fig. 4c); the Fig. S6 legend now reflects the 1,000-set empirical null and the decoupleR ULM formulation.\n> v1.4 (2026-09-10): the table index was extended", 1),
]

with open(P, encoding="utf-8", newline="") as fh:
    t = fh.read()
for old, new, exp in PAIRS:
    n = t.count(old)
    if n != exp:
        raise SystemExit(f"[FAIL] found {n}, expected {exp}\n  >> {old[:150]}")
    t = t.replace(old, new)
with open(P, "w", encoding="utf-8", newline="") as fh:
    fh.write(t)
print("OK  Additional_File_1_Figure_Legends.md", len(PAIRS), "pairs")
