#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""_tmp_patch_batch1b.py — 定点改写：ch1 ch2 ch6 ch7 ch8 ch9（批次1 收尾 + 批次2 部分）"""
import os
import sys

MS = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                  "submission", "manuscript")


def patch(fn, pairs):
    p = os.path.join(MS, fn)
    with open(p, encoding="utf-8", newline="") as fh:
        t = fh.read()
    for old, new, exp in pairs:
        n = t.count(old)
        if n != exp:
            raise SystemExit(f"[FAIL] {fn}: found {n}, expected {exp}\n  >> {old[:130]}")
        t = t.replace(old, new)
    with open(p, "w", encoding="utf-8", newline="") as fh:
        fh.write(t)
    print(f"OK  {fn}  ({len(pairs)} pairs)")


# ============================== ch1 ==============================
patch("01_Title_Abstract_Keywords.md", [
    ("Across 25 gene-by-cell-context colocalisation tests against IgAN, no shared causal structure was supported (maximum PPH4 = 0.031).",
     "Of 25 planned gene-by-cell-context colocalisation tests against IgAN, 14 were computable; none supported a shared causal variant (maximum PPH4 = 0.031).", 1),
    ("East Asian enrichment was strongest for the C1GALT1 lead (frequency 0.522 vs 0.204), and the intervention window mapped to reversible acquired epigenetic control (C1GALT1C1 score 6.0).",
     "East Asian differentiation was strongest at the C1GALT1 lead rs10238682 (0.522 vs 0.204), whereas the discordant GALNT12 lead rs7856182 was rarer in East Asians (0.025 vs 0.154). The intervention window mapped to reversible acquired epigenetic control (C1GALT1C1 score 6.0).", 1),
    ("5. **Upstream regulatory layer (L7, added 2026-09-10)** — an annotation strip beneath the epigenetic window, \"upstream circuit (L7): cytokine \u2192 DNMT1 \u2192 Sp/KLF\", with an arrow feeding into that window; this reports the single-cell regulon result of Results R8.",
     "5. **Upstream regulatory layer** — an annotation strip beneath the epigenetic window, \"upstream circuit: cytokine \u2192 DNMT1 \u2192 Sp/KLF\", with an arrow feeding into that window; this reports the single-cell regulon result on the cytokine\u2013DNMT1\u2013Sp/KLF circuit at the C1GALT1C1 promoter.", 1),
    ("> Revision in v1.4 (2026-09-10): Plan A integration",
     "> Revision in v1.5 (2026-09-10, reviewer-perspective clean-up): internal analysis codes (L0\u2013L7, R1\u2013R8) removed from the Abstract and Graphical-abstract text; the Abstract colocalisation sentence now reads \"25 planned \u2026 14 were computable\"; the East Asian sentence now names the leads (rs10238682 and the discordant rs7856182). Word count re-measured below.\n> Revision in v1.4 (2026-09-10): Plan A integration", 1),
])

# ============================== ch2 ==============================
patch("02_Background.md", [
    ("its systemic output\u2014circulating Gd-IgA1\u2014is a heritable quantitative trait [6]",
     "its systemic output\u2014circulating Gd-IgA1\u2014is a quantitative trait influenced by common genetic variation [6, 7, 8]", 1),
    ("Genetic studies of this axis have so far operated at two disconnected levels.",
     "Genetic studies of this axis have so far proceeded at two disconnected levels.", 1),
    ("> v1.4 (2026-09-10): Plan A integration of evidence layer L7",
     "> v1.5 (2026-09-10, reviewer-perspective clean-up): \"operated at two disconnected levels\" \u2192 \"proceeded at\"; the heritability framing of circulating Gd-IgA1 now reads \"a quantitative trait influenced by common genetic variation [6, 7, 8]\" (the earlier \"heritable\" implied family/twin-sib evidence that reference [6], an association study, does not provide).\n> v1.4 (2026-09-10): Plan A integration of evidence layer L7", 1),
])

# ============================== ch6 ==============================
patch("06_Conclusions_Abbreviations.md", [
    ("| 5-azaC | 5-azacytidine |\n", "", 1),
    ("Internal layer labels (L0\u2013L7, R1\u2013R8, M1\u2013M11) are removed from the submitted text.",
     "Internal analysis labels (the layer and result codes used during the analysis) are not used in the running text of the submitted manuscript.", 1),
    ("> Reference numbers: none in this chapter.",
     "> v1.4 (2026-09-10, reviewer-perspective clean-up): the orphan abbreviation \"5-azaC\" (never used in the text, which spells out 5-azacytidine throughout) was deleted; the note on internal labels was reworded so that it no longer reproduces the codes it declares unused.\n> Reference numbers: none in this chapter.", 1),
])

# ============================== ch7 ==============================
patch("07_TitlePage_Declarations.md", [
    ("including the single-cell atlas and upstream-regulator inference of evidence layer L7)",
     "including the single-cell atlas and upstream-regulator inference)", 1),
    ("This study received no specific grant from any funding agency in the public, commercial or not-for-profit sectors.",
     "This study received no specific grant from any funding agency in the public, commercial or not-for-profit sectors. No funding body had any role in the design of the study, the collection, analysis or interpretation of data, or the writing of the manuscript.", 1),
    ("YW conceived and supervised the study. JB and YW designed the analyses. JB performed the bioinformatics and statistical analyses. TH prepared the figures. MX, NY and WP contributed to data curation and to interpretation of the results. BW and WP provided clinical interpretation. NY and YW wrote and critically revised the manuscript. All authors read and approved the final manuscript.",
     "YW conceived and supervised the study and obtained the resources. JB and BW contributed equally as co-first authors: JB and YW designed the analyses, and JB performed the bioinformatics and statistical analyses and the single-cell reanalysis. BW provided the clinical interpretation of the axis and supervised the clinical framing of the Gd-IgA1\u2013IgAN relationship. TH prepared the figures. MX, NY and WP contributed to data curation and to interpretation of the results. WP provided clinical interpretation. NY and YW wrote and critically revised the manuscript. All authors read and approved the final manuscript.", 1),
])

# ============================== ch8 ==============================
patch("08_Figure_Legends_Tables.md", [
    ("(evidence layers L3\u2013L4: colocalisation, 14 of 25 computable pairs all null, maximum PPH4 = 0.031; allelic-score null)",
     "(colocalisation, 14 of 25 computable pairs all null, maximum PPH4 = 0.031; allelic-score null)", 1),
    ("(evidence layer L6: *C1GALT1C1* promoter, 5-azacytidine-reversible hypermethylation)",
     "(*C1GALT1C1* promoter, 5-azacytidine-reversible hypermethylation)", 1),
    ("the seventh row is the upstream regulatory layer (L7: single-cell regulon inference placing the axis under a ",
     "the seventh row is the upstream regulatory layer (single-cell regulon inference placing the axis under a ", 1),
    ("the independent East Asian regulatory cross-reference of R6 (ImmuNexUT; evidence layer L5)",
     "the independent East Asian regulatory cross-reference (ImmuNexUT)", 1),
    ("> v1.5 (2026-09-10, supplementary-table suffixes)",
     "> v1.6 (2026-09-10, reviewer-perspective clean-up): internal analysis codes removed from the legends (the Fig. 1a reference to \"evidence layers L3\u2013L4\" and \"evidence layer L6\", the Fig. 1c phrase \"(L7: \u2026)\", and the Fig. 5 cross-reference \"R6 (ImmuNexUT; evidence layer L5)\"). Fig. 1c retains its \"(L0\u2013L7)\" index because the panel itself prints the eight numbered layer rows.\n> v1.5 (2026-09-10, supplementary-table suffixes)", 1),
])

# ============================== ch9 ==============================
patch("09_References.md", [
    ("Sci Rep. 2025;16(1):1700.", "Sci Rep. 2026;16(1):1700.", 1),
    ("(2025;16(1):1700, DOI 10.1038/s41598-025-31213-9, PMID 41466033)",
     "(2026;16(1):1700, DOI 10.1038/s41598-025-31213-9, PMID 41466033; published online 30 December 2025, assigned to the 2026 volume)", 1),
    ("> v1.3 (2026-09-10) \u2014 two integrity corrections applied after a full reference audit.",
     "> v1.5 (2026-09-10, reviewer-perspective clean-up): entry [27] Kim et al. volume year corrected from 2025 to **2026** (Sci Rep 16(1):1700; published online 30 December 2025), verified against the publisher page; the manuscript title now reads \"causal-inference framework\". No other bibliographic change.\n> v1.3 (2026-09-10) \u2014 two integrity corrections applied after a full reference audit.", 1),
])

print("batch-1b done")
