# -*- coding: utf-8 -*-
# _tmp_patch_w4_version_notes.py — 2026-09-11 全维度审查修复批次（版本注 + 词数声明）
# 向 01/03/04/05/06/08/09 头部注块末尾追加本轮版本注（行尾自适应 CRLF/LF）。
# 词数均为 _tmp_measure_wordcount_2026-09-11.py 实测值（口径见测量脚本头注）。
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os

M = _paths.at("阶段4.5_复现包/submission/manuscript")

NOTES = {
"01_Title_Abstract_Keywords.md": (
 '> Revision in v1.9 (2026-09-10, third reviewer-perspective pass): the opening Abstract sentence was restructured to remove a dangling participle ("Produced by …, whether genetic variation … remains unresolved" → "Gd-IgA1, the central pathogenic molecule …, is produced by …; whether … remains unresolved"); the word count is unchanged at 348 (limit 350). No data, values or conclusions changed.',
 '> v2.0 (2026-09-11, full-coverage audit): three wording corrections — the graphical-abstract window line now reads "(5-aza-class DNMT inhibitors · EGCG)" (was "5-azacytidine · EGCG"), matching the PubMed-verified split attribution of the demethylation evidence ([20] used 5-aza-2′-deoxycytidine; [19] used 5-azacytidine without a methylation assay); the Abstract eQTL phrase now reads "naive B-cell-enriched" (was "naive B-cell-specific"); and "contextualization" was corrected to the British "contextualisation". Word count re-measured on the structured Abstract: unchanged at **348** (limit 350; labels included, headings excluded).'
),
"03_Methods.md": (
 '> v1.10 (2026-09-10, third reviewer-perspective pass): the direct-look-up distance wording now reads "approximately two orders of magnitude above the genome-wide significance threshold (P = 9.4 × 10⁻⁶ versus 5 × 10⁻⁸)" (was "three orders of magnitude from genome-wide significance"), matching the Results wording for the same value; IDG, the Pharos target-development levels and NBDC are now expanded at first use. No data, values or conclusions changed.',
 '> v1.11 (2026-09-11, full-coverage audit): (i) reference-level corrections — the demethylation literature anchor is now split by drug ([20] 5-aza-2′-deoxycytidine for the IL-4-driven C1GALT1C1 promoter hypermethylation reversal; [19] 5-azacytidine for IL-17 expression restoration), the D = 3 anchor is reworded accordingly, and the EGCG sentence names DNMT1 only (HDAC1 removed; [21] shows no HDAC1 binding); (ii) data corrections — the cis-mQTL query now names the two interpretable autosomal C1GALT1 leads (was "four lead variants"), the lead-variant look-up null result reads "smallest adjusted q = 0.12" (was 0.28, inconsistent with the stated 13-test BH), the look-up is restricted to the pre-specified primary contexts with the antigen-experienced contexts and the X-chromosomal records explicitly excluded, and the Table 1 GoDMC sample-size range is corrected to 24,988–28,181 per CpG (was 25,095–28,181; Table S5a actuals); (iii) wording — the 10x chemistry is given as "5′ gene-expression" without a kit version (three places, including Table 1), "reanalysed/reanalysis" spelling, NK expanded at first use, decoupleR cited as [24] and JASPAR as [25] (reference renumber, see 09_References.md). Word counts re-measured on the chapter body (raw whitespace tokenisation; Markdown separators, table rows and the editorial header block excluded; Table 1 excluded in full — caption, footnote and body rows): **3,635** including the sixteen subsection headings, **3,530** excluding headings (under the previous rows-only-excluded convention the figures are 3,671/3,566, i.e. a net +45 words from this audit; the declared figures going forward use the strict full-exclusion convention).'
),
"04_Results.md": (
 '> v1.9 (2026-09-10, third reviewer-perspective pass): the methylation section now cites Table 5 in the running text (journal style requires every main-text table to be cited in the text, and Table 5 had been uncited); the B-cell-state heterogeneity test is now attributed explicitly to rs13226913; Sp/KLF is expanded at first use; the direct-look-up distance wording is aligned with Methods. No data, values or conclusions changed.',
 '> v2.0 (2026-09-11, full-coverage audit): (i) the PPH1 sentence is reframed to the measured groups — "posterior mass concentrated on PPH1 in the comparisons with the strongest cis-eQTL evidence (PPH1 = 0.79–0.91; Table 4)" (the weaker group spans 0.18–0.29); (ii) a dangling "(Fig. 3a)" call removed; (iii) the direct look-up is described as covering the pre-specified primary contexts; (iv) the Fung et al. paragraph gains the Wang 2021 candidate-level tension sentence (the earlier Gd-IgA1 GWAS reported nominal susceptibility associations, including GALNT12 × C1GALT1 interaction P = 6.6 × 10⁻³ [8], unreplicated in the cross-ancestry meta-analysis); (v) the demethylation-evidence sentence is split by drug (5-aza-2′-deoxycytidine [20] for the IL-4 methylation reversal; 5-azacytidine [19] for IL-17 expression restoration) and the clinical-precedent citation widened to [19, 20]; (vi) [25]→[26], [26]→[27] per the reference renumber. Word counts re-measured on the chapter body (separators, table rows, the planned-table appendix and the editorial header block excluded): **3,292** including the eight subsection headings, **3,175** excluding headings (superseding 3,211/3,094).'
),
"05_Discussion_Limitations.md": (
 '> v1.8 (2026-09-10, author decision on the closing paragraph): the final paragraph no longer states that the three functional tests are "under way in our laboratory" or that they "will be reported separately"; the experimental predictions are now presented as what the framework makes directly testable and as the direction for the next studies of the axis, and the paragraph closes on the boundary-and-direction summary of the study\'s contribution. The change removes a claim about the authors\' own in-progress work that could not be substantiated; no data, values or conclusions changed.',
 '> v1.9 (2026-09-11, full-coverage audit): the chaperone-paragraph demethylation evidence is split by drug ([20] 5-aza-2′-deoxycytidine for the IL-4-driven C1GALT1C1 promoter hypermethylation reversal; [19] 5-azacytidine for IL-17 expression restoration) with the clinical-precedent citation widened to [19, 20]; the sibeprenlimab precedent is now cited to the phase-2 publication [29] (was [2]); the single-cell chaperone detection figure is corrected to 4.1% of all B cells (was 4.2%; Table S7 weighted (737 + 135 + 859)/42,259 = 4.097%); Kim et al. → [28] and Fung et al. → [27] per the reference renumber. Word counts re-measured on the chapter body (separators and the editorial header block excluded): Discussion **1,618** including the section heading / **1,616** excluding (superseding 1,600/1,598); Limitations **605** / **603** (unchanged).'
),
"06_Conclusions_Abbreviations.md": (
 '> v1.6 (2026-09-10, third reviewer-perspective pass): the abbreviation entry "PPH3/PPH4" now reads "PPH0–PPH4" (posterior probability of the five colocalisation hypotheses H0–H4), matching the notation defined in Results; the closing note no longer claims that all gene symbols are expanded at first appearance. No other changes.',
 '> v1.7 (2026-09-11, full-coverage audit): the abbreviation list gains **NK | natural killer**, used by the single-cell atlas marker gate in Methods (the gate text now expands NK at first use there as well). Conclusions text unchanged (148 words including / 146 excluding the section heading).'
),
"08_Figure_Legends_Tables.md": (
 '> v1.9 (2026-09-10, figure-publishability audit): Fig. 1c legend ordinal corrected — the upstream regulatory layer is the **eighth** (final) row of the panel, not the seventh, matching the eight numbered layer rows the panel prints. Companion figure re-renders from the same audit: **Fig. 2 v3.5, Fig. 3 v3.4, Fig. 4 v3.5, Fig. 6 v3.4** (British "grey"/"colocalisation" spellings; Fig. 4c legend moved clear of the PPH4 = 0.05 threshold line; Fig. 6b axis-label order aligned to the plotted column order) and supplementary **Fig. S2** (panel title de-coded; delivered file renamed `Figure_S2_5x5_eQTL_heatmap`) and **Fig. S3** (TIFF regenerated from the corrected render). No data, values or conclusions changed.',
 '> v2.1 (2026-09-11, full-coverage audit): legend-integrity fixes — Fig. 1a window label now reads "DNMT-inhibitor–reversible hypermethylation" (was "5-azacytidine-reversible"); Fig. 3b PPH1 wording quotes the measured groups (strongest-evidence comparisons 0.79–0.91 versus weaker 0.18–0.29; Table 4); Fig. 3c states "seven signal-level tests across the four comparisons with official SuSiE credible sets" (Table S2b: 1 + 2 + 1 + 3 credible sets); the Fig. 6b level-3 entry is corrected to 5-aza-2′-deoxycytidine [20] with the IL-17 5-azacytidine clause [19]; Fig. 6d reads "reversed by 5-aza-2′-deoxycytidine [20]". Companion figure re-renders: **Fig. 1 v3.8, Fig. 6 v3.5, GA v3.4** (in-panel drug labels class-levelled to "5-aza DNMTi" / "5-aza class" / "decitabine"; layout self-checks pass with zero overlaps). No data, values or conclusions changed.'
),
"09_References.md": (
 '> v1.6 (2026-09-11, draft-to-package review): entry [21] Fang 2003 article title corrected to "…reactivates methylation-silenced genes in cancer cell lines." (verified against PubMed PMID 14633667; Cancer Res 2003;63(22):7563–7570 — the journal page gives "cell lines", not "cells"; the absence of a registered DOI for this AACR 2003 article was re-confirmed). The v1.4 audit note above, which mis-typed the Fang entry as "[20]", now reads "[21]". Entry count and numbering unchanged (27).',
 '> v2.0 (2026-09-11, full-coverage audit): one renumber step and two new entries — (i) new entry [24] Badia-i-Mompel 2022 (decoupleR; Bioinform Adv 2(1):vbac016) is introduced in Methods (upstream-regulator inference) ahead of JASPAR, so the former [24]–[27] shift to [25]–[28]; (ii) new entry [29] Mathur 2024 (sibeprenlimab ENVISION phase 2; N Engl J Med 390(1):20–31) is introduced in the Discussion. Entry count 27 → 29; the strictly ascending first-appearance sequence was re-verified across the manuscript. Both new entries were verified against the publisher records (DOIs 10.1093/bioadv/vbac016 and 10.1056/NEJMoa2305635). [27] Fung 2025 remains a medRxiv preprint (re-check immediately before submission).'
),
}

def patch(path, anchor, note):
    with open(path, "r", encoding="utf-8", newline="") as f:
        text = f.read()
    assert text.count(anchor) == 1, f"anchor not unique in {os.path.basename(path)}"
    eol = "\r\n" if "\r\n" in text else "\n"
    text = text.replace(anchor, anchor + eol + note)
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write(text)
    print(f"OK  {os.path.basename(path)}  (version note appended, eol={eol!r})")

if __name__ == "__main__":
    for name, (anchor, note) in NOTES.items():
        patch(M(name), anchor, note)
    print("ALL VERSION NOTES APPENDED")
