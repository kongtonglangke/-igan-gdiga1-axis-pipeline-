# -*- coding: utf-8 -*-
# _tmp_patch_w_full_audit_text.py — 2026-09-11 全维度审查修复批次（文本层）
# 逐条断言、幂等（重复运行会因 old 不存在而报错，属预期）。
# 覆盖：P0 药物归属（5-aza-2′-deoxycytidine/decitabine vs 5-azacytidine）、
# P0 mQTL four→two、P1 批次（PPH1 措辞、Fig3a 悬垂引用、BH q=0.12、
# 10 对口径、Abstract enriched、Wang 2021 张力句、decoupleR/Mathur 引文、
# NK 展开、拼写、图注、索引、cover letter）+ 参考文献重编号 24–27→25–28 + 新 [24]/[29]。
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os, sys

BASE = _paths.at("阶段4.5_复现包/submission")
M = BASE("manuscript")
A1 = BASE("additional_file_1", "Additional_File_1_Figure_Legends.md")

EDITS = {
os.path.join(M, "01_Title_Abstract_Keywords.md"): [
 ("cell-type-specific eQTL contextualization of four published Gd-IgA1 lead variants",
  "cell-type-specific eQTL contextualisation of four published Gd-IgA1 lead variants"),
 ("behaved as a naive B-cell-specific cis-eQTL (P = 9.36 × 10⁻⁶)",
  "behaved as a naive B-cell-enriched cis-eQTL (P = 9.36 × 10⁻⁶)"),
 ("the acquired, reversible epigenetic window at the C1GALT1C1 promoter (5-azacytidine · EGCG).",
  "the acquired, reversible epigenetic window at the C1GALT1C1 promoter (5-aza-class DNMT inhibitors · EGCG)."),
],
os.path.join(M, "03_Methods.md"): [
 ("cis-mQTL associations for the four lead variants with interpretable autosomes were queried from GoDMC (GRCh37) [13]",
  "cis-mQTL associations for the two interpretable autosomal C1GALT1 lead variants (rs13226913 and rs10238682) were queried from GoDMC (GRCh37) [13]"),
 ("the top cis-eQTL variant (lowest nominal p_beta) in every context in which the gene was measurable (10 gene–context pairs) was looked up regionally",
  "the top cis-eQTL variant (lowest nominal p_beta) in every pre-specified primary context in which the gene was measurable—the naive-B-cell and the two bulk-tissue contexts (10 gene–context pairs)—was looked up regionally; the antigen-experienced B-cell contexts were not part of the pre-specified look-up, and the X-chromosomal C1GALT1C1 records are untestable in the meta-analysis"),
 ("the lead-variant look-ups remained null (smallest adjusted q = 0.28)",
  "the lead-variant look-ups remained null (smallest adjusted q = 0.12)"),
 ("reversal of C1GALT1C1 promoter hypermethylation by the demethylating agent 5-azacytidine)",
  "reversal of C1GALT1C1 promoter hypermethylation by the demethylating agent 5-aza-2′-deoxycytidine)"),
 ("(IL-4 and IL-17 reduce C1GALT1/C1GALT1C1 expression through CpG-island hypermethylation, reversed by 5-azacytidine) [19, 20]",
  "(IL-4 reduces C1GALT1 and C1GALT1C1 expression through promoter CpG-island hypermethylation, reversed by the demethylating agent 5-aza-2′-deoxycytidine [20]; IL-17 likewise down-regulates both genes, with expression restored by 5-azacytidine [19])"),
 ("the direct binding and inhibition of DNA methyltransferase 1 and histone deacetylase 1 by the green-tea catechin epigallocatechin-3-gallate (EGCG) [21]",
  "the direct binding and inhibition of DNA methyltransferase 1 by the green-tea catechin epigallocatechin-3-gallate (EGCG) [21]"),
 ("we re-analysed the IgAN peripheral-blood mononuclear-cell (PBMC) single-cell RNA-seq dataset of Kim and colleagues (GSE285335; 26 donors: 6 late-stage IgAN, 11 early-stage IgAN and 9 healthy controls; 10x Genomics 5′ v1; 282,463 single cells)",
  "we reanalysed the IgAN peripheral-blood mononuclear-cell (PBMC) single-cell RNA-seq dataset of Kim and colleagues (GSE285335; 26 donors: 6 late-stage IgAN, 11 early-stage IgAN and 9 healthy controls; 10x Genomics 5′ gene-expression assay; 282,463 single cells)"),
 ("a 10x Genomics 5′ v1 dataset comprising 26 donors",
  "a 10x Genomics 5′ gene-expression dataset comprising 26 donors"),
 ("The output of this re-analysis is presented as exploratory",
  "The output of this reanalysis is presented as exploratory"),
 ("(NK: NKG7 or GNLY; T: CD3D;",
  "(natural killer (NK): NKG7 or GNLY; T: CD3D;"),
 ("computed with decoupleR's univariate linear model (ULM), whose per-cell test statistic",
  "computed with decoupleR's univariate linear model (ULM) [24], whose per-cell test statistic"),
 ("The 879 non-redundant JASPAR 2024 CORE vertebrate position-frequency matrices [24]",
  "The 879 non-redundant JASPAR 2024 CORE vertebrate position-frequency matrices [25]"),
 ("| JASPAR 2024 CORE | jaspar.genereg.net | [24] |",
  "| JASPAR 2024 CORE | jaspar.genereg.net | [25] |"),
 ("Reference numbers are final: [12]–[24] follow the whole-manuscript order of first appearance (renumbered 2026-09-10); [25]–[27] extend the sequence in Results/Discussion.",
  "Reference numbers are final: [12]–[25] follow the whole-manuscript order of first appearance (renumbered 2026-09-11); [26]–[28] extend the sequence in Results/Discussion."),
],
os.path.join(M, "04_Results.md"): [
 ("and posterior mass concentrated on PPH1 in the tissue eQTLs.",
  "and posterior mass concentrated on PPH1 in the comparisons with the strongest cis-eQTL evidence (PPH1 = 0.79–0.91; Table 4)."),
 ("maps to the five axis-gene regions (Fig. 3a).",
  "maps to the five axis-gene regions."),
 ("serum total IgA levels [25] maps",
  "serum total IgA levels [26] maps"),
 ("pending re-analysis when full-effect-size data become available",
  "pending reanalysis when full-effect-size data become available"),
 ("across all measurable contexts (10 gene–context pairs)",
  "across the pre-specified primary contexts (10 gene–context pairs)"),
 ("Fung et al. [26], which found, despite adequate power, no genetic evidence that Gd-IgA1-raising C1GALT1 alleles affect IgAN risk.",
  "Fung et al. [27], which found, despite adequate power, no genetic evidence that Gd-IgA1-raising C1GALT1 alleles affect IgAN risk. The earlier Gd-IgA1 GWAS had itself reported nominal candidate-level associations of the axis loci with IgAN susceptibility in a 2,352-case Chinese panel, including a GALNT12 × C1GALT1 interaction with disease risk (P = 6.6 × 10⁻³) [8]; none of these signals replicated in the much larger cross-ancestry meta-analysis used here, and we therefore regard them as unreplicated candidate-level observations."),
 ("together with hypermethylation of its promoter CpG island that is reversed by the demethylating agent 5-azacytidine [20]. C1GALT1 scored 4.5",
  "together with IL-4-driven hypermethylation of its promoter CpG island that is reversed by the demethylating agent 5-aza-2′-deoxycytidine [20]; in the IL-17 setting, 5-azacytidine likewise restored C1GALT1 and C1GALT1C1 expression [19]. C1GALT1 scored 4.5"),
 ("demethylating agents of the 5-azacytidine class provide clinical precedent [20]. These results frame",
  "demethylating agents of the 5-azacytidine class provide clinical precedent [19, 20]. These results frame"),
 ("Reference numbers are final (renumbered 2026-09-10): Results introduces [25] (serum total IgA loci) and [26] (Fung et al.)",
  "Reference numbers are final (renumbered 2026-09-11): Results introduces [26] (serum total IgA loci) and [27] (Fung et al.)"),
],
os.path.join(M, "05_Discussion_Limitations.md"): [
 ("282,463 single cells) [27] identifies 42,259 B cells",
  "282,463 single cells) [28] identifies 42,259 B cells"),
 ("cytokine-driven down-regulation of the chaperone—mediated by IL-17 and IL-4—together with hypermethylation of its promoter CpG island that is reversed by 5-azacytidine [19, 20].",
  "cytokine-driven down-regulation of the chaperone—mediated by IL-17 [19] and IL-4 [20]—together with IL-4-driven hypermethylation of its promoter CpG island that is reversed by the demethylating agent 5-aza-2′-deoxycytidine [20]; in the IL-17 setting, 5-azacytidine likewise restored C1GALT1 and C1GALT1C1 expression [19]."),
 ("demethylating agents of the 5-azacytidine class provide clinical precedent [20]; class-I",
  "demethylating agents of the 5-azacytidine class provide clinical precedent [19, 20]; class-I"),
 ("is already being tested in IgAN [2],",
  "is already being tested in IgAN [29],"),
 ("(Fung et al. [26]) is a preprint",
  "(Fung et al. [27]) is a preprint"),
 ("Reference numbers are final (renumbered 2026-09-10): Discussion introduces [27] (Kim et al.); all other citations reuse earlier numbers.",
  "Reference numbers are final (renumbered 2026-09-11): Discussion introduces [28] (Kim et al.) and [29] (Mathur et al.); all other citations reuse earlier numbers."),
],
os.path.join(M, "06_Conclusions_Abbreviations.md"): [
 ("| NFE | non-Finnish European |\n| PBMC |",
  "| NFE | non-Finnish European |\n| NK | natural killer |\n| PBMC |"),
],
os.path.join(M, "08_Figure_Legends_Tables.md"): [
 ("in the tissue eQTL comparisons the posterior mass is concentrated on PPH1 (an eQTL signal without a disease signal), which reaches 0.79–0.91 in the whole-blood and kidney-cortex contexts, whereas in several B-cell contexts the posterior is more evenly distributed (PPH1 0.18–0.28);",
  "the posterior mass is concentrated on PPH1 (an eQTL signal without a disease signal) in the comparisons with the strongest cis-eQTL evidence (PPH1 = 0.79–0.91), whereas comparisons with weaker eQTL evidence distribute the posterior more evenly (PPH1 = 0.18–0.29);"),
 ("Signal-level re-testing of the four combined SuSiE credible sets (region-level coloc ABF PPH4, blue, versus SuSiE signal-level PPH4, orange).",
  "Signal-level re-testing across the four comparisons with official SuSiE credible sets (seven signal-level tests; region-level coloc ABF PPH4, blue, versus SuSiE signal-level PPH4, orange)."),
 ("The only level-3 entry is the reversal of C1GALT1C1 promoter hypermethylation by 5-azacytidine in an IgAN-relevant B-cell context.",
  "The only level-3 entry is the reversal of C1GALT1C1 promoter hypermethylation by the demethylating agent 5-aza-2′-deoxycytidine in an IgAN-relevant B-cell context [20]; in an IL-17 context, 5-azacytidine likewise restored C1GALT1 and C1GALT1C1 expression [19]."),
 ("(IL-4-driven hypermethylation of the C1GALT1C1 promoter, reversed by 5-azacytidine [20])",
  "(IL-4-driven hypermethylation of the C1GALT1C1 promoter, reversed by 5-aza-2′-deoxycytidine [20])"),
],
A1: [
 ("Public 10x Genomics 5′ v1 PBMC scRNA-seq data of 26 IgAN donors (6 late-stage, 11 early-stage, 9 healthy controls; 282,463 single cells; Kim et al. [27]; GEO: GSE285335)",
  "Public 10x Genomics 5′ gene-expression PBMC scRNA-seq data of 26 donors (17 with IgAN: 6 late-stage, 11 early-stage; 9 healthy controls; 282,463 single cells; Kim et al. [28]; GEO: GSE285335)"),
 ("aggregated to donor × B-cell-state pseudobulk units (103 units) and related",
  "aggregated to donor × B-cell-state pseudobulk units (103 units across four B-cell states and 26 donors) and related"),
 ("rs10238682: β = +0.186, P = 1.2 × 10⁻⁴",
  "rs10238682: β = +0.185, P = 1.16 × 10⁻⁴"),
 ("attenuated roughly two- to five-fold in antigen-experienced states",
  "attenuated two- to twelve-fold in antigen-experienced states"),
 ("Cochran's Q (p = 0.046 and 0.047; I² ≈ 67%)",
  "Cochran's Q (P = 0.046 and 0.047; I² ≈ 67%)"),
 ("meta-regression (p = 0.018 for both variants)",
  "meta-regression (P = 0.018 for both variants)"),
 ("did not reach significance (p = 0.11)",
  "did not reach significance (P = 0.11)"),
 ("| Table S1 | Full 5-gene × 5-context cis-eQTL matrix (OneK1K naive, memory and intermediate B cells; GTEx v10 whole blood and kidney cortex): −log₁₀(P), per-allele β, standard error, allele frequency and the per-cell call. Corresponds to Fig. 2, Table 3 and Additional file 1: Fig. S2. |",
  "| Table S1 | Full 5-gene × 5-context cis-eQTL matrix (OneK1K naive, memory and intermediate B cells; GTEx v10 whole blood and kidney cortex): the top cis-variant per gene–context with chromosome, position, per-allele β, nominal and permutation P and the number of individuals. Corresponds to Fig. 2, Table 3 and Additional file 1: Fig. S2. |"),
 ("| Table S2b | Computability verdict and robustness summary for the 14 computable colocalisation comparisons: credible-set size, ABF posteriors, SuSiE-versus-ABF posterior difference and the robust-no-colocalisation flag. Corresponds to Fig. 3 and Table 4. |",
  "| Table S2b | SuSiE-versus-ABF robustness summary for the four colocalisation comparisons with official credible sets: credible-set size, ABF posteriors, SuSiE-versus-ABF posterior difference and the robust-no-colocalisation flag (the computability verdicts for all 25 planned comparisons are in Table S2a). Corresponds to Fig. 3 and Table 4. |"),
 ("| Table S4a | Per-lead IgAN lookup: the two C1GALT1 leads (rs13226913, rs10238682) and their cis-eQTL variants looked up in the cross-ancestry IgAN meta-analysis (GCST90018866), with IgAN β, SE, effect-allele frequency and P by ancestry. Corresponds to the direct disease-level lookup (Results). |",
  "| Table S4a | Per-context top-variant IgAN lookup: the top cis-eQTL variant of each axis gene in each pre-specified primary context (10 gene–context pairs) tested for exact presence in the cross-ancestry IgAN meta-analysis (GCST90018866), with IgAN β, SE, effect-allele frequency and P where the variant is found (no top variant was present in the release; the per-lead readout is reported in the main text and the window-level minima in Table S4b). Corresponds to the direct disease-level lookup (Results). |"),
 ("| Table S10 | CollecTRI regulator-by-B-cell-state associations with the composite axis score (ρ, nominal P, FDR q, target-matched random-set 95th percentile and empirical P; 950 combinations). Corresponds to Fig. 1c and Additional file 1: Fig. S6. |",
  "| Table S10 | CollecTRI regulator-by-B-cell-state associations with each axis-gene target and the composite axis score (Spearman ρ, nominal P, Benjamini–Hochberg q and number of panel targets; 950 regulator–state combinations × six targets = 5,700 rows). The target-matched random-set null and the empirical P values are summarised in the Fig. S6 legend and in the main text. Corresponds to Fig. 1c and Additional file 1: Fig. S6. |"),
],
BASE("cover_letter_draft.md"): [
 ("including 5-azacytidine-reversible promoter hypermethylation",
  "including DNMT-inhibitor–reversible promoter hypermethylation"),
 ("**naive B-cell-specific cis-eQTL**",
  "**naive B-cell-enriched cis-eQTL**"),
 ("Mendelian randomisation provides",
  "Mendelian randomization provides"),
],
}

def patch_file(path, pairs):
    with open(path, "r", encoding="utf-8", newline="") as f:
        text = f.read()
    for i, (old, new) in enumerate(pairs):
        n = text.count(old)
        assert n == 1, f"{'MISS' if n==0 else 'DUP'} {os.path.basename(path)} #{i}: {old[:70]!r} (count={n})"
        text = text.replace(old, new)
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write(text)
    print(f"OK  {os.path.basename(path)}  ({len(pairs)} edits)")

def renumber_references(path):
    with open(path, "r", encoding="utf-8", newline="") as f:
        text = f.read()
    # 降序重编号 27→28, 26→27, 25→26, 24→25（仅条目行首）
    for old_n, new_n in [(27, 28), (26, 27), (25, 26), (24, 25)]:
        head = f"\n{old_n}. "
        assert text.count(head) == 1, f"entry {old_n} head count != 1"
        text = text.replace(head, f"\n{new_n}. ")
    decoupler = ("24. Badia-i-Mompel P, Vélez Santiago J, Braunger J, Geiss C, Dimitrov D, Müller-Dott S, et al. "
                 "decoupleR: ensemble of computational methods to infer biological activities from omics data. "
                 "Bioinform Adv. 2022;2(1):vbac016. https://doi.org/10.1093/bioadv/vbac016\n")
    anchor = "\n25. Rauluseviciute"
    assert text.count(anchor) == 1
    text = text.replace(anchor, "\n" + decoupler + anchor)
    mathur = ("29. Mathur M, Barratt J, Chacko B, Chan TM, Kooienga L, Oh KH, et al. "
              "A Phase 2 Trial of Sibeprenlimab in Patients with IgA Nephropathy. "
              "N Engl J Med. 2024;390(1):20–31. https://doi.org/10.1056/NEJMoa2305635")
    assert "\n28. Kim G" in text
    if not text.endswith("\n"):
        text += "\n"
    text += mathur + "\n"
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write(text)
    print("OK  09_References.md  (renumber 24–27→25–28, new [24] decoupleR, new [29] Mathur)")

if __name__ == "__main__":
    for path, pairs in EDITS.items():
        patch_file(path, pairs)
    renumber_references(os.path.join(M, "09_References.md"))
    print("ALL TEXT EDITS APPLIED")
