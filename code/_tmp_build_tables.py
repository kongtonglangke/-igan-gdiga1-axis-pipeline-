# -*- coding: utf-8 -*-
"""Build 章节8_Figure_legends_and_Tables v1.0 EN (step 4.4 integration) — r2."""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import pandas as pd, os

ROOT = _paths.LEGACY
OUT = ROOT("阶段4/初稿/章节8_Figure_legends_and_Tables_v1.0_纯英文版.md")

# ---- anchors: values quoted in the approved Results text keep their text form ----
TEXT_P = {
    6.290992e-11: "6.3×10⁻¹¹", 8.505434e-10: "8.5×10⁻¹⁰", 5.581684e-04: "5.6×10⁻⁴",
    4.619283e-06: "4.6×10⁻⁶", 5.054155e-10: "5.1×10⁻¹⁰", 3.4e-05: "3.4×10⁻⁵",
    0.004294: "0.004", 0.003296: "0.003", 0.002552: "0.003", 0.001737: "0.002",
    0.986257: "0.986", 0.09496: "0.095",
    5.575410e-12: "5.6×10⁻¹²", 1.9062e-05: "1.9×10⁻⁵", 9.358880e-06: "9.36×10⁻⁶",
    5.6044e-02: "0.056", 0.400064: "0.40", 1.162100e-04: "1.16×10⁻⁴",
    2.99947e-11: "3.0×10⁻¹¹", 0.204596: "0.20", 0.231335: "0.23",
    0.31418: "0.31", 0.515319: "0.52",
    1.203e-09: "1.2×10⁻⁹", 2.379e-09: "2.4×10⁻⁹", 0.01088: "0.011", 0.004042: "0.004",
    0.2842: "0.284", 0.222: "0.222", 0.1004: "0.100",
    0.063916: "0.064",
}
SUP = {"0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴",
       "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹"}

def fmt_p(p):
    p = float(p)
    for k, v in TEXT_P.items():
        if abs(p - k) <= max(1e-12, abs(k) * 1e-3):
            return v
    if p >= 0.01:
        return f"{p:.3g}"
    if p >= 0.001:
        return f"{p:.4f}"
    m, e = f"{p:.1e}".split("e")
    e = int(e)
    es = "".join(SUP[c] for c in str(abs(e)))
    return f"{m}×10{'⁻' if e < 0 else ''}{es}"

def fmt_fc(x):
    x = float(x)
    if abs(x) < 0.005:
        return "0.00"
    return f"{x:+.2f}"

lines = []
A = lines.append

A("# Figure legends and Tables (integration appendix)")
A("")
A("> Manuscript: Cell-type-specific genetic and epigenetic control of the galactose-deficient IgA1 axis in IgA nephropathy: a multi-omics causal framework with East Asian perspectives")
A("> Format: Genome Medicine (BMC). This appendix holds the figure legends (to be placed after the Reference list in the main manuscript file) and the main-text Tables 2–7 (Table 1 is embedded in Methods).")
A("> Version: v1.0 (2026-09-06) — step 4.4 integration draft; values transcribed programmatically from the analysis TSVs; P values shown in the text-quoted form where the number is stated in Results (v1.1).")
A("")

A("## Figure legends")
A("")
legends = [
"""**Fig. 1 | The O-glycosylation biosynthesis axis is attenuated in peripheral blood but directionally separated in kidney.** **(a)** Conceptual framework of the study. B/plasma cells synthesise IgA1 through the IgA1 O-glycosylation biosynthesis axis (five enzymes: *C1GALT1*, *C1GALT1C1*/*COSMC*, *GALNT2*, *GALNT12*, *ST6GALNAC2*); the galactose-deficient output (Gd-IgA1) is a quantitative trait with a strong genetic anchor in Han Chinese patients (genome-wide significant lead rs10238682, P = 1.2 × 10⁻⁹), whereas mesangial deposition of galactose-deficient IgA1 and progression to IgAN require downstream hits (mucosal immunity, complement). The framework separates two layers: germline genetics of the axis, which do not measurably drive IgAN risk (evidence layers L3–L4: colocalization, 14 of 25 computable pairs all null, maximum PPH4 = 0.031; allelic-score null), and acquired epigenetic control of the axis, the actionable window of the study (evidence layer L6: *C1GALT1C1* promoter, 5-azacytidine-reversible hypermethylation). **(b)** Directional expression of the five axis genes across three disease transcriptome datasets: log2 fold change (IgAN versus control) in GSE73953 peripheral blood mononuclear cells (IgAN n = 15 versus two pooled healthy controls), GSE115857 kidney biopsy bulk tissue (n = 55 vs 7; DASL platform) and GSE93798 microdissected glomeruli (n = 20 vs 22). Asterisks mark |logFC| > 1 with P < 0.05. All five genes are down-regulated in PBMCs (logFC −1.73 to −4.43) and show a directionally different pattern in kidney (e.g. *GALNT2* and *ST6GALNAC2* up-regulated in kidney bulk; *C1GALT1* down-regulated in glomeruli). Because the PBMC series pools only two healthy controls, the PBMC row is a directional signal. **(c)** Evidence-chain summary of the layer conclusions (L0–L6) with the key measured anchors and the figure/table in which each is reported; Q denotes the pre-specified research question of each block (Q1–Q4) as defined in the Results. PBMC, peripheral blood mononuclear cell; DASL, cDNA-mediated annealing, selection, extension and ligation; PPH4, posterior probability of a shared causal variant in coloc.""",
"""**Fig. 2 | C1GALT1 cis-regulation is enriched in naive B cells: cell-type-specific eQTL context of the axis.** **(a)** cis-eQTL evidence for the five axis genes (rows) across three OneK1K B-lineage contexts (naive, memory, intermediate; n = 876/726/618) and two GTEx v10 tissues (whole blood n = 853; kidney cortex n = 106), coloured by −log₁₀(P) (cap = 8) with per-allele effect β annotated where computable. Hatched grey cells mark contexts in which the gene is not testable: *C1GALT1C1* is X-chromosomal and absent from the eQTL Catalogue release 8 ("X-linked"), whereas *GALNT12* and *ST6GALNAC2* are filtered from the B-lineage contexts by low expression ("B low expr."). **(b)** cis-eQTL replication of the four published Gd-IgA1 lead variants (rs13226913 and rs10238682, *C1GALT1*; rs7856182, *GALNT12*; rs5910940, *C1GALT1C1*) for their axis gene across the same five contexts, plotted as −log₁₀(P) (permutation P where available, nominal otherwise); the red dashed line marks the genome-wide significance threshold (5 × 10⁻⁸) and the red box highlights the rs13226913–naive-B-cell signal (P = 9.36 × 10⁻⁶). Non-testable cells are hatched as in **a**. **(c)** Per-allele cis-eQTL effect (β ± SE) of the two *C1GALT1* leads along the B-cell differentiation gradient (naive → memory → intermediate). Both leads show a monotonic decrease from naive B cells (rs13226913: 0.17 ± 0.04, P = 9.36 × 10⁻⁶; rs10238682: 0.19 ± 0.05, P = 1.16 × 10⁻⁴) to intermediate B cells (0.035 ± 0.042 and 0.016 ± 0.052, respectively), placing the strongest cis-regulation at the naive stage. eQTL, expression quantitative trait locus.""",
"""**Fig. 3 | The genetic signal of the axis is trait-specific for Gd-IgA1 and shows no shared structure with IgAN susceptibility.** **(a)** Association strength (−log₁₀ P) of the four published Gd-IgA1 lead variants with serum Gd-IgA1 levels in Han Chinese patients (Wang et al. 2021 [9], red) and with IgAN susceptibility in the cross-ancestry meta-analysis (Sakaue et al. 2021 [10], blue). Vertical dashed lines mark P = 5 × 10⁻⁸ (genome-wide significance) and P = 0.05 (nominal). rs5910940 maps to the X chromosome, which is absent from the IgAN meta-analysis (hatched, "chrX not in GWAS"). None of the leads reaches nominal significance for IgAN (P = 0.10–0.28). **(b)** Colocalization posteriors (PPH0–PPH4; coloc ABF, p1 = p2 = 1 × 10⁻⁴, p12 = 1 × 10⁻⁵) for the 14 computable gene-by-context comparisons of the 25 planned (11 not computable: *C1GALT1C1* in all contexts and *GALNT12*/*ST6GALNAC2* in the B-cell contexts). PPH4 never exceeds 0.5 (maximum 0.0308, ST6GALNAC2 × whole blood; red dashed line) and PPH1 (eQTL signal only) dominates (0.79–0.91); the allelic score of the three autosomal leads against IgAN is null (β = 0.0039, SE = 0.025, P = 0.876). **(c)** Signal-level re-testing of the four combined SuSiE credible sets (region-level coloc ABF PPH4, blue, versus SuSiE signal-level PPH4, orange). All signal-level PPH4 ≤ 0.031 and none reaches the conventional PPH4 = 0.8 threshold. Gd-IgA1, galactose-deficient IgA1; PPH, posterior probability of hypothesis H (coloc: H1 eQTL only; H3 distinct signals; H4 shared causal variant); GWAS, genome-wide association study.""",
"""**Fig. 4 | Gene-body CpG methylation marks transcriptional coupling without mediating IgAN risk.** **(a)** Landscape of the five significant cis-mQTL CpGs (P < 5 × 10⁻⁸, GoDMC whole-blood mQTLs) across the *C1GALT1* gene body (chr7:7,196,565–7,288,282, hg19). Bars give the cis-mQTL effect per rs13226913 effect allele (sign relative to the eQTL-raising allele; GoDMC whole blood): four intragenic CpGs at +25.7, +28.3, +65.0 and +72.2 kb from the transcription start site (TSS) and one in the 3′-flanking region (+94.9 kb); no significant CpG lies in the promoter (hatched box, TSS ± 1.5 kb). **(b)** Mendelian randomization of CpG methylation on IgAN risk (Wald ratio using the lead variant rs13226913 or the top cis-sentinel variant of each CpG, and inverse-variance weighting where ≥ 3 independent instruments were available). All 11 estimates are null (P > 0.05; smallest P = 0.064, clumped IVW for cg04827551, β = −0.057); bars denote 95% confidence intervals. **(c)** Colocalization (PPH4) of each CpG mQTL with IgAN susceptibility (GCST90018866) and with the *C1GALT1* cis-eQTL in GTEx v10 whole blood, kidney cortex and OneK1K naive B cells. mQTL–IgAN PPH4 ≤ 0.014 for all five CpGs (red), whereas the mQTL–eQTL comparison assigns cg17994788 a shared-variant posterior in whole blood (PPH4 = 0.62), kidney cortex (0.97) and naive B cells (0.79), and cg16101574 in naive B cells (0.82); the remaining CpG–eQTL comparisons are dominated by PPH3 (distinct signals). Dashed line: PPH4 = 0.05. **(d)** Stepwise mediation framework for the germline pathway (lead → methylation → expression → Gd-IgA1 → IgAN). Steps ①–③ are supported (five cis-mQTL CpGs; PPH4 = 0.79–0.82 in naive B cells at step ②; positive lead → expression cis-eQTL at step ③); the mediation proportion at step ④ is not estimable because the Han-Chinese Gd-IgA1 release is P-value-only, and step ⑤ (CpG → IgAN) is null. mQTL, methylation quantitative trait locus; PPH, posterior probability of hypothesis (coloc); MR, Mendelian randomization; TSS, transcription start site.""",
"""**Fig. 5 | The axis shows the strongest and most specific genetic signal in East Asian populations.** **(a)** Effect-allele frequencies of the four published Gd-IgA1 lead variants in gnomAD v4 East Asian (EAS) and non-Finnish European (NFE) populations, with the effect-allele frequency of the cross-ancestry IgAN meta-analysis (GCST90018866) where available. The *C1GALT1* leads are markedly more frequent in EAS than in NFE (rs13226913-C: 0.926 vs 0.420; rs10238682-G: 0.522 vs 0.204), whereas the *GALNT12* lead rs7856182-T shows the opposite pattern (0.025 vs 0.154; red box, "REVERSED"); the X-chromosomal *C1GALT1C1* lead rs5910940-A has comparable frequencies (0.546 vs 0.511) and no meta-analysis frequency because chromosome X is absent from the IgAN meta-analysis (hatched, "chrX NA"). **(b)** Association strength of the four leads for the Han-Chinese Gd-IgA1 phenotype (Wang et al. 2021 [9]) versus IgAN susceptibility in the cross-ancestry meta-analysis (Sakaue et al. 2021 [10]), plotted as −log₁₀(P). The two genome-wide significant Gd-IgA1 leads (rs10238682, rs7856182; red) are strong for the production phenotype but not for disease susceptibility; the triangle marks the X-chromosomal lead rs5910940 (no IgAN association available). Dashed lines mark P = 5 × 10⁻⁸. EAS, East Asian; NFE, non-Finnish European.""",
"""**Fig. 6 | The actionable window lies in acquired epigenetic control of the biosynthesis axis.** **(a)** Intervenability scores of the five axis genes (0–8), decomposed into genetic anchor (G), drug handle (D), functional evidence (F) and risk penalty (−P) according to the prespecified scoring framework (Methods). *C1GALT1C1* receives the highest score (6.0), followed by *C1GALT1* (4.5); *GALNT12*, *GALNT2* and *ST6GALNAC2* each score 2.5. **(b)** Evidence levels (0–3) for candidate interventions per gene across five mechanism classes: nucleoside DNA methyltransferase (DNMT) inhibitors (5-azacytidine class, "5-Aza-CR"), the non-nucleoside DNMT inhibitor (−)-epigallocatechin-3-gallate (EGCG), class I HDAC inhibitors, research-grade HDAC inhibitors and sialyltransferase inhibitors. The only level-3 entry is the reversal of *C1GALT1C1* promoter hypermethylation by 5-azacytidine in an IgAN-relevant B-cell context. **(c)** Direct druggability of the axis members (Pharos tier) versus their intervenability score. Four of the five genes are Tbio/Tdark (no approved drug; red shading) and only *GALNT2* is Tchem (three research-grade ligands, no approved drug); enzyme-active-site inhibition is therefore not a realistic lever. **(d)** Summary of the actionable window: direct enzyme inhibition is not supported by the annotation, whereas the transcriptional/epigenetic layer is reachable by demethylating agents (5-azacytidine class), EGCG and class I HDAC inhibitors, with the reversible *C1GALT1C1* promoter CpG-island window (reversal of IL-4/IL-17-driven hypermethylation by 5-azacytidine [20]) as the most concrete entry point. Scores are ranking aids and make no causal claim. HDAC, histone deacetylase.""",
]
for lg in legends:
    A(lg); A("")
A("_The graphical abstract is a supporting file and carries its own caption in the submission system._")
A("")

ctx_names = {
    "OneK1K_B_naive": "OneK1K B naive", "OneK1K_B_memory": "OneK1K B memory",
    "OneK1K_B_intermediate": "OneK1K B intermediate", "GTEx_v10_blood": "GTEx v10 whole blood",
    "GTEx_v10_kidney_cortex": "GTEx v10 kidney cortex",
}

# ------------------------------------------------------------------ Table 2
A("## Table 2. Differential expression of the five axis genes across disease transcriptome datasets")
A("")
A("| Dataset (tissue, cases/controls) | Gene | logFC | P |")
A("|---|---|---|---|")
t2_sets = [
    ("GSE73953 (PBMC, 15/2 pooled)", ROOT("/阶段1/M1_GEO表达锚点/GSE73953_IgAN_vs_HC_best.csv")),
    ("GSE115857 (kidney bulk, 55/7)", ROOT("/阶段1/M1_GEO表达锚点/GSE115857_IgAN_vs_Ctrl_LD_best.csv")),
    ("GSE93798 (glomeruli, 20/22)", ROOT("/阶段1/M1_GEO表达锚点/GSE93798_IgAN_vs_Ctrl_best.csv")),
]
for label, p in t2_sets:
    df = pd.read_csv(p, sep=",")
    for _, r in df.iterrows():
        A(f"| {label} | *{r['gene']}* | {fmt_fc(r['logFC'])} | {fmt_p(r['pvalue'])} |")
A("")
A("logFC, log2 fold change (IgAN versus control); PBMC, peripheral blood mononuclear cells. P values are nominal from Welch's two-sample t-test on the pre-specified five-gene panel (Methods). The PBMC series pools two healthy controls and supports a directional inference only; the DASL platform of GSE115857 has a narrow dynamic range. Corresponds to R1 and Fig. 1b.")
A("")

# ------------------------------------------------------------------ Table 3
A("## Table 3. cis-eQTL replication of the published Gd-IgA1 lead variants across contexts")
A("")
A("| Lead variant (gene) | Context | cis-eQTL P (axis gene) | β (SE) |")
A("|---|---|---|---|")
m = pd.read_csv(ROOT("/阶段2/主证据1_lead_eQTL/lead_eQTL_matrix.tsv"), sep="\t")
lead_order = {"rs13226913": 0, "rs10238682": 1, "rs7856182": 2, "rs5910940": 3}
ctx_rank = {"OneK1K_B_naive": 0, "OneK1K_B_memory": 1, "OneK1K_B_intermediate": 2,
            "GTEx_v10_blood": 3, "GTEx_v10_kidney_cortex": 4}
m["lead_o"] = m["lead_rsid"].map(lead_order)
m["ctx_o"] = m["dataset"].map(ctx_rank)
mm = m.sort_values(["lead_o", "ctx_o"])
for _, r in mm.iterrows():
    gene = r["lead_gene"]
    tested = str(r["axis_gene_tested"]).strip().lower() == "yes"
    if tested:
        A(f"| {r['lead_rsid']} (*{gene}*) | {ctx_names[r['dataset']]} | {fmt_p(float(r['axis_min_p']))} | {r['axis_beta']:+.3f} ({r['axis_se']:.3f}) |")
    else:
        if gene == "C1GALT1C1":
            note = "Not testable (X chromosome absent from eQTL Catalogue r8)"
        else:
            note = "Not testable (gene filtered from B-lineage eQTL by low expression)"
        A(f"| {r['lead_rsid']} (*{gene}*) | {ctx_names[r['dataset']]} | — | {note} |")
A("")
A("β is the per-allele cis-eQTL effect on the axis gene (SE in parentheses); P is the nominal p_beta (permutation p_perm where available). Corresponds to R2 and Fig. 2b. The full 5-gene × 5-context matrix, including *GALNT12*/*ST6GALNAC2* bulk-tissue rows, is given in Additional file 1: Table S2.")
A("")

# ------------------------------------------------------------------ Table 4
A("## Table 4. Colocalization of axis cis-eQTLs with IgAN susceptibility")
A("")
A("| Gene | Context | n_SNP | PPH1 | PPH4 |")
A("|---|---|---|---|---|")
c = pd.read_csv(ROOT("/阶段2/主证据2_共享结构/coloc_5genes_x_contexts.tsv"), sep="\t")
ok = c[(c["n_SNP"].astype(str).str.isdigit()) & (c["n_SNP"].astype(int) > 0)]
gene_order = ["C1GALT1", "GALNT2", "GALNT12", "ST6GALNAC2"]
for g in gene_order:
    for ck in ["GTEx_v10_blood", "GTEx_v10_kidney_cortex", "OneK1K_B_naive", "OneK1K_B_memory", "OneK1K_B_intermediate"]:
        row = ok[(ok["gene"] == g) & (ok["context"] == ck)]
        if not row.empty:
            r = row.iloc[0]
            A(f"| *{g}* | {ctx_names[ck]} | {int(r['n_SNP'])} | {r['PPH1']:.3f} | {r['PPH4']:.3f} |")
A("")
A("Posteriors from coloc ABF (p1 = p2 = 1 × 10⁻⁴, p12 = 1 × 10⁻⁵) for the 14 computable gene-by-context comparisons; 11 of the 25 planned comparisons are not computable (*C1GALT1C1* in all contexts because the IgAN meta-analysis lacks chromosome X; *GALNT12* and *ST6GALNAC2* in the three B-cell contexts because their cis-eQTL regions are not catalogued). No PPH4 reaches 0.5; PPH4 is the posterior of a shared causal variant (colocalization) and PPH1 of an eQTL signal without a disease signal. SuSiE signal-level re-testing of the four comparisons with official credible sets (seven signal-level tests) reproduced the null (maximum PPH4 = 0.031; Additional file 1: Table S2). An unweighted allelic score of the three autosomal leads against IgAN was null (β = 0.0039, SE = 0.025, P = 0.876; reported as a direction-consistency summary, not an instrument-variable analysis). Corresponds to R3 and Fig. 3b, c.")
A("")

# ------------------------------------------------------------------ Table 5
A("## Table 5. Significant cis-mQTL CpGs of C1GALT1: position, cis-mQTL, Mendelian randomization and colocalization")
A("")
sig = pd.read_csv(ROOT("/阶段3/M4_中介/sigCpG_pos.tsv"), sep="\t")
dirc = pd.read_csv(ROOT("/阶段3/M4_中介/M4_方向一致性_lead.tsv"), sep="\t")
mr = pd.read_csv(ROOT("/阶段3/M4_中介/M4_MR_CpG_x_IgAN.tsv"), sep="\t")
coloc1 = pd.read_csv(ROOT("/阶段3/M4_中介/M4_coloc_mQTL_x_IgAN.tsv"), sep="\t")
coloc2 = pd.read_csv(ROOT("/阶段3/M4_中介/M4_coloc_mQTL_x_eQTL.tsv"), sep="\t")
TSS = 7196565
A("| CpG | hg19 position | kb from TSS | Region | cis-mQTL P | mQTL β (rs13226913)† | MR β (P)‡ | Coloc PPH4 (mQTL × IgAN) | Coloc PPH4 (mQTL × B-naive eQTL) |")
A("|---|---|---|---|---|---|---|---|---|")
for _, r in sig.sort_values("pos37").iterrows():
    cpg = r["cpg"]
    kb = (r["pos37"] - TSS) / 1000.0
    region = "Gene body" if r["pos37"] <= 7288282 else "3′ flanking"
    dr = dirc[(dirc["lead"] == "rs13226913") & (dirc["cpg"] == cpg) & (dirc["context"] == "GTEx_v10_blood")]
    b13 = float(dr.iloc[0]["mqtl_beta_on_eQTLalt"]) if not dr.empty else float("nan")
    P13 = fmt_p(float(dr.iloc[0]["mqtl_P"])) if not dr.empty else "—"
    sub = dirc[dirc["cpg"] == cpg]
    minp = float(sub["mqtl_P"].min())
    # smallest-P MR estimate among: Wald_lead(rs13226913) / Wald_topcis(first) / any IVW_clumped*
    s = mr[mr["cpg"] == cpg]
    cands = []
    t = s[(s["method"] == "Wald_lead") & (s["instrument"] == "rs13226913")]
    if not t.empty:
        cands.append(t.iloc[0])
    t = s[s["method"] == "Wald_topcis"]
    if not t.empty:
        cands.append(t.iloc[0])
    t = s[s["method"].str.startswith("IVW_clumped")]
    t = t[t["mr_p"].notna()]
    if not t.empty:
        cands.append(t.iloc[0])
    if cands:
        best = min(cands, key=lambda x: x["mr_p"])
        mrlab = f"{best['mr_beta']:+.3f} ({fmt_p(float(best['mr_p']))})"
    else:
        mrlab = "—"
    p4_igan = float(coloc1[coloc1["cpg"] == cpg].iloc[0]["PPH4"])
    p4_nb = float(coloc2[(coloc2["cpg"] == cpg) & (coloc2["context"] == "OneK1K_B_naive")].iloc[0]["PPH4"])
    A(f"| {cpg} | {int(r['pos37'])} | +{kb:.1f} | {region} | {fmt_p(minp)} | {fmt_fc(b13)} ({P13}) | {mrlab} | {p4_igan:.3f} | {p4_nb:.3f} |")
A("")
A("cis-mQTL P is the minimum whole-blood mQTL P across the two *C1GALT1* leads (GoDMC). †Per-allele methylation effect of rs13226913 with the sign oriented to the eQTL-raising allele (GoDMC whole blood); P in parentheses. ‡Mendelian-randomization estimate of CpG methylation on IgAN risk with the smallest P per CpG (Wald ratio on rs13226913, Wald ratio on the top cis-sentinel variant, or clumped inverse-variance weighting; method in Additional file 1: Table S4); P in parentheses. In the mQTL–eQTL colocalization, PPH4 for cg17994788 was 0.62 (whole blood) and 0.97 (kidney cortex); all other CpG–eQTL comparisons were dominated by PPH3 (distinct causal variants). Corresponds to R4 and Fig. 4a–c.")
A("")

# ------------------------------------------------------------------ Table 6
A("## Table 6. Allele frequencies and association P values of the Gd-IgA1 lead variants")
A("")
f5 = pd.read_csv(ROOT("/阶段3/M5_东亚/M5_gnomAD_EASEUR_频率_20260903.tsv"), sep="\t")
A("| Lead variant (gene) | Effect allele | gnomAD v4 EAS AF | gnomAD v4 NFE AF | IgAN meta EAF | Han Gd-IgA1 P | IgAN meta P |")
A("|---|---|---|---|---|---|---|")
for _, r in f5.iterrows():
    eaf = f"{r['IgAN_meta_EAF']:.3f}" if pd.notna(r["IgAN_meta_EAF"]) else "NA (chrX)"
    metaP = fmt_p(float(r["IgAN_meta_P"])) if pd.notna(r["IgAN_meta_P"]) else "NA (chrX)"
    A(f"| {r['lead_rsid']} (*{r['lead_gene']}*) | {r['alt']} | {r['gnomAD_EAS_AF']:.3f} | {r['gnomAD_EUR_NFE_AF']:.3f} | {eaf} | {fmt_p(float(r['Wang2021_Han_P']))} | {metaP} |")
A("")
A("AF, allele frequency; EAS, East Asian; NFE, non-Finnish European; EAF, effect-allele frequency. Gd-IgA1 P values are from the Han-Chinese Gd-IgA1 GWAS (Wang et al. 2021 [9]); IgAN meta P values are from the cross-ancestry meta-analysis (Sakaue et al. 2021 [10]; GCST90018866), which does not include chromosome X (rs5910940). Corresponds to R6 and Fig. 5a.")
A("")

# ------------------------------------------------------------------ Table 7
A("## Table 7. Intervenability scoring of the five axis genes")
A("")
s7 = pd.read_csv(ROOT("/阶段3/M6_干预/M6_可干预性评分.tsv"), sep="\t")
pharos = {"C1GALT1": "Tbio (0 approved)", "C1GALT1C1": "Tbio (0 approved)",
          "GALNT12": "Tdark (0 ligands)", "GALNT2": "Tchem (3 research ligands, 0 approved)",
          "ST6GALNAC2": "Tbio (0 approved)"}
A("| Gene | Pharos tier | G (0–3) | D (0–3) | F (0–2) | −P (0–2) | Total (0–8) |")
A("|---|---|---|---|---|---|---|")
for _, r in s7.iterrows():
    pv = float(r['评分_风险扣分P(0-2)'])
    Pcol = f"−{pv:.1f}" if pv > 0 else "0.0"
    A(f"| *{r['基因']}* | {pharos[r['基因']]} | {r['评分_遗传锚点G(0-3)']:.1f} | {r['评分_药物手柄D(0-3)']:.1f} | {r['评分_功能证据F(0-2)']:.1f} | {Pcol} | {r['可干预性总分(0-8)']:.1f} |")
A("")
A("Intervention score = G + D + F − P, where G anchors the genetic/regulatory evidence to the Gd-IgA1 production phenotype, D the availability of an epigenetic drug handle, F the functional and in vivo pharmacological evidence in IgAN-relevant cells or models, and P a penalty for specificity, safety and off-target risks (criteria prespecified in Methods). Scores rank and visualise candidate intervention windows and make no causal claim. Corresponds to R7 and Fig. 6a.")
A("")

with open(OUT, "w", encoding="utf-8") as f:
    f.write("\n".join(lines))
print("WROTE", OUT)
