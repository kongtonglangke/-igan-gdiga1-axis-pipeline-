#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""_tmp_patch_batch3.py — L7 零分布 200→1000 的数字同步 + cover letter + data availability"""
import os
import sys

PKG = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MS = os.path.join(PKG, "submission", "manuscript")
SUB = os.path.join(PKG, "submission")
AF = os.path.join(SUB, "additional_file_1")


def patch(path, pairs):
    with open(path, encoding="utf-8", newline="") as fh:
        t = fh.read()
    for old, new, exp in pairs:
        n = t.count(old)
        if n != exp:
            raise SystemExit(f"[FAIL] {path}: found {n}, expected {exp}\n  >> {old[:140]}")
        t = t.replace(old, new)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write(t)
    print(f"OK  {os.path.basename(path)}  ({len(pairs)} pairs)")


# ---------- ch3 ----------
patch(os.path.join(MS, "03_Methods.md"), [
    ("significance was assessed against an empirical null of 200 random gene sets of matched target number drawn from the expressed panel",
     "significance was assessed against an empirical null of 1,000 random gene sets of matched target number drawn from the expressed panel", 1),
])

# ---------- ch4 ----------
patch(os.path.join(MS, "04_Results.md"), [
    ("significance was judged against an empirical null distribution of 200 target-number-matched random sets rather than by \u03c1 or by false-discovery rate alone",
     "significance was judged against an empirical null distribution of 1,000 target-number-matched random sets rather than by \u03c1 or by false-discovery rate alone", 1),
    ("reached Spearman \u03c1 of 0.73\u20130.78 against the axis score",
     "reached a median Spearman \u03c1 of 0.73\u20130.78 against the axis score (95th percentile 0.84\u20130.88)", 1),
    ("the DNA methyltransferase DNMT1 (plasma cells; \u03c1 = 0.85, empirical P = 0.010), and the CpG-centred Sp/KLF family (KLF2 and SP1 in naive B cells; empirical P = 0.005 and 0.015, respectively).",
     "the DNA methyltransferase DNMT1 (plasma cells; \u03c1 = 0.85, empirical P = 0.025), and the CpG-centred Sp/KLF family (KLF2 and SP1 in naive B cells; empirical P = 0.006 and 0.019, respectively).", 1),
    ("left nine factors\u2014KLF2, SP1, EGR1, HIF1A, MAZ, NFIC, SMAD3, ETV6 and TFAP2A\u2014whose motifs",
     "left nine factors\u2014KLF2, SP1, EGR1, ETV6, HIF1A, MAZ, MYC, NFIC and SMAD3\u2014whose motifs", 1),
    ("and the two single-cell B-cell partitions (3.2% versus 4.2% C1GALT1C1 detection) are now explicitly distinguished at the point of use.",
     "and the two single-cell B-cell partitions (3.2% versus 4.2% C1GALT1C1 detection) are now explicitly distinguished at the point of use. The L7 empirical null was re-run at 1,000 target-matched random sets (was 200), which resolves the empirical P values off the 1/200 floor; the pass set is unchanged (36 of 950 combinations) but the intersected factor list is now KLF2, SP1, EGR1, ETV6, HIF1A, MAZ, MYC, NFIC and SMAD3 (TFAP2A no longer passes at the higher resolution).", 1),
])

# ---------- ch8 ----------
patch(os.path.join(MS, "08_Figure_Legends_Tables.md"), [
    ("> v1.6 (2026-09-10, reviewer-perspective clean-up)",
     "> v1.7 (2026-09-10, L7 null re-run at 1,000 sets): Fig. S6 legend/summary index now reflect the higher-resolution empirical null; legends that quote the L7 numbers were checked for stale values.\n> v1.6 (2026-09-10, reviewer-perspective clean-up)", 1),
])

# ---------- Additional_File_1 legends ----------
patch(os.path.join(AF, "Additional_File_1_Figure_Legends.md"), [
    ("exceed the empirical null (200 target-number-matched random gene sets; empirical P < 0.05)",
     "exceed the empirical null (1,000 target-number-matched random gene sets; empirical P < 0.05)", 1),
])

# ---------- cover letter ----------
patch(os.path.join(SUB, "cover_letter_draft.md"), [
    ("> **\u72b6\u6001**\uff1adraft v1.0\uff082026-09-06\uff09\u2014 \u5f85\u4f5c\u8005\u7ec8\u6821\u540e\u5b9a\u7a3f\u7b7e\u540d\u3002\u65b9\u62ec\u53f7 `[ ]` \u4e3a\u5f85\u586b\u9879\u3002",
     "> **\u72b6\u6001**\uff1adraft v1.1\uff082026-09-10\uff09\u2014 \u5df2\u6309\u5ba1\u7a3f\u89c6\u89d2\u6e05\u6d17\uff08\u4e03\u5c42\u8bc1\u636e\u5c42\uff1b25\u219214 \u53e3\u5f84\uff1b\u8865 GSE285335 \u4e0e ImmuNexUT\uff1bcausal-inference \u63aa\u8f9e\uff09\uff0c\u5f85\u4f5c\u8005\u7ec8\u6821\u540e\u5b9a\u7a3f\u7b7e\u540d\u3002\u65b9\u62ec\u53f7 `[ ]` \u4e3a\u5f85\u586b\u9879\u3002", 1),
    ("whole-blood mQTLs (GoDMC), and population allele frequencies (gnomAD v4), and interrogate the five-gene biosynthesis axis (*C1GALT1, C1GALT1C1, GALNT2, GALNT12, ST6GALNAC2*) across six sequential evidence layers.",
     "whole-blood mQTLs (GoDMC), population allele frequencies (gnomAD v4), IgAN peripheral-blood single-cell transcriptomes (GSE285335), and an independent East Asian immune-cell-type eQTL panel (ImmuNexUT), and interrogate the five-gene biosynthesis axis (*C1GALT1, C1GALT1C1, GALNT2, GALNT12, ST6GALNAC2*) across seven sequential evidence layers.", 1),
    ("(iii) across 25 gene-by-context colocalisation tests, the genetic signal of the axis shows **no shared causal structure with IgAN susceptibility** (maximum PPH4 = 0.031), whereas the same leads reach genome-wide significance for Gd-IgA1 quantity in Han Chinese patients;",
     "(iii) of the 25 planned gene-by-context colocalisation tests, 14 were computable and none supported a **shared causal variant with IgAN susceptibility** (maximum PPH4 = 0.031), whereas the same leads reach genome-wide significance for Gd-IgA1 quantity in Han Chinese patients;", 1),
    ("and (vi) the actionable window maps to **reversible acquired epigenetic control** (*C1GALT1C1*, including 5-azacytidine-reversible promoter hypermethylation), rather than to direct enzyme inhibition of the axis.",
     "and (vi) the actionable window maps to **reversible acquired epigenetic control** (*C1GALT1C1*, including 5-azacytidine-reversible promoter hypermethylation), rather than to direct enzyme inhibition of the axis; and (vii) single-cell regulon and motif inference places that window inside a concrete upstream circuit, in which cytokine signalling (IL-4/IL-17 via STAT6/STAT3), DNMT1 and methylation-sensitive Sp/KLF occupancy converge on the *C1GALT1C1* promoter.", 1),
    ("it resolves, within one causal framework, the cell type and the ancestry context",
     "it resolves, within one causal-inference framework, the cell type and the ancestry context", 1),
    ("*Cover letter draft v1.0 \u00b7 \u9636\u6bb54.5 \u590d\u73b0\u5305 \u00b7 \u6295\u7a3f\u524d\u8bf7\u5168\u4f53\u4f5c\u8005\u6821\u9605\u5e76\u786e\u8ba4\u901a\u8baf\u4f5c\u8005\u7f72\u540d\u4e0e\u8054\u7cfb\u65b9\u5f0f*",
     "*Cover letter draft v1.1 (2026-09-10) \u00b7 \u9636\u6bb54.5 \u590d\u73b0\u5305 \u00b7 \u6295\u7a3f\u524d\u8bf7\u5168\u4f53\u4f5c\u8005\u6821\u9605\u5e76\u786e\u8ba4\u901a\u8baf\u4f5c\u8005\u7f72\u540d\u4e0e\u8054\u7cfb\u65b9\u5f0f*", 1),
])

# ---------- data availability ----------
patch(os.path.join(SUB, "data_availability_draft.md"), [
    ("- GSE93798 \u2014 microdissected glomeruli transcriptome (20/22), Affymetrix (GPL22945).",
     "- GSE93798 \u2014 microdissected glomeruli transcriptome (20/22), Affymetrix (GPL22945).\n- GSE285335 \u2014 IgAN PBMC single-cell RNA-seq (26 donors: 6 late-stage and 11 early-stage IgAN patients and 9 healthy controls; 10x Genomics 5\u2032 v1; 282,463 cells), used for the cell-level cross-check of the B-cell-state signal and for the single-cell atlas and upstream-regulator inference.", 1),
    ("**GWAS summary statistics (GWAS Catalog):**",
     "**East Asian regulatory cross-reference (NBDC Human Database, Japan):**\n- ImmuNexUT (E-GEAD-398) \u2014 immune-cell-type-specific conditional cis-eQTL summary statistics (416 individuals; 28 purified peripheral-blood subsets; GRCh38; FDR < 0.05 conditional signals).\n\n**GWAS summary statistics (GWAS Catalog):**", 1),
    ("- Pharos / IDG druggability annotations (pharos.nih.gov; Tclin/Tchem/Tbio/Tdark tiers).",
     "- Pharos / IDG druggability annotations (pharos.nih.gov; Tclin/Tchem/Tbio/Tdark tiers).\n- CollecTRI regulon resource (curated transcription-factor\u2013target network; saezlab/CollecTRI) \u2014 used for the single-cell upstream-regulator inference.\n- JASPAR 2024 CORE non-redundant vertebrate motif matrices (jaspar.genereg.net) \u2014 used for the methylation-sensitive motif analysis.", 1),
    ("*v1.0 (2026-09-06) \u00b7 reproducibility package \u00b7 consistent line-by-line with ch7 Declarations and ch3 Table 1*",
     "*v1.1 (2026-09-10) \u00b7 reproducibility package \u00b7 four resources that were present in ch7 Declarations but missing from this明细 were added (GSE285335; ImmuNexUT E-GEAD-398; CollecTRI; JASPAR 2024 CORE), so that this file and ch7 Declarations are now consistent line-by-line with each other and with ch3 Table 1.*", 1),
])

# ---------- residual scan ----------
print("\n== residual scan ==")
bad = 0
for dp, dn, fn in os.walk(SUB):
    for f in fn:
        if f.endswith(".md") and "Manuscript_GM_fulltext" not in f:
            p = os.path.join(dp, f)
            t = open(p, encoding="utf-8", newline="").read()
            for pat in ["TFAP2A", "empirical P = 0.010", "200 target", "200 random", "six sequential evidence",
                        "(R1", "(R2", "(R3", "(R4", "(R5", "(R6", "(R7", "(R8", "evidence layer L"]:
                c = t.count(pat)
                if c:
                    print(f"   {pat:34s} x{c}  {os.path.relpath(p, SUB)}")
                    bad += c
print("residual total:", bad)
