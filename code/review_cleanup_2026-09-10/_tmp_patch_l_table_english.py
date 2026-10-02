# -*- coding: utf-8 -*-
"""Deliverable-cleanliness pass for the supplementary tables (submission/additional_file_1).

Fixes:
  1) filenames: drop internal module tokens (_M4_/_M6_) from the five affected part-files;
  2) Table S6a / S6b: translate the Chinese headers and content to English (values unchanged);
  3) Table S2a note column, S4b exact_hit column, S5c feature column: translate the Chinese notes;
  4) Table S4a: drop the internal reference from the column name eQTL_beta(来自M2a) -> eQTL_beta.

No numeric value, comparison, verdict or conclusion is altered by this pass.
"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import io, os

BASE = _paths.at("阶段4.5_复现包/submission/additional_file_1")

RENAMES = {
    "Table_S5b_M4_MR_CpG_x_IgAN.tsv": "Table_S5b_MR_CpG_x_IgAN.tsv",
    "Table_S5c_M4_coloc_mQTL_x_IgAN.tsv": "Table_S5c_coloc_mQTL_x_IgAN.tsv",
    "Table_S5d_M4_coloc_mQTL_x_eQTL.tsv": "Table_S5d_coloc_mQTL_x_eQTL.tsv",
    "Table_S6a_M6_intervenability_scores.tsv": "Table_S6a_intervenability_scores.tsv",
    "Table_S6b_M6_drug_directory.tsv": "Table_S6b_drug_directory.tsv",
}

# ---- Table S6a (headers only; values unchanged) -------------------------------
S6A = """gene\tG(0-3)\tD(0-3)\tF(0-2)\tP(0-2)\ttotal(0-8)
C1GALT1C1\t2.0\t3.0\t2.0\t1.0\t6.0
C1GALT1\t3.0\t1.5\t1.0\t1.0\t4.5
GALNT12\t2.0\t0.5\t0.0\t0.0\t2.5
GALNT2\t1.5\t1.0\t0.0\t0.0\t2.5
ST6GALNAC2\t1.5\t1.0\t0.0\t0.0\t2.5
"""

# ---- Table S6b (full English translation; values unchanged) -------------------
S6B = """Drug/molecule\tClass\tActual target enzyme (epigenetic)\tEpiFactors role\tChemical type / toxicity\tApproval status (public record)\tB/plasma-cell-relevant evidence\tPotential use in the IgAN axis (hypothesis)\tEvidence source\tRisks / notes
EGCG (epigallocatechin-3-gallate)\tNatural polyphenol (green tea)\tDNMT1/DNMT3B (direct binding inhibition, Ki ~ 6.9 uM, non-nucleoside); HDAC1/2 (reversible inhibition); activity against HATs also reported\tUpstream inhibition of the writer (DNMT) and eraser (HDAC) arms; a natural epigenetic modulator\tNon-nucleoside; non-cytotoxic mechanism\tNot marketed as a prescription drug (topical preparations Veregen/Polyphenon E for warts; dietary/research grade)\tImmunomodulatory effect on B lymphocytes (review level); reactivates methylation-silenced genes in vitro\tCandidate natural "transcriptional restoration" molecule: hypothesised to up-regulate B/plasma-cell C1GALT1/C1GALT1C1 transcription -> lower Gd-IgA1; the most safety-friendly preclinical handle\tFang 2003 (Ki 6.89 uM); Negri 2018 Nutrients review; HeLa: direct binding to DNMT1/DNMT3B/HDAC1; Modern Immunology 2018 review\tConcentration-dependent: pro-oxidant/cytotoxic at high concentration; low oral bioavailability; multi-target promiscuity
Azacitidine (5-azacytidine, 5-aza-CR; Vidaza)\tNucleoside DNMT inhibitor\tDNMT1 (incorporated into RNA/DNA, then trapped by DNMT -> degradation -> passive demethylation)\twriter inhibitor\tNucleoside (cytotoxic background)\tFDA approved 2004 (MDS and other myeloid neoplasms)\tReverses IL-4/IL-17-induced IgA1 hypogalactosylation in vitro (IgAN B cells/cell lines), restoring C1GALT1C1 expression - direct functional evidence in an IgAN context (citation via the T-cell review)\tMechanistic probe / proof-of-concept drug: validates the causal chain "demethylation restores C1GALT1C1/C1GALT1 -> lowers Gd-IgA1"\tClin Exp Nephrol 2019 (review of T cells in IgAN) and the IL-4/IL-17 primary studies cited therein\tSystemic myelosuppression/hepatotoxicity/teratogenicity; unsuitable as a direct therapeutic; for mechanistic validation of the window only
Decitabine (5-aza-2'-deoxycytidine, 5-aza-CdR; Dacogen)\tNucleoside DNMT inhibitor\tDNMT1 (DNA incorporation, trapping and depletion of DNMT)\twriter inhibitor\tNucleoside (cytotoxic background)\tFDA approved 2006 (MDS/AML)\tSame as 5-aza-CR (same mechanistic family); direct IgAN B-cell evidence weaker than the 5-aza-CR literature\tAs for 5-aza-CR; an alternative tool compound for the demethylation window\tPublic pharmacology records; the IgAN application is an extrapolated hypothesis\tSame nucleoside-class toxicity; long-term demethylation carries a secondary-malignancy risk
RG108\tNon-nucleoside DNMT1 inhibitor (research use)\tDNMT1 (direct active-site binding, non-incorporating)\twriter inhibitor\tNon-nucleoside\tNot approved (laboratory tool compound)\tNo direct IgAN/B-cell evidence\tAlternative low-toxicity demethylation probe (if nucleoside-class toxicity must be avoided)\tPublic pharmacology records\tWeak in vitro activity; limited cellular uptake
Vorinostat (SAHA; Zolinza)\tHydroxamic-acid HDAC inhibitor (broad-spectrum class I/II)\tHDAC1/2/3/6, etc. (Zn2+-dependent)\teraser inhibitor\tSmall molecule (broad-spectrum epigenetic)\tFDA approved 2006 (cutaneous T-cell lymphoma, CTCL)\tMouse: impairs the primary antibody response while preserving memory B cells (Waibel 2015 Nat Commun) - supports that the B-cell response window is modifiable\tIf C1GALT1C1/C1GALT1 silencing has a deacetylation component, class I HDACi may reactivate it (hypothesis; requires functional testing); also represents the "plasma-cell/antibody-output modulation" direction\tWaibel et al. Nat Commun 6:7838 (2015)\tMyelosuppression/diarrhoea/cardiac QT risk; broad-spectrum epigenetic off-target effects
Romidepsin (Istodax)\tCyclic-peptide HDAC inhibitor (class I selective: HDAC1/2)\tHDAC1/HDAC2 (class I)\teraser inhibitor\tSmall molecule (class I selective)\tFDA approved 2009 CTCL / 2011 PTCL\tB/plasma-cell neoplasms are sensitive; no direct IgAN evidence\tclass I (HDAC1/2) selectivity - preferable if the axis epigenetic window is mediated by HDAC1/2 (hypothesis)\tPublic pharmacology records; MM review (HDACi in plasma-cell neoplasms)\tQT/myelotoxicity; intravenous administration
Panobinostat (LBH589; Farydak)\tHydroxamic-acid pan-HDAC inhibitor (class I/II/IV, greater class I potency)\tHDAC1/2/3/6, etc.\teraser inhibitor\tSmall molecule (broad-spectrum)\tFDA approved 2015 (multiple myeloma, in combination with bortezomib and dexamethasone)\tMarkedly reduces autoreactive plasma cells, autoantibodies and nephritis (MRL/lpr lupus mice, Waibel 2015); clinically approved for multiple myeloma (a plasma-cell neoplasm)\tThe strongest animal-level concept evidence for the "plasma-cell output/autoantibody reduction" direction; in principle lowers the source of Gd-IgA1 production (not the transcriptional-restoration direction)\tWaibel et al. Nat Commun 6:7838 (2015); FDA 2015\tMarked myelosuppression/gastrointestinal toxicity; broad-spectrum off-target effects; does not support direct use in IgAN
Entinostat (MS-275)\tBenzamide HDAC inhibitor (class I selective: HDAC1/2/3)\tHDAC1/2/3 (class I)\teraser inhibitor\tSmall molecule (class I selective, oral)\tNot approved (previously entered multi-phase trials in solid and haematological tumours)\tNo direct IgAN evidence; the class I spectrum fits the C1GALT1C1 epigenetic-window hypothesis\tAlternative class I-selective oral probe (if approval prospects exist)\tPublic pharmacology records\tSafety data from clinical development are incomplete
Chidamide (Epidaza)\tBenzamide HDAC inhibitor (class I selective HDAC1/2/3/10)\tHDAC1/2/3/10 (class I)\teraser inhibitor\tSmall molecule (class I selective, oral)\tChina NMPA approved 2015 (peripheral T-cell lymphoma, PTCL); a representative of domestically accessible epigenetic drugs in China\tNo direct IgAN evidence\tRepresentative of "class I HDACi accessibility" in the Chinese context; hypothesis as above\tNMPA public records\tHaematological toxicity; off-label use requires ethical and regulatory assessment
Belinostat (Beleodaq)\tHydroxamic-acid HDAC inhibitor (class I/II)\tHDAC1/2/3/6, etc.\teraser inhibitor\tSmall molecule (broad-spectrum)\tFDA approved 2014 (PTCL)\tNo direct IgAN evidence\tAdditional same-class option (intravenous); rationale as for vorinostat\tPublic pharmacology records\tSame toxicity spectrum as the hydroxamic acids
Trichostatin A (TSA)\tHydroxamic-acid pan-HDAC inhibitor (tool compound)\tHDAC1/2/3/6, etc. (class I/II)\teraser inhibitor\tSmall molecule (research tool)\tNot approved (classic in vitro tool compound)\tReference for mechanistic studies in the HDACi field (in vitro)\tIn vitro "HDAC window" probe\tPublic records\tUnstable in vivo; in vitro only
"""

# ---- in-place substitutions for the remaining tables --------------------------
SUBS = {
    "Table_S2a_coloc_25comparisons.tsv": [
        ("IgAN GWAS 不含该染色体, coloc 不可行",
         "IgAN GWAS does not include this chromosome; colocalisation not feasible"),
        ("该语境无此基因 cis-eQTL; 且 IgAN GWAS 不含 X 染色体, coloc 物理不可行",
         "no cis-eQTL for this gene in this context; the IgAN GWAS also lacks chromosome X, so colocalisation is physically impossible"),
        ("该语境无此基因 cis-eQTL (region_cache 缺/未收录)",
         "no cis-eQTL for this gene in this context (absent / not catalogued in region_cache)"),
    ],
    "Table_S4a_cisQTL_topvar_IgAN_lookup.tsv": [
        ("eQTL_beta(来自M2a)", "eQTL_beta"),
    ],
    "Table_S4b_cisQTL_topvar_IgAN_locus_lookup.tsv": [
        ("NO(参考面板SNP差异)", "NO (reference-panel SNP difference)"),
    ],
    "Table_S5c_coloc_mQTL_x_IgAN.tsv": [
        ("基因体内(intragenic)", "intragenic"),
        ("3'下游<5kb", "3' downstream < 5 kb"),
    ],
}


def read(p):
    return io.open(p, encoding="utf-8").read()


def write(p, s):
    with io.open(p, "w", encoding="utf-8", newline="") as fh:
        fh.write(s)


def main():
    log = []

    # (1) rename
    for old, new in RENAMES.items():
        op, np_ = BASE(old), BASE(new)
        if os.path.exists(op):
            os.replace(op, np_)
            log.append(f"RENAME {old} -> {new}")
        else:
            log.append(f"SKIP   {old} (absent)")

    # (2) full rewrites
    write(BASE("Table_S6a_intervenability_scores.tsv"), S6A)
    log.append("WRITE  Table_S6a_intervenability_scores.tsv (English headers)")
    write(BASE("Table_S6b_drug_directory.tsv"), S6B)
    log.append("WRITE  Table_S6b_drug_directory.tsv (full English translation)")

    # (3) in-place substitutions
    for name, repls in SUBS.items():
        p = BASE(name)
        if not os.path.exists(p):
            log.append(f"MISS   {name}")
            continue
        s = read(p)
        n = 0
        for a, b in repls:
            c = s.count(a)
            s = s.replace(a, b)
            n += c
        write(p, s)
        log.append(f"SUBS   {name}: {n} replacement(s)")

    print("\n".join(log))


if __name__ == "__main__":
    main()
