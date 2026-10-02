# -*- coding: utf-8 -*-
"""
_tmp_patch_p_round3_fixes.py — 第三轮审稿视角核查（2026-09-10 深夜）的整改补丁

本脚本对 `submission/manuscript/` 的定稿章节 + `submission/cover_letter_draft.md`
执行 10 处修正，每处均内置「命中次数 == 1」断言，命中数不符即抛错中止（防止误改/漏改）。
不做任何数值、图件、结论改动。

修正清单
  F1  ch4 Results（甲基化节）补 Table 5 正文引用 —— 主文表 5 此前未被任何正文句引用，
      BMC/Genome Medicine 要求正文表按序在文内被引用。
  F2  ch3 Methods「three orders of magnitude」→「two orders of magnitude」——
      同一事实（P = 9.4e-6 vs 全基因组阈 5e-8，比值 188 ≈ 2.3 个数量级）在 ch4 写作
      "two orders"，两处口径互相矛盾。
  F3  ch5 Discussion「attenuate along the naive-to-memory transition」→「attenuated in the
      antigen-experienced states (with partial recovery in memory)」——与已修正的 Results
      正文/Fig. 2c 图注（"not a monotonic decline / not a gradient"）冲突。
  F4  ch1 Abstract 首两句悬垂分词修正（"Produced by …, whether …"—分词主语错位）——
      改为同位语结构，词数净变 0（348 不变，限 350）。
  F5  ch5 Discussion 单细胞表达「essentially flat … (P = 2.9e-75)」——效应量与显著性口径
      自相矛盾，改为给出实测均值并说明大细胞数下的小位移亦显著。
  F6  ch4 Results 异质性检验补明所测 lead（rs13226913；Table S11 中两个 lead 各有其值）。
  F7  ch6 缩写表「PPH3/PPH4」→「PPH0–PPH4」，与 Results 正文定义的 PPH0–PPH4 记法一致。
  F8  ch4 Results 首次出现的 Sp/KLF 补全称；ch3 Methods 补 IDG 全称、Pharos 四级全称、
      NBDC 全称（此前均未展开）。
  F9  ch6 缩写表注「gene symbols … are expanded at first appearance」表述与实况不符
      （多数 TF 基因符号未展开），改为准确表述。
  F10 cover letter：七层 → 八层，并补入遗漏的第 8 层「疾病层面直接查证」，其余条目顺延；
      疾病层面阴性证据清单口径写准（coloc + 直接回查 + 等位评分 + CpG 甲基化 MR）。

用法（包根目录）：python reproducibility/code/review_cleanup_2026-09-10/_tmp_patch_p_round3_fixes.py
"""
import os
import sys

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
M = os.path.join(BASE, "submission", "manuscript")
SUB = os.path.join(BASE, "submission")

EDITS = [
    # ---------------- F1: Table 5 cited in the running text ----------------
    ("ch4", os.path.join(M, "04_Results.md"),
     "none lies within the promoter (TSS ± 1.5 kb) (Fig. 4a).",
     "none lies within the promoter (TSS ± 1.5 kb) (Fig. 4a; Table 5)."),

    # ---------------- F2: order-of-magnitude口径统一 ----------------
    ("ch3", os.path.join(M, "03_Methods.md"),
     "a rare-variant window scan that remains approximately three orders of magnitude from genome-wide significance and does not involve any axis lead variant.",
     "a rare-variant window scan that remains approximately two orders of magnitude above the genome-wide significance threshold (P = 9.4 × 10⁻⁶ versus 5 × 10⁻⁸) and does not involve any axis lead variant."),
    ("ch4", os.path.join(M, "04_Results.md"),
     "still approximately two orders of magnitude from significance (complete lookup",
     "still approximately two orders of magnitude above the genome-wide significance threshold (complete lookup"),

    # ---------------- F3: Discussion 与 Fig. 2c/Results 口径一致 ----------------
    ("ch5", os.path.join(M, "05_Discussion_Limitations.md"),
     "act as cis-eQTLs whose effects are concentrated in naive B cells and attenuate along the naive-to-memory transition (Fig. 2)",
     "act as cis-eQTLs whose effects are concentrated in naive B cells and are attenuated in the antigen-experienced states, with partial recovery in memory (Fig. 2)"),

    # ---------------- F4: Abstract 悬垂分词（词数净变 0） ----------------
    ("ch1", os.path.join(M, "01_Title_Abstract_Keywords.md"),
     "Galactose-deficient IgA1 (Gd-IgA1) is the central pathogenic molecule of the multi-hit model of IgA nephropathy (IgAN). Produced by the O-glycosylation machinery of B cells and plasma cells, whether genetic variation acting through this axis governs Gd-IgA1 levels, IgAN susceptibility, or both, and in which cell types, remains unresolved.",
     "Galactose-deficient IgA1 (Gd-IgA1), the central pathogenic molecule of the multi-hit model of IgA nephropathy (IgAN), is produced by the O-glycosylation machinery of B cells and plasma cells; whether genetic variation acting through this axis governs Gd-IgA1 levels, IgAN susceptibility, or both, and in which cell types, remains unresolved."),

    # ---------------- F5: 单细胞表达量口径 ----------------
    ("ch5", os.path.join(M, "05_Discussion_Limitations.md"),
     "whereas the expression level among expressing cells is essentially flat for C1GALT1 and shows a modest decline for C1GALT1C1 (Mann–Whitney P = 2.9 × 10⁻⁷⁵ and P = 1.5 × 10⁻²⁵, respectively;",
     "whereas the expression level among expressing cells changes little in magnitude for C1GALT1 (mean non-zero log1p CPM 1.18 in naive cells versus 1.05 in plasma cells) and declines modestly for C1GALT1C1 (1.02 versus 0.85); with this cell number both displacements are formally significant (Mann–Whitney P = 2.9 × 10⁻⁷⁵ and P = 1.5 × 10⁻²⁵, respectively;"),

    # ---------------- F6: 异质性检验补明 lead ----------------
    ("ch4", os.path.join(M, "04_Results.md"),
     "A formal test on the uniformly recomputed per-subset effects confirmed that the effect is gated by B-cell state (Cochran's Q = 6.15",
     "A formal test on the uniformly recomputed per-subset effects confirmed that the effect of rs13226913 is gated by B-cell state (Cochran's Q = 6.15"),

    # ---------------- F7: 缩写表 PPH 记法 ----------------
    ("ch6", os.path.join(M, "06_Conclusions_Abbreviations.md"),
     "| PPH3/PPH4 | posterior probability of colocalisation hypotheses H3/H4 |",
     "| PPH0–PPH4 | posterior probability of the five colocalisation hypotheses H0–H4 |"),

    # ---------------- F8: 未展开缩写补全称 ----------------
    ("ch4", os.path.join(M, "04_Results.md"),
     "and the CpG-centred Sp/KLF family (KLF2 and SP1 in naive B cells;",
     "and the CpG-centred specificity protein/Krüppel-like factor (Sp/KLF) family (KLF2 and SP1 in naive B cells;"),
    ("ch3", os.path.join(M, "03_Methods.md"),
     "and druggability annotation from Pharos/IDG [15].",
     "and druggability annotation from Pharos/IDG (Illuminating the Druggable Genome) [15]."),
    ("ch3", os.path.join(M, "03_Methods.md"),
     "Direct druggability of the five axis genes was annotated from Pharos/IDG (Tclin/Tchem/Tbio/Tdark) [15].",
     "Direct druggability of the five axis genes was annotated from Pharos/IDG target development levels (Tclin, approved drug target; Tchem, potent small-molecule ligand; Tbio, characterised biology; Tdark, understudied) [15]."),
    ("ch3", os.path.join(M, "03_Methods.md"),
     "were obtained from the NBDC Human Database (accession E-GEAD-398; GRCh38)",
     "were obtained from the National Bioscience Database Center (NBDC) Human Database (accession E-GEAD-398; GRCh38)"),

    # ---------------- F9: 缩写表注表述准确性 ----------------
    ("ch6", os.path.join(M, "06_Conclusions_Abbreviations.md"),
     "*Notes: Gene symbols (C1GALT1, C1GALT1C1, GALNT2, GALNT12, ST6GALNAC2, COSMC) follow HGNC nomenclature and are expanded at first appearance in the text; they are not repeated here.",
     "*Notes: Gene symbols (C1GALT1, C1GALT1C1, GALNT2, GALNT12, ST6GALNAC2, COSMC) follow HGNC nomenclature; the five axis genes are expanded at first appearance in the text and are therefore not repeated here, and all other gene symbols are used in their standard HGNC form."),

    # ---------------- F10: cover letter 七层→八层 ----------------
    ("cover-header", os.path.join(SUB, "cover_letter_draft.md"),
     "> **状态**：draft v1.1（2026-09-10）— 已按审稿视角清洗（七层证据层；25→14 口径；补 GSE285335 与 ImmuNexUT；causal-inference 措辞），待作者终校后定稿签名。方括号 `[ ]` 为待填项。",
     "> **状态**：draft v1.2（2026-09-10）— 已按审稿视角清洗（**八层**证据层，与正文 Abstract/Results/Fig. 1c 一致；25→14 口径；补 GSE285335 与 ImmuNexUT；causal-inference 措辞），待作者终校后定稿签名。方括号 `[ ]` 为待填项。"),
    ("cover-layers", os.path.join(SUB, "cover_letter_draft.md"),
     "across seven sequential evidence layers",
     "across eight sequential evidence layers"),
    ("cover-enum", os.path.join(SUB, "cover_letter_draft.md"),
     "(iv) gene-body CpG methylation marks transcription-coupled regulation of *C1GALT1* **without mediating IgAN risk**; (v) the axis signal is strongest and most specific in East Asian populations; (vi) the actionable window maps to **reversible acquired epigenetic control** (*C1GALT1C1*, including 5-azacytidine-reversible promoter hypermethylation), rather than to direct enzyme inhibition of the axis; and (vii) single-cell regulon and motif inference places that window inside a concrete upstream circuit",
     "(iv) gene-body CpG methylation marks transcription-coupled regulation of *C1GALT1* **without mediating IgAN risk**; (v) a direct variant- and locus-level look-up of the axis leads and cis-eQTL top variants in the IgAN meta-analysis is likewise null; (vi) the axis signal is strongest and most specific in East Asian populations, cross-referenced against an independent Japanese immune-cell eQTL panel; (vii) the actionable window maps to **reversible acquired epigenetic control** (*C1GALT1C1*, including 5-azacytidine-reversible promoter hypermethylation), rather than to direct enzyme inhibition of the axis; and (viii) single-cell regulon and motif inference places that window inside a concrete upstream circuit"),
    ("cover-null", os.path.join(SUB, "cover_letter_draft.md"),
     "the joint null for disease susceptibility across colocalisation, Mendelian randomization and allelic scoring provides a principled negative result",
     "the joint null for disease susceptibility across colocalisation, direct variant look-up, allelic scoring and CpG-methylation Mendelian randomisation provides a principled negative result"),
    ("cover-footer", os.path.join(SUB, "cover_letter_draft.md"),
     "*Cover letter draft v1.1 (2026-09-10) · 阶段4.5 复现包 · 投稿前请全体作者校阅并确认通讯作者署名与联系方式*",
     "*Cover letter draft v1.2 (2026-09-10) · 阶段4.5 复现包 · 投稿前请全体作者校阅并确认通讯作者署名与联系方式*"),
]


def main():
    failures = []
    for tag, path, old, new in EDITS:
        t = open(path, encoding="utf-8").read()
        n = t.count(old)
        if n != 1:
            failures.append((tag, os.path.basename(path), n, old[:70]))
            continue
        open(path, "w", encoding="utf-8").write(t.replace(old, new, 1))
        print(f"  OK   {tag:14s} {os.path.basename(path):38s} 1 hit")
    if failures:
        print("\n!!! FAILED (hit count != 1):")
        for f in failures:
            print("   ", f)
        sys.exit(1)
    print("\nall edits applied.")


if __name__ == "__main__":
    main()
