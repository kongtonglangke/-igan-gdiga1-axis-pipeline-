# -*- coding: utf-8 -*-
# _tmp_patch_w2_full_audit_resume.py — 2026-09-11 全维度审查修复批次（续跑）
# 前段 _tmp_patch_w_full_audit_text.py 已完成 01/03/04/05；本脚本只覆盖剩余文件：
#   06（CRLF 行尾，NK 缩写行）、08（图注 5 处，含 Fig 1a 补漏 + 斜体感知的 Fig 6b/6d）、
#   Additional_File_1（11 处，斜体 *et al.* / *C1GALT1* 感知）、cover_letter（3 处）、
#   09 重编号 + 新 [24] decoupleR / [29] Mathur、Table_S6b 药物归属 P0 修正（2 行 × 2 列）。
# 逐条断言（count==1），失败即停，不产生半写状态。
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os

BASE = _paths.at("阶段4.5_复现包/submission")
M = BASE("manuscript")
A1 = BASE("additional_file_1", "Additional_File_1_Figure_Legends.md")

EDITS = {
# ---- 06：CRLF 文件，old/new 均用 \r\n ----
os.path.join(M, "06_Conclusions_Abbreviations.md"): [
 ("| NFE | non-Finnish European |\r\n| PBMC |",
  "| NFE | non-Finnish European |\r\n| NK | natural killer |\r\n| PBMC |"),
],
# ---- 08：图注（LF 文件，单行编辑）----
os.path.join(M, "08_Figure_Legends_Tables.md"): [
 # Fig 1a（补漏：与 Fig 6d 同源的药物归属错误）
 ("(*C1GALT1C1* promoter, 5-azacytidine-reversible hypermethylation)",
  "(*C1GALT1C1* promoter, DNMT-inhibitor–reversible hypermethylation)"),
 # Fig 3b PPH1 措辞（Table 4 实数：强 0.79–0.91 / 弱 0.177–0.289）
 ("in the tissue eQTL comparisons the posterior mass is concentrated on PPH1 (an eQTL signal without a disease signal), which reaches 0.79–0.91 in the whole-blood and kidney-cortex contexts, whereas in several B-cell contexts the posterior is more evenly distributed (PPH1 0.18–0.28);",
  "the posterior mass is concentrated on PPH1 (an eQTL signal without a disease signal) in the comparisons with the strongest cis-eQTL evidence (PPH1 = 0.79–0.91), whereas comparisons with weaker eQTL evidence distribute the posterior more evenly (PPH1 = 0.18–0.29);"),
 # Fig 3c 七次信号级检验（S2b：4 比较、n_cs 1+2+1+3=7）
 ("Signal-level re-testing of the four combined SuSiE credible sets (region-level coloc ABF PPH4, blue, versus SuSiE signal-level PPH4, orange).",
  "Signal-level re-testing across the four comparisons with official SuSiE credible sets (seven signal-level tests; region-level coloc ABF PPH4, blue, versus SuSiE signal-level PPH4, orange)."),
 # Fig 6b P0（文件内 C1GALT1C1 带斜体星号）
 ("The only level-3 entry is the reversal of *C1GALT1C1* promoter hypermethylation by 5-azacytidine in an IgAN-relevant B-cell context.",
  "The only level-3 entry is the reversal of *C1GALT1C1* promoter hypermethylation by the demethylating agent 5-aza-2′-deoxycytidine in an IgAN-relevant B-cell context [20]; in an IL-17 context, 5-azacytidine likewise restored C1GALT1 and C1GALT1C1 expression [19]."),
 # Fig 6d P0（带斜体星号）
 ("(IL-4-driven hypermethylation of the *C1GALT1C1* promoter, reversed by 5-azacytidine [20])",
  "(IL-4-driven hypermethylation of the *C1GALT1C1* promoter, reversed by 5-aza-2′-deoxycytidine [20])"),
],
# ---- Additional file 1（CRLF 文件，单行编辑）----
A1: [
 # Fig S4：26 donors 中 9 为健康对照（原“26 IgAN donors”错误）；[27]→[28]；去化学版本号
 ("Public 10x Genomics 5′ v1 PBMC scRNA-seq data of 26 IgAN donors (6 late-stage, 11 early-stage, 9 healthy controls; 282,463 single cells; Kim *et al.* [27]; GEO: GSE285335)",
  "Public 10x Genomics 5′ gene-expression PBMC scRNA-seq data of 26 donors (17 with IgAN: 6 late-stage, 11 early-stage; 9 healthy controls; 282,463 single cells; Kim *et al.* [28]; GEO: GSE285335)"),
 # Fig S6：103 units 口径（26 donors × 4 B-cell states − 1）
 ("aggregated to donor × B-cell-state pseudobulk units (103 units) and related",
  "aggregated to donor × B-cell-state pseudobulk units (103 units across four B-cell states and 26 donors) and related"),
 # Fig S7：Table 3 实测 β=+0.185 / P=1.16e-4（原 0.186/1.2e-4 精度不符）
 ("rs10238682: β = +0.186, P = 1.2 × 10⁻⁴",
  "rs10238682: β = +0.185, P = 1.16 × 10⁻⁴"),
 # Fig S7：衰减倍数实测 2.1–11.6×（原 two- to five-fold 低估）
 ("attenuated roughly two- to five-fold in antigen-experienced states",
  "attenuated two- to twelve-fold in antigen-experienced states"),
 # Fig S7：p→P 体例统一 ×3
 ("Cochran's Q (p = 0.046 and 0.047; I² ≈ 67%)",
  "Cochran's Q (P = 0.046 and 0.047; I² ≈ 67%)"),
 ("meta-regression (p = 0.018 for both variants)",
  "meta-regression (P = 0.018 for both variants)"),
 ("did not reach significance (p = 0.11)",
  "did not reach significance (P = 0.11)"),
 # 索引 S1：按 S1 实际列（top variant/chr/pos/β/nominal+permutation P/n）重写
 ("| Table S1 | Full 5-gene × 5-context cis-eQTL matrix (OneK1K naive, memory and intermediate B cells; GTEx v10 whole blood and kidney cortex): −log₁₀(P), per-allele β, standard error, allele frequency and the per-cell call. Corresponds to Fig. 2, Table 3 and Additional file 1: Fig. S2. |",
  "| Table S1 | Full 5-gene × 5-context cis-eQTL matrix (OneK1K naive, memory and intermediate B cells; GTEx v10 whole blood and kidney cortex): the top cis-variant per gene–context with chromosome, position, per-allele β, nominal and permutation P and the number of individuals. Corresponds to Fig. 2, Table 3 and Additional file 1: Fig. S2. |"),
 # 索引 S2b：4 行 SuSiE-vs-ABF（可算性判定在 S2a）
 ("| Table S2b | Computability verdict and robustness summary for the 14 computable colocalisation comparisons: credible-set size, ABF posteriors, SuSiE-versus-ABF posterior difference and the robust-no-colocalisation flag. Corresponds to Fig. 3 and Table 4. |",
  "| Table S2b | SuSiE-versus-ABF robustness summary for the four colocalisation comparisons with official credible sets: credible-set size, ABF posteriors, SuSiE-versus-ABF posterior difference and the robust-no-colocalisation flag (the computability verdicts for all 25 planned comparisons are in Table S2a). Corresponds to Fig. 3 and Table 4. |"),
 # 索引 S4a：10 行 top-variant 精确命中查询（全 NOT_FOUND；带斜体 *C1GALT1*）
 ("| Table S4a | Per-lead IgAN lookup: the two *C1GALT1* leads (rs13226913, rs10238682) and their cis-eQTL variants looked up in the cross-ancestry IgAN meta-analysis (GCST90018866), with IgAN β, SE, effect-allele frequency and P by ancestry. Corresponds to the direct disease-level lookup (Results). |",
  "| Table S4a | Per-context top-variant IgAN lookup: the top cis-eQTL variant of each axis gene in each pre-specified primary context (10 gene–context pairs) tested for exact presence in the cross-ancestry IgAN meta-analysis (GCST90018866), with IgAN β, SE, effect-allele frequency and P where the variant is found (no top variant was present in the release; the per-lead readout is reported in the main text and the window-level minima in Table S4b). Corresponds to the direct disease-level lookup (Results). |"),
 # 索引 S10：5,700 行 = 950 组合 × 6 目标
 ("| Table S10 | CollecTRI regulator-by-B-cell-state associations with the composite axis score (ρ, nominal P, FDR q, target-matched random-set 95th percentile and empirical P; 950 combinations). Corresponds to Fig. 1c and Additional file 1: Fig. S6. |",
  "| Table S10 | CollecTRI regulator-by-B-cell-state associations with each axis-gene target and the composite axis score (Spearman ρ, nominal P, Benjamini–Hochberg q and number of panel targets; 950 regulator–state combinations × six targets = 5,700 rows). The target-matched random-set null and the empirical P values are summarised in the Fig. S6 legend and in the main text. Corresponds to Fig. 1c and Additional file 1: Fig. S6. |"),
],
# ---- cover letter（CRLF 文件，单行编辑）----
BASE("cover_letter_draft.md"): [
 ("including 5-azacytidine-reversible promoter hypermethylation",
  "including DNMT-inhibitor–reversible promoter hypermethylation"),
 ("**naive B-cell-specific cis-eQTL**",
  "**naive B-cell-enriched cis-eQTL**"),
 ("Mendelian randomisation provides",
  "Mendelian randomization provides"),
],
# ---- Table S6b 药物目录：P0 药物归属修正（TSV 纯文本替换）----
BASE("additional_file_1", "Table_S6b_drug_directory.tsv"): [
 # 5-azacytidine 行：IL-4 证据实为 decitabine（Sun 2015），本药仅 IL-17（Lin 2018，DAKIKI）
 ("Reverses IL-4/IL-17-induced IgA1 hypogalactosylation in vitro (IgAN B cells/cell lines), restoring C1GALT1C1 expression — direct functional evidence in an IgAN context (citation via the T-cell review)",
  "Reverses IL-17-induced IgA1 hypogalactosylation in vitro (DAKIKI IgA1-producing cells), restoring C1GALT1 and C1GALT1C1 expression — direct functional evidence in an IgAN context; the parallel IL-4-driven methylation-reversal evidence was obtained with decitabine (see next row), not with 5-azacytidine"),
 ("Clin Exp Nephrol 2019 (review of T cells in IgAN) and the IL-4/IL-17 primary studies cited therein",
  "Lin 2018 (Ren Fail 40(1):60–67; IL-17 with 5-azacytidine in DAKIKI cells); Clin Exp Nephrol 2019 review of T cells in IgAN"),
 # Decitabine 行：直接 IgAN B 细胞甲基化证据在本药（原称“证据更弱”错误）
 ("Same as 5-aza-CR (same mechanistic family); direct IgAN B-cell evidence weaker than the 5-aza-CR literature",
  "Reverses IL-4-driven C1GALT1C1 (Cosmc) promoter hypermethylation in IgAN-patient B cells, restoring chaperone expression — the primary methylation-readout evidence for the C1GALT1C1 window in an IgAN context (Sun 2015)"),
 ("Public pharmacology records; the IgAN application is an extrapolated hypothesis",
  "Sun 2015 (PLoS One 10(2):e0112305; IL-4 with 5-aza-2′-deoxycytidine in paediatric IgAN B cells)"),
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
    print("ALL RESUME EDITS APPLIED")
