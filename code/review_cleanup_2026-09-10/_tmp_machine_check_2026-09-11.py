# -*- coding: utf-8 -*-
# _tmp_machine_check_2026-09-11.py — 全维度审查修复后机检（母本 v1.12 + 附件 + cover letter）
# 检查项：Abstract 词数上限、引文首现序 [1..29] 严格递增、正文中文 0、内部编号 0、
# 陈旧串 0（药物误归属/q=0.28/25,095/4.2%/5′ v1/26 IgAN donors/specific/GCST90018888 等）、
# 正文表均被引用、关键新句已入母本。
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os, re, sys

BASE = _paths.at("阶段4.5_复现包/submission")
MS = BASE("manuscript", "Manuscript_GM_fulltext.md")
A1 = BASE("additional_file_1", "Additional_File_1_Figure_Legends.md")
CL = BASE("cover_letter_draft.md")

fails = []
def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    if not ok:
        fails.append(name)

ft = open(MS, encoding="utf-8").read()
# 去头部装配注（<!-- --> 块）与 References 节，供引文序检查
body = re.sub(r"<!--.*?-->", "", ft, flags=re.S)

# --- 1. Abstract 词数 ---
lines = body.splitlines()
def sec(lines, pat):
    out, on = [], False
    for ln in lines:
        s = ln.strip()
        if re.match(pat, s):
            if on: break
            on = True
        elif on and s.startswith("## "):
            break
        if on: out.append(ln)
    return out
ab = sec(lines, r"^## Abstract\b")
n = sum(len(l.strip().split()) for l in ab if l.strip() and not l.strip().startswith("#"))
check("Abstract words == 348 (limit 350)", n == 348, f"measured {n}")

# --- 2. 引文首现序 [1..29] 严格递增（正文，不含 References 节）---
main_part = body.split("## References")[0]
first_pos = {}
for m in re.finditer(r"\[(\d{1,2})(?:[,\s–\-\]\d]*)", main_part):
    # 取方括号组内所有编号
    grp = m.group(0)
    for nm in re.findall(r"\d{1,2}", grp):
        k = int(nm)
        if 1 <= k <= 29 and k not in first_pos:
            first_pos[k] = m.start()
seq = [first_pos.get(k, -1) for k in range(1, 30)]
missing = [k for k in range(1, 30) if first_pos.get(k, -1) == -1]
# 非递减：同组并列首现（如 [6, 7, 8] 同括号组同位置）合法；真倒挂（后号先现）判 FAIL
ascending = all(seq[i] <= seq[i+1] for i in range(28)) and not missing
check("citations [1..29] first-appearance non-decreasing (ties allowed within one bracket group)", ascending,
      f"missing={missing}" if missing else ("order broken" if not ascending else "1..29 OK"))

# --- 3. 正文中文 0 ---
cjk = re.findall(r"[\u4e00-\u9fff]", main_part)
check("no CJK in main text", len(cjk) == 0, f"{len(cjk)} chars")

# --- 4. 内部编号（白名单：Fig 1c 图注的 (L0–L7)）---
tmp = main_part.replace("(L0–L7)", "")
codes = re.findall(r"\b(?:L[0-7]|R[1-8])\b", tmp)
check("no internal analysis codes in running text", len(codes) == 0, f"{codes[:6]}")

# --- 5. 陈旧串 0（母本）---
STALE_MASTER = [
 "5-azacytidine-reversible",
 "reversed by 5-azacytidine [20]",
 "hypermethylation by 5-azacytidine in an IgAN-relevant",
 "smallest adjusted q = 0.28",
 "25,095–28,181",
 "4.2% of all B cells",
 "10x Genomics 5′ v1",
 "naive B-cell-specific",
 "cis-mQTL associations for the four lead variants",
 "GCST90018888",
 "5-Aza-CR",
 "(PPH1 0.18–0.28)",
 "in the tissue eQTL comparisons the posterior mass",
 "the four combined SuSiE credible sets",
 "is already being tested in IgAN [2],",
 "serum total IgA levels [25]",
 "Fung et al. [26]",
 "282,463 single cells) [27]",
 "(Fung et al. [26]) is a preprint",
 "histone deacetylase 1 by the green-tea catechin",
 "10 gene–context pairs) was looked up regionally.\n",
 "| +0.20 (0.025) |",
 "| +0.14 (0.012) |",
 "| +0.14 (0.013) |",
 "| −0.34 (0.023) |",
 "| +0.42 (0.010) |",
 "gives a higher estimate (4.2%)",
 "(GoDMC whole blood); standard error in parentheses. ‡Mendelian-randomization estimate",
]
bad = [s for s in STALE_MASTER if s in main_part]
check("stale strings in master == 0", not bad, f"{bad}")

# --- 6. 正文表均被正文引用 ---
for t in ["Table 2", "Table 3", "Table 4", "Table 5", "Table 6", "Table 7"]:
    cnt = len(re.findall(re.escape(t) + r"(?![\d.])", main_part.split("## Figure legends")[0]))
    check(f"{t} cited in running text", cnt >= 1, f"{cnt} occurrence(s)")

# --- 7. 关键新句已入母本 ---
NEW_MASTER = [
 "5-aza-2′-deoxycytidine [20]; in the IL-17 setting, 5-azacytidine likewise restored",
 "two interpretable autosomal C1GALT1 lead variants (rs13226913 and rs10238682)",
 "smallest adjusted q = 0.12",
 "24,988–28,181 per CpG",
 "unreplicated candidate-level observations",
 "(PPH1 = 0.79–0.91; Table 4)",
 "detected in 4.1% of all B cells",
 "is already being tested in IgAN [29]",
 "naive B-cell-enriched cis-eQTL",
 "seven signal-level tests",
 "reversed by 5-aza-2′-deoxycytidine [20])",
 "| NK | natural killer |",
 "24. Badia-i-Mompel P,",
 "29. Mathur M,",
 "10x Genomics 5′ gene-expression",
 "decoupleR's univariate linear model (ULM) [24]",
 "position-frequency matrices [25]",
 "| +0.20 (0.066) |",
 "| +0.14 (0.014) |",
 "| +0.14 (0.019) |",
 "| −0.34 (0.057) |",
 "| +0.42 (0.011) |",
 "random-effects meta-analysis estimates, with the random-effects standard error in parentheses",
 "gives a higher estimate (4.1%)",
 "Assembled file v1.13",
]
missing_new = [s for s in NEW_MASTER if s not in ft]
check("new corrected sentences present in master", not missing_new, f"{missing_new}")

# --- 8. Additional file 1 陈旧串 ---
a1 = open(A1, encoding="utf-8").read()
STALE_A1 = ["10x Genomics 5′ v1", "26 IgAN donors", "Kim *et al.* [27]",
            "β = +0.186, P = 1.2 × 10⁻⁴", "two- to five-fold",
            "(p = 0.046", "(p = 0.018", "(p = 0.11)",
            "standard error, allele frequency and the per-cell call",
            "the 14 computable colocalisation comparisons: credible-set size",
            "Per-lead IgAN lookup",
            "The target-matched random-set null and the empirical P values are summarised in the Fig. S6 legend",
            "(ρ, nominal P, FDR q, target-matched random-set 95th percentile and empirical P; 950 combinations)",
            "CpG identifier, chromosome, hg19 position, lead-variant hits and annotation"]
bad_a1 = [s for s in STALE_A1 if s in a1]
check("stale strings in Additional file 1 == 0", not bad_a1, f"{bad_a1}")
NEW_A1 = ["26 donors (17 with IgAN: 6 late-stage, 11 early-stage; 9 healthy controls",
          "Kim *et al.* [28]", "β = +0.185, P = 1.16 × 10⁻⁴", "two- to twelve-fold",
          "103 units across four B-cell states and 26 donors",
          "Per-context top-variant IgAN lookup", "5,700 rows",
          "axis_composite_null_q95", "axis_composite_emp_p"]
missing_a1 = [s for s in NEW_A1 if s not in a1]
check("new corrected strings present in Additional file 1", not missing_a1, f"{missing_a1}")
# A1 中文 0
cjk_a1 = re.findall(r"[\u4e00-\u9fff]", a1)
check("no CJK in Additional file 1", len(cjk_a1) == 0, f"{len(cjk_a1)}")

# --- 9. cover letter ---
cl = open(CL, encoding="utf-8").read()
bad_cl = [s for s in ["5-azacytidine-reversible", "naive B-cell-specific", "Mendelian randomisation"] if s in cl]
check("stale strings in cover letter == 0", not bad_cl, f"{bad_cl}")
ok_cl = all(s in cl for s in ["DNMT-inhibitor–reversible promoter hypermethylation",
                              "naive B-cell-enriched cis-eQTL", "Mendelian randomization"])
check("new strings present in cover letter", ok_cl)

# --- 10. Table S10 结构（P1-2 合并后）---
import csv
S10 = BASE("additional_file_1", "Table_S10_regulon_axis_associations.tsv")
with open(S10, encoding="utf-8", newline="") as f:
    rdr = csv.DictReader(f, delimiter="\t")
    s10 = list(rdr)
    s10cols = list(rdr.fieldnames)
check("S10 has 5,700 rows", len(s10) == 5700, f"{len(s10)}")
check("S10 carries axis_composite_null_q95/emp_p",
      s10cols[-2:] == ["axis_composite_null_q95", "axis_composite_emp_p"], f"{s10cols}")
comp = [r for r in s10 if r["target"] == "AxisComposite"]
check("S10 composite rows annotated (950)", len(comp) == 950 and
      all(r["axis_composite_emp_p"] for r in comp), f"{len(comp)}")
per_gene_empty = all(not r["axis_composite_emp_p"] for r in s10 if r["target"] != "AxisComposite")
check("S10 per-gene rows leave composite fields empty", per_gene_empty)
dnmt1 = [r for r in comp if r["b_subset"] == "B_plasma" and r["tf"] == "DNMT1"]
check("S10 DNMT1 plasma emp_p == 0.025", len(dnmt1) == 1 and dnmt1[0]["axis_composite_emp_p"] == "0.025")

print()
if fails:
    print("TOTAL FAIL:", len(fails), fails)
    sys.exit(1)
print("ALL MACHINE CHECKS PASSED")
