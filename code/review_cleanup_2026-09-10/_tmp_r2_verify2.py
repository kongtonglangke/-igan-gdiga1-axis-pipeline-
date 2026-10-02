# -*- coding: utf-8 -*-
"""Verify2: word counts, md5, manifest, misc residual checks."""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os, re, json, csv, hashlib, glob

ROOT = _paths.LEGACY
MS = ROOT(r"阶段4.5_复现包\submission\manuscript")
SUB = ROOT(r"阶段4.5_复现包\submission\additional_file_1")
REPRO = ROOT(r"阶段4.5_复现包\reproducibility")

# ---------- 2. WORD COUNTS ----------
def wc_section(text, keep_headings):
    """raw whitespace tokenisation; drop '>' header lines, '---' separators, table rows."""
    n = 0
    for ln in text.splitlines():
        s = ln.strip()
        if not s or s.startswith(">") or set(s) <= set("-: ") or s.startswith("|"):
            continue
        if s.startswith("#"):
            if keep_headings:
                n += len(s.lstrip("#").strip().split())
            continue
        n += len(s.split())
    return n

def get_body(fp):
    with open(fp, encoding="utf-8") as f:
        t = f.read()
    # body starts after the editorial header: first '---' line
    parts = re.split(r"(?m)^---\s*$", t)
    return parts[-1] if len(parts) > 1 else t

def sec(text, start_pat, end_pats):
    lines = text.splitlines()
    out, on = [], False
    for ln in lines:
        if re.match(start_pat, ln):
            on = True; out.append(ln); continue
        if on and any(re.match(p, ln) for p in end_pats):
            break
        if on: out.append(ln)
    return "\n".join(out)

print("### WORD COUNTS (declared: Abstract 348; Methods 3626/3521; Results 3211/3094; Disc 1600/1598; Lim 605/603; Concl 148/146) ###")
# Abstract: ## Abstract section of 01
b1 = get_body(os.path.join(MS, "01_Title_Abstract_Keywords.md"))
ab = sec(b1, r"^## Abstract\b", [r"^## "])
print("Abstract incl/excl:", wc_section(ab, True), wc_section(ab, False), "(decl 348; v1.9 note: labels incl, headings excl)")

# Methods 03: exclude Table 1 block (caption line **Table 1.** + table rows already excluded + note paragraph)
b3 = get_body(os.path.join(MS, "03_Methods.md"))
m = re.search(r"(?m)^## Methods\s*$", b3)
b3m = b3[m.start():] if m else b3
# remove Table 1: from '**Table 1.**' to the line starting '\* GEO series'
t1s = b3m.find("**Table 1.**")
t1e = b3m.find("GEO accession numbers throughout")
if t1s != -1 and t1e != -1:
    t1e = b3m.find("\n", t1e)
    b3m = b3m[:t1s] + b3m[t1e:]
print("Methods incl/excl (Table1 removed):", wc_section(b3m, True), wc_section(b3m, False))

# Results 04
b4 = get_body(os.path.join(MS, "04_Results.md"))
m = re.search(r"(?m)^## Results\s*$", b4)
b4r = b4[m.start():] if m else b4
print("Results incl/excl:", wc_section(b4r, True), wc_section(b4r, False))

# Discussion & Limitations 05
b5 = get_body(os.path.join(MS, "05_Discussion_Limitations.md"))
d = sec(b5, r"^## Discussion\b", [r"^## "])
l = sec(b5, r"^## Limitations\b", [r"^## "])
print("Discussion incl/excl:", wc_section(d, True), wc_section(d, False))
print("Limitations incl/excl:", wc_section(l, True), wc_section(l, False))

# Conclusions 06
b6 = get_body(os.path.join(MS, "06_Conclusions_Abbreviations.md"))
c = sec(b6, r"^## Conclusions\b", [r"^## "])
print("Conclusions incl/excl:", wc_section(c, True), wc_section(c, False))

# ---------- 3. MD5 ----------
print("\n### MD5 fig_scripts vs 阶段4/主图/scripts ###")
d1 = ROOT(r"阶段4\主图\scripts")
d2 = os.path.join(REPRO, r"code\fig_scripts")
def md5(p):
    with open(p, "rb") as f: return hashlib.md5(f.read()).hexdigest()
f1 = {f for f in os.listdir(d1) if f.endswith(".py") and not f.startswith("_tmp")}
f2 = {f for f in os.listdir(d2) if f.endswith(".py")}
for f in sorted(f1 | f2):
    p1, p2 = os.path.join(d1, f), os.path.join(d2, f)
    if f not in f1: print("ONLY in repro:", f)
    elif f not in f2: print("ONLY in stage4:", f)
    elif md5(p1) != md5(p2): print("DIVERGE:", f, md5(p1)[:8], md5(p2)[:8])
    else: print("same:", f)

# ---------- 4. MANIFEST ----------
print("\n### data manifest vs tree ###")
with open(os.path.join(REPRO, "README.md"), encoding="utf-8") as f:
    rd = f.read()
print("README total lines:", len(rd.splitlines()))
m = re.search(r"## 数据目录\n(.*?)\n## ", rd, re.S)
man = m.group(1) if m else ""
print("manifest paragraph:", man[:600])
actual = sorted(os.listdir(os.path.join(REPRO, "data")))
print("actual data subdirs:", actual)
for d in actual:
    print(d, "in-manifest-text:", d in rd)

# ---------- misc ----------
print("\n### S5b mr_se (broadcast check) ###")
with open(os.path.join(SUB, "Table_S5b_MR_CpG_x_IgAN.tsv"), encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["method"]=="Wald_lead":
            z = float(r["mr_beta"])/float(r["mr_se"]) if r["mr_se"] else None
            print(r["cpg"], r["instrument"], "b=",r["mr_beta"][:7], "se=",r["mr_se"][:7], "z=", None if z is None else round(z,4), "p=", r["mr_p"][:20])

print("\n### emp_p<0.05 count ###")
with open(os.path.join(REPRO, r"results\sc_regulatory\regulon_axis_composite.tsv"), encoding="utf-8") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))
np_ = sum(1 for r in rows if r["emp_p"] and float(r["emp_p"]) < 0.05)
print("rows:", len(rows), "emp_p<0.05:", np_, "(claim 36/950)")

print("\n### GoDMC sample sizes (Table 1: 25,095-28,181) ###")
with open(os.path.join(SUB, "Table_S5a_sigCpG_5CpG_pos.tsv"), encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        info = r["info"]
        m2 = re.search(r"'samplesize': (\d+)", info)
        print(r["cpg"], "n=", m2.group(1) if m2 else "?")

print("\n### Table 7 Pharos tiers (M6) ###")
for fn in [r"M6_可干预性评分.tsv", r"M6_靶点药物证据表.tsv"]:
    p = os.path.join(REPRO, "data", "M6_intervention", fn)
    if os.path.exists(p):
        with open(p, encoding="utf-8", errors="replace") as f:
            print(fn, "->", f.read()[:900])

print("\n### C1 step1 summary (ImmuNexUT 416/28) ###")
with open(os.path.join(REPRO, r"results\C1_eas_eqtl\step1_summary.json"), encoding="utf-8") as f:
    print(f.read()[:800])

print("\n### step3 summary (JASPAR 879) ###")
with open(os.path.join(REPRO, r"results\sc_regulatory\step3_summary.json"), encoding="utf-8") as f:
    print(f.read()[:900])

print("\n### GCST90018866 accession consistency ###")
for fp in glob.glob(os.path.join(MS, "*.md")) + glob.glob(os.path.join(SUB, "*.md")):
    with open(fp, encoding="utf-8", errors="replace") as f:
        t = f.read()
    for acc in set(re.findall(r"GCST\d+", t)):
        if acc not in ("GCST90018866","GCST90011884"):
            print(os.path.basename(fp), "unexpected acc:", acc)
print("done acc scan")

print("\n### DV group split 6/11/9 ###")
p = os.path.join(REPRO, r"results\directional_validation\bcell_group_comparison.tsv")
with open(p, encoding="utf-8", errors="replace") as f:
    print(f.read()[:700])
