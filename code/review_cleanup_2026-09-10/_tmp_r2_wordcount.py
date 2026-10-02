# -*- coding: utf-8 -*-
"""Word count v2: body = everything after FIRST '---' line; raw whitespace
tokenisation; drop '>' notes, '---' separators, table rows.
incl = count heading lines fully (incl. ## markers as tokens); excl = skip headings."""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os, re

MS = _paths.at("阶段4.5_复现包/submission/manuscript")

def body_after_first_sep(fp):
    with open(fp, encoding="utf-8") as f:
        lines = f.read().splitlines()
    for i, ln in enumerate(lines):
        if ln.strip() == "---":
            return lines[i+1:]
    return lines

def wc(lines, keep_headings):
    n = 0
    for ln in lines:
        s = ln.strip()
        if not s or s.startswith(">") or s == "---" or s.startswith("|"):
            continue
        if s.startswith("#"):
            if keep_headings:
                n += len(s.split())   # '## Methods' -> 2 tokens (marker + word)
            continue
        n += len(s.split())
    return n

def sec(lines, start_pat, end_prefix="## "):
    out, on = [], False
    for ln in lines:
        s = ln.strip()
        if re.match(start_pat, s):
            on = True
        elif on and s.startswith(end_prefix):
            break
        if on:
            out.append(ln)
    return out

decl = {"Abstract": (348, 348), "Methods": (3626, 3521), "Results": (3211, 3094),
        "Discussion": (1600, 1598), "Limitations": (605, 603), "Conclusions": (148, 146)}

# Abstract (## Abstract section of 01)
b = body_after_first_sep(MS("01_Title_Abstract_Keywords.md"))
ab = sec(b, r"^## Abstract\b")
print("Abstract incl/excl:", wc(ab, True), wc(ab, False), "decl 348 (labels incl, headings excl)")

# Methods (exclude Table 1 block incl. its note paragraph)
b = body_after_first_sep(MS("03_Methods.md"))
mm = sec(b, r"^## Methods\b")
txt = "\n".join(mm)
t1s = txt.find("**Table 1.**")
t1e = txt.find("GEO accession numbers throughout")
if t1s != -1 and t1e != -1:
    t1e = txt.find("\n", t1e)
    txt = txt[:t1s] + txt[t1e:]
mm2 = txt.splitlines()
print("Methods incl/excl:", wc(mm2, True), wc(mm2, False), "decl 3626/3521")

# Results
b = body_after_first_sep(MS("04_Results.md"))
rr = sec(b, r"^## Results\b")
print("Results incl/excl:", wc(rr, True), wc(rr, False), "decl 3211/3094")

# Discussion / Limitations
b = body_after_first_sep(MS("05_Discussion_Limitations.md"))
d = sec(b, r"^## Discussion\b")
l = sec(b, r"^## Limitations\b")
print("Discussion incl/excl:", wc(d, True), wc(d, False), "decl 1600/1598")
print("Limitations incl/excl:", wc(l, True), wc(l, False), "decl 605/603")

# Conclusions
b = body_after_first_sep(MS("06_Conclusions_Abbreviations.md"))
c = sec(b, r"^## Conclusions\b")
print("Conclusions incl/excl:", wc(c, True), wc(c, False), "decl 148/146")
