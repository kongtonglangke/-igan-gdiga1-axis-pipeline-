# -*- coding: utf-8 -*-
"""Word count v3: measure on assembled fulltext Manuscript_GM_fulltext.md.
Same rules: raw whitespace tokens; drop '>' notes, '---' separators, table rows;
Methods excludes Table 1 (caption+rows+note); Abstract = its section only."""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os, re

MS = _paths.at("阶段4.5_复现包/submission/manuscript")
ft = MS("Manuscript_GM_fulltext.md")
with open(ft, encoding="utf-8") as f:
    lines = f.read().splitlines()

def wc(lines, keep_headings):
    n = 0
    for ln in lines:
        s = ln.strip()
        if not s or s.startswith(">") or s == "---" or s.startswith("|"):
            continue
        if s.startswith("#"):
            if keep_headings:
                n += len(s.split())
            continue
        n += len(s.split())
    return n

def sec(lines, start_pat):
    out, on = [], False
    for ln in lines:
        s = ln.strip()
        if re.match(start_pat, s):
            if on:  # hit a second occurrence -> stop
                break
            on = True
        elif on and s.startswith("## "):
            break
        if on:
            out.append(ln)
    return out

# list headings to know structure
print("HEADINGS:")
for ln in lines:
    if ln.strip().startswith("## "):
        print("  ", ln.strip()[:80])

ab = sec(lines, r"^## Abstract\b")
print("Abstract incl/excl:", wc(ab, True), wc(ab, False), "decl 348")

mm = sec(lines, r"^## Methods\b")
txt = "\n".join(mm)
t1s = txt.find("**Table 1.**")
t1e = txt.find("GEO accession numbers throughout")
if t1s != -1 and t1e != -1:
    t1e = txt.find("\n", t1e)
    txt = txt[:t1s] + txt[t1e:]
print("Methods incl/excl:", wc(txt.splitlines(), True), wc(txt.splitlines(), False), "decl 3626/3521")
# also Methods WITH table-1 note kept variant (only rows auto-dropped)
print("Methods(raw, no T1 removal) incl/excl:", wc(mm, True), wc(mm, False))

rr = sec(lines, r"^## Results\b")
print("Results incl/excl:", wc(rr, True), wc(rr, False), "decl 3211/3094")
d = sec(lines, r"^## Discussion\b")
l = sec(lines, r"^## Limitations\b")
print("Discussion incl/excl:", wc(d, True), wc(d, False), "decl 1600/1598")
print("Limitations incl/excl:", wc(l, True), wc(l, False), "decl 605/603")
c = sec(lines, r"^## Conclusions\b")
print("Conclusions incl/excl:", wc(c, True), wc(c, False), "decl 148/146")

# chapter-file 01 Abstract section headings detail
with open(MS("01_Title_Abstract_Keywords.md"), encoding="utf-8") as f:
    t1 = f.read().splitlines()
started = False
for i, ln in enumerate(t1):
    if ln.strip() == "---" and not started:
        started = i
    if started and ln.strip().startswith("##"):
        print("01 heading @", i, ln.strip()[:70])
