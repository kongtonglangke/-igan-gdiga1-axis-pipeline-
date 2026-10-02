# -*- coding: utf-8 -*-
# _tmp_measure_wordcount_2026-09-11.py — 全维度审查修复后词数实测（章节文件直测）
# 口径（与既往声明一致并明确化）：原始空白分词；剔除 '>' 头部注、'---' 分隔符、
# 表格行（'|' 起始）；incl = 含小节标题，excl = 不含；Abstract 仅 ## Abstract 节；
# Methods 整段剔除 Table 1（caption "**Table 1.**" 起至 "GEO accession numbers throughout"
# 所在行止，含表体与表注）；Results 止于下一 "## "；Discussion/Limitations/Conclusions 各节独立。
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os, re

MS = _paths.at("阶段4.5_复现包/submission/manuscript")

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
            if on:
                break
            on = True
        elif on and s.startswith("## "):
            break
        if on:
            out.append(ln)
    return out

def load(name):
    with open(MS(name), encoding="utf-8") as f:
        return f.read().splitlines()

c01 = load("01_Title_Abstract_Keywords.md")
ab = sec(c01, r"^## Abstract\b")
# 与装配脚本一致：剔除 "*Word count" 词数声明行（_build_submission_fulltext.py:88）
ab = [l for l in ab if not l.strip().startswith("*Word count")]
print("Abstract        incl/excl: %d / %d   (decl 348; GM limit 350)" % (wc(ab, True), wc(ab, False)))

c03 = load("03_Methods.md")
mm = sec(c03, r"^## Methods\b")
txt = "\n".join(mm)
t1s = txt.find("**Table 1.**")
t1e = txt.find("GEO accession numbers throughout")
assert t1s != -1 and t1e != -1, "Table 1 anchors missing in 03"
t1e = txt.find("\n", t1e)
m_not1 = (txt[:t1s] + txt[t1e:]).splitlines()
print("Methods (excl Table 1 caption+rows+note) incl/excl: %d / %d" % (wc(m_not1, True), wc(m_not1, False)))
print("Methods (Table 1 rows only auto-dropped)          incl/excl: %d / %d" % (wc(mm, True), wc(mm, False)))

c04 = load("04_Results.md")
rr = sec(c04, r"^## Results\b")
print("Results         incl/excl: %d / %d   (decl 3211/3094)" % (wc(rr, True), wc(rr, False)))

c05 = load("05_Discussion_Limitations.md")
d = sec(c05, r"^## Discussion\b")
l = sec(c05, r"^## Limitations\b")
print("Discussion      incl/excl: %d / %d   (decl 1600/1598)" % (wc(d, True), wc(d, False)))
print("Limitations     incl/excl: %d / %d   (decl 605/603)" % (wc(l, True), wc(l, False)))

c06 = load("06_Conclusions_Abbreviations.md")
c = sec(c06, r"^## Conclusions\b")
print("Conclusions     incl/excl: %d / %d   (decl 148/146)" % (wc(c, True), wc(c, False)))

# 正文合计（Methods 取整段剔除 Table 1 口径）
tot_i = wc(ab, True) + wc(m_not1, True) + wc(rr, True) + wc(d, True) + wc(l, True) + wc(c, True)
tot_e = wc(ab, False) + wc(m_not1, False) + wc(rr, False) + wc(d, False) + wc(l, False) + wc(c, False)
print("Main-text total (Abstract+Methods excl-T1+Results+Discussion+Limitations+Conclusions)")
print("                incl/excl: %d / %d" % (tot_i, tot_e))
