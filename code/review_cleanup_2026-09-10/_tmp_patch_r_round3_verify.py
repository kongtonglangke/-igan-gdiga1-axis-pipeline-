# -*- coding: utf-8 -*-
"""
_tmp_patch_r_round3_verify.py — 第三轮整改后的机检复核（修订口径 v2）

口径修正（v1 的误报已定位）：
  · 引用首现序只在 ch1–ch8 正文上计算，**排除 ch9 参考文献表**（期刊卷期号如 (7993)/(6509)
    会被括号正则误当作引用）；匹配的括号组必须是**纯引用列表**（逗号分隔的 1–2 位整数，
    不含连接号），从而排除 G(0–3)/D(0–3)/F(0–2)/P(0–2)/(0–10) 这类评分区间。
  · Abstract 词数排除 `*Word count: …*` 编辑注行。
  · 主文表只校验 Table 1–7（本文无 Table 8）。
  · 陈旧串检查前先剥掉母本的 HTML 装配说明块（其中引述了旧串名）。
"""
import os
import re

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
M = os.path.join(BASE, "submission", "manuscript")
FILES = ["01_Title_Abstract_Keywords.md", "02_Background.md", "03_Methods.md", "04_Results.md",
         "05_Discussion_Limitations.md", "06_Conclusions_Abbreviations.md", "07_TitlePage_Declarations.md",
         "08_Figure_Legends_Tables.md"]
MOTHER = os.path.join(M, "Manuscript_GM_fulltext.md")


def body(f):
    out = []
    for l in open(os.path.join(M, f), encoding="utf-8").read().split("\n"):
        if l.startswith(">") or re.match(r'^\s*(-{3,}|\|)', l) or l.startswith("Figure calls in text:"):
            continue
        out.append(l)
    return "\n".join(out)


def main():
    ok = True
    allb = "\n".join(body(f) for f in FILES)

    # A. Abstract
    ch1 = open(os.path.join(M, "01_Title_Abstract_Keywords.md"), encoding="utf-8").read()
    ab = ch1[ch1.find("## Abstract"):ch1.find("## Keywords")]
    paras = [p for p in ab.split("\n")
             if p.strip() and not p.startswith("#") and not p.startswith("*Word count")
             and not re.fullmatch(r'-{3,}', p.strip())]
    labelled = " ".join(paras)
    nolab = re.sub(r'\*\*(Background|Methods|Results|Conclusions):\*\*', '', labelled)
    wl, wn = len(labelled.split()), len(nolab.split())
    print(f"A. Abstract words: {wl} (labels incl) / {wn} (labels excl); limit 350 -> {'PASS' if wl <= 350 else 'FAIL'}")
    ok &= wl <= 350

    # B. main-text tables 1-7 cited in ascending order
    mainb = "\n".join(body(f) for f in FILES if f != "08_Figure_Legends_Tables.md")
    tabs = [(n, len(re.findall(r'\bTable\s+%d\b' % n, mainb))) for n in range(1, 8)]
    for n, c in tabs:
        print(f"B. Table {n}: {c}x cited" + ("" if c else "  <-- UNCITED"))
    b_ok = all(c >= 1 for _, c in tabs)
    print(f"B. every main-text table (1-7) cited -> {'PASS' if b_ok else 'FAIL'}")
    ok &= b_ok

    # C. reference first-appearance order (body text only; pure citation lists only)
    seq = []
    for m in re.finditer(r'[\[\(](\d{1,2}(?:,\s*\d{1,2})*)[\]\)]', allb):
        seq.extend(int(x) for x in re.findall(r'\d{1,2}', m.group(1)))
    first = {}
    for k, n in enumerate(seq):
        first.setdefault(n, k)
    fo = sorted(first, key=lambda k: first[k])
    c_ok = fo == list(range(1, 28))
    print(f"C. citation first-appearance order -> {'PASS (1..27 strictly ascending)' if c_ok else 'FAIL ' + str(fo)}")
    ok &= c_ok

    # D. Chinese characters
    zh = re.findall(r'[\u4e00-\u9fff]', allb)
    print(f"D. Chinese characters in submission text: {len(zh)} -> {'PASS' if not zh else 'FAIL'}")
    ok &= not zh

    # E. internal codes in mother (HTML assembly block stripped)
    mt = re.sub(r'<!--.*?-->', '', open(MOTHER, encoding="utf-8").read(), flags=re.S)
    h = {'evidence layer L': len(re.findall(r'evidence layer L\d', mt)),
         '(Rn)': len(re.findall(r'\(R[1-8]\)', mt)),
         'bare Rn': len(re.findall(r'(?<![\w.])R[1-8](?![\w])', mt)),
         'bare Mn': len(re.findall(r'(?<![\w.])M[1-6](?![\w])', mt)),
         '(L0-L7)': len(re.findall(r'\(L0[–-]L7\)', mt))}
    e_ok = h['evidence layer L'] == 0 and h['(Rn)'] == 0 and h['bare Rn'] == 0 and h['bare Mn'] == 0 and h['(L0-L7)'] == 1
    print(f"E. internal codes {h} -> {'PASS' if e_ok else 'FAIL'}")
    ok &= e_ok

    # F. stale / new strings (mother treated with the assembly block stripped)
    stale = ["three orders of magnitude", "attenuate along the naive-to-memory transition",
             "essentially flat for C1GALT1", "seven sequential evidence layers", "PPH3/PPH4",
             "colocalization", "naïve", "six sequential evidence"]
    for k in stale:
        n = allb.count(k) + mt.count(k)
        print(f"F. stale {k!r}: {n} -> {'PASS' if n == 0 else 'FAIL'}")
        ok &= n == 0
    checks = [("two orders of magnitude above the genome-wide significance threshold", 2),
              ("with partial recovery in memory", 1),
              ("PPH0–PPH4", 2),
              ("specificity protein/Krüppel-like factor (Sp/KLF)", 1),
              ("Illuminating the Druggable Genome", 1),
              ("National Bioscience Database Center (NBDC)", 1),
              ("effect of rs13226913 is gated by B-cell state", 1),
              ("(Fig. 4a; Table 5)", 1)]
    for k, v in checks:
        n = allb.count(k)
        print(f"F. new {k!r}: {n} (expect {v}) -> {'PASS' if n == v else 'FAIL'}")
        ok &= n == v

    print(f"\nOVERALL: {'ALL CHECKS PASS' if ok else 'SOME CHECKS FAILED'}")


if __name__ == "__main__":
    main()
