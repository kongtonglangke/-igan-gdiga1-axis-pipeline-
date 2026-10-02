# -*- coding: utf-8 -*-
"""
_tmp_patch_u_round3b_header_verify.py
第三轮作者裁定办结（② 删「实验室在研实验」前瞻段）——装配说明升 v1.10 + 母本重建 + 全量机检。

改动：
  · reproducibility/code/_build_submission_fulltext.py 头部装配说明 v1.9 → v1.10，
    新增 "Author-decision follow-up 2026-09-10 (v1.10)" 说明行。
  · 重建 submission/manuscript/Manuscript_GM_fulltext.md。
  · 机检 A–H（沿用第三轮口径，新增本轮陈旧串/新串）。

用法：python reproducibility/code/review_cleanup_2026-09-10/_tmp_patch_u_round3b_header_verify.py
"""
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))  # 阶段4.5_复现包
BUILD = os.path.join(BASE, "reproducibility", "code", "_build_submission_fulltext.py")
M = os.path.join(BASE, "submission", "manuscript")
MOTHER = os.path.join(M, "Manuscript_GM_fulltext.md")
FILES = ["01_Title_Abstract_Keywords.md", "02_Background.md", "03_Methods.md", "04_Results.md",
         "05_Discussion_Limitations.md", "06_Conclusions_Abbreviations.md",
         "07_TitlePage_Declarations.md", "08_Figure_Legends_Tables.md"]

NEW_NOTE = (
    "Author-decision follow-up 2026-09-10 (v1.10): the closing paragraph of the Discussion no longer "
    "states that the three functional tests are \"under way in our laboratory\" or that they \"will be "
    "reported separately\"; the predictions are now presented as what the framework makes directly "
    "testable, and the paragraph closes on the boundary-and-direction summary of the study's "
    "contribution (ch5 v1.8). The single-cell case-control computation stays out of the manuscript "
    "(author decision: not to include it). No data, values, figures or conclusions changed."
)


def patch_header():
    s = open(BUILD, encoding="utf-8").read()
    old_v = "Assembled file v1.9 (2026-09-10), machine-assembled from the pure-English chapter files"
    new_v = "Assembled file v1.10 (2026-09-10), machine-assembled from the pure-English chapter files"
    if s.count(old_v) == 1:
        s = s.replace(old_v, new_v)
        print("header: v1.9 -> v1.10")
    elif s.count(new_v) == 1:
        print("header: already v1.10 (idempotent skip)")
    else:
        print("FAIL: assembler version string hit count", s.count(old_v), s.count(new_v))
        sys.exit(1)

    if "Author-decision follow-up 2026-09-10 (v1.10)" not in s:
        lines = s.split("\n")
        idx = [k for k, l in enumerate(lines)
               if l.startswith("Third reviewer-perspective pass 2026-09-10 (v1.9):")]
        if len(idx) != 1:
            print("FAIL: v1.9 note line hit count =", len(idx)); sys.exit(1)
        lines.insert(idx[0] + 1, NEW_NOTE)
        s = "\n".join(lines)
        print("header: v1.10 note inserted")
    else:
        print("header: v1.10 note already present (idempotent skip)")

    open(BUILD, "w", encoding="utf-8", newline="\n").write(s)


def rebuild():
    r = subprocess.run([sys.executable, BUILD], cwd=BASE, capture_output=True, text=True)
    print((r.stdout or r.stderr).strip())
    if r.returncode != 0:
        sys.exit(1)


def body(f):
    out = []
    for l in open(os.path.join(M, f), encoding="utf-8").read().split("\n"):
        if l.startswith(">") or re.match(r'^\s*(-{3,}|\|)', l) or l.startswith("Figure calls in text:"):
            continue
        out.append(l)
    return "\n".join(out)


def verify():
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
    print(f"A. Abstract words: {wl} incl labels / {wn} excl; limit 350 -> {'PASS' if wl <= 350 else 'FAIL'}")
    ok &= wl <= 350

    # B. main-text tables 1-7 cited
    mainb = "\n".join(body(f) for f in FILES if f != "08_Figure_Legends_Tables.md")
    tabs = [(n, len(re.findall(r'\bTable\s+%d\b' % n, mainb))) for n in range(1, 8)]
    b_ok = all(c >= 1 for _, c in tabs)
    print("B. table citations:", {n: c for n, c in tabs}, "->", "PASS" if b_ok else "FAIL")
    ok &= b_ok

    # C. reference first-appearance order
    seq = []
    for m in re.finditer(r'[\[\(](\d{1,2}(?:,\s*\d{1,2})*)[\]\)]', allb):
        seq.extend(int(x) for x in re.findall(r'\d{1,2}', m.group(1)))
    first = {}
    for k, n in enumerate(seq):
        first.setdefault(n, k)
    fo = sorted(first, key=lambda k: first[k])
    c_ok = fo == list(range(1, 28))
    print(f"C. citation first-appearance order -> {'PASS (1..27)' if c_ok else 'FAIL ' + str(fo)}")
    ok &= c_ok

    # D. Chinese characters
    zh = re.findall(r'[\u4e00-\u9fff]', allb)
    print(f"D. Chinese chars in submission text: {len(zh)} -> {'PASS' if not zh else 'FAIL'}")
    ok &= not zh

    # E. internal codes in mother
    mt = re.sub(r'<!--.*?-->', '', open(MOTHER, encoding="utf-8").read(), flags=re.S)
    h = {'evidence layer L': len(re.findall(r'evidence layer L\d', mt)),
         '(Rn)': len(re.findall(r'\(R[1-8]\)', mt)),
         'bare Rn': len(re.findall(r'(?<![\w.])R[1-8](?![\w])', mt)),
         'bare Mn': len(re.findall(r'(?<![\w.])M[1-6](?![\w])', mt)),
         '(L0-L7)': len(re.findall(r'\(L0[–-]L7\)', mt))}
    e_ok = (h['evidence layer L'] == 0 and h['(Rn)'] == 0 and h['bare Rn'] == 0
            and h['bare Mn'] == 0 and h['(L0-L7)'] == 1)
    print(f"E. internal codes {h} -> {'PASS' if e_ok else 'FAIL'}")
    ok &= e_ok

    # F. stale strings
    stale = ["three orders of magnitude", "attenuate along the naive-to-memory transition",
             "essentially flat for C1GALT1", "seven sequential evidence layers", "PPH3/PPH4",
             "colocalization", "naïve", "six sequential evidence",
             "under way in our laboratory", "will be reported separately",
             "explicit tests for the orthogonal"]
    for k in stale:
        n = allb.count(k) + mt.count(k)
        print(f"F. stale {k!r}: {n} -> {'PASS' if n == 0 else 'FAIL'}")
        ok &= n == 0

    # G. new strings present
    checks = [("two orders of magnitude above the genome-wide significance threshold", 2),
              ("with partial recovery in memory", 1),
              ("PPH0–PPH4", 2),
              ("specificity protein/Krüppel-like factor (Sp/KLF)", 1),
              ("Illuminating the Druggable Genome", 1),
              ("National Bioscience Database Center (NBDC)", 1),
              ("effect of rs13226913 is gated by B-cell state", 1),
              ("(Fig. 4a; Table 5)", 1),
              ("What this framework contributes is therefore a boundary and a direction", 1),
              ("restores galactosylation capacity in B cells from patients with IgAN", 1)]
    for k, v in checks:
        n = allb.count(k)
        print(f"G. new {k!r}: {n} (expect {v}) -> {'PASS' if n == v else 'FAIL'}")
        ok &= n == v

    print(f"\nOVERALL: {'ALL CHECKS PASS' if ok else 'SOME CHECKS FAILED'}")
    return ok


def main():
    patch_header()
    rebuild()
    if not verify():
        sys.exit(1)


if __name__ == "__main__":
    main()
