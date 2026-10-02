# -*- coding: utf-8 -*-
"""Post-edit machine checks for the submission package (2026-09-10 v1.8 pass).

  A. Abstract word count (ch1 and the assembled mother file must agree).
  B. Reference first-appearance order in the mother file must be strictly 1..27.
  C. No bare internal analysis codes (L0-L7 / R1-R8) in running text.
  D. No CJK characters anywhere in the manuscript chapters.
  E. No leftover "seven ... evidence layer(s)" phrasing.
"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import io, re, glob, os

PKG = _paths.at("阶段4.5_复现包")
MAN = PKG("submission", "manuscript")


def rd(p):
    return io.open(p, encoding="utf-8").read()


def abstract_count(txt):
    m = re.search(r"## Abstract\s*(.*?)(?=\n## )", txt, re.S)
    blk = m.group(1)
    tot = 0
    for lab in ("Background:", "Methods:", "Results:", "Conclusions:"):
        mm = re.search(re.escape(lab) + r"(.*?)(?=\n\n|\Z)", blk, re.S)
        if mm:
            tot += len(mm.group(1).replace("**", "").split())
    return tot + 4  # the four labels


def strip_assembly(txt):
    return re.sub(r"<!--.*?-->", "", txt, flags=re.S)


def main():
    mother = rd(os.path.join(MAN, "Manuscript_GM_fulltext.md"))
    ch1 = rd(os.path.join(MAN, "01_Title_Abstract_Keywords.md"))

    print("A. Abstract word count (labels included, limit 350)")
    print(f"     ch1            = {abstract_count(ch1)}")
    print(f"     mother file    = {abstract_count(mother)}")

    body = strip_assembly(mother)
    print("B. Reference first-appearance order")
    nums = []
    for m in re.finditer(r"\[(\d+(?:\s*[,\u2013\u2014-]\s*\d+)*)\]", body):
        for part in re.split(r"[,,\u2013\u2014-]", m.group(1)):
            part = part.strip()
            if part.isdigit():
                nums.append(int(part))
    first_order, seen = [], set()
    for n in nums:
        if n not in seen:
            seen.add(n)
            first_order.append(n)
    expect = list(range(1, max(first_order) + 1)) if first_order else []
    print(f"     distinct cited = {len(first_order)}  (max = {max(first_order) if first_order else 0})")
    print(f"     first-appearance order = {first_order}")
    print(f"     strictly increasing 1..N : {first_order == expect}")
    gaps = sorted(set(expect) - seen) if expect else []
    print(f"     missing numbers         : {gaps if gaps else 'none'}")

    print("C. Bare internal codes (L0-L7 / R1-R8) in the assembled body")
    hits = []
    for i, ln in enumerate(body.splitlines(), 1):
        if "(L0\u2013L7)" in ln or "(L0-L7)" in ln:      # intentional in-panel index
            continue
        if re.search(r"\b(L[0-7]|R[1-8])\b", ln):
            hits.append((i, ln[:110]))
    for i, t in hits[:20]:
        print(f"     line {i}: {t}")
    print(f"     total = {len(hits)}")

    print("D. CJK characters in the manuscript chapters")
    tot = 0
    for f in sorted(glob.glob(os.path.join(MAN, "0*.md"))):
        n = len(re.findall(r"[\u4e00-\u9fff]", rd(f)))
        tot += n
        if n:
            print(f"     {os.path.basename(f)}: {n}")
    print(f"     total = {tot}")

    print("E. Leftover 'seven ... evidence layer' phrasing")
    for f in sorted(glob.glob(os.path.join(MAN, "0*.md"))):
        for i, ln in enumerate(rd(f).splitlines(), 1):
            if ln.lstrip().startswith(">"):
                continue
            if re.search(r"seven[^.]{0,40}evidence layer", ln, re.I):
                print(f"     {os.path.basename(f)}:{i}: {ln[:110]}")
    print("     (no output above = clean)")


if __name__ == "__main__":
    main()
