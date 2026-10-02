# -*- coding: utf-8 -*-
"""Count the structured Abstract of Chapter 1 under the project convention:
whitespace tokenisation of the four labelled Abstract paragraphs only
('Background:', 'Methods:', 'Results:', 'Conclusions:' labels INCLUDED;
section headings and italic notes EXCLUDED).
"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import io, os, re, sys

CH1 = _paths.legacy("阶段4.5_复现包/submission/manuscript/01_Title_Abstract_Keywords.md")


def extract_abstract(path):
    txt = io.open(path, encoding="utf-8").read()
    lines = txt.splitlines()
    out, on = [], False
    for ln in lines:
        if ln.strip().startswith("## Abstract"):
            on = True
            continue
        if on and ln.strip().startswith("## "):
            break
        if on:
            out.append(ln)
    return "\n".join(out)


def count(path):
    block = extract_abstract(path)
    total = 0
    per = {}
    for label in ("Background:", "Methods:", "Results:", "Conclusions:"):
        # find the paragraph starting with the label
        m = re.search(re.escape(label) + r"(.*?)(?=\n\n|\Z)", block, re.S)
        if not m:
            per[label] = None
            continue
        para = m.group(1)
        para = para.replace("**", "").strip()
        n = len(para.split())
        per[label] = n
        total += n
    return per, total


if __name__ == "__main__":
    p = sys.argv[1] if len(sys.argv) > 1 else CH1
    per, total = count(p)
    for k, v in per.items():
        print(f"  {k:14s} {v}")
    print(f"  TOTAL (labels included, headings excluded) = {total}   limit 350 -> {'OK' if total <= 350 else 'OVER'}")
