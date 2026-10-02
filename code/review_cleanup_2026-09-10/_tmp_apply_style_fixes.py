#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_tmp_apply_style_fixes.py — 一次性体例批改（2026-09-10 review v1.0 → v1.1）
  1) 标题/正文：multi-omics causal framework -> multi-omics causal-inference framework
                layered causal framework      -> layered causal-inference framework
  2) 拼写统一：naïve -> naive
  3) 拼写统一：colocalization/Colocalization -> colocalisation/Colocalisation
只作用于 submission/ 下的纯英文稿件章节 + cover letter + data availability；
09_References.md 仅替换稿件标题行（文献条目内的 "colocalisation" 保持原样）。
运行后逐条打印替换计数，任何一项为 0 即报错（防静默未生效）。
"""
import os
import sys

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SUB = os.path.join(BASE, "submission")
MS = os.path.join(SUB, "manuscript")

CH = [f"{i:02d}_" + n for i, n in [
    (1, "Title_Abstract_Keywords.md"),
    (2, "Background.md"),
    (3, "Methods.md"),
    (4, "Results.md"),
    (5, "Discussion_Limitations.md"),
    (6, "Conclusions_Abbreviations.md"),
    (7, "TitlePage_Declarations.md"),
    (8, "Figure_Legends_Tables.md"),
]]

CAUSAL = [
    ("multi-omics causal framework", "multi-omics causal-inference framework"),
    ("layered causal framework", "layered causal-inference framework"),
]
SPELL = [("na\u00efve", "naive")]
COLOC = [("Colocalization", "Colocalisation"), ("colocalization", "colocalisation")]


def apply(path, pairs, tag, require=True):
    with open(path, encoding="utf-8", newline="") as fh:   # newline="" 保留原始行尾
        t = fh.read()
    n_tot = 0
    for a, b in pairs:
        n = t.count(a)
        if n:
            t = t.replace(a, b)
            n_tot += n
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write(t)
    if n_tot == 0 and require:
        raise SystemExit(f"[FAIL] {tag}: 0 replacement in {path}")
    print(f"  {tag:34s} {os.path.basename(path):38s} {n_tot}")
    return n_tot


grand = 0
print("== causal-inference（标题 + 正文）==")
for f in CH:
    grand += apply(os.path.join(MS, f), CAUSAL, "causal-inference", require=False)
apply(os.path.join(MS, "09_References.md"), [CAUSAL[0]], "causal-inference(title-line)")
apply(os.path.join(SUB, "cover_letter_draft.md"), CAUSAL, "causal-inference(cov)")
apply(os.path.join(SUB, "data_availability_draft.md"), CAUSAL, "causal-inference(da)", require=False)

print("== naive（去分音符）==")
for f in CH:
    grand += apply(os.path.join(MS, f), SPELL, "naive", require=False)
apply(os.path.join(MS, "09_References.md"), SPELL, "naive(refs)", require=False)
apply(os.path.join(SUB, "cover_letter_draft.md"), SPELL, "naive(cov)", require=False)
apply(os.path.join(SUB, "data_availability_draft.md"), SPELL, "naive(da)", require=False)

print("== colocalisation（-ise 统一）==")
for f in CH:
    grand += apply(os.path.join(MS, f), COLOC, "colocalisation", require=False)
apply(os.path.join(SUB, "cover_letter_draft.md"), COLOC, "colocalisation(cov)")
apply(os.path.join(SUB, "data_availability_draft.md"), COLOC, "colocalisation(da)", require=False)

print("total replacements:", grand)
