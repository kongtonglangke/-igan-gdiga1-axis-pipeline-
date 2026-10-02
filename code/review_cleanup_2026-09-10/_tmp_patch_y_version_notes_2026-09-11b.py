# -*- coding: utf-8 -*-
# Patch Y (2026-09-11, second batch): version notes for the P1-1 / R8-4.1% / P1-2 fixes.
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
from pathlib import Path

BASE = Path(_paths.legacy("阶段4.5_复现包/submission/manuscript"))

def add_note(fname, anchor_startswith, note):
    p = BASE / fname
    raw = p.read_bytes()
    eol = "\r\n" if b"\r\n" in raw else "\n"
    t = raw.decode("utf-8")
    lines = t.split("\r\n" if eol == "\r\n" else "\n")
    hits = [i for i, ln in enumerate(lines) if ln.startswith(anchor_startswith)]
    assert len(hits) == 1, (fname, anchor_startswith, hits)
    lines.insert(hits[0] + 1, note)
    p.write_bytes((eol.join(lines)).encode("utf-8"))
    print(f"{fname}: note inserted after line {hits[0]+1} (eol={eol!r})")

add_note(
    "04_Results.md",
    "> v2.0 (2026-09-11, full-coverage audit):",
    "> v2.1 (2026-09-11, audit fix batch 2): the B-cell detection-rate cross-check in the upstream-regulator subsection now quotes the marker-gated estimate as 4.1% (was 4.2%), matching Table S7 — the weighted detection rate across the gated B-cell states is (737 + 135 + 859) / 42,259 = 4.097%. No other text changes.",
)

add_note(
    "08_Figure_Legends_Tables.md",
    "> v2.1 (2026-09-11, full-coverage audit):",
    "> v2.2 (2026-09-11, audit fix batch 2): Table 5 mQTL standard errors corrected to the random-effects values (cg19603390 0.066, cg19473623 0.014, cg17994788 0.019, cg04827551 0.057, cg16101574 0.011) — the printed SEs had been the fixed-effect (ARE) values while the β column and the cis-mQTL P column are random-effects (MRE); the corrected SEs match Additional file 1: Table S5b (se_mqtl) exactly. Footnote † now states \u201crandom-effects meta-analysis estimates, with the random-effects standard error in parentheses\u201d. Point estimates, P values and conclusions unchanged.",
)
print("ALL PATCH-Y EDITS OK")
