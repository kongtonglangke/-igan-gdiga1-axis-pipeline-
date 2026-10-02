# -*- coding: utf-8 -*-
"""Retag the three Figure S1 PNGs from 199 dpi to 300 dpi.

Rationale: every other supplementary PNG carries ~300 dpi metadata (S2-S7: 299 dpi)
while the S1 PNGs alone are tagged 199 dpi, which implies a 53-cm-wide figure and
would trip the journal's raster-resolution check. The pixel data are unchanged; only
the pHYs density chunk is rewritten. Pixel bytes are verified identical before/after.
"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import io, hashlib, os
from PIL import Image

BASE = _paths.at("阶段4.5_复现包/submission/additional_file_1")
FILES = [
    "Figure_S1_GSE73953_5genes_boxplot.png",
    "Figure_S1_GSE115857_5genes_boxplot.png",
    "Figure_S1_GSE93798_5genes_boxplot.png",
]


def px_hash(p):
    with Image.open(p) as im:
        return hashlib.md5(im.convert("RGBA").tobytes()).hexdigest()


for name in FILES:
    p = BASE(name)
    before_dpi = Image.open(p).info.get("dpi")
    h_before = px_hash(p)
    with Image.open(p) as im:
        im.load()
        im.info.pop("dpi", None)
        im.save(p, format="PNG", dpi=(300, 300), optimize=True)
    after_dpi = Image.open(p).info.get("dpi")
    h_after = px_hash(p)
    print(f"{name}\n   dpi {before_dpi} -> {after_dpi} | pixels identical: {h_before == h_after}")
    assert h_before == h_after, "pixel data changed!"
print("done")
