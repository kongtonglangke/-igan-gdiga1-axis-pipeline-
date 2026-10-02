# -*- coding: utf-8 -*-
"""
600 dpi 重渲染包装器（用于补充图 TIFF 交付件）
用法： python _render_600.py <目标脚本.py> <输出目录>
原理：monkeypatch Figure.savefig —— 把 dpi 强制为 600，并把输出文件名改写为
      <stem>_600.png（不覆盖 300 dpi 的投稿 PNG/PDF）。
"""
import os, sys, runpy
import matplotlib
matplotlib.use("Agg")
from matplotlib.figure import Figure

TARGET = sys.argv[1]
OUTDIR = sys.argv[2] if len(sys.argv) > 2 else None
_orig = Figure.savefig


def _patched(self, fname, *a, **kw):
    f = str(fname)
    root, ext = os.path.splitext(f)
    if OUTDIR and not os.path.isabs(f):
        root = os.path.join(OUTDIR, os.path.basename(root))
    kw["dpi"] = 600
    return _orig(self, root + "_600" + ext, *a, **kw)


Figure.savefig = _patched
os.makedirs(OUTDIR, exist_ok=True) if OUTDIR else None
runpy.run_path(TARGET, run_name="__main__")
print("[600dpi] done ->", OUTDIR)
