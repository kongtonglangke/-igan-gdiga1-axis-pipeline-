# -*- coding: utf-8 -*-
"""_tmp_patch_f_ga_delL7.py — 作者决策③：GA 图内 "(L7)" 字样去编号（v3.2 → v3.3）

两处副本同步：
  A = 阶段4/主图/scripts/ga_render.py（canonical，OUT_DIR=阶段4/主图/out）
  B = 阶段4.5_复现包/reproducibility/code/fig_scripts/ga_render.py（交付副本）
改动：
  1) 渲染字符串 "upstream circuit (L7)\n..." -> "upstream circuit\n..."
  2) 默认 OUT_STEM "GA_v3.2" -> "GA_v3.3"
  3) docstring 追加 v3.3 记录
"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import io, os, sys

ROOT = _paths.LEGACY
FILES = [
    ROOT(r"阶段4\主图\scripts\ga_render.py"),
    ROOT(r"阶段4.5_复现包\reproducibility\code\fig_scripts\ga_render.py"),
]

OLD_TEXT = '"upstream circuit (L7)\\ncytokine \\u2192 DNMT1 \\u2192 Sp/KLF",'
NEW_TEXT = '"upstream circuit\\ncytokine \\u2192 DNMT1 \\u2192 Sp/KLF",'

OLD_STEM = 'os.environ.get("OUT_STEM", "GA_v3.2")'
NEW_STEM = 'os.environ.get("OUT_STEM", "GA_v3.3")'

DOC_ANCHOR = ('     对应正文 Results R8 / Fig. 1c 第七行 L7；图注草案同步补此句。\n"""')
DOC_NEW = ('     对应正文 Results R8 / Fig. 1c 第七行 L7；图注草案同步补此句。\n'
           'v3.3 2026-09-10 作者决策③「去编号」：图内注释条文字由\n'
           '     "upstream circuit (L7)" 改为 "upstream circuit"（去掉内部层号，\n'
           '     与投稿正文全面去内部编号的口径一致）；几何、配色、箭头均不变。\n'
           '     其余元素与 v3.2 完全一致。v3.0/v3.1/v3.2 留档不删。\n"""')

for p in FILES:
    t = io.open(p, encoding="utf-8", newline="").read()
    n = 0
    if OLD_TEXT in t:
        t = t.replace(OLD_TEXT, NEW_TEXT); n += 1
        print("  OK  text-string  <- %s" % os.path.basename(os.path.dirname(p)))
    else:
        print("  MISS text-string:", p)
    if OLD_STEM in t:
        t = t.replace(OLD_STEM, NEW_STEM)
        print("  OK  OUT_STEM")
    else:
        print("  MISS OUT_STEM:", p)
    if "v3.3 2026-09-10 作者决策③" not in t and DOC_ANCHOR in t:
        t = t.replace(DOC_ANCHOR, DOC_NEW)
        print("  OK  docstring v3.3")
    else:
        print("  SKIP docstring (already or anchor miss)")
    io.open(p, "w", encoding="utf-8", newline="").write(t)
    print("  saved:", p, "\n")
print("done")
