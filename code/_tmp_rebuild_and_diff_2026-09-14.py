# -*- coding: utf-8 -*-
"""备份母本 -> 重建 -> 逐行 diff（只看是否只改了预期位置）"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import difflib
import os
import shutil

BASE = _paths.at("阶段4.5_复现包")
OUT = BASE("submission", "manuscript", "Manuscript_GM_fulltext.md")
TMP = BASE("reproducibility", "code", "_tmp_pre_author_edit.md")

shutil.copyfile(OUT, TMP)
old = open(OUT, encoding="utf-8").read().split("\n")

os.system(
    '"C:/Users/user/.workbuddy/binaries/python/versions/3.13.12/python.exe" "'
    + BASE("reproducibility", "code", "_build_submission_fulltext.py")
    + '"'
)

new = open(OUT, encoding="utf-8").read().split("\n")
d = list(difflib.unified_diff(old, new, "before", "after", lineterm="", n=1))
print("diff hunk lines:", len(d))
for l in d:
    print(l[:220])
print("old chars:", sum(len(x) + 1 for x in old), "new chars:", sum(len(x) + 1 for x in new))
os.remove(TMP)
