# -*- coding: utf-8 -*-
"""_tmp_patch_h_additional_ise.py — 补充材料图注索引英式拼写统一（2026-09-10）

发现：`Additional_File_1_Figure_Legends.md` 的补充材料表索引中仍有 6 处美式
"colocalization"（v1.5 的 -ise 统一只覆盖了 ch1–ch9，漏了这个直接交付件）。
另修 GA 图注 v03 变更说明里对美式拼写的字面引用。
"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import io, re, os

PKG = _paths.at("阶段4.5_复现包")
P_ADD = PKG("submission", "additional_file_1",
                     "Additional_File_1_Figure_Legends.md")
P_GA = _paths.legacy("阶段4/主图/data_check/ga_caption_v03.md")

# 1) Additional file 1：-ise 统一（大小写各一）
t = io.open(P_ADD, encoding="utf-8", newline="").read()
n1 = t.count("colocalization"); n2 = t.count("Colocalization")
t = t.replace("colocalization", "colocalisation").replace("Colocalization", "Colocalisation")
io.open(P_ADD, "w", encoding="utf-8", newline="").write(t)
print(f"[1] Additional file 1: colocalization {n1} -> colocalisation；Colocalization {n2} -> Colocalisation")
rest = len(re.findall(r"[Cc]olocalization", t))
print(f"    残留美式拼写 = {rest}（应为 0）")

# 2) GA caption v03：变更说明改写，避免字面出现美式拼写
t = io.open(P_GA, encoding="utf-8", newline="").read()
old = "美式 colocalization 统一为英式 colocalisation，与全稿一致"
new = "美式 -ization 拼写统一为英式 -isation，与全稿一致"
if old in t:
    t = t.replace(old, new); print("[2] GA caption v03 变更说明已改写")
else:
    print("[2] GA caption v03 目标串未命中")
io.open(P_GA, "w", encoding="utf-8", newline="").write(t)
print("    GA v03 残留 colocaliz* =", len(re.findall(r"[Cc]olocaliz\w*", t)), "| 残留 L7 =", len(re.findall(r"L7", t)))
print("done")
