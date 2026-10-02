# -*- coding: utf-8 -*-
"""
图表审查整改 · 第二批：补充材料图件交付同步
2026-09-10

背景：第一批(A1-A7)只修了脚本源码；其中 S2/S3 的标题修正未重新渲染到交付件，
      且 S2 交付文件名仍含内部模块代码 "M2a"。本批把修正后的渲染落盘到投稿目录，
      并把补充材料 TIFF 统一到"与 PNG 同像素、600 dpi 元数据、RGBA、LZW"的既有约定。

动作：
  C1 S2: 用修正后的渲染(eQTL_5genes_heatmap.png/pdf)覆盖交付件，文件名去 "M2a"
         Figure_S2_M2a_5x5_eQTL_heatmap.*  ->  Figure_S2_5x5_eQTL_heatmap.*
  C2 S3: 由当前 300 dpi PNG 重新生成交付 TIFF（旧 TIFF 为改标题前的版本）
  C3 S4/S7: TIFF dpi 元数据 300 -> 600（与 S1/S2/S5/S6 统一，仅元数据不变像素）
  C4 清理：旧 S2 文件、_600 临时件、空的 _hi600/
"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import io, os, shutil
from PIL import Image

ROOT = _paths.LEGACY
AF   = ROOT(r"阶段4.5_复现包\submission\additional_file_1")
S2SRC= ROOT(r"阶段1\M2a_eQTL语境\figures")
CLEAN= ROOT(r"阶段4.5_复现包\reproducibility\code\review_cleanup_2026-09-10")

def save_tiff_from_png(png_path, tiff_path, dpi=600):
    im = Image.open(png_path).convert("RGBA")
    im.save(tiff_path, format="TIFF", compression="tiff_lzw", dpi=(dpi, dpi))
    chk = Image.open(tiff_path)
    print(f"   OK tiff {os.path.basename(tiff_path)} size={chk.size} dpi={chk.info.get('dpi')} mode={chk.mode} comp={chk.info.get('compression')}")

print("== C1  S2  交付件刷新 + 去内部代码 M2a ==")
s2_new = os.path.join(AF, "Figure_S2_5x5_eQTL_heatmap")
shutil.copyfile(os.path.join(S2SRC, "eQTL_5genes_heatmap.png"), s2_new + ".png")
shutil.copyfile(os.path.join(S2SRC, "eQTL_5genes_heatmap.pdf"), s2_new + ".pdf")
print("   OK png/pdf -> Figure_S2_5x5_eQTL_heatmap")
save_tiff_from_png(s2_new + ".png", s2_new + ".tiff", 600)
for ext in ("png", "pdf", "tiff"):
    old = os.path.join(AF, f"Figure_S2_M2a_5x5_eQTL_heatmap.{ext}")
    if os.path.exists(old):
        os.remove(old); print(f"   removed old Figure_S2_M2a_5x5_eQTL_heatmap.{ext}")

print("== C2  S3  交付 TIFF 重生成（旧版为改标题前） ==")
s3_png = os.path.join(AF, "Figure_S3_lead_eQTL_full_matrix_heat.png")
save_tiff_from_png(s3_png, os.path.join(AF, "Figure_S3_lead_eQTL_full_matrix_heat.tiff"), 600)

print("== C3  S4/S7 TIFF dpi 统一为 600（仅元数据） ==")
for name in ("Figure_S4_Bcell_C1GALT1_trajectory", "Figure_S7_Bcell_eQTL_heterogeneity"):
    t = os.path.join(AF, name + ".tiff")
    im = Image.open(t).convert("RGBA")
    im.save(t, format="TIFF", compression="tiff_lzw", dpi=(600, 600))
    chk = Image.open(t)
    print(f"   OK {name}.tiff size={chk.size} dpi={chk.info.get('dpi')}")

print("== C4  清理临时件 ==")
for p in [
    os.path.join(S2SRC, "eQTL_5genes_heatmap_600.png"),
    os.path.join(S2SRC, "eQTL_5genes_heatmap_600.pdf"),
    os.path.join(S2SRC, "eQTL_5genes_forest_600.png"),
    os.path.join(S2SRC, "eQTL_5genes_forest_600.pdf"),
    os.path.join(AF, "Figure_S3_lead_eQTL_full_matrix_heat_600.png"),
    os.path.join(AF, "Figure_S3_lead_eQTL_full_matrix_heat_600.pdf"),
]:
    if os.path.exists(p):
        os.remove(p); print("   removed", os.path.relpath(p, ROOT))
hi = os.path.join(CLEAN, "_hi600")
if os.path.isdir(hi):
    try:
        os.rmdir(hi); print("   removed empty _hi600/")
    except OSError as e:
        print("   _hi600 not empty:", e)

print("\nDONE")
