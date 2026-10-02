# -*- coding: utf-8 -*-
"""
图表审查整改 · 第一批：脚本级字符串修正 + 双副本同步
2026-09-10

修补项（逐条计数校验）：
  A1 fig_scripts/fig2_render.py   : "gray = no cis-eQTL record" -> "grey = ..."
  A2 fig_scripts/fig3_render.py   : y 轴 "colocalization" -> "colocalisation"
  A3 fig_scripts/fig4_render.py   : panel c 图例移出受压区（右下 -> 左上空白区）
  A4 fig_scripts/fig6_render.py   : x 轴括号枚举顺序对齐列序（vorinostat/EGCG 颠倒）
  A5 阶段1/M2a_eQTL语境/plot_eqtl_5genes.py : 标题去内部代码 "M2a" + gray->grey
  A6 reproducibility/code/M2a_eqtl/plot_eqtl_5genes.py : 同上
  A7 reproducibility/code/additional_file_1/figure_s3_...py : 标题去内部标签 "primary evidence 1"
  B1 阶段4/主图/scripts/{fig3,fig4,fig6}_render.py <- 由 fig_scripts 覆盖（消除副本分叉）
"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import io, os, shutil, sys

ROOT = _paths.LEGACY
FS   = ROOT(r"阶段4.5_复现包\reproducibility\code\fig_scripts")
MSTD = ROOT(r"阶段4\主图\scripts")

def patch(path, pairs):
    t = io.open(path, encoding="utf-8", newline="").read()
    for old, new, exp, lab in pairs:
        n = t.count(old)
        if n != exp:
            print(f"   MISS [{lab}] count={n} expect={exp}")
            continue
        t = t.replace(old, new)
        print(f"   OK   [{lab}] x{n}")
    io.open(path, "w", encoding="utf-8", newline="").write(t)

print("== A1 fig2_render.py ==")
patch(os.path.join(FS, "fig2_render.py"), [
    ("gray = no cis-eQTL record", "grey = no cis-eQTL record", 1, "fig2 gray->grey"),
])

print("== A2 fig3_render.py ==")
patch(os.path.join(FS, "fig3_render.py"), [
    ('"PPH4 (colocalization posterior)"', '"PPH4 (colocalisation posterior)"', 1, "fig3c ylabel -ise"),
])

print("== A3 fig4_render.py ==")
patch(os.path.join(FS, "fig4_render.py"), [
    ('ax.legend(loc="lower right", bbox_to_anchor=(0.995, 0.02),',
     'ax.legend(loc="upper left", bbox_to_anchor=(0.02, 0.98),', 1, "fig4c legend 移左上"),
])

print("== A4 fig6_render.py ==")
patch(os.path.join(FS, "fig6_render.py"), [
    ('"Interventions (5-azacytidine / vorinostat / EGCG / research / sialyl-Ti)"',
     '"Interventions (5-azacytidine / EGCG / vorinostat / research / sialyl-Ti)"', 1, "fig6b xlabel 顺序"),
])

S2_OLD = ('"M2a  cis-eQTL significance [-log10(p_beta)] for O-glycosylation genes\\n'
          '(red depth = stronger eQTL; gray = no record in dataset)"')
S2_NEW = ('"cis-eQTL significance [\u2212log10(p_beta)] for O-glycosylation genes\\n'
          '(red depth = stronger eQTL; grey = no record in dataset)"')
print("== A5/A6 S2 plot_eqtl_5genes.py (两副本) ==")
for p in [ROOT(r"阶段1\M2a_eQTL语境\plot_eqtl_5genes.py"),
          ROOT(r"阶段4.5_复现包\reproducibility\code\M2a_eqtl\plot_eqtl_5genes.py")]:
    print(" --", p)
    patch(p, [(S2_OLD, S2_NEW, 1, "S2 title 去 M2a")])

print("== A7 S3 脚本 ==")
patch(ROOT(r"阶段4.5_复现包\reproducibility\code\additional_file_1\figure_s3_lead_eQTL_full_matrix_heat.py"), [
    ('"cis-eQTL attribution of the four Gd-IgA1 lead SNPs across contexts\\n'
     '(primary evidence 1: B-cell/plasma-cell mechanism)"',
     '"cis-eQTL attribution of the four Gd-IgA1 lead SNPs across contexts\\n'
     '(B-cell / plasma-cell mechanism)"', 1, "S3 title 去内部标签"),
])

print("== B1 同步 fig3/4/6 副本 (fig_scripts -> 阶段4/主图/scripts) ==")
for f in ["fig3_render.py", "fig4_render.py", "fig6_render.py"]:
    src, dst = os.path.join(FS, f), os.path.join(MSTD, f)
    shutil.copyfile(src, dst)
    same = io.open(src, "rb").read() == io.open(dst, "rb").read()
    print(f"   OK   {f} 同步 -> {same}")

print("\nDONE")
