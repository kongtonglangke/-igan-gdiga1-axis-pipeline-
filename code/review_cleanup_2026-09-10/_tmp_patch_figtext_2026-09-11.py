# -*- coding: utf-8 -*-
# _tmp_patch_figtext_2026-09-11.py — 2026-09-11 全维度审查修复批次（图件文字层）
# P0 药物归属延伸至图内标签：正文/图注已将「5-azacytidine 逆转 C1GALT1C1 启动子甲基化」
# 更正为「5-aza-2′-deoxycytidine (decitabine) [20]；5-azacytidine 仅限 IL-17 表达恢复 [19]」，
# 图内标签同步为类药物名（5-aza DNMTi / 5-aza class / decitabine），避免单药误归属。
# 同步修改两处脚本副本（fig_scripts 为准 + 阶段4/主图/scripts），改后强制 md5 一致。
# 注：脚本内 5-azacytidine 字样的「历史变更注释」行保留不改（属日志性质）。
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os, hashlib, shutil

ROOT = _paths.LEGACY
FIGSCRIPTS = ROOT(r"阶段4.5_复现包\reproducibility\code\fig_scripts")
STAGE4 = ROOT(r"阶段4\主图\scripts")

EDITS = {
"fig1_render.py": [
 (r'"ACTIONABLE WINDOW\n(C1GALT1C1 promoter; 5-azacytidine)"',
  r'"ACTIONABLE WINDOW\n(C1GALT1C1 promoter; 5-aza DNMTi)"'),
 (r'"C1GALT1C1=6.0 > C1GALT1=4.5; 5-azacytidine [20]"',
  r'"C1GALT1C1=6.0 > C1GALT1=4.5; 5-aza DNMTi [19, 20]"'),
],
"fig6_render.py": [
 ('("DNMTi",        "5-azacytidine"),',
  '("DNMTi",        "5-aza class"),'),
 ('"Interventions (5-azacytidine / EGCG / vorinostat / research / sialyl-Ti)"',
  '"Interventions (5-aza class / EGCG / vorinostat / research / sialyl-Ti)"'),
 ('"DNMTi (5-azacytidine)",',
  '"DNMTi (5-aza class)",'),
 ('"reversible by 5-azacytidine [20]"],',
  '"reversible by decitabine [20]"],'),
],
"ga_render.py": [
 (r'"5-azacytidine \u00b7 EGCG",',
  r'"5-aza DNMTi \u00b7 EGCG",'),
 (r'"(IL-4/IL-17 induce C1GALT1C1 promoter hypermethylation;\nreversible by DNMTi; X chr caveat applies)"',
  r'"(IL-4 induces C1GALT1C1 promoter hypermethylation;\nreversible by DNMTi; X chr caveat applies)"'),
],
}

def md5(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()

def patch(path, pairs):
    with open(path, "r", encoding="utf-8", newline="") as f:
        text = f.read()
    for i, (old, new) in enumerate(pairs):
        n = text.count(old)
        assert n == 1, f"{'MISS' if n==0 else 'DUP'} {os.path.basename(path)} #{i}: {old[:60]!r} (count={n})"
        text = text.replace(old, new)
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write(text)
    print(f"OK  {os.path.basename(path)} @ fig_scripts ({len(pairs)} edits)")

if __name__ == "__main__":
    for name, pairs in EDITS.items():
        patch(os.path.join(FIGSCRIPTS, name), pairs)
    # 同步到 阶段4/主图/scripts 并验证 md5 一致
    for name in EDITS:
        src = os.path.join(FIGSCRIPTS, name)
        dst = os.path.join(STAGE4, name)
        shutil.copy2(src, dst)
        assert md5(src) == md5(dst), f"md5 mismatch after sync: {name}"
        print(f"SYNC {name} -> 阶段4/主图/scripts (md5 {md5(src)[:8]} OK)")
    print("ALL FIG-SCRIPT EDITS APPLIED")
