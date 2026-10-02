#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""_readme_update_v1.8.py — 包 README 升 v1.8（补漏 8 处裸 Rn + 词数微调）"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os

ROOT = _paths.at("阶段4.5_复现包")


def patch(fn, pairs):
    p = ROOT(fn)
    with open(p, encoding="utf-8", newline="") as fh:
        t = fh.read()
    for old, new, exp in pairs:
        n = t.count(old)
        if n != exp:
            raise SystemExit(f"[FAIL] {fn}: found {n}, expected {exp}\n  >> {old[:150]}")
        t = t.replace(old, new)
    with open(p, "w", encoding="utf-8", newline="") as fh:
        fh.write(t)
    print(f"OK  {fn}  ({len(pairs)} pairs)")


patch("README.md", [
    ("# 阶段 4.5 复现包 / 投稿材料包（v1.7 · 2026-09-10）",
     "# 阶段 4.5 复现包 / 投稿材料包（v1.8 · 2026-09-10）", 1),

    ("**Abstract 346 ≤ 350**。作者终校/建仓回填前。",
     "**Abstract 346 ≤ 350**。**本轮（v1.8）＝补漏二次机检**：用更宽正则复扫发现 v1.1 的"
     "\u201c内部编号残留 0\u201d口径偏窄（只检括号式 `(R2)`、**漏检裸写法**），**补清投稿正文内 8 处裸 `Rn`**"
     "（ch3\u201creported in R2\u201d、ch4\u201cthe marker-gated partition of R2\u201d、ch8 表注\u201cCorresponds to Rn and Fig. X\u201d×6）；"
     "同时修正机检正则口径（须**同时覆盖括号式与裸式**、并**允许分隔符后空格**——后者曾致母本引用序假告警）。"
     "复核：母本裸 `Rn`/`Mn` 均 0、引用 **[1..27] 严格递增**、Abstract 346 ≤ 350。作者终校/建仓回填前。", 1),

    ("    │   └── review_cleanup_2026-09-10/  审稿整改补丁脚本 5 件（可审计，非分析流程）",
     "    │   └── review_cleanup_2026-09-10/  审稿整改补丁脚本 6 件（可审计，非分析流程）", 1),

    ("ch1 **v1.5**（审稿清洗：去编号 + 东亚句点名 rs10238682/rs7856182）",
     "ch1 **v1.6**（审稿清洗：去编号 + 东亚句点名 rs10238682/rs7856182；补漏：字数注去 L7 字样）", 1),

    ("ch3 **v1.8**（16 小节；去编号 + 单细胞 QC 澄清 + **L4 FDR 敏感性** + 共定位先验依据 + decoupleR/IVW 口径） | **3,593 含 / 3,488 不含**",
     "ch3 **v1.9**（16 小节；去编号 + 单细胞 QC 澄清 + **L4 FDR 敏感性** + 共定位先验依据 + decoupleR/IVW 口径；补漏：去裸 R2） | **3,592 含 / 3,487 不含**", 1),

    ("ch4 **v1.7**（去编号 + **Q1–Q4 挂标** + R1 量化 + 双分区检出率澄清；**L7 零分布改 1,000 组重跑**） | **3,203 含 / 3,086 不含**（8 小节）",
     "ch4 **v1.8**（去编号 + **Q1–Q4 挂标** + R1 量化 + 双分区检出率澄清；**L7 零分布改 1,000 组重跑**；补漏：去裸 R2） | **3,201 含 / 3,084 不含**（8 小节）", 1),

    ("ch8 **v1.7**（图注去编号；Fig. 1c 保留 (L0–L7) 图内索引；Fig. S6 图注同步 1,000 组零分布）",
     "ch8 **v1.8**（图注去编号；Fig. 1c 保留 (L0–L7) 图内索引；Fig. S6 图注同步 1,000 组零分布；补漏：6 处表注去裸 Rn）", 1),

    ("**9,879 含 / 9,649 不含**（= BG 762/760 + Methods 3,593/3,488 + Results 3,203/3,086 + Disc 1,568/1,566 + Lim 605/603 + Concl 148/146；2026-09-10 整改后统一脚本实测）",
     "**9,876 含 / 9,646 不含**（= BG 762/760 + Methods 3,592/3,487 + Results 3,201/3,084 + Disc 1,568/1,566 + Lim 605/603 + Concl 148/146；2026-09-10 整改 + 补漏后统一脚本实测）", 1),

    ("附：本稿正文合计 **9,879** 含小节标题",
     "附：本稿正文合计 **9,876** 含小节标题", 1),

    ("v1.5（2026-09-10）：**审稿整改后重建**——全稿去内部编号、去分音符、colocalisation 英式统一、标题 causal→causal-inference；装配说明同步（去编号状态 + 补充表后缀 + 期刊字数结论）。重建后 **[1..27] 严格递增无缺口**，正文内部编号残留仅 Fig. 1c 图例的 (L0–L7)。",
     "v1.5（2026-09-10）：**审稿整改后重建**——全稿去内部编号、去分音符、colocalisation 英式统一、标题 causal→causal-inference；装配说明同步（去编号状态 + 补充表后缀 + 期刊字数结论）。重建后 **[1..27] 严格递增无缺口**，正文内部编号残留仅 Fig. 1c 图例的 (L0–L7)。\n"
     "v1.6（2026-09-10）：**补漏重建**——清除上一轮遗漏的 8 处**裸 `Rn`**（ch3/ch4 正文各 1 处、ch8 表注 6 处）；重建后正文 109,804 chars，母本裸 `Rn`/`Mn` = 0。", 1),

    ("⑬ 整改补丁脚本 5 件归档至 `reproducibility/code/review_cleanup_2026-09-10/`。",
     "⑬ 整改补丁脚本 5 件归档至 `reproducibility/code/review_cleanup_2026-09-10/`。\n"
     "- v1.8（2026-09-10）：**补漏二次机检**。用更宽正则（同时覆盖括号式与裸式）复扫全文，发现 v1.7 的\u201c内部编号残留 0\u201d结论口径偏窄——只检了 `(R2)` 一类**括号式**，**漏检裸写法**。补清投稿正文内 **8 处裸 `Rn`**：ch3\u201creported in R2\u201d→\u201creported above\u201d、ch4\u201cthe marker-gated partition of R2\u201d→\u201cthe marker-gated partition\u201d、ch8 六处表注\u201cCorresponds to Rn and Fig. X\u201d→\u201cCorresponds to Fig. X.\u201d；ch1 字数注\u201cafter the L7 integration\u201d→\u201cafter the upstream-regulator integration\u201d（该行不进母本）。**保留不动**（经核验不在母本、拼装时自动剔除）：ch4 集成附录\u201cFigure calls in text\u201d行、`reproducibility/data/M1_geo … M6_intervention` 数据目录路径、各章 `>` 头部修订注；Fig. 1c 图例 `(L0–L7)` 仍保留（面板自印）。版本：ch1 **v1.6**、ch3 **v1.9**、ch4 **v1.8**、ch8 **v1.8**；补丁脚本 `_tmp_patch_e_residual_Rn.py`（归档第 6 件）；整改记录升 **v1.2**（第七节）。**机检正则口径修正**：须①覆盖括号式+裸式；②允许分隔符后空格（旧正则曾致 `[6, 7, 8]` 漏匹配、引用序**假告警**）。复核：母本 `evidence layer L*`/`(Rn)`/裸 `Rn`/裸 `Mn`/`(Lx)` 全 0、引用 **[1..27] 严格递增无缺口**、`(L0–L7)` 保留 1。词数微调：Methods **3,592/3,487**、Results **3,201/3,084**、合计 **9,876/9,646**；Abstract **346 不变**。", 1),
])

print("done README.md")
