# -*- coding: utf-8 -*-
"""_readme_update_v1.9.py — 包 README v1.8→v1.9、reproducibility README v1.3→v1.4"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import io, os

PKG = _paths.at("阶段4.5_复现包")
P1 = PKG("README.md")
P2 = PKG("reproducibility", "README.md")


def upd(path, pairs, require_all=True):
    t = io.open(path, encoding="utf-8", newline="").read()
    for old, new, label in pairs:
        n = t.count(old)
        if n == 0:
            print(f"  MISS [{label}]  ({n})")
            if require_all:
                continue
            continue
        t = t.replace(old, new)
        print(f"  OK   [{label}] x{n}")
    io.open(path, "w", encoding="utf-8", newline="").write(t)
    return t


print("== 包 README ==")
pairs1 = [
 ("# 阶段 4.5 复现包 / 投稿材料包（v1.8 · 2026-09-10）",
  "# 阶段 4.5 复现包 / 投稿材料包（v1.9 · 2026-09-10）", "标题版本"),

 ("复核：母本裸 `Rn`/`Mn` 均 0、引用 **[1..27] 严格递增**、Abstract 346 ≤ 350。作者终校/建仓回填前。",
  "复核：母本裸 `Rn`/`Mn` 均 0、引用 **[1..27] 严格递增**、Abstract 346 ≤ 350。"
  "**本轮（v1.9）＝作者定夺项办结**（答复 ②否 / ③去 / ④删 / ⑤当日复核；①未答）："
  "**②** 通讯作者 Yan Wang 与参考文献 [8] Wang YN **非同一人**（作者确认“否”）→ 无需在 Competing interests / Acknowledgements 披露，ch7 文本零改动；"
  "**③** GA 图内内部层号 **已去除**——`ga_render.py` 两处副本同步把渲染字符串 `\"upstream circuit (L7)\"` 改为 `\"upstream circuit\"`，默认输出升 **GA v3.3**，重渲染后 `_ga_check.py` → 13 texts / 4 solid boxes / leaks = 0；ch1 GA 版本引用 v3.2→**v3.3**（ch1 → **v1.7**）；GA 图注草案 v02→**v03**（去编号 + -isation 统一）；投稿件 `submission/figures/GA.{png,pdf,tiff}` 已换 v3.3；"
  "**④** `Additional_File_1_Figure_Legends.md` 的 `>` 编辑版本注头块 **已删除**（5 行原样入档 `reproducibility/code/review_cleanup_2026-09-10/Additional_File_1_Figure_Legends_header_notes_archived.md`），该交付件现已无 `>` 头注；"
  "**⑤** [26] Fung **2026-09-10 当日复核仍为 medRxiv 预印本 v1**（posted 2025-10-25，DOI 10.1101/2025.10.25.25338806），未见正式期刊版 → 引用维持 `medRxiv [preprint]. 2025.`；"
  "**顺带修复**：补充材料表索引内 **6 处美式 `colocalization`**（v1.7 的 -isation 统一漏了这个直接交付件）已改为 `colocalisation`/`Colocalisation`。"
  "机检（母本重建 110,371 chars）：全部交付件内部编号/陈旧串 **全 0**、引用 **[1..27] 严格递增**；余留 **① 署名/CRediT 终形态**、**⑥ Funding 待立项**。作者终校/建仓回填前。", "概述追加 v1.9"),

 ("Fig1 v3.6 / Fig2 v3.4 / Fig3 v3.2 / Fig4 v3.2 / Fig5 v3.0 / Fig6 v3.2 / GA v3.2",
  "Fig1 v3.6 / Fig2 v3.4 / Fig3 v3.2 / Fig4 v3.2 / Fig5 v3.0 / Fig6 v3.2 / GA v3.3", "figures 行 GA 版本"),

 ("审稿整改补丁脚本 6 件（可审计，非分析流程）",
  "审稿整改补丁脚本 9 件 + README 升版脚本 4 件 + 图注头块归档件 1 件（可审计，非分析流程）", "目录树脚本件数"),

 ("| Title/Abstract/Keywords | ch1 **v1.6**（审稿清洗：去编号 + 东亚句点名 rs10238682/rs7856182；补漏：字数注去 L7 字样） | Abstract **346**（限 350 ✅） | 全稿头部 |",
  "| Title/Abstract/Keywords | ch1 **v1.7**（审稿清洗：去编号 + 东亚句点名 rs10238682/rs7856182；补漏：字数注去 L7 字样；**GA 版本引用 v3.2→v3.3**） | Abstract **346**（限 350 ✅） | 全稿头部 |", "组件表 ch1"),

 ("""**本轮整改新增的投稿前事项（详见整改记录第六节）**：① 署名次序/CRediT 终形态由作者确认；② 通讯作者 Yan Wang 与 [8] Wang YN 是否同一人待确认；③ GA 图内 "upstream circuit (L7)" 是否同去编号（改一行 + 重渲染）；④ `Additional_File_1_Figure_Legends.md` 头部 `>` 编辑版本注投稿前需删（母本的 `<!-- ASSEMBLY NOTE -->` 同理）；⑤ Funding 若投稿前立项需改资助表述。
4. **[26] Fung 2025 投稿前复核正式发表版**（现为 medRxiv preprint）；
5. Cover letter 定稿签名。""",
  """**本轮整改新增的投稿前事项（整改记录第六 / 八节）**：② 通讯 Yan Wang 与 [8] Wang YN —— 作者已确认**非同一人** → 无需披露 ✅；③ GA 图内 "upstream circuit (L7)" —— 作者决定**去编号**，已重渲染为 **GA v3.3** ✅；④ `Additional_File_1_Figure_Legends.md` 头部 `>` 编辑版本注 —— 作者决定**删**，已删并入档 ✅；**余留**：① 署名次序 / CRediT 终形态待全体作者确认；⑥ Funding 若投稿前立项需改资助表述。（母本 `<!-- ASSEMBLY NOTE -->` 属拼装内部注、不随投稿件提交。）
4. **[26] Fung 2025 —— 2026-09-10 当日复核：仍为 medRxiv 预印本 v1（2025-10-25，DOI 10.1101/2025.10.25.25338806），未见正式期刊版；引用维持 `medRxiv [preprint]. 2025.`**，投稿当日可再查一次；
5. Cover letter 定稿签名。""", "开放项重写"),

 ("v1.8（2026-09-10）：**补漏二次机检**。", "v1.9A_MARKER_v1.8（2026-09-10）：**补漏二次机检**。", "v1.8 锚点临时"),
]

upd(P1, pairs1)

# 在 v1.8 日志段之后插入 v1.9 日志段
t = io.open(P1, encoding="utf-8", newline="").read()
t = t.replace("v1.9A_MARKER_v1.8（2026-09-10）：**补漏二次机检**。",
              "v1.8（2026-09-10）：**补漏二次机检**。", 1)
V19 = ("- v1.9（2026-09-10 傍晚）：**作者定夺项办结 + 补充材料拼写修复**。作者答复 ②否 / ③去 / ④删 / ⑤当日复核（①未答）。"
       "②**Yan Wang ≠ [8] Wang YN**（非同一人）→ 无需披露，ch7 零改动；"
       "③**GA 去内部编号**：`ga_render.py`（`阶段4/主图/scripts/` 与 `reproducibility/code/fig_scripts/` 两处副本）渲染字符串去 `\"(L7)\"`、默认 `OUT_STEM`→`GA_v3.3`，重渲染 + `_ga_check.py` 全 0，投稿件 `submission/figures/GA.*` 换 v3.3（v3.0/v3.1/v3.2 留档），ch1→**v1.7**（GA 版本引用 v3.2→v3.3），GA 图注草案→**v03**；"
       "④**删补充材料 `>` 头注**：`Additional_File_1_Figure_Legends.md` 头块 5 行入档 `review_cleanup_2026-09-10/Additional_File_1_Figure_Legends_header_notes_archived.md` 后从交付件删除；"
       "⑤**[26] Fung 复核**：仍 medRxiv 预印本 v1（2025-10-25），引用不变；"
       "⑥**顺带修**：补充材料表索引 6 处美式 `colocalization`→`colocalisation`/`Colocalisation`；"
       "⑦母本重建 **110,371 chars**，机检全部交付件内部编号/陈旧串 0、引用 [1..27] 严格递增；整改记录升 **v1.3**（第八节）。\n")

anchor = "- v1.6（2026-09-10）：**作者四项决策全部落实**。"
assert anchor in t
t = t.replace(anchor, V19 + anchor, 1)
io.open(P1, "w", encoding="utf-8", newline="").write(t)
print("  OK   [版本记录 v1.9]")

# ---------------------------------------------------------------- repro README
print("\n== reproducibility README ==")
pairs2 = [
 ("| `code/review_cleanup_2026-09-10/` | **审稿整改补丁脚本 6 件**（体例统一、去内部编号、定点改写、L7 数字同步、Additional file 图注修正、**裸 Rn 补漏**）；仅供审计，不参与分析流程 | — |",
  "| `code/review_cleanup_2026-09-10/` | **审稿整改补丁脚本 9 件**（体例统一、去内部编号、定点改写、L7 数字同步、Additional file 图注修正、**裸 Rn 补漏**、**GA 去编号**、**作者决策 ②③④ 落实**、**补充材料 -isation**）+ README 升版脚本 4 件 + `Additional_File_1_Figure_Legends_header_notes_archived.md`（补充材料头块归档）；仅供审计，不参与分析流程 | — |", "目录映射行"),

 ("| `code/fig_scripts/` | **主图 Fig1–Fig6 + GA 冻结件渲染脚本**（`style_common.py` 统一样式；PNG300 + PDF + TIFF600） | 投稿图件 |",
  "| `code/fig_scripts/` | **主图 Fig1–Fig6 + GA 冻结件渲染脚本**（`style_common.py` 统一样式；PNG300 + PDF + TIFF600；**GA 默认输出 v3.3**，2026-09-10 去内部层号） | 投稿图件 |", "fig_scripts 行"),
]
upd(P2, pairs2)

t2 = io.open(P2, encoding="utf-8", newline="").read()
old_f = "*v1.3 · 2026-09-10 · 阶段4.5 复现包* — v1.2→v1.3"
new_f = ("*v1.4 · 2026-09-10（傍晚）· 阶段4.5 复现包* — v1.3→v1.4：`fig_scripts/ga_render.py` 渲染字符串去内部层号 `\"(L7)\"`、默认输出升 **GA_v3.3**（`阶段4/主图/scripts/` 与 `code/fig_scripts/` 两处副本同步；`_ga_check.py` 全 0）；`review_cleanup_2026-09-10/` 补登 `_tmp_patch_f_ga_delL7.py`、`_tmp_patch_g_author_decisions.py`、`_tmp_patch_h_additional_ise.py`、`_tmp_report_v13.py` 及补充材料头块归档件（作者决策 ②③④ 落实）；本文件的目录映射「对应稿件」列仍沿用内部编号作定位索引（投稿正文已不使用）。\n" + old_f)
assert old_f in t2
t2 = t2.replace(old_f, new_f, 1)
io.open(P2, "w", encoding="utf-8", newline="").write(t2)
print("  OK   [repro 版本行 v1.4]")
print("\ndone")
