# -*- coding: utf-8 -*-
# Patch AB (2026-09-11): 由旧版留档件生成五章双语版新文件（ch1 v2.0 / ch3 v1.11 / ch4 v2.1 / ch5 v1.9 / ch6 v1.7）。
# 头部版本/配套名/引用状态/词数口径更新，并插入对应「同步记录」；§1 英文随后由 sync_初稿双语版.py 注入。
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
from pathlib import Path

DIR = Path(_paths.legacy("阶段4/初稿"))

def build(src_name, dst_name, reps, record):
    t = (DIR / src_name).read_text(encoding="utf-8")
    for old, new in reps:
        assert t.count(old) == 1, (dst_name, old[:50], t.count(old))
        t = t.replace(old, new)
    # 在头部 `---` 之后、第一个既有 `## ` 标题之前插入新同步记录
    i = t.find("\n## ")
    assert i != -1, dst_name
    t = t[:i + 1] + record + "\n\n---\n\n" + t[i + 1:]
    (DIR / dst_name).write_text(t, encoding="utf-8", newline="\n")
    print(f"built {dst_name} ({len(t)} chars)")

# ---------- ch1 v2.0 ----------
build(
 "章节1_Title_Abstract_Keywords_v1.9_中英对照注释版.md",
 "章节1_Title_Abstract_Keywords_v2.0_中英对照注释版.md",
 [
  ("# 章节 1（v1.9）Title page + Abstract + Keywords + Graphical abstract — 中英对照注释版",
   "# 章节 1（v2.0）Title page + Abstract + Keywords + Graphical abstract — 中英对照注释版"),
  ("`章节1_Title_Abstract_Keywords_v1.9_纯英文版.md`",
   "`章节1_Title_Abstract_Keywords_v2.0_纯英文版.md`"),
  ("> **版本**：**v1.9（2026-09-10 全量对齐重做）**｜历史版本 v1.0 / v1.1 / v1.2 留档不删",
   "> **版本**：**v2.0（2026-09-11 全维度审查批次 1）**｜历史版本 v1.0 / v1.1 / v1.2 / v1.9 留档不删"),
 ],
 """## v2.0 同步记录（2026-09-11 · 全维度审查批次 1）

- 包内 ch1 升 **v2.0**（版本注记定格；Abstract 词数 **348** 不变、四段结构与结论不变）。
- §1 英文定稿随包内 `01_Title_Abstract_Keywords.md`（v2.0）逐字节注入；中文对照层无需改动。
""",
)

# ---------- ch3 v1.11 ----------
build(
 "章节3_Methods_v1.10_中英对照注释版.md",
 "章节3_Methods_v1.11_中英对照注释版.md",
 [
  ("# 章节 3（v1.10）Methods — 中英对照注释版",
   "# 章节 3（v1.11）Methods — 中英对照注释版"),
  ("`章节3_Methods_v1.10_纯英文版.md`",
   "`章节3_Methods_v1.11_纯英文版.md`"),
  ("> **版本**：**v1.10（2026-09-10 全量对齐重做）**｜历史版本 v1.0–v1.3 留档不删",
   "> **版本**：**v1.11（2026-09-11 全维度审查批次 1）**｜历史版本 v1.0–v1.3 / v1.10 留档不删"),
  ("> **引用状态**：本章引入 [6]–[24] 中属 Methods 独占的条目（[12] 起首次出现），见 §3 引用映射；[25]–[27] 首现于 Results/Discussion",
   "> **引用状态**：本章引入 [6]–[25] 中属 Methods 独占的条目（[12] 起首次出现；[24] decoupleR、[25] JASPAR 均首现于本章），见 §3 引用映射；[26]–[29] 首现于 Results/Discussion"),
  ("> **词数口径**：原始空白分词；**3,626** 含 16 个小节标题 / **3,521** 不含（Markdown 分隔符、表格行、编辑头块、附录均排除）",
   "> **词数口径**：原始空白分词；**3,637** 含 16 个小节标题 / **3,532** 不含（Markdown 分隔符、表格行、编辑头块、附录均排除；严格全 Table 1 剔除口径）"),
 ],
 """## v1.11 同步记录（2026-09-11 · 全维度审查批次 1）

- 包内 ch3 升 **v1.11**：① Table 1 GoDMC 行更正为 **24,988–28,181 per CpG**（对齐 Table S5a 实测样本量）；② 文献锚句恢复首现序 **[19]→[20]**（IL-17/5-azacytidine 在先、IL-4/5-aza-2′-deoxycytidine 在后）；③ GSE285335 化学版本表述去版本号（\"10x Genomics 5′ gene-expression\"）；④ 引文编号同步至 29 条新序（[24] decoupleR 插入）。
- 词数实测 **3,637 / 3,532**（+47；口径不变）。
- §1 英文定稿随包内 `03_Methods.md`（v1.11）逐字节注入。
""",
)

# ---------- ch4 v2.1 ----------
build(
 "章节4_Results_v1.9_中英对照注释版.md",
 "章节4_Results_v2.1_中英对照注释版.md",
 [
  ("# 章节 4（v1.9）Results — 中英对照注释版",
   "# 章节 4（v2.1）Results — 中英对照注释版"),
  ("`章节4_Results_v1.9_纯英文版.md`",
   "`章节4_Results_v2.1_纯英文版.md`"),
  ("> **版本**：**v1.9（2026-09-10 全量对齐重做）**｜历史版本 v1.0–v1.2（R1–R7 七小节）留档不删",
   "> **版本**：**v2.1（2026-09-11 全维度审查批次 1+2）**｜历史版本 v1.0–v1.2（R1–R7 七小节）/ v1.9 留档不删"),
  ("> **引用状态**：本章引入 [25]（血清总 IgA 位点）与 [26]（Fung 等）；[22] ImmuNexUT 在 Methods 已引入。GEO 系列以访问号引用，不挂编号",
   "> **引用状态**：本章引入 [26]（血清总 IgA 位点）、[27]（Fung 等）与 [28]（Kim，GSE285335 归属）；[22] ImmuNexUT 在 Methods 已引入。GEO 系列以访问号引用，不挂编号"),
  ("> **词数口径**：原始空白分词、排除分隔符/表格行/计划表附录 → **3,211** 含 8 个小节标题 / **3,094** 不含",
   "> **词数口径**：原始空白分词、排除分隔符/表格行/计划表附录 → **3,292** 含 8 个小节标题 / **3,175** 不含"),
 ],
 """## v2.1 同步记录（2026-09-11 · 全维度审查批次 1+2）

- 批次 1（ch4 v2.0）：PPH1 句重构为实测两组（0.79–0.91 vs 0.18–0.29；Table 4）；去悬垂 \"(Fig. 3a)\"；Fung 段补 Wang 2021 候选层张力句；去甲基化证据句按药物拆分（5-aza-2′-deoxycytidine [20]＝IL-4 甲基化逆转；5-azacytidine [19]＝IL-17 表达恢复），临床先例引文扩为 [19, 20]；编号同步 29 条新序。
- 批次 2（ch4 v2.1）：R8 双分区检出率 cross-check 由 4.2% 更正为 **4.1%**（Table S7 加权 (737+135+859)/42,259=4.097%）。
- 词数实测 **3,292 / 3,175**；§1 英文定稿随包内 `04_Results.md`（v2.1）逐字节注入。
""",
)

# ---------- ch5 v1.9 ----------
build(
 "章节5_Discussion_Limitations_v1.8_中英对照注释版.md",
 "章节5_Discussion_Limitations_v1.9_中英对照注释版.md",
 [
  ("# 章节 5（v1.8）Discussion + Limitations — 中英对照注释版",
   "# 章节 5（v1.9）Discussion + Limitations — 中英对照注释版"),
  ("`章节5_Discussion_Limitations_v1.8_纯英文版.md`",
   "`章节5_Discussion_Limitations_v1.9_纯英文版.md`"),
  ("> **版本**：**v1.8（2026-09-10 全量对齐重做）**｜历史版本 v1.0 / v1.1 留档不删",
   "> **版本**：**v1.9（2026-09-11 全维度审查批次 1）**｜历史版本 v1.0 / v1.1 / v1.8 留档不删"),
  ("> **引用状态**：本章零新增编号（[27] Kim 于 Discussion 引入，为全稿末位编号）；见 §3",
   "> **引用状态**：本章引入 [29] Mathur（sibeprenlimab Ⅱ期临床先例，为全稿末位编号）；[28] Kim 首现于 Results ⑧；见 §3"),
 ],
 """## v1.9 同步记录（2026-09-11 · 全维度审查批次 1）

- 包内 ch5 升 **v1.9**：版本注记定格（Discussion **1,618/1,616**、Limitations **605/603**）；干预窗临床先例句引文为 [29]（Mathur NEJM 2024）。
- §1 英文定稿随包内 `05_Discussion_Limitations.md`（v1.9）逐字节注入。
""",
)

# ---------- ch6 v1.7 ----------
build(
 "章节6_Conclusions_缩略语_v1.6_中英对照注释版.md",
 "章节6_Conclusions_缩略语_v1.7_中英对照注释版.md",
 [
  ("# 章节 6（v1.6）Conclusions + List of abbreviations — 中英对照注释版",
   "# 章节 6（v1.7）Conclusions + List of abbreviations — 中英对照注释版"),
  ("`章节6_Conclusions_缩略语_v1.6_纯英文版.md`",
   "`章节6_Conclusions_缩略语_v1.7_纯英文版.md`"),
  ("> **版本**：**v1.6（2026-09-10 全量对齐重做）**｜历史版本 v1.0 / v1.1 留档不删",
   "> **版本**：**v1.7（2026-09-11 全维度审查批次 1）**｜历史版本 v1.0 / v1.1 / v1.6 留档不删"),
 ],
 """## v1.7 同步记录（2026-09-11 · 全维度审查批次 1）

- 包内 ch6 升 **v1.7**：缩略语表补全 NK 行（\"| NK | natural killer |\"）。
- §1 英文定稿随包内 `06_Conclusions_Abbreviations.md`（v1.7）逐字节注入。
""",
)

print("ALL PATCH-AB BUILDS OK")
