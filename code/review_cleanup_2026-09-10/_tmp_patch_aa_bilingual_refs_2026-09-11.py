# -*- coding: utf-8 -*-
# Patch AA (2026-09-11): 生成 References_v2.0_中英对照注释版.md（基于 v1.6 留档件）。
# 内容：27→29 条重编号（[24] decoupleR 插入、JASPAR/Liu/Fung/Kim 顺移、[29] Mathur 新增）、
# §2 表修正（[20] Sun 2015 药物改 5-氮杂-2′-脱氧胞苷（地西他滨）、删错误的「无 DOI」、行 83 注 [20]→[21]）、
# §3 映射重排、§4 核查数更新、§5 审定记录追加。§1 英文定稿随后由 sync_初稿双语版.py 注入。
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
from pathlib import Path

SRC = Path(_paths.legacy("阶段4/初稿/References_v1.6_中英对照注释版.md"))
DST = Path(_paths.legacy("阶段4/初稿/References_v2.0_中英对照注释版.md"))

t = SRC.read_text(encoding="utf-8")
n0 = len(t)

def rep(old, new):
    global t
    assert t.count(old) == 1, f"count!=1: {old[:60]!r} -> {t.count(old)}"
    t = t.replace(old, new)

# ---- 头部 ----
rep("# References [1]–[27]（v1.6）— 中英对照注释版", "# References [1]–[29]（v2.0）— 中英对照注释版")
rep("`References_v1.6_纯英文版.md`", "`References_v2.0_纯英文版.md`")
rep("> **版本**：**v1.6（2026-09-11 初稿→包复核）**｜历史版本 v1.0 / v1.1 / v1.2 / v1.5 留档不删",
    "> **版本**：**v2.0（2026-09-11 全维度审查批次 1）**｜历史版本 v1.0 / v1.1 / v1.2 / v1.5 / v1.6 留档不删")
rep("> **编号规则**：**严格按正文首现序 [1] → [27] 递增**（脚本机检通过）",
    "> **编号规则**：**严格按正文首现序 [1] → [29] 递增**（脚本机检通过）")

# ---- v2.0 同步记录（插在 v1.6 记录之后、§1 之前）----
rep("## 1. References 英文定稿（投稿文本）",
"""## v2.0 同步记录（2026-09-11 · 全维度审查批次 1）

| 变更 | v1.6（旧，27 条） | v2.0（本版，29 条） |
|---|---|---|
| **条目数** | 27 | **29**（新增 [24] decoupleR（Badia-i-Mompel 2022, Bioinform Adv 2(1):vbac016）与 [29] Mathur 2024（NEJM 390(1):20–31，sibeprenlimab Ⅱ期）；原 [24]–[27] 顺移为 [25]–[28]） |
| **新增条目来源** | — | L7 调控子活性方法（decoupleR ULM，原稿已用但未引）+ Discussion 干预窗临床先例（抗 APRIL Ⅱ期） |
| **首现序** | [1..27] 严格递增 | **[1..29] 机检非递减通过**（同括号组并列首现合法） |
| **§2 表修正** | [20] Sun 2015 注记「5-氮杂胞苷逆转；无 DOI」 | 更正为「5-氮杂-2′-脱氧胞苷（地西他滨）逆转」（PubMed 核：Sun 2015 用 decitabine；IL-17/5-azacytidine 属 [19] Lin 2018）；删「无 DOI」误注（[20] 有 DOI 10.1371/journal.pone.0112305；无 DOI 者为 [21] Fang 2003） |

- §1 英文定稿随包内 `09_References.md`（v2.0）逐字节注入（`sync_初稿双语版.py`）。

---

## 1. References 英文定稿（投稿文本）""")

# ---- §2 表 ----
rep("> 下表列出每条文献的**标识、类型与在稿中的作用**；DOI 一栏为投稿版实际所附链接（[20] 无 DOI）。",
    "> 下表列出每条文献的**标识、类型与在稿中的作用**；DOI 一栏为投稿版实际所附链接（[21] 无 DOI）。")
rep("| 20 | Sun Q, et al. DNA methylation in Cosmc promoter region and aberrantly glycosylated IgA1 associated with pediatric IgAN. | 2015, PLoS One | 原始研究 | **C1GALT1C1(Cosmc) 启动子甲基化** + 5-氮杂胞苷逆转（C1GALT1C1 得分 6.0 的核心锚点）；**无 DOI** |",
    "| 20 | Sun Q, et al. DNA methylation in Cosmc promoter region and aberrantly glycosylated IgA1 associated with pediatric IgAN. | 2015, PLoS One | 原始研究 | **C1GALT1C1(Cosmc) 启动子甲基化**（IL-4 驱动）+ **5-氮杂-2′-脱氧胞苷（地西他滨）逆转**（C1GALT1C1 得分 6.0 的核心锚点） |")
rep("| 24 | Rauluseviciute I, et al. JASPAR 2024: 20th anniversary of the open-access database of transcription factor binding profiles. | 2024, Nucleic Acids Res | 资源论文 | **JASPAR 2024 CORE** 基序矩阵（甲基化敏感基序分析） |",
    "| 24 | Badia-i-Mompel P, et al. decoupleR: ensemble of computational methods to infer biological activities from omics data. | 2022, Bioinform Adv | 方法论文 | **decoupleR ULM**（L7 调控子活性的 per-cell 估计方法） |\n| 25 | Rauluseviciute I, et al. JASPAR 2024: 20th anniversary of the open-access database of transcription factor binding profiles. | 2024, Nucleic Acids Res | 资源论文 | **JASPAR 2024 CORE** 基序矩阵（甲基化敏感基序分析） |")
rep("| 25 | Liu L, et al.", "| 26 | Liu L, et al.")
rep("| 26 | Fung WWS, et al.", "| 27 | Fung WWS, et al.")
rep("| 27 | Kim G, et al. Chromatin accessibility of circulating CD8⁺ T cells differentiates disease severity in IgAN. | 2026, Sci Rep | 原始研究 | **GSE285335 单细胞数据集**的正式来源（原会议摘要引用已更正） |",
    "| 28 | Kim G, et al. Chromatin accessibility of circulating CD8⁺ T cells differentiates disease severity in IgAN. | 2026, Sci Rep | 原始研究 | **GSE285335 单细胞数据集**的正式来源（原会议摘要引用已更正） |\n| 29 | Mathur M, et al. A Phase 2 Trial of Sibeprenlimab in Patients with IgA Nephropathy. | 2024, N Engl J Med | 原始研究（Ⅱ期临床） | **抗 APRIL（sibeprenlimab）降低 Gd-IgA1 的临床先例**（Discussion 干预窗论证） |")

# ---- §3 映射 ----
rep("""| [23] | Methods ⑬ | 单细胞图谱小节 | Results ⑧ |
| [24] | Methods ⑭ | 基序分析小节 | Results ⑧ |
| [25] | Results ③ | 糖基特异性句 | — |
| [26] | Results ⑤ | 一致性句 | Discussion |
| [27] | Table 1 / Discussion | GSE285335 归属 | Methods ⑫（以数据集名/访问号引用正文） |""",
"""| [23] | Methods ⑬ | 单细胞图谱小节 | Results ⑧ |
| [24] | Methods ⑬ | 调控子活性估计（decoupleR ULM） | Results ⑧ |
| [25] | Methods ⑭ | 基序分析小节 | Results ⑧ |
| [26] | Results ③ | 糖基特异性句 | — |
| [27] | Results ⑤ | 一致性句 | Discussion |
| [28] | Results ⑧ | GSE285335 归属 | Table 1 / Methods ⑫ / Discussion（以数据集名/访问号引用） |
| [29] | Discussion | 临床先例句（sibeprenlimab Ⅱ期） | — |""")
rep("> **首现序自洽**：[12]–[24] 严格按全稿首现序重排（2026-09-10）；[25]–[27] 延伸至 Results/Discussion。**GEO 系列以访问号引用、不挂编号**，故 Table 1 的 GEO 行标注为 \"—*\"。",
    "> **首现序自洽**：[12]–[25] 严格按全稿首现序（2026-09-11 终排）；[26]–[29] 延伸至 Results/Discussion。**GEO 系列以访问号引用、不挂编号**，故 Table 1 的 GEO 行标注为 \"—*\"。")

# ---- §4 核查 ----
rep("- **条目数**：27（[1]–[27]，编号连续无缺号）。", "- **条目数**：29（[1]–[29]，编号连续无缺号）。")
rep("- **首现序**：正文提取的引用序列 **[1..27] 严格递增**，无回跳（脚本 `_ref_audit_v12.py` 升级版核验通过）。",
    "- **首现序**：正文提取的引用序列 **[1..29] 非递减通过**（同括号组并列首现合法；母本机检 `_tmp_machine_check_2026-09-11.py`，2026-09-11）。")
rep("- **DOI 完整性**：26/27 条附 DOI 且经 Crossref 解析成功；[21] Fang 2003 **按事实留空**（无注册 DOI）。",
    "- **DOI 完整性**：28/29 条附 DOI 且经 Crossref 解析成功；[21] Fang 2003 **按事实留空**（无注册 DOI）。")
rep("- **作者串体例**：≤ 6 位全列；≥ 7 位列前 6 位 + \"et al.\"，全 27 条一致。",
    "- **作者串体例**：≤ 6 位全列；≥ 7 位列前 6 位 + \"et al.\"，全 29 条一致。")
rep("- **未追随项**：**[26] 为预印本**——不附卷期页，DOI 为 medRxiv 版本；**投稿当日须复核是否已正式发表**。",
    "- **未追随项**：**[27] 为预印本**——不附卷期页，DOI 为 medRxiv 版本；**投稿当日须复核是否已正式发表**。")

# ---- §5 追加 ----
rep("**引用链演进（2026-09-10）**：",
    "**v2.0（2026-09-11，全维度审查批次 1）**：27 → **29 条**——插入 [24] decoupleR（L7 调控子活性方法，原稿已用未引）、追加 [29] Mathur NEJM 2024（sibeprenlimab Ⅱ期临床先例）；原 [24]–[27] 顺移 [25]–[28]；§2 [20] Sun 2015 注记更正（地西他滨、删「无 DOI」误注）；§3 映射重排（[28] Kim 首现订正为 Results ⑧）。\n\n**引用链演进（2026-09-10）：")

DST.write_text(t, encoding="utf-8", newline="\n")
print(f"bilingual References v2.0 written: {len(t)} chars (src {n0})")
print("ALL PATCH-AA EDITS OK")
