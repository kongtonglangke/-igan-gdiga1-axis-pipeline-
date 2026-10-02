# -*- coding: utf-8 -*-
# Patch Z (2026-09-11): 包 README v2.6 + reproducibility README v2.1（含 F4 sc_regulatory 补录、
# 组件表版本/词数定格、母本 v1.12/v1.13 建造日志、词数口径声明）。
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
from pathlib import Path

PKG = Path(_paths.LEGACY, "README.md")
REP = Path(_paths.repo("README.md"))

def replace_line(lines, startswith, new_line, expect=1):
    hits = [i for i, ln in enumerate(lines) if ln.startswith(startswith)]
    assert len(hits) == expect, (startswith[:40], hits)
    lines[hits[0]] = new_line
    return hits[0]

# ================= 包 README =================
t = PKG.read_text(encoding="utf-8")
lines = t.split("\n")

# 1) 标题版本
assert lines[0].startswith("# 阶段 4.5 复现包 / 投稿材料包（v2.5"), lines[0][:60]
lines[0] = lines[0].replace("（v2.5 · 2026-09-11）", "（v2.6 · 2026-09-11）")

# 2) 组件表
replace_line(lines, "| Title/Abstract/Keywords |",
 "| Title/Abstract/Keywords | ch1 **v2.0**（全维度审查批次 1：版本注记定格；Abstract 348 不变。此前：第三轮悬垂分词修正；审稿清洗去编号 + 东亚句点名 rs10238682/rs7856182；GA 版本引用 v3.3） | Abstract **348**（限 350 ✅） | 全稿头部 |")
replace_line(lines, "| Methods（含 Table 1） |",
 "| Methods（含 Table 1） | ch3 **v1.11**（全维度审查批次 1：Table 1 GoDMC 行更正 24,988–28,181 per CpG；文献锚句恢复首现序 [19]→[20]（IL-17/5-azacytidine 在先）；GSE285335 化学版本表述去版本号；引文编号同步至 29 条新序） | **3,637 含 / 3,532 不含**（严格全 Table 1 剔除口径） | 正文 |")
replace_line(lines, "| Results R1–R8 |",
 "| Results R1–R8 | ch4 **v2.1**（批次 1：PPH1 句重构为实测两组、药物归属按 [19]/[20] 拆分、Fung 段补 Wang 2021 候选层张力句、编号同步；批次 2：R8 双分区检出率 4.2%→**4.1%**（对齐 Table S7 加权 (737+135+859)/42,259=4.097%）） | **3,292 含 / 3,175 不含**（8 小节） | 正文 |")
replace_line(lines, "| Discussion + Limitations |",
 "| Discussion + Limitations | ch5 **v1.9**（全维度审查批次 1 版本注记定格；此前：删「在研正交实验」前瞻段 + 末段「边界与方向」结语；去「单调梯度」表述 + 单细胞表达量口径） | Discussion **1,618/1,616** / Limitations **605/603** | 正文 |")
replace_line(lines, "| Conclusions + 缩略语表 |",
 "| Conclusions + 缩略语表 | ch6 **v1.7**（全维度审查批次 1：缩略语表 NK 行补全；此前八层口径同步 + PPH0–PPH4 记法） | Conclusions **148 含 / 146 不含** | 正文 |")
replace_line(lines, "| Figure legends + Tables 2–7 |",
 "| Figure legends + Tables 2–7 | ch8 **v2.2**（批次 1：Fig. 1a/3b/3c/6b/6d 五处图注修正 + 配套重渲 Fig. 1 v3.8/Fig. 6 v3.5/GA v3.4；批次 2：**Table 5 mQTL SE 更正为随机效应值**（0.066/0.014/0.019/0.057/0.011，与 Table S5b se_mqtl 一致），脚注 † 明记 random-effects 估计） | — | 正文尾部（图注/表） |")
replace_line(lines, "| References |",
 "| References | **v2.0**（**27→29 条**：decoupleR 插入为 [24]（Bioinform Adv 2022;2(1):vbac016），Mathur NEJM 2024 Sibeprenlimab 追加为 [29]；[21] Fang 2003 题名更正 \"cancer cell lines\"；首现序 [1..29] 机检非递减通过） | **29 条** | 正文尾部 |")
replace_line(lines, "| **正文合计（章节2–6，含小节标题）** |",
 "| **正文合计（章节2–6，含小节标题）** | — | **10,062 含 / 9,832 不含**（= BG 762/760 + Methods 3,637/3,532 + Results 3,292/3,175 + Disc 1,618/1,616 + Lim 605/603 + Concl 148/146；2026-09-11 全维度审查后按统一口径脚本实测） | 脚本实测定格 |")

# 3) 词数口径备注行（第 42 行备注后补一条口径说明）
hits = [i for i, ln in enumerate(lines) if ln.startswith("- 备注：**Discussion/Limitations 长度已定案")]
assert len(hits) == 1
lines.insert(hits[0] + 1,
 "- 词数统计口径（2026-09-11 起统一声明）：原始空白分词；排除各章头部 `>` 注记、`---` 分隔符、表格行与 `*Word count*` 行；**Methods 另严格全 Table 1 剔除**（表题+表注+表行）；「含/不含」指含/不含小节标题；正文合计**含 Limitations**。脚本：`reproducibility/code/review_cleanup_2026-09-10/_tmp_measure_wordcount_2026-09-11.py`。")

# 4) 母本建造日志：v1.11 后插 v1.12/v1.13
hits = [i for i, ln in enumerate(lines) if ln.startswith("v1.11（2026-09-11）：**初稿→包复核后重建**")]
assert len(hits) == 1
lines.insert(hits[0] + 1,
 "v1.12（2026-09-11）：**全维度审查批次 1 后重建**——吸收参考文献 27→29 条重编号（ch9 v2.0）、ch3 文献锚句首现序恢复（[19]→[20]）、Table 1 GoDMC 行更正、药物归属按 [19]/[20] 拆分（ch3/ch4/ch8/Table S6b + 图件 Fig. 1 v3.8/Fig. 6 v3.5/GA v3.4）、Table S5a 重排可读列；母本 **116,021** chars（装配说明版本签一度因同文件并行编辑丢失，v1.13 补回）。机检全过。")
lines.insert(hits[0] + 2,
 "v1.13（2026-09-11）：**审查批次 2 后重建**——ch8 **v2.2**（Table 5 mQTL SE 改随机效应值 0.066/0.014/0.019/0.057/0.011 + 脚注明记 random-effects，与 Table S5b 一致）、ch4 **v2.1**（R8 检出率 4.2%→4.1%）；装配说明补登 v1.12/v1.13 两条注记；母本 **117,066** chars。机检全过（Abstract 348、[1..29] 非递减、正文中文 0、内部编号 0、陈旧串 0、S10 结构五项）。")

# 5) 追加 v2.6 版本注（最后一个非空行之后）
tail = [i for i, ln in enumerate(lines) if ln.strip()]
assert lines[tail[-1]].startswith("- v2.5（2026-09-11）"), lines[tail[-1]][:40]
insert_at = tail[-1] + 1
lines.insert(insert_at,
 "- v2.6（2026-09-11）：**全维度审查（四人专家组 + 数据复核 r2）修复全部落地**。两轮修复：**批次 1**——参考文献 27→**29** 条（[24] decoupleR、[29] Mathur NEJM）、Methods 文献锚句首现序恢复、Table 1 GoDMC 行 24,988–28,181、药物归属拆分（5-aza-2′-deoxycytidine [20]＝IL-4 甲基化逆转/Sun 2015；5-azacytidine [19]＝IL-17 表达恢复/Lin 2018，PubMed 核）、GoDMC P 列随机效应口径确认、BH q 经**两轮独立复算 + 疑点方正式撤回**（13 检验族 6.1e-5、lead 0.12 均正确，勿动）、Table S5a 重排 11 可读列（原档留档）、图件 Fig. 1 v3.8/Fig. 6 v3.5/GA v3.4、Table S6b 药物行重写；**批次 2**——**Table 5 mQTL SE 更正为随机效应值**（原混入固定效应 ARE 的 SE；β 与 cis-mQTL P 列本来就是 MRE，S5b se_mqtl 可证）、R8 检出率 4.2%→4.1%、**Table S10 补 axis_composite_null_q95/axis_composite_emp_p 两列**（950 组合行，源 regulon_axis_composite.tsv；逐基因行留空，原档留档）、reproducibility/README 数据目录补 sc_regulatory（F4）。组件表与词数按统一口径实测定格（正文合计 **10,062/9,832**；Abstract **348** 不变）。留痕脚本 `_tmp_patch_{w2,w3,w4,x,y,z}_*.py`、`_tmp_patch_figtext_2026-09-11.py`、`_tmp_machine_check_2026-09-11.py`（全 PASS）、`_tmp_measure_wordcount_2026-09-11.py`。**至今无未裁定遗留项**（开放项仍为作者终校/docx·pdf/建仓 DOI/[26] 当日复核/cover letter 定稿）。")

PKG.write_text("\n".join(lines), encoding="utf-8")
print("package README -> v2.6 OK")

# ================= reproducibility README =================
t2 = REP.read_text(encoding="utf-8")
l2 = t2.split("\n")

# F4: data/ 目录清单补 sc_regulatory
hits = [i for i, ln in enumerate(l2) if ln.startswith("M1_geo（表达与差异）")]
assert len(hits) == 1
old = "M6_intervention（评分/药物）。清单与文件来源"
new = "M6_intervention（评分/药物）、sc_regulatory（单细胞 atlas/上游调控子/axis 关联与经验零分布）。清单与文件来源"
assert l2[hits[0]].count(old) == 1
l2[hits[0]] = l2[hits[0]].replace(old, new)

# v2.1 版本注插到 *v2.0 之前（新者在前）
hits = [i for i, ln in enumerate(l2) if ln.startswith("*v2.0 · 2026-09-11")]
assert len(hits) == 1
l2.insert(hits[0],
 "*v2.1 · 2026-09-11 · 阶段4.5 复现包* — v2.0→v2.1：**全维度审查修复的包内同步**。① `data/` 目录清单补登 `sc_regulatory/`（F4 遗漏更正；该目录一直在包内，仅清单漏列）；② 留痕脚本新增批次 1–2：`_tmp_patch_w2_full_audit_resume.py`、`_tmp_patch_w3_r2_findings.py`（含 Table_S5a 重排，原档 `Table_S5a_sigCpG_5CpG_pos.before_reformat_2026-09-11.tsv` 留档）、`_tmp_patch_w4_version_notes.py`、`_tmp_patch_figtext_2026-09-11.py`（Fig. 1 v3.8/Fig. 6 v3.5/GA v3.4 重渲）、`_tmp_patch_x_s10_nullq95_2026-09-11.py`（Table S10 补 axis_composite_null_q95/emp_p 两列，原档 `.before_nullq95_2026-09-11.tsv` 留档）、`_tmp_patch_y_version_notes_2026-09-11b.py`、`_tmp_patch_z_readmes_2026-09-11.py`；③ `_build_submission_fulltext.py` 装配说明升 **v1.13**（补登 v1.12 注记），母本 **117,066** chars；④ `_tmp_machine_check_2026-09-11.py` 扩充 S10 结构五项检查，全 PASS；⑤ 词数口径统一声明见包 README（Methods 严格全 Table 1 剔除；正文合计 10,062/9,832）。")

REP.write_text("\n".join(l2), encoding="utf-8")
print("reproducibility README -> v2.1 OK (F4 sc_regulatory listed)")
print("ALL PATCH-Z EDITS OK")
