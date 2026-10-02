# -*- coding: utf-8 -*-
"""_tmp_report_v13.py — 整改记录升 v1.3：作者定夺项办结 + 新发现入档"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import io, os

P = _paths.legacy("阶段4.5_复现包/审稿视角全文检查报告_v1.1_整改记录_2026-09-10.md")
t = io.open(P, encoding="utf-8", newline="").read()

OLD_SEC6 = """## 六、仍需作者/用户定夺（本轮未擅自处理）

1. **署名次序与 CRediT 的最终形态**（P1-7）——本轮的改写只在"不虚构事实"的前提下使贡献描述与并列一作自洽；**若作者希望第一署名者的贡献更实，需作者提供真实分工**。
2. **通讯作者 Yan Wang 与参考文献 [8] Wang YN 是否同一人**（P1-7）——若为同一人须在 Competing interests 或 Acknowledgements 明示；文件内此问仍未决。
3. **GA 图内的 "upstream circuit (L7)" 字样**是否同去编号——正文已无编号，Fig. 1c 保留 L0–L7 行索引属"面板自印图层"，GA 的单个 "(L7)" 与之呼应但属图件本体（本轮按"图片不在范围"未动）；如要统一，改 `ga_render.py` 一行字符串 + 重渲染即可。
4. **Additional file 1 与各章文件头部的 `>` 编辑版本注**：母本已在拼装时剔除，但 `Additional_File_1_Figure_Legends.md` 是直接交付件，其 `>` 头块投稿前需删（与母本 `<!-- ASSEMBLY NOTE -->` 同理）。
5. **[26] Fung 仍是 medRxiv 预印本**：投稿当日必须再查一次是否已正式发表（原报告 P0 延伸项，未变）。
6. **Funding 表述**：若 2027 省基金在投稿前立项，须改为基金资助字样。"""

NEW_SEC6 = """## 六、作者/用户定夺项（09-10 傍晚作者已答复，逐条办结见第八节）

| # | 事项 | 作者答复 | 状态 |
|---|---|---|---|
| ① | 署名次序与 CRediT 的最终形态（P1-7） | **未答复** | ⏳ 待全体作者确认；本轮的改写只在"不虚构事实"前提下使贡献描述与并列一作自洽，若希望第一署名者贡献更实，需作者提供真实分工 |
| ② | 通讯作者 Yan Wang 与参考文献 [8] Wang YN 是否同一人 | **否** | ✅ 办结（非同一人 → 无需披露） |
| ③ | GA 图内 "upstream circuit (L7)" 是否同去编号 | **去** | ✅ 办结（GA v3.2→**v3.3**，两处脚本副本同步 + 重渲染） |
| ④ | Additional file 1 直接交付件的 `>` 编辑版本注头块 | **删** | ✅ 办结（入档后从交付件删除） |
| ⑤ | [26] Fung medRxiv → 正式版复核 | 当日复核 | ✅ 办结（09-10 复核：仍为 medRxiv 预印本 v1） |
| ⑥ | Funding 表述（2027 省基金立项后） | — | ⏳ 立项后改 |"""

assert OLD_SEC6 in t, "section 六 未命中"
t = t.replace(OLD_SEC6, NEW_SEC6, 1)
print("OK 六、替换")

OLD_FOOT = "*v1.2 · 2026-09-10 · 在 v1.1 基础上补漏 8 处裸 `Rn` 残留并修正机检正则口径；图片本体按作者要求未动。*"
NEW_FOOT = OLD_FOOT + """

---

## 八、作者定夺项办结 + 新发现（v1.3，2026-09-10 傍晚）

作者答复：**②否、③去、④删、⑤当日复核**；①未答复（仍待）。

| # | 作者答复 | 本轮处置 | 涉及文件 / 版本 |
|---|---|---|---|
| ② | **否**（Yan Wang ≠ [8] Wang YN） | 非同一人 → 无共同利益需披露；ch7 `Competing interests` 维持 "no competing interests"，**未改文本**，仅本条入档备查 | ch7（未改） |
| ③ | **去** | `ga_render.py` 渲染字符串 `"upstream circuit (L7)"` → `"upstream circuit"`（`阶段4/主图/scripts/` 与 `reproducibility/code/fig_scripts/` **两处副本同步**）；默认 `OUT_STEM` → `GA_v3.3`；重渲染 + `_ga_check.py` → **13 texts / 4 solid boxes / leaks = 0**；ch1 的 GA 版本引用 v3.2→**v3.3**（ch1 → **v1.7**）；GA 图注草案 v02→**v03**；投稿件 `submission/figures/GA.{png,pdf,tiff}` 已换为 v3.3（v3.0/v3.1/v3.2 留档不删） | `ga_render.py`（×2）、`GA_v3.3.{png,pdf,tiff}`、`ga_caption_v03.md`、ch1 v1.7 |
| ④ | **删** | 头块 5 行**原样入档** `reproducibility/code/review_cleanup_2026-09-10/Additional_File_1_Figure_Legends_header_notes_archived.md`，随后从交付件删除；交付件现首行=标题、次行=`---`，`>` 头注 **0** | `Additional_File_1_Figure_Legends.md`、归档件 |
| ⑤ | 当日复核 | 2026-09-10 复核（medRxiv / Scilit / Altmetric）：[26] Fung 仍为 **medRxiv 预印本 v1（posted 2025-10-25；DOI 10.1101/2025.10.25.25338806）**，未见正式期刊版 → 引用维持 `medRxiv [preprint]. 2025.`，**未改文本** | ch9 [26]（未改） |

**决策④执行中机检新命中并修复**：`Additional_File_1_Figure_Legends.md` 的**补充材料表索引**里仍有 **6 处美式 `colocalization`**——v1.5 的 -isation 统一只覆盖了 ch1–ch9 正文，**漏了这个直接交付件**；已全数改为 `colocalisation` / `Colocalisation`（残留 0）。该件现已无版本注头，版本信息改记于本报告与包 README。

**余留未办结**：① 署名次序与 CRediT 最终形态（待全体作者确认）；⑥ 2027 省基金若立项后改 Funding 表述。

**终检（09-10 傍晚，母本重建 110,371 chars）**：全部交付件（剔 `>` 头注 / HTML 块 / ch4 集成附录）——`evidence layer L*` / `(Rn)` / 裸 `Rn` / 裸 `Mn` / `(Lx)` / 陈旧串（causal framework、colocalization、naïve、six sequential、TFAP2A）**全 0**；`(L0–L7)` 保留 1（Fig. 1c 面板自印）；引用首现序 **[1..27] 严格递增**。

*v1.3 · 2026-09-10（傍晚）· 在 v1.2 基础上办结作者决策 ②③④⑤ 并修复补充材料索引 6 处美式拼写。*"""

assert OLD_FOOT in t, "footer 未命中"
t = t.replace(OLD_FOOT, NEW_FOOT, 1)
print("OK 八、追加")

io.open(P, "w", encoding="utf-8", newline="").write(t)
print("saved, len:", len(t))
