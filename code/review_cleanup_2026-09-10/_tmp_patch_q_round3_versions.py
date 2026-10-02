# -*- coding: utf-8 -*-
"""
_tmp_patch_q_round3_versions.py — 第三轮整改的版本注落位

在 ch1/ch3/ch4/ch5/ch6 头部修订注块末尾各追加一条 v-next 记录（与既有风格一致），
并把母本装配说明升到 v1.9。不做任何正文改动。
"""
import os
import sys

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
M = os.path.join(BASE, "submission", "manuscript")

NOTES = {
    "01_Title_Abstract_Keywords.md":
        "> Revision in v1.9 (2026-09-10, third reviewer-perspective pass): the opening Abstract sentence was "
        "restructured to remove a dangling participle (\"Produced by …, whether genetic variation … remains "
        "unresolved\" → \"Gd-IgA1, the central pathogenic molecule …, is produced by …; whether … remains "
        "unresolved\"); the word count is unchanged at 348 (limit 350). No data, values or conclusions changed.",
    "03_Methods.md":
        "> v1.10 (2026-09-10, third reviewer-perspective pass): the direct-look-up distance wording now reads "
        "\"approximately two orders of magnitude above the genome-wide significance threshold (P = 9.4 × 10⁻⁶ "
        "versus 5 × 10⁻⁸)\" (was \"three orders of magnitude from genome-wide significance\"), matching the "
        "Results wording for the same value; IDG, the Pharos target-development levels and NBDC are now expanded "
        "at first use. No data, values or conclusions changed.",
    "04_Results.md":
        "> v1.9 (2026-09-10, third reviewer-perspective pass): the methylation section now cites Table 5 in the "
        "running text (journal style requires every main-text table to be cited in the text, and Table 5 had been "
        "uncited); the B-cell-state heterogeneity test is now attributed explicitly to rs13226913; Sp/KLF is "
        "expanded at first use; the direct-look-up distance wording is aligned with Methods. No data, values or "
        "conclusions changed.",
    "05_Discussion_Limitations.md":
        "> v1.7 (2026-09-10, third reviewer-perspective pass): the cell-type paragraph no longer describes the "
        "C1GALT1 cis-eQTL effect as attenuating \"along the naive-to-memory transition\" — wording that implied the "
        "monotonic gradient the Results and the Fig. 2c legend explicitly exclude — and now reads \"attenuated in "
        "the antigen-experienced states, with partial recovery in memory\"; the single-cell per-cell expression "
        "read-out is now reported with its measured magnitudes instead of \"essentially flat\" placed beside a "
        "highly significant P value. No data, values or conclusions changed.",
    "06_Conclusions_Abbreviations.md":
        "> v1.6 (2026-09-10, third reviewer-perspective pass): the abbreviation entry \"PPH3/PPH4\" now reads "
        "\"PPH0–PPH4\" (posterior probability of the five colocalisation hypotheses H0–H4), matching the notation "
        "defined in Results; the closing note no longer claims that all gene symbols are expanded at first "
        "appearance. No other changes.",
}

OLD_ASSEMBLY = "Assembled file v1.8 (2026-09-10), machine-assembled from the pure-English chapter files in this folder; no text added or altered by assembly."
NEW_ASSEMBLY = ("Assembled file v1.9 (2026-09-10), machine-assembled from the pure-English chapter files in this folder; "
                "no text added or altered by assembly.\n"
                "Third reviewer-perspective pass 2026-09-10 (v1.9): ten text-integrity fixes - Table 5 now cited in the "
                "running text; \"three orders of magnitude\" aligned to \"two orders\" for P = 9.4e-6 (ch3 vs ch4); the "
                "Discussion no longer implies a monotonic naive-to-memory gradient (aligned with Fig. 2c/Results); the "
                "Abstract dangling participle removed (word count unchanged at 348/350); the single-cell expression "
                "read-out reported with magnitudes; the heterogeneity test attributed to rs13226913; PPH0-PPH4 notation "
                "in the abbreviation list; Sp/KLF, IDG, the Pharos T-levels and NBDC expanded at first use; the "
                "abbreviation-note gene-symbol claim corrected; and the cover letter raised from seven to eight evidence "
                "layers. No data, values, figures or conclusions changed.")


def insert_note(path, note):
    lines = open(path, encoding="utf-8").read().split("\n")
    last = max(i for i, l in enumerate(lines) if l.startswith(">"))
    assert not any(l.strip() == note for l in lines), f"note already present in {path}"
    lines.insert(last + 1, note)
    open(path, "w", encoding="utf-8").write("\n".join(lines))
    print(f"  OK   version note added -> {os.path.basename(path)}")


def main():
    for fn, note in NOTES.items():
        insert_note(os.path.join(M, fn), note)
    p = os.path.join(BASE, "reproducibility", "code", "_build_submission_fulltext.py")
    t = open(p, encoding="utf-8").read()
    assert t.count(OLD_ASSEMBLY) == 1, "assembly header line not found exactly once"
    open(p, "w", encoding="utf-8").write(t.replace(OLD_ASSEMBLY, NEW_ASSEMBLY, 1))
    print("  OK   assembler header -> v1.9")
    print("done.")


if __name__ == "__main__":
    main()
