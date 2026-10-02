# -*- coding: utf-8 -*-
# _tmp_patch_w3_r2_findings.py — 2026-09-11 全维度审查修复批次（r2 数据核查 FAIL 项）
# F2: 03 Table 1 GoDMC 行样本量范围 25,095–28,181 → 24,988–28,181（Table_S5a 实测：
#     25560/25095/24988/28181/24988，最小值 24,988）
# F3: 05 Discussion "detected in 4.2% of all B cells" → 4.1%
#     （Table_S7 加权 (737+135+859)/42,259 = 4.097% ≈ 4.1%；A1 S4 图注 "≈ 4%" 不变）
# 另: Table_S5a info 列 Python dict dump → 可读列（原档留痕至 review_cleanup_2026-09-10），
#     A1 索引 S5a 条目同步新口径。
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os, shutil, csv

BASE = _paths.at("阶段4.5_复现包/submission")
M = BASE("manuscript")
A1D = BASE("additional_file_1")
A1 = os.path.join(A1D, "Additional_File_1_Figure_Legends.md")
S5A = os.path.join(A1D, "Table_S5a_sigCpG_5CpG_pos.tsv")
ARCH = _paths.legacy("阶段4.5_复现包/reproducibility/code/review_cleanup_2026-09-10/Table_S5a_sigCpG_5CpG_pos.before_reformat_2026-09-11.tsv")

EDITS = {
os.path.join(M, "03_Methods.md"): [
 ("| GoDMC | Whole-blood methylation QTLs | 25,095–28,181 per CpG | GRCh37 | godmc.org.uk (API v0.1) | [13] |",
  "| GoDMC | Whole-blood methylation QTLs | 24,988–28,181 per CpG | GRCh37 | godmc.org.uk (API v0.1) | [13] |"),
],
os.path.join(M, "05_Discussion_Limitations.md"): [
 ("(C1GALT1C1 is detected in 4.2% of all B cells)",
  "(C1GALT1C1 is detected in 4.1% of all B cells)"),
],
A1: [
 ("| Table S5a | The five significant cis-mQTL CpGs of *C1GALT1*: CpG identifier, chromosome, hg19 position, lead-variant hits and annotation. Corresponds to Fig. 4 and Table 5. |",
  "| Table S5a | The five significant cis-mQTL CpGs of *C1GALT1*: CpG identifier, chromosome, hg19 position, lead-variant hits, probe type, weighted mean and SD of the methylation β value, mQTL sample size, QC flags (Zhou; TwinsUK) and association class. Corresponds to Fig. 4 and Table 5. |"),
],
}

def patch_file(path, pairs):
    with open(path, "r", encoding="utf-8", newline="") as f:
        text = f.read()
    for i, (old, new) in enumerate(pairs):
        n = text.count(old)
        assert n == 1, f"{'MISS' if n==0 else 'DUP'} {os.path.basename(path)} #{i}: {old[:70]!r} (count={n})"
        text = text.replace(old, new)
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write(text)
    print(f"OK  {os.path.basename(path)}  ({len(pairs)} edits)")

def reformat_s5a():
    # 留痕原档
    shutil.copy2(S5A, ARCH)
    with open(S5A, "r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert len(rows) == 5 and set(rows[0].keys()) == {"cpg", "chr", "pos37", "lead_hits", "info"}
    # 从 dict 字面量提取字段（纯手工文件，无生成脚本；字段名固定）
    import ast
    out_rows = []
    for r in rows:
        d = ast.literal_eval(r["info"])
        assert d["name"] == r["cpg"] and str(d["pos"]) == r["pos37"] and str(d["chr"]) == r["chr"]
        out_rows.append({
            "cpg": r["cpg"], "chr": r["chr"], "pos37": r["pos37"], "lead_hits": r["lead_hits"],
            "probe_type": d["probetype"],
            "weighted_mean_methylation_beta": repr(d["weighted_mean"]),
            "weighted_sd": repr(d["weighted_sd"]),
            "mQTL_sample_size": str(d["samplesize"]),
            "qc_zhou": str(d["qc_zhou"]), "qc_twinsuk": str(d["qc_twinsuk"]),
            "assoc_class": d["assoc_class"],
        })
    ns = sorted(int(r["mQTL_sample_size"]) for r in out_rows)
    assert ns[0] == 24988 and ns[-1] == 28181, f"sample-size range changed: {ns}"
    cols = ["cpg", "chr", "pos37", "lead_hits", "probe_type",
            "weighted_mean_methylation_beta", "weighted_sd", "mQTL_sample_size",
            "qc_zhou", "qc_twinsuk", "assoc_class"]
    with open(S5A, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out_rows)
    print(f"OK  Table_S5a reformatted (5 rows, 11 cols; original archived to {os.path.basename(ARCH)})")

if __name__ == "__main__":
    for path, pairs in EDITS.items():
        patch_file(path, pairs)
    reformat_s5a()
    print("ALL R2-FINDING FIXES APPLIED")
