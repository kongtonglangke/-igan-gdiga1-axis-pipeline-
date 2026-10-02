# -*- coding: utf-8 -*-
"""阶段2 主证据3 第③部分：cis-eQTL 顶变异反查 IgAN 信号（稳健性阴性说明）

输入：
  - 00_rawdata/GCST90018866/GCST90018866.h.tsv.gz（Sakaue 2021 IgAN，GRCh38，23.5M 行）
  - M2a 输出的 top cis-eQTL 变异表（GRCh38）
输出：阶段2/主证据3_机制轴不驱动疾病易感/cisEQL_topvar_IgAN_lookup.tsv
逻辑：对 M2a 各显著 top variant（剔除 X 染色体），精确 chr+pos 匹配 IgAN GWAS，
      报告 P/beta/se/EA/OA/EAF → 若全部无信号（P>0.05 或低效应）即支持
      "Gd-IgA1 机制 cis-eQTL 不驱动 IgAN 疾病易感"的稳健性说明。
"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import gzip, os, time

IGAN = _paths.legacy("00_rawdata/GCST90018866/GCST90018866.h.tsv.gz")
OUT_DIR = _paths.at("阶段2/主证据3_机制轴不驱动疾病易感")
os.makedirs(OUT_DIR, exist_ok=True)
OUT = OUT_DIR("cisEQL_topvar_IgAN_lookup.tsv")

# (gene, dataset, chr, pos)  —— 来自 M2a 显著/有趣 top variant（GRCh38）
VARIANTS = [
    ("C1GALT1",  "OneK1K_B_naive",        "7",  7183473),
    ("GALNT2",   "OneK1K_B_naive",        "1",  230279547),
    ("C1GALT1",  "GTEx_v10_blood",        "7",  7223408),
    ("C1GALT1",  "GTEx_v10_kidney_cortex","7",  7182144),
    ("GALNT12",  "GTEx_v10_blood",        "9",  99691091),
    ("GALNT12",  "GTEx_v10_kidney_cortex","9",  99315366),
    ("ST6GALNAC2","GTEx_v10_blood",       "17", 76575845),
    ("ST6GALNAC2","GTEx_v10_kidney_cortex","17",76566488),
    ("GALNT2",   "GTEx_v10_blood",        "1",  230529311),
    ("GALNT2",   "GTEx_v10_kidney_cortex","1",  231049260),
]
key2meta = {(c, p): (g, d) for g, d, c, p in VARIANTS}
WANT = {(c, p) for _, _, c, p in VARIANTS}

t0 = time.time()
hits = {}   # (chr,pos) -> row
with gzip.open(IGAN, "rt") as f:
    header = f.readline().rstrip("\n").split("\t")
    hd = {c: j for j, c in enumerate(header)}
    i_chr, i_pos, i_ea, i_oa = hd["chromosome"], hd["base_pair_location"], hd["effect_allele"], hd["other_allele"]
    i_beta, i_se, i_eaf, i_p = hd["beta"], hd["standard_error"], hd["effect_allele_frequency"], hd["p_value"]
    for n, line in enumerate(f, 1):
        if n % 5_000_000 == 0:
            print(f"  ... {n:,} rows ({time.time()-t0:.0f}s)", flush=True)
        p = line.rstrip("\n").split("\t")
        k = (p[i_chr], p[i_pos])
        if k in WANT and k not in hits:
            hits[k] = (p[i_ea], p[i_oa], p[i_beta], p[i_se], p[i_eaf], p[i_p])

lines = ["gene\tdataset\tchr\tpos\teQTL_beta(来自M2a)\tIgAN_EA\tIgAN_OA\tIgAN_beta\tIgAN_SE\tIgAN_EAF\tIgAN_P"]
# 补 eQTL beta（来自 M2a 表，便于对照方向）
eqtl_beta = {
    ("C1GALT1", "7", 7183473): 0.1811, ("GALNT2", "1", 230279547): 0.2841,
    ("C1GALT1", "7", 7223408): -0.5073, ("C1GALT1", "7", 7182144): 0.9444,
    ("GALNT12", "9", 99691091): 0.5277, ("GALNT12", "9", 99315366): 0.4449,
    ("ST6GALNAC2", "17", 76575845): 0.2111, ("ST6GALNAC2", "17", 76566488): -0.5707,
    ("GALNT2", "1", 230529311): 0.0960, ("GALNT2", "1", 231049260): -0.3746,
}
for gene, ds, ch, pos in VARIANTS:
    eb = eqtl_beta.get((gene, ch, pos), "?")
    if (ch, pos) in hits:
        ea, oa, b, se, eaf, p = hits[(ch, pos)]
        lines.append(f"{gene}\t{ds}\t{ch}\t{pos}\t{eb}\t{ea}\t{oa}\t{b}\t{se}\t{eaf}\t{p}")
    else:
        lines.append(f"{gene}\t{ds}\t{ch}\t{pos}\t{eb}\tNA\tNA\tNA\tNA\tNA\tNOT_FOUND")

with open(OUT, "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")
print(f"[OUT] {OUT}  ({time.time()-t0:.0f}s)")
print("=== 结果 ===")
for l in lines[1:]:
    print("  " + l)
