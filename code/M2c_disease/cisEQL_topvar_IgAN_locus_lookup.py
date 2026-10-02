# -*- coding: utf-8 -*-
"""阶段2 主证据3 第③部分 v2：cis-eQTL 顶变异 → IgAN 区域信号检验（稳健性阴性说明）

v2 修正（2026-09-02）：
  - v1 用精确坐标匹配 → 全 NOT_FOUND。根因：Sakaue（BBJ 日本人群参考面板）与
    OneK1K/GTEx（TopMed/1KG）存在 SNP 级 imputation 差异，精确命中本不应预期。
  - 正确的检验单元 = cis 区域（±500kb，与 eQTL cis 定义一致）：
    对每个 top cis-eQTL variant，报告其区域在 IgAN GWAS 中最显著 SNP 的 P/beta/se。
    区域最小 P 全 >5e-8（或接近基因组显著但需 LD 判断）→ 支持
    "Gd-IgA1 机制 cis-eQTL 区域在 IgAN 疾病易感中无信号"阴性结论。
  - 另附精确坐标命中与否（供审阅诚实记录）。

输入：00_rawdata/GCST90018866/GCST90018866.h.tsv.gz（Sakaue 2021 IgAN, GRCh38, 23.5M 行）
输出：cisEQL_topvar_IgAN_locus_lookup.tsv
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
OUT = OUT_DIR("cisEQL_topvar_IgAN_locus_lookup.tsv")
FLANK = 500_000  # cis 区域 ±500kb

# (gene, dataset, chr, pos)  —— M2a 各语境 top cis-eQTL variant（GRCh38）
VARIANTS = [
    ("C1GALT1",  "OneK1K_B_naive",         "7",  7183473),
    ("GALNT2",   "OneK1K_B_naive",         "1",  230279547),
    ("C1GALT1",  "GTEx_v10_blood",         "7",  7223408),
    ("C1GALT1",  "GTEx_v10_kidney_cortex", "7",  7182144),
    ("GALNT12",  "GTEx_v10_blood",         "9",  99691091),
    ("GALNT12",  "GTEx_v10_kidney_cortex", "9",  99315366),
    ("ST6GALNAC2","GTEx_v10_blood",        "17", 76575845),
    ("ST6GALNAC2","GTEx_v10_kidney_cortex","17", 76566488),
    ("GALNT2",   "GTEx_v10_blood",         "1",  230529311),
    ("GALNT2",   "GTEx_v10_kidney_cortex", "1",  231049260),
]
# 合并重叠窗口：同 chr 且区间重叠的归一组
import bisect
intervals = []  # (lo, hi, meta_list)
for gene, ds, ch, pos in VARIANTS:
    lo, hi = pos - FLANK, pos + FLANK
    merged = False
    for iv in intervals:
        if iv[0] == ch and not (hi < iv[1][0] or lo > iv[1][1]):
            iv[1][0] = min(iv[1][0], lo); iv[1][1] = max(iv[1][1], hi)
            iv[2].append((gene, ds, pos))
            merged = True
            break
    if not merged:
        intervals.append([ch, [lo, hi], [(gene, ds, pos)]])

print(f"合并后窗口数: {len(intervals)}")
best = {id(k): None for k in intervals}  # 窗口 -> 区域内最显著 SNP
exact = {}  # (gene,ds,chr,pos) -> 精确命中行

t0 = time.time()
with gzip.open(IGAN, "rt") as f:
    header = f.readline().rstrip("\n").split("\t")
    hd = {c: j for j, c in enumerate(header)}
    i_chr, i_pos = hd["chromosome"], hd["base_pair_location"]
    i_ea, i_oa, i_rsid = hd["effect_allele"], hd["other_allele"], hd["rsid"]
    i_beta, i_se, i_eaf, i_p = hd["beta"], hd["standard_error"], hd["effect_allele_frequency"], hd["p_value"]
    for n, line in enumerate(f, 1):
        if n % 5_000_000 == 0:
            print(f"  ... {n:,} rows ({time.time()-t0:.0f}s)", flush=True)
        p = line.rstrip("\n").split("\t")
        ch, pos = p[i_chr], int(p[i_pos])
        for iv in intervals:
            if iv[0] != ch or not (iv[1][0] <= pos <= iv[1][1]):
                continue
            try:
                pv = float(p[i_p])
            except Exception:
                continue
            cur = best[id(iv)]
            if cur is None or pv < cur[0]:
                best[id(iv)] = (pv, ch, pos, p[i_rsid], p[i_ea], p[i_oa], p[i_beta], p[i_se], p[i_eaf])
            # 精确命中记录
            for gene, ds, vpos in iv[2]:
                if pos == vpos and (gene, ds, ch, vpos) not in exact:
                    exact[(gene, ds, ch, vpos)] = (p[i_rsid], p[i_ea], p[i_oa], p[i_beta], p[i_se], p[i_eaf], p[i_p])

# 输出
eqtl_beta = {
    ("C1GALT1", "7", 7183473): 0.1811, ("GALNT2", "1", 230279547): 0.2841,
    ("C1GALT1", "7", 7223408): -0.5073, ("C1GALT1", "7", 7182144): 0.9444,
    ("GALNT12", "9", 99691091): 0.5277, ("GALNT12", "9", 99315366): 0.4449,
    ("ST6GALNAC2", "17", 76575845): 0.2111, ("ST6GALNAC2", "17", 76566488): -0.5707,
    ("GALNT2", "1", 230529311): 0.0960, ("GALNT2", "1", 231049260): -0.3746,
}
lines = ["gene\tdataset\tcis_eQTL_variant(chr:pos)\teQTL_beta\texact_hit\tlocus_minP_SNP\tlocus_minP_P\tlocus_minP_beta\tlocus_minP_SE\tlocus_minP_EAF"]
for gene, ds, ch, pos in VARIANTS:
    key = (gene, ch, pos)
    eb = eqtl_beta[key]
    ex = exact.get(key)
    ex_str = "YES" if ex else "NO(参考面板SNP差异)"
    # 找含该 pos 的窗口
    win = None
    for iv in intervals:
        if iv[0] == ch and iv[1][0] <= pos <= iv[1][1]:
            win = iv; break
    b = best[id(win)]
    if b:
        lines.append(f"{gene}\t{ds}\t{ch}:{pos}\t{eb}\t{ex_str}\t{b[1]}:{b[2]}({b[3]})\t{b[0]:.3g}\t{b[6]}\t{b[7]}\t{b[8]}")
    else:
        lines.append(f"{gene}\t{ds}\t{ch}:{pos}\t{eb}\t{ex_str}\tNA\tNA\tNA\tNA\tNA")

with open(OUT, "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")
print(f"[OUT] {OUT} ({time.time()-t0:.0f}s)")
print("=== 结果（区域最显著 SNP）===")
for l in lines[1:]:
    print("  " + l)
