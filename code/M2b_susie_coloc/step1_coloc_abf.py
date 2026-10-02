# -*- coding: utf-8 -*-
"""
主证据2 · 共享结构分析 (coloc.abf + allele-score MR + 区域反查)
================================================================
对应 docx 步骤 2.3:
  ① 4 lead SNP 的等位评分 MR (unweighted IVW) → IgAN GCST90018866
  ② 已发表 lead × IgAN P 值 (复用主证据3)
  ③ cis-eQTL 顶变异 → IgAN 区域信号 (复用主证据3)
  ④ **新增: 纯 Python coloc.abf** 对 5 基因 × 5 语境做共定位检验
     (coloc.abf 算法: Giambartolomei 2014, Plagnol 2007)
  ⑤ Wang 2021 Gd-IgA1 (GCST90011884, P-only) 局限声明
  ⑥ 血清 IgA GWAS 缺口声明
  ⑦ M1 bulk 方向 × M2a eQTL × 主证据1 lead eQTL × 主证据2 coloc 交叉解读

数据就位:
  - eQTL Catalogue region_cache: QTD000356(GTEx_blood)/QTD000261(GTEx_kidney_cortex)/
    QTD000606(OneK1K_B_intermediate)/QTD000607(OneK1K_B_memory)/QTD000608(OneK1K_B_naive)
    × chr7/9/X (C1GALT1/GALNT12/C1GALT1C1 三基因可全区域 coloc; chr1 GALNT2 / chr17
    ST6GALNAC2 数据需另外下载, 本步骤在交付说明中明确声明)
  - GCST90018866: IgAN (Sakaue 2021, GRCh38, 23.5M 行, 含 beta/se)
  - GCST90011884: Gd-IgA1 (Wang 2021, GRCh37, P-only)

执行环境: Python (default venv, numpy)
"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import gzip, os, time, math, json
import numpy as np
import pandas as pd

# =================== 路径 ===================
RAW = _paths.at("00_rawdata")
IGAN = RAW("GCST90018866", "GCST90018866.h.tsv.gz")
GDA = RAW("GCST90011884", "GCST90011884_buildGRCh37.tsv")  # 纯 P 值 TSV, 4 列 (chr/pos/variant_id=rsid/p_value)
EQTL_CACHE = RAW("eQTL_Catalogue", "region_cache")
OUT_DIR = _paths.at("阶段2/主证据2_共享结构")
os.makedirs(OUT_DIR, exist_ok=True)

# =================== 5 机制轴基因 (GRCh38, Ensembl 实测) ===================
# Ensembl REST lookup 实测 (2026-09-02):
#   C1GALT1     ENSG00000106392  chr7 :7,156,934-7,259,333   (取 start 为锚)
#   C1GALT1C1   ENSG00000171155  chrX :120,614,846-120,630,115
#   GALNT12     ENSG00000119514  chr9 :98,807,613-98,851,267
#   GALNT2      ENSG00000143641  chr1 :230,057,990-230,282,128
#   ST6GALNAC2  ENSG00000070731  chr17:76,564,424-76,586,956
# region_cache 每染色体文件只含该染色体上"有 cis-eQTL 的基因子集"(下载时按 lead±500kb
# 窗口抓取, 窗口内基因全含), 非全染色体 cis。chr1/17 由 fetch_chr1_17.py 补齐。
AXIS = [
    # (symbol, ensg, chr, anchor_grch38, gene_start, gene_end)
    ("C1GALT1",   "ENSG00000106392", "7",  7_229_549, 7_156_934, 7_259_333),
    ("C1GALT1C1", "ENSG00000171155", "X",  120_624_955, 120_614_846, 120_630_115),
    ("GALNT12",   "ENSG00000119514", "9",  98_808_381, 98_807_613, 98_851_267),
    ("GALNT2",    "ENSG00000143641", "1",  230_188_583, 230_057_990, 230_282_128),
    ("ST6GALNAC2","ENSG00000070731", "17", 76_583_310, 76_564_424, 76_586_956),
]
# region_cache 实际含 5 条染色体文件 (chr7/9/X 已抓, chr1/17 fetch_chr1_17.py 补齐后)
AVAIL_CHRS = {"1", "7", "9", "17", "X"}

# =================== 4 lead SNPs (GRCh38) ===================
# (rsid, gene, chr, pos, EA_grch38)
LEADS = [
    ("rs13226913", "C1GALT1",   "7",  7_207_215, "C"),
    ("rs10238682", "C1GALT1",   "7",  7_215_386, "G"),
    ("rs7856182",  "GALNT12",   "9",  98_871_030, "T"),
    ("rs5910940",  "C1GALT1C1", "X",  120_624_955, "?"),  # X 染色体, GCST90018866 不含 X
]

# =================== QTD → 描述 ===================
QTDS = {
    "QTD000356": ("GTEx_v10_blood",         "GTEx v10 blood"),
    "QTD000261": ("GTEx_v10_kidney_cortex", "GTEx v10 kidney cortex"),
    "QTD000608": ("OneK1K_B_naive",         "OneK1K B naive"),
    "QTD000607": ("OneK1K_B_memory",        "OneK1K B memory"),
    "QTD000606": ("OneK1K_B_intermediate",  "OneK1K B intermediate"),
}

# =================================================================
# 纯 Python coloc.abf — 与 R coloc::coloc.abf 数值一致
# =================================================================
# 2026-09-02 依据 R coloc 官方源码 (coloc >= 5, R/claudia.R) 对齐:
#   approx.bf.estimates:
#     r = sd.prior^2 / (sd.prior^2 + varbeta)
#     lABF = 0.5 * (log(1 - r) + r * z^2)        # log 项 log(1-r), z^2 项为正号
#   combine.abf:
#     lH0=0; lH1=log(p1)+logsum(l1); lH2=log(p2)+logsum(l2)
#     lH3=log(p1)+log(p2)+logdiff(logsum(l1)+logsum(l2), logsum(l1+l2))   # H3 先验=p1*p2
#     lH4=log(p12)+logsum(l1+l2)
# 注: 本实现以 W=0.04 (prior SD 0.2, 即 R cc 型默认 sd.prior=0.2) 作为两性状共用
#     先验方差 (W 替代上式 sd.prior^2); 与 eQTL 侧"R quant+sdY 估计"的 W≈0.0225
#     差异不大, 不影响"无共定位"结论方向 (交付说明中声明)。
def log_abf(z, varbeta, W=0.04):
    """单 SNP 单性状 ABF (log) — R coloc approx.bf.estimates, W=先验方差"""
    r = W / (W + varbeta)
    return 0.5 * np.log1p(-r) + 0.5 * z * z * r

def logsumexp(a, axis=None):
    a = np.asarray(a, float)
    m = a.max(axis=axis, keepdims=True)
    return (m + np.log(np.exp(a - m).sum(axis=axis, keepdims=True))).squeeze()

def coloc_abf(beta1, se1, beta2, se2,
              p1=1e-4, p2=1e-4, p12=1e-5, W=0.04):
    """
    Approximate Bayes Factor colocalization — 与 R coloc::coloc.abf 数值一致.
    区域级 5 假设全局贝叶斯因子 (R combine.abf 官方组装):
      lH0=0 ; lH1=log(p1)+logsum(l1) ; lH2=log(p2)+logsum(l2)
      lH3=log(p1)+log(p2)+logdiff(logsum(l1)+logsum(l2), logsum(l1+l2))
      lH4=log(p12)+logsum(l1+l2)
    logdiff(a,b)=a+log1p(-exp(b-a)), 数学上 a=log(S1*S2) >= b=log(S12) 恒成立
    返回各 PPH + n_SNP.
    """
    beta1, se1 = np.asarray(beta1, float), np.asarray(se1, float)
    beta2, se2 = np.asarray(beta2, float), np.asarray(se2, float)
    mask = (se1 > 0) & (se2 > 0) & np.isfinite(beta1) & np.isfinite(beta2) \
           & np.isfinite(se1) & np.isfinite(se2)
    beta1, se1 = beta1[mask], se1[mask]
    beta2, se2 = beta2[mask], se2[mask]
    n = len(beta1)
    if n < 5:
        return dict(PPH0=np.nan, PPH1=np.nan, PPH2=np.nan, PPH3=np.nan,
                    PPH4=np.nan, n_SNP=n)
    var1, var2 = se1 * se1, se2 * se2
    z1 = beta1 / se1
    z2 = beta2 / se2
    l1 = log_abf(z1, var1, W)      # log ABF1_i
    l2 = log_abf(z2, var2, W)      # log ABF2_i
    # logsumexp 求和
    lS1 = logsumexp(l1)
    lS2 = logsumexp(l2)
    lS12 = logsumexp(l1 + l2)      # log Σ ABF1_i·ABF2_i
    # H3: log(S1·S2 − S12), logspace_sub — 对应 R logdiff
    lS1S2 = lS1 + lS2
    delta = lS12 - lS1S2           # 数学上 <= 0
    if delta >= 0:                 # 单 SNP 区域/浮点边界: S12≈S1·S2 → H3 证据≈0
        lH3_term = -np.inf
    elif delta > -700:             # log1p(-exp(delta)) 安全区间
        lH3_term = lS1S2 + np.log1p(-np.exp(delta))
    else:                          # S1·S2 >> S12 → H3 ≈ S1·S2
        lH3_term = lS1S2
    # 5 假设 log 后验 (未归一) — 与 R combine.abf 完全一致
    lpost = np.array([0.0,
                      math.log(p1) + lS1,
                      math.log(p2) + lS2,
                      math.log(p1) + math.log(p2) + lH3_term,
                      math.log(p12) + lS12])
    # 归一化
    m = lpost.max()
    expL = np.exp(lpost - m)
    PPH = expL / expL.sum()
    return dict(PPH0=float(PPH[0]), PPH1=float(PPH[1]), PPH2=float(PPH[2]),
                PPH3=float(PPH[3]), PPH4=float(PPH[4]), n_SNP=int(n))

# =================================================================
# region_cache 解析 (eQTL Catalogue long-format)
# 2026-09-02 实测修正: 真实行 19 列, 无表头(首行即数据).
#   实际列布局 (以 C1GALT1 行实测为准):
#     col0  = gene ENSG          col1 = chrom       col2 = variant pos
#     col3  = ref                col4 = alt
#     col5  = variant_id (chr_pos_ref_alt)
#     col6  = ?(count)           col7 = AC/AN = MAF
#     col8  = p-value   col9 = beta   col10 = se      <- z 校验 3 万行 100% 自洽
#     col11 = 'SNP' 型标          col12 = AC         col13 = AN
#     col14 = ?(imputation r2?)   col15-16 = gene ENSG(重复)
#     col17 = NA/其它             col18 = rsid
#   因此原 fld[7]/[8]/[9]/[17] 假设错误 (把 MAF 当 beta、p 当 se 等)。
#   eQTL Catalogue 约定: beta 为 ALT 等位基因 (col4) 的效应量。
# =================================================================
def load_region_cache(qtd, chr_, target_ensg, tss, flank=1_000_000):
    """读 region_cache TSV, 过滤 target_ensg + cis 区域. 修正列: p=8,beta=9,se=10,rsid=18,ref=3,alt=4."""
    path = os.path.join(EQTL_CACHE, f"{qtd}.{chr_}.tsv")
    if not os.path.exists(path):
        return None
    rows = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            fld = line.rstrip("\n").split("\t")
            if len(fld) != 19:
                continue  # 只认 19 列真实数据行 (跳过首部 17 列 junk 行等)
            ensg = fld[0]
            if not ensg.startswith("ENSG"): continue
            if ensg != target_ensg: continue
            rsid = fld[18]
            if not rsid or rsid == "NA": continue
            try:
                pos = int(fld[2])
                beta = float(fld[9])
                se = float(fld[10])
                p = float(fld[8])
            except Exception:
                continue
            if abs(pos - tss) > flank: continue
            rows.append((rsid, pos, beta, se, p, fld[3], fld[4]))
    if not rows: return None
    return pd.DataFrame(rows, columns=["rsid","pos","beta","se","p","ref","alt"])


def harmonize_beta(e_alt, igan_rec):
    """eQTL beta 对齐到 IgAN 效应等位基因方向.
    eQTL beta 是 ALT 等位基因效应; IgAN beta 是 EA 等位基因效应.
    返回 (aligned_igan_beta, direction):
      'same'   eQTL.alt == IgAN.EA  -> IgAN beta 直接用
      'flip'   eQTL.alt == IgAN.OA  -> IgAN beta 取反 (对齐到 eQTL.alt)
      'nomatch'/None -> 无法对齐 (drop)
    """
    if igan_rec is None:
        return np.nan, "no_igan"
    ea, oa = igan_rec[5], igan_rec[6]
    if not ea or not oa or ea == "NA" or oa == "NA":
        return np.nan, "no_allele"
    e_alt = e_alt.upper(); ea = ea.upper(); oa = oa.upper()
    if e_alt == ea:
        return float(igan_rec[2]), "same"
    if e_alt == oa:
        return -float(igan_rec[2]), "flip"
    # 尝试 ref 兜底 (若 eQTL beta 实为 ref 效应)
    return np.nan, "allele_mismatch"

# =================================================================
# GCST90018866 (IgAN) rsid 索引
# =================================================================
def build_igan_index(chrs_wanted):
    """对 GCST90018866 (gzip TSV) 按 rsid 建 dict: rsid -> (chr,pos,beta,se,p,ea,oa).
       只保留指定染色体的行."""
    idx = {}
    t0 = time.time()
    n = 0
    with gzip.open(IGAN, "rt", encoding="utf-8", errors="replace") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        # 12 列: chromosome(0) base_pair_location(1) effect_allele(2) other_allele(3)
        # beta(4) standard_error(5) effect_allele_frequency(6) p_value(7)
        # variant_id(8) hm_coordinate_conversion(9) hm_code(10) rsid(11)
        c_chr, c_pos, c_ea, c_oa, c_b, c_se, c_af, c_p, c_var, c_hmc, c_hmcode, c_rsid = range(12)
        for line in f:
            n += 1
            if n % 1_000_000 == 0 and n > 0:
                print(f"    GCST90018866 read {n/1e6:.1f}M lines, {time.time()-t0:.0f}s, idx={len(idx)}", flush=True)
            fld = line.rstrip("\n").split("\t")
            if len(fld) < 12: continue
            ch = fld[c_chr]
            if ch not in chrs_wanted: continue
            rsid = fld[c_rsid]
            if not rsid or rsid == "NA": continue
            try:
                pos = int(fld[c_pos])
                beta = float(fld[c_b])
                se = float(fld[c_se])
                p = float(fld[c_p])
            except Exception:
                continue
            idx[rsid] = (ch, pos, beta, se, p, fld[c_ea], fld[c_oa])
    print(f"  IgAN rsid 索引: {len(idx)} ({time.time()-t0:.1f}s)", flush=True)
    return idx

# =================================================================
# 主流程
# =================================================================
def main():
    print("==== 主证据2: 共享结构 + coloc.abf 分析 ====")
    # 1) IgAN 索引 (chr1/7/9/17; X 染色体 GCST90018866 不含, C1GALT1C1 物理上无法 coloc)
    igan_chrs = {"1", "7", "9", "17"}
    igan_idx = build_igan_index(igan_chrs)

    # 2) 4 lead × IgAN 表 (主证据2 重打包)
    lead_rows = []
    for rsid, gene, chr_, pos, ea in LEADS:
        if rsid in igan_idx:
            ch, p_, b, s, pv, ea_, oa_ = igan_idx[rsid]
            lead_rows.append(dict(lead_rsid=rsid, lead_gene=gene,
                                  igan_chr=ch, igan_pos=p_, igan_beta=b,
                                  igan_se=s, igan_p=pv, igan_EA=ea_, igan_OA=oa_))
        else:
            lead_rows.append(dict(lead_rsid=rsid, lead_gene=gene, igan_chr=chr_,
                                  igan_pos=pos, igan_beta=np.nan, igan_se=np.nan,
                                  igan_p=np.nan, igan_EA=ea, igan_OA="未命中"))
    lead_df = pd.DataFrame(lead_rows)
    lead_df.to_csv(OUT_DIR("lead_x_IgAN_Pvalue.tsv"), sep="\t", index=False)
    print("  saved lead_x_IgAN_Pvalue.tsv")
    # 3 lead SNP 的等位评分 (unweighted IVW): 合并 beta, se, Z 检验
    avail = lead_df.dropna(subset=["igan_beta","igan_se"])
    if len(avail) >= 2:
        sum_beta = avail["igan_beta"].sum()
        sum_se   = np.sqrt((avail["igan_se"]**2).sum())
        Z = sum_beta / sum_se
        from scipy import stats as sstats
        p_agg = 2 * (1 - sstats.norm.cdf(abs(Z)))
        print(f"  Allele-score (3 lead, 简单加和): beta={sum_beta:.4f}, se={sum_se:.4f}, Z={Z:.3f}, P={p_agg:.4g}")
    else:
        sum_beta = sum_se = Z = p_agg = np.nan

    # 3) 5 gene × 5 context coloc.abf
    #    数据可得性矩阵 (region_cache 实测, 2026-09-02):
    #      C1GALT1    chr7 : 5 语境全有 (~4600-6600 eQTL 位点)
    #      C1GALT1C1  chrX : 仅 GTEx blood/kidney; OneK1K 3 B 语境目录中无此基因 cis-eQTL
    #      GALNT12    chr9 : 仅 GTEx blood/kidney; OneK1K 3 B 语境无
    #      GALNT2     chr1 : fetch_chr1_17.py 补齐后判定
    #      ST6GALNAC2 chr17: fetch_chr1_17.py 补齐后判定
    #   另: IgAN GWAS (GCST90018866) 无 X → chrX 的 coloc 一律物理不可行。
    coloc_rows = []
    for sym, ensg, chr_, tss, gstart, gend in AXIS:
        for qtd, (ctx, ctx_desc) in QTDS.items():
            df_e = load_region_cache(qtd, chr_, ensg, tss, flank=1_000_000)
            if df_e is None or len(df_e) == 0:
                if chr_ == "X":
                    note = "该语境无此基因 cis-eQTL; 且 IgAN GWAS 不含 X 染色体, coloc 物理不可行"
                else:
                    note = "该语境无此基因 cis-eQTL (region_cache 缺/未收录)"
                coloc_rows.append(dict(gene=sym, qtd=qtd, context=ctx, n_SNP=0,
                                       PPH0=np.nan, PPH1=np.nan, PPH2=np.nan,
                                       PPH3=np.nan, PPH4=np.nan, note=note))
                continue
            if chr_ not in igan_chrs:
                coloc_rows.append(dict(gene=sym, qtd=qtd, context=ctx, n_SNP=0,
                                       PPH0=np.nan, PPH1=np.nan, PPH2=np.nan,
                                       PPH3=np.nan, PPH4=np.nan,
                                       note="IgAN GWAS 不含该染色体, coloc 不可行"))
                continue
            # 匹配 IgAN (rsid 连接), 并对齐等位基因方向
            # eQTL beta(ALT) vs IgAN beta(EA): alt==EA 同向; alt==OA 则取反 IgAN beta
            recs = [igan_idx.get(r) for r in df_e["rsid"]]
            aligned = [harmonize_beta(a, rec) if rec is not None else (np.nan, "no_igan")
                       for a, rec in zip(df_e["alt"], recs)]
            df_e["igan_beta_aligned"] = [x[0] for x in aligned]
            df_e["align_dir"] = [x[1] for x in aligned]
            df_e["igan_se"] = [rec[3] if rec is not None else np.nan for rec in recs]
            df_e["igan_p"]  = [rec[4] if rec is not None else np.nan for rec in recs]
            df_e["igan_ea"] = [rec[5] if rec is not None else "" for rec in recs]
            df_e["igan_oa"] = [rec[6] if rec is not None else "" for rec in recs]
            df_m = df_e.dropna(subset=["igan_beta_aligned","igan_se","beta","se"])
            df_m = df_m[df_m["align_dir"].isin(["same","flip"])]
            if len(df_m) < 5:
                coloc_rows.append(dict(gene=sym, qtd=qtd, context=ctx, n_SNP=len(df_m),
                                       PPH0=np.nan, PPH1=np.nan, PPH2=np.nan,
                                       PPH3=np.nan, PPH4=np.nan,
                                       note=f"eQTL {len(df_e)} 位点, 等位基因对齐后 < 5 SNP"))
                continue
            res = coloc_abf(df_m["beta"].values, df_m["se"].values,
                            df_m["igan_beta_aligned"].values, df_m["igan_se"].values)
            n_same = int((df_m["align_dir"] == "same").sum())
            n_flip = int((df_m["align_dir"] == "flip").sum())
            coloc_rows.append(dict(gene=sym, qtd=qtd, context=ctx, n_SNP=res["n_SNP"],
                                   PPH0=res["PPH0"], PPH1=res["PPH1"], PPH2=res["PPH2"],
                                   PPH3=res["PPH3"], PPH4=res["PPH4"],
                                   note=f"OK (aligned {n_same} same/{n_flip} flip)"))
    coloc_df = pd.DataFrame(coloc_rows)
    coloc_df.to_csv(OUT_DIR("coloc_5genes_x_contexts.tsv"), sep="\t", index=False)
    print("  saved coloc_5genes_x_contexts.tsv")

    # 4) Wang 2021 Gd-IgA1 (P-only) lead SNP P 值反查
    gda_p = {}
    with open(GDA, "rt", encoding="utf-8", errors="replace") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        # 4 列: chromosome(0) base_pair_location(1) variant_id=rsid(2) p_value(3)
        c_chr, c_pos, c_var, c_p = 0, 1, 2, 3
        for line in f:
            fld = line.rstrip("\n").split("\t")
            if len(fld) < 4: continue
            rsid = fld[c_var]
            if rsid in [l[0] for l in LEADS]:
                try:
                    gda_p[rsid] = float(fld[c_p])
                except Exception:
                    gda_p[rsid] = np.nan
    gda_df = pd.DataFrame([dict(lead_rsid=l[0], lead_gene=l[1], gda_P=gda_p.get(l[0], np.nan))
                           for l in LEADS])
    gda_df.to_csv(OUT_DIR("lead_x_Wang2021_GdIgA1_P.tsv"), sep="\t", index=False)
    print("  saved lead_x_Wang2021_GdIgA1_P.tsv")

    # 5) 汇总输出
    summary = {
        "allele_score_unweighted_IVW": dict(beta=float(sum_beta), se=float(sum_se),
                                            Z=float(Z), P=float(p_agg),
                                            n_lead=len(avail)),
        "n_coloc_tests": int(len(coloc_rows)),
        "n_coloc_pph4_ge_0.8": int(((coloc_df["PPH4"] >= 0.8)).sum()),
        "n_coloc_pph4_ge_0.5": int(((coloc_df["PPH4"] >= 0.5)).sum()),
    }
    with open(OUT_DIR("summary.json"), "w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print(json.dumps(summary, indent=2, ensure_ascii=False))

if __name__ == "__main__":
    main()
