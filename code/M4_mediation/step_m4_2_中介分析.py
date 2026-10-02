#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M4-2 中介分析主流程 (公开数据可计算范围)
====================================================
输入:
  - cpg_full_assoc_*.tsv          GoDMC API 全量关联 (5 显著 CpG, GRCh37, 含 a1/a2/clumped)
  - lead_x_mQTL_sigCpG_raw.tsv    lead 级关联
  - GCST90018866.h.tsv.gz         IgAN 全量汇总 (GRCh38, rsid)
  - eQTL region_cache chr7        GTEx blood/kidney + OneK1K B naive/mem/inter (GRCh38, rsid)
输出 (阶段3/M4_中介/):
  M4_输入层_QC.tsv                CpG 注释: 基因特征/TSS 距离/工具变量/lead 统计
  M4_MR_CpG_x_IgAN.tsv            Wald(lead 与 top cis sentinel) + IVW(clumped) → IgAN
  M4_coloc_mQTL_x_IgAN.tsv        per CpG mQTL×IgAN coloc.abf
  M4_coloc_mQTL_x_eQTL.tsv        per CpG × 5 eQTL 语境 coloc.abf
  M4_方向一致性_lead.tsv          lead SNP → CpG 甲基化 vs → C1GALT1 表达 方向对比
方法学声明:
  - GoDMC assoc_meta 仅收录显著 SNP–CpG 对 → coloc 为"单信号同源性指示", 交付时如实声明
  - GCST90011884 (Gd-IgA1) 仅 P 值 → 正式中介比例无法计算, 输出缺口声明
"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import gzip, os, json, math, time
import numpy as np
import pandas as pd

RAW = _paths.at("00_rawdata")
IGAN = RAW("GCST90018866", "GCST90018866.h.tsv.gz")
EQTL_CACHE = RAW("eQTL_Catalogue", "region_cache")
OUT = _paths.at("阶段3/M4_中介")
os.makedirs(OUT, exist_ok=True)

# ---------------- 常量 ----------------
C1GALT1_ENSG = "ENSG00000106392"
# C1GALT1 hg19 (mygene.info 2026-09-03 实测): chr7:7,196,565–7,288,282, 正链 (+)
GENE_HG19 = dict(start=7_196_565, end=7_288_282, tss=7_196_565)
CIS_FLANK = 1_000_000
P_THRESH = 5e-8

# 5 显著 CpG (GoDMC API, GRCh37)
SIG_CPGS = ["cg19603390", "cg19473623", "cg17994788", "cg04827551", "cg16101574"]
CPG_POS = {"cg19603390": 7_222_222, "cg19473623": 7_224_869,
           "cg17994788": 7_261_594, "cg04827551": 7_268_805,
           "cg16101574": 7_291_514}
CPG_N = {"cg19603390": 24988, "cg19473623": 24988, "cg17994788": 28181,
         "cg04827551": 25560, "cg16101574": 25095}
LEADS = {"rs13226913": "C1GALT1", "rs10238682": "C1GALT1",
         "rs7856182": "GALNT12", "rs5910940": "C1GALT1C1"}
# eQTL 语境 (QTD → 描述); region_cache 仅 chr7 (C1GALT1 所在)
QTDS = [("QTD000356", "GTEx_v10_blood"),
        ("QTD000261", "GTEx_v10_kidney_cortex"),
        ("QTD000608", "OneK1K_B_naive"),
        ("QTD000607", "OneK1K_B_memory"),
        ("QTD000606", "OneK1K_B_intermediate")]

# ---------------- coloc.abf (与阶段2 step1 完全一致) ----------------
def log_abf(z, varbeta, W=0.04):
    r = W / (W + varbeta)
    return 0.5 * np.log1p(-r) + 0.5 * z * z * r

def logsumexp(a, axis=None):
    a = np.asarray(a, float)
    m = a.max(axis=axis, keepdims=True)
    return (m + np.log(np.exp(a - m).sum(axis=axis, keepdims=True))).squeeze()

def coloc_abf(beta1, se1, beta2, se2, p1=1e-4, p2=1e-4, p12=1e-5, W=0.04):
    beta1, se1 = np.asarray(beta1, float), np.asarray(se1, float)
    beta2, se2 = np.asarray(beta2, float), np.asarray(se2, float)
    mask = (se1 > 0) & (se2 > 0) & np.isfinite(beta1) & np.isfinite(beta2) \
           & np.isfinite(se1) & np.isfinite(se2)
    beta1, se1, beta2, se2 = beta1[mask], se1[mask], beta2[mask], se2[mask]
    n = len(beta1)
    if n < 5:
        return dict(PPH0=np.nan, PPH1=np.nan, PPH2=np.nan, PPH3=np.nan,
                    PPH4=np.nan, n_SNP=n)
    var1, var2 = se1 ** 2, se2 ** 2
    z1, z2 = beta1 / se1, beta2 / se2
    l1, l2 = log_abf(z1, var1, W), log_abf(z2, var2, W)
    lS1, lS2 = logsumexp(l1), logsumexp(l2)
    lS12 = logsumexp(l1 + l2)
    lS1S2 = lS1 + lS2
    delta = lS12 - lS1S2
    if delta >= 0:
        lH3_term = -np.inf
    elif delta > -700:
        lH3_term = lS1S2 + np.log1p(-np.exp(delta))
    else:
        lH3_term = lS1S2
    lpost = np.array([0.0, math.log(p1) + lS1, math.log(p2) + lS2,
                      math.log(p1) + math.log(p2) + lH3_term, math.log(p12) + lS12])
    m = lpost.max()
    expL = np.exp(lpost - m)
    PPH = expL / expL.sum()
    return dict(PPH0=float(PPH[0]), PPH1=float(PPH[1]), PPH2=float(PPH[2]),
                PPH3=float(PPH[3]), PPH4=float(PPH[4]), n_SNP=int(n))

# ---------------- 加载 ----------------
def load_cpg_full(cpg):
    """读取 cpg_full_assoc_*.tsv → DataFrame; 解析 snp 'chr7:pos:SNP'"""
    p = OUT("cpg_full_assoc_%s.tsv" % cpg)
    df = pd.read_csv(p, sep="\t", dtype=str)
    for col in ["rsid", "cpg", "cistrans", "clumped", "a1", "a2"]:
        df[col] = df[col].fillna("")
    def _pos(s):
        try:
            return int(str(s).split(":")[1])
        except Exception:
            return np.nan
    df["pos37"] = df["snp"].map(_pos)
    for col in ["beta_a1", "se_mre", "pval_mre", "beta_are_a1", "se_are", "pval_are",
                "freq_a1", "samplesize"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df

def load_eqtl_gene(qtd):
    """region_cache: target gene 的 cis-eQTL (beta=ALT). 列布局沿用阶段2实测: p=8,beta=9,se=10,rsid=18,alt=4"""
    p = os.path.join(EQTL_CACHE, "%s.7.tsv" % qtd)
    if not os.path.exists(p):
        return None
    rows = []
    with open(p, encoding="utf-8", errors="replace") as f:
        for line in f:
            fld = line.rstrip("\n").split("\t")
            if len(fld) != 19:
                continue
            if fld[0] != C1GALT1_ENSG:
                continue
            rsid = fld[18]
            if not rsid or rsid == "NA":
                continue
            try:
                rows.append((rsid, int(fld[2]), float(fld[8]), float(fld[9]),
                             float(fld[10]), fld[3], fld[4]))
            except Exception:
                continue
    if not rows:
        return None
    return pd.DataFrame(rows, columns=["rsid", "pos38", "p", "beta", "se", "ref", "alt"])

def build_igan_chr7_index():
    """GCST90018866 → 仅 chr7 rsid 索引: rsid -> (beta,se,p,ea,oa)"""
    idx = {}
    t0 = time.time()
    with gzip.open(IGAN, "rt", encoding="utf-8", errors="replace") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        # 12 列: 0 chr 1 pos 2 EA 3 OA 4 beta 5 se 6 EAF 7 P 8 var_id 9 hm 10 hmcode 11 rsid
        for line in f:
            fld = line.rstrip("\n").split("\t")
            if len(fld) < 12:
                continue
            if fld[0] != "7":
                continue
            rsid = fld[11]
            if not rsid or rsid == "NA":
                continue
            try:
                beta, se = float(fld[4]), float(fld[5])
                p = float(fld[7])
            except Exception:
                continue
            idx[rsid] = (beta, se, p, fld[2], fld[3])
    print("  IgAN chr7 rsid 索引: %d (%ds)" % (len(idx), time.time() - t0), flush=True)
    return idx

def align_to_allele(a1, a2, beta_a1, ref_allele):
    """把 mQTL beta_a1 对齐到 ref_allele 方向. ref_allele 可=EA 或 eQTL alt.
    返回 (aligned_beta, ok)"""
    a1u, a2u = a1.upper(), a2.upper()
    ra = ref_allele.upper()
    if a1u == ra:
        return beta_a1, True
    if a2u == ra:
        return -beta_a1, True
    return np.nan, False

# ================= 主流程 =================
def main():
    print("==== M4-2 中介分析 (公开数据可计算范围) ====", flush=True)
    igan = build_igan_chr7_index()

    # 预载 eQTL
    eqtl = {qtd: load_eqtl_gene(qtd) for qtd, _ in QTDS}
    for qtd, desc in QTDS:
        d = eqtl[qtd]
        print("  eQTL %s(%s): %d 位点" % (qtd, desc, 0 if d is None else len(d)), flush=True)

    # lead 级 mQTL 原始 (取 lead 行作为 IV 信息)
    lead_raw = pd.read_csv(OUT("lead_x_mQTL_sigCpG_raw.tsv"), sep="\t", dtype=str)

    # 输入层 QC 行
    qc_rows, mr_rows, coloc_igan_rows, coloc_eqtl_rows, dir_rows = [], [], [], [], []

    for cpg in SIG_CPGS:
        pos = CPG_POS[cpg]
        df = load_cpg_full(cpg)
        df["cis"] = df["pos37"].sub(pos).abs() <= CIS_FLANK
        dfm = df[(df["rsid"] != "") & df["cis"].fillna(False)].copy()

        # ---- 基因特征注释 ----
        s, e, tss = GENE_HG19["start"], GENE_HG19["end"], GENE_HG19["tss"]
        if pos < s - 1500:
            feat = "5'上游>1.5kb"; dtss = pos - tss
        elif pos <= s + 1500:
            feat = "启动子TSS±1.5kb"; dtss = pos - tss
        elif pos <= e:
            feat = "基因体内(intragenic)"; dtss = pos - tss
        elif pos <= e + 5000:
            feat = "3'下游<5kb"; dtss = pos - tss
        else:
            feat = "3'下游>5kb"; dtss = pos - tss

        # 显著 IV (cis, P<5e-8)
        iv = dfm[dfm["pval_mre"].fillna(1) < P_THRESH]
        iv_clumped = iv[iv["clumped"].astype(str) == "True"]
        lead_rows_here = dfm[dfm["rsid"].isin(LEADS)]

        # ---- lead 级工具: 取该 CpG 的 lead mQTL 行 ----
        lead_iv = {}
        for _, r in lead_rows_here.iterrows():
            lead_iv.setdefault(r["rsid"], r)

        # 与 IgAN 合并
        dfm["igan_beta"], dfm["igan_se"], dfm["igan_p"] = np.nan, np.nan, np.nan
        dfm["igan_ea"], dfm["igan_oa"] = "", ""
        recs = [igan.get(r) for r in dfm["rsid"]]
        dfm["igan_beta"] = [x[0] if x else np.nan for x in recs]
        dfm["igan_se"] = [x[1] if x else np.nan for x in recs]
        dfm["igan_p"] = [x[2] if x else np.nan for x in recs]
        dfm["igan_ea"] = [x[3] if x else "" for x in recs]
        dfm["igan_oa"] = [x[4] if x else "" for x in recs]

        # mQTL beta 对齐到 IgAN EA 方向 → coloc 输入
        al = dfm.apply(lambda r: align_to_allele(r["a1"], r["a2"], r["beta_a1"],
                                                 r["igan_ea"] if r["igan_ea"] else ""), axis=1)
        dfm["bm_EA"] = [x[0] for x in al]
        dfm["align_ok"] = [x[1] for x in al]
        dfc = dfm[(dfm["align_ok"]) & dfm["igan_se"].notna()].copy()
        dfc = dfc[dfc["igan_se"] > 0]

        # ---- coloc mQTL × IgAN ----
        if len(dfc) >= 5:
            res = coloc_abf(dfc["bm_EA"].values, dfc["se_mre"].values,
                            dfc["igan_beta"].values, dfc["igan_se"].values)
            top_m = dfm.loc[dfm["pval_mre"].fillna(1).idxmin()]
            top_i = dfc.loc[dfc["igan_p"].fillna(1).idxmin()]
            coloc_igan_rows.append(dict(cpg=cpg, cpg_pos37=pos, feature=feat, n_shared=res["n_SNP"],
                                        mqtl_top=top_m["rsid"], mqtl_top_P=top_m["pval_mre"],
                                        igan_top=top_i["rsid"], igan_top_P=top_i["igan_p"],
                                        PPH0=res["PPH0"], PPH1=res["PPH1"], PPH2=res["PPH2"],
                                        PPH3=res["PPH3"], PPH4=res["PPH4"]))
        else:
            coloc_igan_rows.append(dict(cpg=cpg, cpg_pos37=pos, feature=feat, n_shared=len(dfc),
                                        mqtl_top="", mqtl_top_P=np.nan, igan_top="", igan_top_P=np.nan,
                                        PPH0=np.nan, PPH1=np.nan, PPH2=np.nan, PPH3=np.nan, PPH4=np.nan))

        # ---- MR: CpG → IgAN ----
        # (i) Wald, 工具=lead SNP (存在时)
        for rs, gene in LEADS.items():
            if rs not in lead_iv:
                continue
            r = lead_iv[rs]
            igan_rec = igan.get(rs)
            if igan_rec is None:
                continue
            b_m, ok = align_to_allele(r["a1"], r["a2"], float(r["beta_a1"]), igan_rec[3])
            if not ok:
                continue
            b_y, se_y = igan_rec[0], igan_rec[1]
            mr_beta = b_y / b_m
            mr_se = se_y / abs(b_m)
            mr_p = 2 * (1 - _norm_cdf(abs(mr_beta / mr_se)))
            mr_rows.append(dict(cpg=cpg, method="Wald_lead", instrument=rs, gene=gene,
                                b_mqtl=b_m, se_mqtl=float(r["se_mre"]),
                                b_igan=b_y, se_igan=se_y, mr_beta=mr_beta, mr_se=mr_se, mr_p=mr_p))
        # (ii) Wald, 工具=cis top sentinel
        if len(dfm) > 0:
            tm = dfm.loc[dfm["pval_mre"].fillna(1).idxmin()]
            igan_rec = igan.get(tm["rsid"])
            if igan_rec is not None:
                b_m, ok = align_to_allele(tm["a1"], tm["a2"], float(tm["beta_a1"]), igan_rec[3])
                if ok:
                    b_y, se_y = igan_rec[0], igan_rec[1]
                    mr_beta = b_y / b_m; mr_se = se_y / abs(b_m)
                    mr_p = 2 * (1 - _norm_cdf(abs(mr_beta / mr_se)))
                    mr_rows.append(dict(cpg=cpg, method="Wald_topcis", instrument=tm["rsid"],
                                        gene="", b_mqtl=b_m, se_mqtl=float(tm["se_mre"]),
                                        b_igan=b_y, se_igan=se_y, mr_beta=mr_beta, mr_se=mr_se, mr_p=mr_p))
        # (iii) IVW over clumped cis IVs (仅当 ≥3 且都与 IgAN 对齐)
        ivc = iv_clumped[(iv_clumped["rsid"] != "")]
        ivals = []
        for _, r in ivc.iterrows():
            igr = igan.get(r["rsid"])
            if igr is None:
                continue
            b_m, ok = align_to_allele(r["a1"], r["a2"], float(r["beta_a1"]), igr[3])
            if ok and abs(b_m) > 1e-8:
                ivals.append((b_m, igr[0], igr[1]))
        if len(ivals) >= 3:
            bms = np.array([x[0] for x in ivals]); bys = np.array([x[1] for x in ivals])
            sys_ = np.array([x[2] for x in ivals])
            w = (bms ** 2) / (sys_ ** 2)
            ivw = float((w * (bys / bms)).sum() / w.sum())
            ivw_se = float(np.sqrt(1 / w.sum()))
            ivw_p = 2 * (1 - _norm_cdf(abs(ivw / ivw_se)))
            mr_rows.append(dict(cpg=cpg, method="IVW_clumped_%d" % len(ivals),
                                instrument=";".join([r["rsid"] for _, r in ivc.iterrows()][:20]),
                                gene="", b_mqtl=np.nan, se_mqtl=np.nan, b_igan=np.nan,
                                se_igan=np.nan, mr_beta=ivw, mr_se=ivw_se, mr_p=ivw_p))
        else:
            mr_rows.append(dict(cpg=cpg, method="IVW_clumped", instrument="",
                                gene="", b_mqtl=np.nan, se_mqtl=np.nan, b_igan=np.nan,
                                se_igan=np.nan, mr_beta=np.nan, mr_se=np.nan, mr_p=np.nan))

        # ---- coloc mQTL × eQTL (per dataset) ----
        for qtd, desc in QTDS:
            de = eqtl[qtd]
            if de is None or len(de) == 0:
                coloc_eqtl_rows.append(dict(cpg=cpg, qtd=qtd, context=desc, n_shared=0,
                                            mqtl_top="", eqtl_top="", eqtl_top_P=np.nan,
                                            PPH0=np.nan, PPH1=np.nan, PPH2=np.nan,
                                            PPH3=np.nan, PPH4=np.nan, note="无 C1GALT1 cis-eQTL 数据"))
                continue
            # mQTL beta 对齐到 eQTL alt 方向 (eQTL 侧列加前缀避免与 GoDMC 同名列冲突)
            de2 = de[["rsid", "beta", "se", "p", "alt"]].rename(columns={
                "beta": "e_beta", "se": "e_se", "p": "e_p", "alt": "e_alt"})
            m2 = dfm.merge(de2, on="rsid", how="inner")
            if len(m2) == 0:
                coloc_eqtl_rows.append(dict(cpg=cpg, qtd=qtd, context=desc, n_shared=0,
                                            mqtl_top="", eqtl_top="", eqtl_top_P=np.nan,
                                            PPH0=np.nan, PPH1=np.nan, PPH2=np.nan,
                                            PPH3=np.nan, PPH4=np.nan, note="mQTL∩eQTL 无共享 rsid"))
                continue
            al2 = m2.apply(lambda r: align_to_allele(r["a1"], r["a2"], r["beta_a1"], r["e_alt"]), axis=1)
            m2["bm_alt"] = [x[0] for x in al2]
            m2["ok2"] = [x[1] for x in al2]
            m2c = m2[(m2["ok2"]) & (m2["e_se"] > 0) & m2["se_mre"].notna() & (m2["se_mre"] > 0)]
            if len(m2c) >= 5:
                res = coloc_abf(m2c["bm_alt"].values, m2c["se_mre"].values,
                                m2c["e_beta"].values, m2c["e_se"].values)
                etop = m2c.loc[m2c["e_p"].fillna(1).idxmin()]
                coloc_eqtl_rows.append(dict(cpg=cpg, qtd=qtd, context=desc,
                                            n_shared=res["n_SNP"], mqtl_top="",
                                            eqtl_top=etop["rsid"], eqtl_top_P=etop["e_p"],
                                            PPH0=res["PPH0"], PPH1=res["PPH1"], PPH2=res["PPH2"],
                                            PPH3=res["PPH3"], PPH4=res["PPH4"],
                                            note="OK"))
            else:
                coloc_eqtl_rows.append(dict(cpg=cpg, qtd=qtd, context=desc, n_shared=len(m2c),
                                            mqtl_top="", eqtl_top="", eqtl_top_P=np.nan,
                                            PPH0=np.nan, PPH1=np.nan, PPH2=np.nan,
                                            PPH3=np.nan, PPH4=np.nan, note="对齐后 <5 SNP"))

        # ---- 输入层 QC / 方向一致性 (lead → CpG vs lead → C1GALT1 表达) ----
        n_iv_cis = int(len(iv))
        n_iv_cl = int(len(iv_clumped))
        top = dfm.loc[dfm["pval_mre"].fillna(1).idxmin()] if len(dfm) else None
        for rs, gene in LEADS.items():
            if rs not in lead_iv:
                continue
            r = lead_iv[rs]
            igr = igan.get(rs)
            # eQTL lead 行 (C1GALT1, 该 rsid)
            for qtd, desc in QTDS:
                de = eqtl[qtd]
                if de is None:
                    continue
                sub = de[de["rsid"] == rs]
                if len(sub) == 0:
                    continue
                e = sub.iloc[0]
                # mQTL 对齐到 eQTL alt
                bm_e, ok = align_to_allele(r["a1"], r["a2"], float(r["beta_a1"]), e["alt"])
                if not ok:
                    continue
                dir_rows.append(dict(lead=rs, gene=gene, cpg=cpg, context=desc,
                                     mqtl_beta_on_eQTLalt=bm_e, mqtl_P=float(r["pval_mre"]),
                                     eqtl_beta=e["beta"], eqtl_se=e["se"], eqtl_P=e["p"],
                                     sign_consistent=(bm_e * e["beta"] > 0)))
        qc_rows.append(dict(cpg=cpg, chr="7", pos37=pos, feature=feat, dist_TSS=dtss,
                            gene="C1GALT1", n_full_assoc=len(df), n_cis=len(dfm),
                            n_sig_cis_IV= n_iv_cis, n_clumped_IV=n_iv_cl,
                            top_cis_rsid="" if top is None else top["rsid"],
                            top_cis_P=np.nan if top is None else top["pval_mre"],
                            top_cis_beta=np.nan if top is None else top["beta_a1"],
                            lead_rsids=",".join(sorted(lead_iv.keys())),
                            sample_n=CPG_N[cpg], note=""))
        print("  [%s] %s feat=%s n_cis=%d clumpedIV=%d" % (cpg, feat, feat, len(dfm), n_iv_cl), flush=True)

    # ---------------- 落盘 ----------------
    pd.DataFrame(qc_rows).to_csv(OUT("M4_输入层_QC.tsv"), sep="\t", index=False)
    pd.DataFrame(mr_rows).to_csv(OUT("M4_MR_CpG_x_IgAN.tsv"), sep="\t", index=False)
    pd.DataFrame(coloc_igan_rows).to_csv(OUT("M4_coloc_mQTL_x_IgAN.tsv"), sep="\t", index=False)
    pd.DataFrame(coloc_eqtl_rows).to_csv(OUT("M4_coloc_mQTL_x_eQTL.tsv"), sep="\t", index=False)
    pd.DataFrame(dir_rows).to_csv(OUT("M4_方向一致性_lead.tsv"), sep="\t", index=False)
    print("  saved 5 files", flush=True)
    print("  [缺口声明] GCST90011884 (Gd-IgA1) P-only → 正式中介比例需全量效应量 (Kiryluk dbGaP / 作者)", flush=True)
    print("[DONE] M4-2", flush=True)

def _norm_cdf(x):
    from scipy import stats
    return stats.norm.cdf(x)

if __name__ == "__main__":
    main()
