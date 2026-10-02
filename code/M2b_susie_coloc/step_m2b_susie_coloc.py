# -*- coding: utf-8 -*-
"""
M2b · SuSiE fine-mapping 级共定位复核 (coloc.susie 语义, 纯 Python)
==================================================================
目的: 对主证据2 (coloc.abf, 单因果变异假设) 的 14 组区域共定位结果,
追加 SuSiE fine-mapping 复核 —— 验证放松"每性状区域仅 1 个因果变异"
假设后, PP.H3 / PP.H4 后验是否稳健 (仍低 → "无共定位"结论稳健)。

方法 (对照 R coloc::coloc.susie / coloc.bf_bf):
  R coloc.susie 的逐对检验数学 = 把 SuSiE 每个独立信号当作一个
  "候选因果分布" (信号级 LBF 向量), 两性状各取一个信号, 在共同 SNP
  上做 combine.abf (H0-H4 贝叶斯因子组装)。本脚本:

  eQTL 侧 (真实 SuSiE fine-map, eQTL Catalogue 官方 credible_sets):
    对每个 cs_id (=1 个独立信号), 其 CS 成员即 SuSiE 后验支持的
    因果候选变异集。逐 SNP 信号级证据量级用 credible_sets 官方发布的
    边际 z/se 重建 Wakefield ABF (log_abf, W=0.04) —— 这正是 SuSiE
    lbf_variable 的输入统计量; 对解析良好的独立信号, 信号级 lbf ≈
    边际 ABF (逐信号加性常数在 H1/H3/H4 与 H0 之外一致, 故量级可比
    coloc.abf)。CS 之外 SNP 的 l1 置 -1e3 (≈0, 即该信号不允许因果
    变异落于其 CS 之外 → SuSiE 定位信息生效)。
    → 语义: "若 eQTL 真信号位于该信号 CS 内, 与 IgAN 是否共享因果?"
    (注: 不用 log(pip) 作 l1 —— pip 是条件概率, 丢失"信号 vs 无效"
    证据量级, 会把强 eQTL 的后验质量误归 H0; 本版定位用 CS, 量级用
    边际 ABF, 与主证据2 coloc.abf 同公式/同先验/同 IgAN 背景 → 干净
    隔离"单因果 vs 信号级"这一变量。)
  IgAN 侧 (无 SuSiE, 单信号):
    与主证据2 coloc.abf 完全相同 —— 区域(±1Mb)共享 SNP 的边际
    Wakefield ABF (R approx.bf.estimates, W=0.04), 即 coloc.abf 的
    单因果假设背景。设计上**仅改变 eQTL 侧模型**(区域边际 ABF →
    信号级 SuSiE pip), IgAN 背景/SNP 全集/先验全同 → 干净隔离
    "单因果假设"这一变量的影响。
  先验: p1=p2=1e-4, p12=1e-5 (与主证据2 完全一致, 隔离模型差异);
        另附 p12=5e-6 (R coloc.susie 默认) 敏感性列。

数据 (2026-09-02 下载, gzip 校验通过):
  00_rawdata/eQTL_Catalogue/susie/{QTD}.credible_sets.tsv.gz
  仅 4 个 基因×语境 组合有官方 SuSiE CS (其余 0 CS, 未达 eQTL Cat
  全转录组 SuSiE 显著门槛 → 无信号可供复核):
    C1GALT1  × GTEx_v10_blood        (QTD000356): L1+L2 两个独立信号
    ST6GALNAC2× GTEx_v10_blood       (QTD000356): L1+L2+L3 三个独立信号
    C1GALT1  × GTEx_v10_kidney_cortex(QTD000261): L1 单信号
    GALNT2   × OneK1K_B_naive        (QTD000608): L1 单信号
  版本一致性: susie QTD 编号与主证据2 region_cache 同为 GTEx_v10/
  OneK1K 语境 → 无 r6/r8 混用问题 (CS lead 与区域 top 信号一致)。

输出:
  M2b_susie_coloc_recheck.tsv  逐 (信号×IgAN) 配对复核表 (coloc.susie 语义)
  M2b_susie_coloc_verdict.tsv  按基因×语境汇总的稳健性结论
  figures/M2b_susie_coloc_PPH.png/pdf  PP.H 对比图 (coloc.abf vs coloc.susie)
执行: python step_m2b_susie_coloc.py
"""
import gzip, os, sys, time, math, json
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import step1_coloc_abf as s1          # 复用 log_abf / region_cache / igan 索引
from step1_coloc_abf import log_abf, load_region_cache, build_igan_index, AXIS

SUSIE = os.path.join(s1.RAW, "eQTL_Catalogue", "susie")
OUT_DIR = HERE
FIG_DIR = os.path.join(OUT_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# 语境 → (qtd, 染色体) — 4 个有 CS 的组合 (实测 2026-09-02)
COMBOS = [
    # (gene, ensg, qtd, chr, tss_anchor, cs_id_list, context描述)
    ("C1GALT1",    "ENSG00000106392", "QTD000356", "7",  7_229_549),
    ("ST6GALNAC2", "ENSG00000070731", "QTD000356", "17", 76_583_310),
    ("C1GALT1",    "ENSG00000106392", "QTD000261", "7",  7_229_549),
    ("GALNT2",     "ENSG00000143641", "QTD000608", "1",  230_188_583),
]
QTD_DESC = {
    "QTD000356": "GTEx_v10_blood",
    "QTD000261": "GTEx_v10_kidney_cortex",
    "QTD000608": "OneK1K_B_naive",
}

# =================================================================
# combine.abf 核心 (与 R coloc 一致; 接受任意 log-BF 向量 l1/l2)
# =================================================================
def logsumexp(a):
    a = np.asarray(a, float)
    m = a.max()
    return m + math.log(np.exp(a - m).sum())

def combine_abf_pp(l1, l2, p1=1e-4, p2=1e-4, p12=1e-5):
    """给定两性状的逐 SNP log-BF 向量 (同长), 输出 H0-H4 后验.
    数学与 s1.coloc_abf 内嵌的 R combine.abf 组装完全一致."""
    l1 = np.asarray(l1, float); l2 = np.asarray(l2, float)
    n = len(l1)
    if n < 5:
        return dict(PPH0=np.nan, PPH1=np.nan, PPH2=np.nan, PPH3=np.nan,
                    PPH4=np.nan, n_SNP=n)
    lS1 = logsumexp(l1)
    lS2 = logsumexp(l2)
    lS12 = logsumexp(l1 + l2)                 # log Σ ABF1·ABF2
    lS1S2 = lS1 + lS2
    delta = lS12 - lS1S2                       # ≤0
    if delta >= 0:
        lH3_term = -np.inf
    elif delta > -700:
        lH3_term = lS1S2 + np.log1p(-np.exp(delta))
    else:
        lH3_term = lS1S2
    lpost = np.array([0.0,
                      math.log(p1) + lS1,
                      math.log(p2) + lS2,
                      math.log(p1) + math.log(p2) + lH3_term,
                      math.log(p12) + lS12])
    m = lpost.max()
    PPH = np.exp(lpost - m); PPH = PPH / PPH.sum()
    return dict(PPH0=float(PPH[0]), PPH1=float(PPH[1]), PPH2=float(PPH[2]),
                PPH3=float(PPH[3]), PPH4=float(PPH[4]), n_SNP=int(n))


# =================================================================
# credible_sets 解析 → 每 CS 的 {variant_key: pip}
# =================================================================
def parse_variant(vid):
    """'chr7_7223408_C_G' → (7, 7223408, 'C', 'G'); 兼容 indel (T_TTA/TAC_T)."""
    p = vid.split("_")
    return p[0].replace("chr", ""), int(p[1]), p[2].upper(), p[3].upper()

def load_cs_pip(qtd, ensg):
    """读 credible_sets, 返回 {cs_id: {'rows':[dict], 'pip':{key:pip}}}.
    key=(pos,ref,alt)."""
    fn = os.path.join(SUSIE, f"{qtd}.credible_sets.tsv.gz")
    out = {}
    with gzip.open(fn, "rt", encoding="utf-8", errors="replace") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            fld = line.rstrip("\n").split("\t")
            if len(fld) != len(hdr): continue
            rec = dict(zip(hdr, fld))
            if rec["gene_id"] != ensg: continue
            cs = rec["cs_id"]
            d = out.setdefault(cs, {"rows": [], "pip": {}})
            d["rows"].append(rec)
    for cs, d in out.items():
        for rec in d["rows"]:
            ch, pos, ref, alt = parse_variant(rec["variant"])
            d["pip"][(pos, ref, alt)] = float(rec["pip"])
    return out


# =================================================================
# 主流程
# =================================================================
def main():
    t0 = time.time()
    print("==== M2b: SuSiE fine-mapping 级共定位复核 ====", flush=True)

    # 1) IgAN rsid 索引 (chr1/7/17)
    igan_idx = build_igan_index({"1", "7", "17"})

    # 2) 主证据2 coloc.abf 基线 (用于对照)
    base = pd.read_csv(os.path.join(OUT_DIR, "coloc_5genes_x_contexts.tsv"),
                       sep="\t")
    base = base.set_index(["gene", "qtd"])

    rows = []
    for gene, ensg, qtd, chr_, tss in COMBOS:
        ctx = QTD_DESC[qtd]
        # --- eQTL 区域共享 SNP 集 (与主证据2 完全一致) ---
        df_e = load_region_cache(qtd, chr_, ensg, tss, flank=1_000_000)
        if df_e is None or len(df_e) == 0:
            print(f"  [skip] {gene}×{ctx}: 无 region_cache", flush=True)
            continue
        recs = [igan_idx.get(r) for r in df_e["rsid"]]
        aligned = [s1.harmonize_beta(a, rec) if rec is not None else (np.nan, "no_igan")
                   for a, rec in zip(df_e["alt"], recs)]
        df_e["igan_beta_aligned"] = [x[0] for x in aligned]
        df_e["align_dir"] = [x[1] for x in aligned]
        df_e["igan_se"] = [rec[3] if rec is not None else np.nan for rec in recs]
        df_e["igan_p"]  = [rec[4] if rec is not None else np.nan for rec in recs]
        df_m = df_e.dropna(subset=["igan_beta_aligned", "igan_se", "beta", "se"])
        df_m = df_m[df_m["align_dir"].isin(["same", "flip"])].copy()
        df_m = df_m[(df_m["igan_se"] > 0) & (df_m["se"] > 0)]
        df_m = df_m.reset_index(drop=True)
        print(f"  [{gene}×{ctx}] 区域共享 SNP 集: {len(df_m)}", flush=True)

        # --- SuSiE CS ---
        css = load_cs_pip(qtd, ensg)
        if not css:
            rows.append(dict(gene=gene, context=ctx, qtd=qtd, cs_id="(无CS)",
                             note="eQTL 侧无 SuSiE CS, 无信号可复核"))
            continue

        # 区域 top IgAN 变异 (描述 hit2 用)
        i_order = np.argsort(df_m["igan_p"].values)
        top_ig = df_m.iloc[i_order[0]]

        for cs, d in sorted(css.items()):
            cs_rows = sorted(d["rows"], key=lambda r: -float(r["pip"]))
            lead = cs_rows[0]
            # 逐共享 SNP 赋 l1: CS 成员 = 边际 Wakefield ABF (credible_sets 的
            # z/se, 信号级证据量级), 其余 -1e3 (≈0 → 该信号因果不允许出 CS)
            cs_zse = {}
            for rec in d["rows"]:
                ch, pos, ref, alt = parse_variant(rec["variant"])
                cs_zse[(int(pos), ref, alt)] = (float(rec["z"]), float(rec["se"]))
            l1 = np.full(len(df_m), -1e3)
            cs_hit = np.zeros(len(df_m), bool)
            for k, (pos, ref, alt) in enumerate(zip(df_m["pos"].astype(int),
                                                    df_m["ref"].str.upper(),
                                                    df_m["alt"].str.upper())):
                key = (int(pos), ref, alt)
                if key in cs_zse:
                    zz, ssee = cs_zse[key]
                    if ssee > 0:
                        l1[k] = log_abf(zz, ssee * ssee, W=0.04)
                        cs_hit[k] = True
            n_cs_in_shared = int(cs_hit.sum())
            if n_cs_in_shared == 0:
                rows.append(dict(gene=gene, ensg=ensg, qtd=qtd, context=ctx,
                                 cs_id=cs, cs_rank=cs.rsplit("_", 1)[-1],
                                 cs_size=len(cs_rows),
                                 note="CS 成员无一对上 IgAN 区域共享 SNP 集, 无法复核"))
                continue
            key_set = set()
            for pos, ref, alt in zip(df_m["pos"].astype(int),
                                     df_m["ref"].str.upper(), df_m["alt"].str.upper()):
                key_set.add((int(pos), ref, alt))
            cs_pip_total = float(sum(d["pip"].values()))
            cs_pip_shared = float(sum(p for key, p in d["pip"].items()
                                      if (int(key[0]), key[1].upper(), key[2].upper())
                                      in key_set))
            cs_pip_frac = cs_pip_shared / cs_pip_total if cs_pip_total > 0 else np.nan

            # IgAN 侧 l2 (Wakefield ABF, W=0.04; z² 符号无关, 无需再对齐)
            se2 = df_m["igan_se"].values
            z2 = df_m["igan_beta_aligned"].values / se2
            l2 = log_abf(z2, se2 ** 2, W=0.04)

            # 逐对 coloc (p12=1e-5 主比较; p12=5e-6 R coloc.susie 默认敏感性)
            pp = combine_abf_pp(l1, l2, p1=1e-4, p2=1e-4, p12=1e-5)
            pp5 = combine_abf_pp(l1, l2, p1=1e-4, p2=1e-4, p12=5e-6)

            # H4 驱动变异 = argmax(l1+l2) (coloc 的 SNP 级 PP.H4 主贡献者)
            w = l1 + l2
            drv = int(np.argmax(w))
            # CS 内 IgAN 证据: 在 CS 成员中挑 IgAN p 最小 / l2 最大
            cs_l2 = np.where(cs_hit, l2, -np.inf)
            cs_best = int(np.argmax(cs_l2)) if cs_hit.any() else -1

            b_key = (gene, qtd)
            base_row = base.loc[b_key] if b_key in base.index else None
            abf = (dict(PPH1=float(base_row["PPH1"]), PPH3=float(base_row["PPH3"]),
                        PPH4=float(base_row["PPH4"]), n_SNP=int(base_row["n_SNP"]))
                   if base_row is not None else {})

            rows.append(dict(
                gene=gene, ensg=ensg, qtd=qtd, context=ctx,
                cs_id=cs, cs_rank=cs.rsplit("_", 1)[-1],
                cs_size=len(cs_rows), cs_min_r2=float(lead["cs_min_r2"]),
                lead_variant=lead["variant"], lead_rsid=lead["rsid"],
                lead_pip=float(lead["pip"]), lead_eqtl_p=float(lead["pvalue"]),
                n_region_shared=len(df_m),
                n_cs_in_shared=n_cs_in_shared, cs_pip_in_shared=cs_pip_shared,
                cs_pip_frac=round(cs_pip_frac, 4),
                top_igan_rsid=str(top_ig["rsid"]), top_igan_p=float(top_ig["igan_p"]),
                cs_best_igan_rsid=(str(df_m.iloc[cs_best]["rsid"]) if cs_best >= 0 else "NA"),
                cs_best_igan_p=(float(df_m.iloc[cs_best]["igan_p"]) if cs_best >= 0 else np.nan),
                driver_variant=(f"chr{chr_}_{df_m.iloc[drv]['pos']}_{df_m.iloc[drv]['ref']}_{df_m.iloc[drv]['alt']}"
                                if df_m.iloc[drv]["rsid"] in igan_idx else "NA"),
                driver_rsid=str(df_m.iloc[drv]["rsid"]),
                driver_is_cs=bool(cs_hit[drv]),
                driver_igan_p=float(df_m.iloc[drv]["igan_p"]),
                PPH0=pp["PPH0"], PPH1=pp["PPH1"], PPH2=pp["PPH2"],
                PPH3=pp["PPH3"], PPH4=pp["PPH4"],
                PPH4_p12_5e6=pp5["PPH4"],
                abf_PPH1=abf.get("PPH1", np.nan), abf_PPH3=abf.get("PPH3", np.nan),
                abf_PPH4=abf.get("PPH4", np.nan), abf_n_SNP=abf.get("n_SNP", np.nan),
                note="OK"))
            print(f"    {cs}: size={len(cs_rows)} shared={n_cs_in_shared}/"
                  f"{len(df_m)} pip_frac={cs_pip_frac:.3f} "
                  f"PPH4={pp['PPH4']:.4g} (abf {abf.get('PPH4',np.nan):.4g}) "
                  f"PPH3={pp['PPH3']:.4g} (abf {abf.get('PPH3',np.nan):.4g})",
                  flush=True)

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(OUT_DIR, "M2b_susie_coloc_recheck.tsv"),
               sep="\t", index=False)
    print(f"\n  saved M2b_susie_coloc_recheck.tsv  ({len(out)} 行, {time.time()-t0:.0f}s)",
          flush=True)

    # ---- 3) 稳健性结论表 (按 基因×语境 汇总) ----
    def bestmax(x):
        x = pd.to_numeric(x, errors="coerce")
        return x.max() if len(x) else np.nan
    g = out.groupby(["gene", "qtd", "context"], as_index=False).agg(
        n_cs=("cs_id", "size"),
        max_PPH3_pair=("PPH3", bestmax),
        max_PPH4_pair=("PPH4", bestmax),
        max_PPH4_p12_5e6=("PPH4_p12_5e6", bestmax),
        abf_PPH1=("abf_PPH1", "first"),
        abf_PPH3=("abf_PPH3", "first"),
        abf_PPH4=("abf_PPH4", "first"),
    )
    g["robust_no_coloc"] = (g["max_PPH4_pair"] < 0.5) & (g["max_PPH4_pair"] == g["max_PPH4_pair"])
    g["H4_delta_susie_vs_abf"] = g["max_PPH4_pair"] - g["abf_PPH4"]
    g.to_csv(os.path.join(OUT_DIR, "M2b_susie_coloc_verdict.tsv"),
             sep="\t", index=False)
    print("  saved M2b_susie_coloc_verdict.tsv")
    print(g.to_string(index=False))

    # ---- 4) 汇总 json ----
    summary = dict(
        n_combos_with_cs=len(g),
        n_pair_tests=int(len(out)),
        all_pairs_PPH4_below_0_05=bool((out["PPH4"] < 0.05).all()),
        all_pairs_PPH4_below_0_5=bool((out["PPH4"] < 0.5).all()),
        max_PPH4_across_pairs=float(out["PPH4"].max()),
        max_PPH3_across_pairs=float(out["PPH3"].max()),
    )
    with open(os.path.join(OUT_DIR, "summary_m2b.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print(json.dumps(summary, indent=2, ensure_ascii=False))

    # ---- 5) 对比图 ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        _plot(out, g)
    except Exception as e:
        print(f"  [warn] 绘图失败: {e}", flush=True)
    print(f"  done in {time.time()-t0:.0f}s", flush=True)


def _short_ctx(ctx):
    m = {"GTEx_v10_blood": "GTEx blood", "GTEx_v10_kidney_cortex": "GTEx kidney",
         "OneK1K_B_naive": "OneK1K B naive", "OneK1K_B_memory": "OneK1K B memory",
         "OneK1K_B_intermediate": "OneK1K B intermed"}
    return m.get(ctx, ctx)


def _plot(out, g):
    """Supplementary figures (English, SCI-ready):
    A. per-signal PP.H4 (coloc.susie re-check) vs region-wide coloc.abf PP.H4 (log scale)
    B. heatmap of PP.H4 across gene-context x test (abf / per-signal SuSiE)
    """
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 8, "axes.spines.top": False,
                         "axes.spines.right": False})
    CTX = _short_ctx

    # ---- A. PP.H4 comparison (log scale, rotated x-labels) ----
    labels, susie_h4, abf_h4 = [], [], []
    for _, r in out.iterrows():
        labels.append(f"{r['gene']}\n{CTX(r['context'])} {r['cs_rank']}")
        susie_h4.append(float(r["PPH4"])); abf_h4.append(float(r["abf_PPH4"]))
    susie_h4 = np.maximum(susie_h4, 1e-9); abf_h4 = np.maximum(abf_h4, 1e-9)
    x = np.arange(len(labels)); w = 0.36
    fig, ax = plt.subplots(figsize=(7.6, 3.9))
    ax.bar(x - w/2, abf_h4, w, label="coloc.abf (region, single causal)", color="#4C72B0")
    ax.bar(x + w/2, susie_h4, w, label="SuSiE signal-level re-check (M2b)", color="#DD8452")
    for xi, v in zip(x - w/2, abf_h4):
        ax.text(xi, v * 1.25, f"{v:.1e}", ha="center", va="bottom", fontsize=6.2)
    for xi, v in zip(x + w/2, susie_h4):
        ax.text(xi, v * 1.25, f"{v:.1e}", ha="center", va="bottom", fontsize=6.2)
    ax.axhline(0.05, ls="--", lw=0.8, color="grey")
    ax.text(len(labels) - 0.4, 0.05 * 1.6, "PP.H4 = 0.05", fontsize=6.5, color="grey")
    ax.set_yscale("log")
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=7, rotation=25, ha="right")
    ax.set_ylabel("PP.H4 (colocalization posterior)")
    ax.set_ylim(3e-7, 1.0)
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    ax.set_title("SuSiE re-check: per-signal PP.H4 vs coloc.abf region PP.H4", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "M2b_susie_coloc_PPH4.png"), dpi=300)
    fig.savefig(os.path.join(FIG_DIR, "M2b_susie_coloc_PPH4.pdf"))
    plt.close(fig)
    print("  saved figures/M2b_susie_coloc_PPH4.png/pdf")

    # ---- B. heatmap: PP.H4 by gene-context x test type (abf / SuSiE per-signal) ----
    import matplotlib.colors as mcolors
    combos = list(set(zip(out["gene"], out["context"])))  # unique combos
    combos = sorted(combos, key=lambda x: (x[0], x[1]))
    # col header: coloc.abf, then per-signal SuSiE
    test_cols = ["coloc.abf"] + [f"SuSiE {r}" for r in
                                  sorted({c.rsplit("_", 1)[-1] for c in out["cs_id"]})]
    M = np.full((len(combos), len(test_cols)), np.nan)
    annot = np.empty(M.shape, dtype=object)
    for i, (g_, ctx_) in enumerate(combos):
        for j, t in enumerate(test_cols):
            if t == "coloc.abf":
                row = out[(out["gene"] == g_) & (out["context"] == ctx_)].iloc[0]
                v = float(row["abf_PPH4"])
            else:
                rnk = t.split()[-1]
                sub = out[(out["gene"] == g_) & (out["context"] == ctx_) &
                          (out["cs_id"].str.endswith("_" + rnk))]
                if len(sub) == 0:
                    continue
                v = float(sub["PPH4"].iloc[0])
            M[i, j] = v; annot[i, j] = f"{v:.1e}" if v < 0.01 else f"{v:.3f}"
    # colormap: low = light, high = dark red, viridis-ish but on white background
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "wh2red", ["#FFFFFF", "#FEE5D9", "#FB6A4A", "#A50F15"])
    vmax = 0.05
    fig, ax = plt.subplots(figsize=(6.4, 0.55 * len(combos) + 1.0))
    im = ax.imshow(M, cmap=cmap, vmin=0, vmax=vmax, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, annot[i, j], ha="center", va="center", fontsize=7.5,
                    color="white" if M[i, j] and M[i, j] > 0.025 else "black")
    ax.set_xticks(range(len(test_cols))); ax.set_xticklabels(test_cols, fontsize=7.5)
    row_lbls = [f"{g_} · {CTX(ctx_)}" for g_, ctx_ in combos]
    ax.set_yticks(range(len(combos))); ax.set_yticklabels(row_lbls, fontsize=7.5)
    ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.04)
    cbar.set_label("PP.H4", fontsize=7.5)
    cbar.ax.tick_params(labelsize=6.5)
    ax.set_title("PP.H4 heatmap: coloc.abf baseline vs SuSiE per-signal re-check",
                 fontsize=9, pad=18)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "M2b_susie_coloc_PPH_heatmap.png"), dpi=300)
    fig.savefig(os.path.join(FIG_DIR, "M2b_susie_coloc_PPH_heatmap.pdf"))
    plt.close(fig)
    print("  saved figures/M2b_susie_coloc_PPH_heatmap.png/pdf")


if __name__ == "__main__":
    if "--plot-only" in sys.argv:
        # 快速重绘: 只读已生成的 TSV, 跳过 IgAN 索引重建
        out = pd.read_csv(os.path.join(OUT_DIR, "M2b_susie_coloc_recheck.tsv"), sep="\t")
        g = pd.read_csv(os.path.join(OUT_DIR, "M2b_susie_coloc_verdict.tsv"), sep="\t")
        import matplotlib
        matplotlib.use("Agg")
        _plot(out, g)
        print("  --plot-only done")
    else:
        main()
