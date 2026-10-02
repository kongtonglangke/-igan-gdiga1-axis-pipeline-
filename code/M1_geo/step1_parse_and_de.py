# -*- coding: utf-8 -*-
"""
M1 · 步骤1.1 GEO 表达锚点分析 (阶段1/M1_GEO表达锚点)
=====================================================
对应《...报告》步骤1.1 (M1)：
  GSE73953 (PBMC, IgAN15 vs MN8 vs HC-pooled2)      [主]
  GSE115857(肾活检bulk, IgAN55 vs 对照7+Living donor等) [辅]
  GSE93798 (肾小球, IgAN20 vs Control22)            [辅]
5 基因: C1GALT1 / C1GALT1C1 / ST6GALNAC2 / GALNT2 / GALNT12
分析: IgAN vs 对照 差异表达 (limma 语义)
执行环境: Python (default venv, numpy/scipy/pandas/matplotlib)
"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import gzip, os, re, json
import numpy as np
import pandas as pd

RAW = _paths.at("00_rawdata/GEO_M1")
OUT = _paths.at("阶段1/M1_GEO表达锚点")
os.makedirs(OUT, exist_ok=True)

AXIS = {  # symbol -> (entrez, ensg)
    "C1GALT1":   (56913, "ENSG00000106392"),
    "C1GALT1C1": (29071, "ENSG00000171155"),
    "ST6GALNAC2":(10610, "ENSG00000070731"),
    "GALNT2":    (2590,  "ENSG00000143641"),
    "GALNT12":   (79695, "ENSG00000119514"),
}

# ---------------- series matrix 解析 ----------------
def parse_series_matrix(path):
    """返回 (meta_dict, expr_df)。expr_df: rows=probe, cols=GSM."""
    txt = gzip.open(path, "rt", encoding="utf-8", errors="replace").read()
    meta = {}
    for line in txt.split("\n"):
        if line.startswith("!") and "\t" in line:
            key, rest = line.split("\t", 1)
            vals = [v.strip('"') for v in rest.split("\t") if v.strip('"')]
            meta.setdefault(key, vals)
    body = txt.split("!series_matrix_table_begin")[1].split("!series_matrix_table_end")[0]
    lines = [l for l in body.split("\n") if l.strip()]
    # SOFT 表格段第一行即表头 (ID_REF + GSM...), 无 '!' 前缀
    cols = [c.strip('"') for c in lines[0].split("\t")]
    data = []
    for l in lines[1:]:
        f = [c.strip('"') for c in l.split("\t")]
        if len(f) == len(cols):
            data.append(f)
    df = pd.DataFrame(data, columns=cols).set_index(cols[0])
    df = df.apply(pd.to_numeric, errors="coerce")
    return meta, df

# ---------------- 分组构造 ----------------
def gse73953_groups(meta):
    titles = meta.get("!Sample_title", [])
    accs   = meta.get("!Sample_geo_accession", [])
    groups = {}
    for a, t in zip(accs, titles):
        tl = t.lower()
        if "iga nephrop" in tl or "igan" in tl:
            groups[a] = "IgAN"
        elif "membranous" in tl:
            groups[a] = "MN"
        elif "healthy" in tl or "control" in tl:
            groups[a] = "HC"
        else:
            groups[a] = "OTHER"
    return groups

def gse115857_groups(meta):
    accs = meta.get("!Sample_geo_accession", [])
    srcs = meta.get("!Sample_source_name_ch1", [])
    groups = {}
    for a, s in zip(accs, srcs):
        sl = s.lower()
        if "igan" in sl:
            groups[a] = "IgAN"
        elif "living donor" in sl:
            groups[a] = "Ctrl_LD"
        elif "membranous" in sl:
            groups[a] = "MN"
        elif "minimal change" in sl:
            groups[a] = "MCD"
        elif "focal" in sl:
            groups[a] = "FSGS"
        else:
            groups[a] = "OTHER"
    return groups

def gse93798_groups(meta):
    accs = meta.get("!Sample_geo_accession", [])
    srcs = meta.get("!Sample_source_name_ch1", [])
    groups = {}
    for a, s in zip(accs, srcs):
        sl = s.lower()
        if "control" in sl or "con" == sl.strip():
            groups[a] = "Ctrl"
        elif "igan" in sl:
            groups[a] = "IgAN"
        else:
            groups[a] = "OTHER"
    return groups

# ---------------- 探针 → 基因映射 ----------------
def map_gpl4133(annot_path):
    """GPL4133 annot: ID(numeric) Gene symbol ... -> {probe: symbol}"""
    lines = gzip.open(annot_path, "rt", encoding="utf-8", errors="replace").read().split("\n")
    start = None
    for i, l in enumerate(lines):
        if l.startswith("!platform_table_begin"):
            start = i; break
    hdr = lines[start+1].split("\t")
    id_i, sym_i = hdr.index("ID"), hdr.index("Gene symbol")
    m = {}
    for l in lines[start+2:]:
        if not l.strip() or l.startswith("!"):
            continue
        f = l.split("\t")
        if len(f) <= max(id_i, sym_i):
            continue
        m[f[id_i]] = f[sym_i]
    return m

def build_probe_gene(gse, expr_df, meta, annot_map=None):
    """返回 Series: probe -> (gene_symbol) 只保留能映射到机制基因的探针"""
    probes = list(expr_df.index)
    if gse == "GSE93798":
        # BrainArray custom CDF: probe = {entrez}_at
        out = {}
        for p in probes:
            mm = re.match(r"^(\d+)_at$", p)
            if mm:
                e = int(mm.group(1))
                for sym, (eg, _) in AXIS.items():
                    if e == eg:
                        out[p] = sym
        return out
    # GSE115857: ILMN 探针 -> Entrez -> symbol (illuminaHumanv4.db, GPL14951)
    if gse == "GSE115857":
        out = {}
        mpath = RAW("illuminaHumanv4_probe2gene.tsv")
        if os.path.exists(mpath):
            tab = pd.read_csv(mpath, sep="\t", dtype={"entrez": str})
            eg2sym = {str(eg): sym for sym, (eg, _) in AXIS.items()}
            for _, row in tab.iterrows():
                eg = str(row["entrez"]).split(".")[0]
                if eg in eg2sym:
                    out[row["probe_id"]] = eg2sym[eg]
        return out
    # GSE73953: numeric feature ID -> annot (GPL4133)
    out = {}
    if annot_map is not None:
        for p in probes:
            s = annot_map.get(p, "")
            if s in AXIS:
                out[p] = s
    return out

# ---------------- limma 语义差异检验 (两组, Welch/经验贝叶斯近似) ----------------
from scipy import stats as sstats

def de_two_group(x_grp, y_grp):
    """x_grp: 病例 (IgAN) 值向量; y_grp: 对照。返回 (log2FC=mean_x-mean_y, t, p)。"""
    x, y = np.asarray(x_grp, float), np.asarray(y_grp, float)
    if len(x) < 2 or len(y) < 2:
        return np.nan, np.nan, np.nan
    # Welch t-test (两组等方差不设限, 小样本稳健)
    t, p = sstats.ttest_ind(x, y, equal_var=False)
    fc = x.mean() - y.mean()
    return fc, t, p

# ---------------- 全流程 ----------------
def run_one(gse, path, group_fn, annot_path=None, compare=("IgAN", "Ctrl"), log2=False):
    meta, df = parse_series_matrix(path)
    if log2:
        df = np.log2(df.clip(lower=1e-3))  # 原始信号 -> log2 表达
    groups = group_fn(meta)
    annot_map = map_gpl4133(annot_path) if annot_path else None
    probe_gene = build_probe_gene(gse, df, meta, annot_map)
    if not probe_gene:
        print(gse, "WARN no probe mapped"); return None
    # 只保留能映射的探针
    keep = [p for p in df.index if p in probe_gene]
    sub = df.loc[keep]
    rows = []
    for p in keep:
        sym = probe_gene[p]
        ca = [c for c in sub.columns if groups.get(c) == compare[0]]
        co = [c for c in sub.columns if groups.get(c) == compare[1]]
        vals_ca = sub.loc[p, ca].astype(float).values
        vals_co = sub.loc[p, co].astype(float).values
        fc, t, pv = de_two_group(vals_ca, vals_co)
        rows.append(dict(gene=sym, probe=p, n_case=len(ca), n_ctrl=len(co),
                         mean_case=np.nanmean(vals_ca), mean_ctrl=np.nanmean(vals_co),
                         logFC=fc, t=t, pvalue=pv))
    out = pd.DataFrame(rows)
    # 多探针取最显著
    best = out.loc[out.groupby("gene")["pvalue"].idxmin()].reset_index(drop=True)
    return {"meta": meta, "groups": groups, "all_probe": out, "best_gene": best,
            "case": compare[0], "ctrl": compare[1], "expr": df}

def save_outputs(r, gse):
    """统一输出: best/allprobe/groups/5genes_expr(log2 尺度, 多探针均值)."""
    case, ctrl = r["case"], r["ctrl"]
    tag = f"{case}_vs_{ctrl}"
    r["best_gene"].to_csv(OUT(f"{gse}_{tag}_best.csv"), index=False)
    r["all_probe"].to_csv(OUT(f"{gse}_{tag}_allprobe.csv"), index=False)
    pd.Series(r["groups"]).to_csv(OUT(f"{gse}_groups.tsv"), sep="\t")
    keep = list(r["all_probe"]["probe"])
    sub = r["expr"].loc[keep]
    gmap = r["all_probe"].set_index("probe")["gene"]
    sub = sub.groupby(gmap).mean()
    sub.to_csv(OUT(f"{gse}_5genes_expr.csv"))
    print(f"  saved: {gse}_{tag}_best.csv / _allprobe.csv / {gse}_groups.tsv / {gse}_5genes_expr.csv")

def main():
    print("== GSE73953 (PBMC, IgAN vs HC-pooled) ==")
    r1 = run_one("GSE73953", RAW("GSE73953_series_matrix.txt.gz"),
                 gse73953_groups, annot_path=RAW("GPL4133.annot.gz"),
                 compare=("IgAN", "HC"), log2=True)
    if r1:
        print("  probe mapped:", len(r1["all_probe"]), "| genes:", sorted(r1["best_gene"]["gene"]))
        print(r1["best_gene"].to_string(index=False))
        save_outputs(r1, "GSE73953")

    print("\n== GSE115857 (肾活检 bulk, IgAN vs Living donor) ==")
    r2 = run_one("GSE115857", RAW("GSE115857_series_matrix.txt.gz"),
                 gse115857_groups, annot_path=None, compare=("IgAN", "Ctrl_LD"))
    if r2:
        print("  probe mapped:", len(r2["all_probe"]), "| genes:", sorted(r2["best_gene"]["gene"]))
        print(r2["best_gene"].to_string(index=False))
        save_outputs(r2, "GSE115857")
    else:
        print("  WARN 无探针映射 — 检查 illuminaHumanv4_probe2gene.tsv 与 ILMN 探针交集")

    print("\n== GSE93798 (肾小球, IgAN vs Control) ==")
    r3 = run_one("GSE93798", RAW("GSE93798_series_matrix.txt.gz"),
                 gse93798_groups, annot_path=None, compare=("IgAN", "Ctrl"))
    if r3:
        print("  probe mapped:", len(r3["all_probe"]), "| genes:", sorted(r3["best_gene"]["gene"]))
        print(r3["best_gene"].to_string(index=False))
        save_outputs(r3, "GSE93798")

if __name__ == "__main__":
    main()
