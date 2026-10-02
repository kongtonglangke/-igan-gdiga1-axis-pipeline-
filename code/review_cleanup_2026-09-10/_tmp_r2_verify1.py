# -*- coding: utf-8 -*-
"""Verify1: table cell-by-cell vs underlying + key numbers."""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os, json, csv

ROOT = _paths.LEGACY
SUB = ROOT(r"阶段4.5_复现包\submission\additional_file_1")
DATA = ROOT(r"阶段4.5_复现包\reproducibility\data")
RES = ROOT(r"阶段4.5_复现包\reproducibility\results")
MS = ROOT(r"阶段4.5_复现包\submission\manuscript")

def load_tsv(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))

fails = []
def chk(label, ms_val, und_val, tol=None):
    ok = False
    if tol is not None:
        try: ok = abs(float(ms_val) - float(und_val)) <= tol
        except: ok = str(ms_val) == str(und_val)
    else:
        ok = str(ms_val).strip() == str(und_val).strip()
    print(("PASS" if ok else "FAIL"), label, "| ms:", ms_val, "| und:", und_val)
    if not ok:
        fails.append((label, ms_val, und_val))

print("### Table 4 vs S2a (14 rows: n_SNP, PPH1, PPH4) ###")
s2a = {(r["gene"], r["context"]): r for r in load_tsv(os.path.join(SUB, "Table_S2a_coloc_25comparisons.tsv"))}
t4 = [
 ("C1GALT1","GTEx_v10_blood",5205,0.865,0.008),("C1GALT1","GTEx_v10_kidney_cortex",4186,0.909,0.011),
 ("C1GALT1","OneK1K_B_naive",4049,0.792,0.011),("C1GALT1","OneK1K_B_memory",4098,0.256,0.006),
 ("C1GALT1","OneK1K_B_intermediate",4072,0.177,0.005),("GALNT2","GTEx_v10_blood",4465,0.213,0.007),
 ("GALNT2","GTEx_v10_kidney_cortex",3499,0.222,0.005),("GALNT2","OneK1K_B_naive",4340,0.891,0.006),
 ("GALNT2","OneK1K_B_memory",4324,0.233,0.009),("GALNT2","OneK1K_B_intermediate",4320,0.276,0.007),
 ("GALNT12","GTEx_v10_blood",3616,0.196,0.005),("GALNT12","GTEx_v10_kidney_cortex",2776,0.239,0.004),
 ("ST6GALNAC2","GTEx_v10_blood",4129,0.850,0.031),("ST6GALNAC2","GTEx_v10_kidney_cortex",3273,0.289,0.007),
]
for g,c,n,p1,p4 in t4:
    r = s2a.get((g,c))
    if not r:
        print("FAIL missing", g, c); fails.append((f"T4 {g}/{c}","present","MISSING in S2a")); continue
    chk(f"T4 {g} {c} n_SNP", n, r["n_SNP"])
    chk(f"T4 {g} {c} PPH1", p1, round(float(r["PPH1"]),3), 0.0006)
    chk(f"T4 {g} {c} PPH4", p4, round(float(r["PPH4"]),3), 0.0006)
# S2a computable count & max PPH4
comp = [r for r in s2a.values() if r.get("PPH4")]
print("S2a total rows:", len(s2a), "computable(PPH4 non-empty):", len(comp),
      "max PPH4:", max(float(r["PPH4"]) for r in comp))

print("\n### Table 5 vs S5a/S5b/S5c/S5d ###")
s5a = {r["cpg"]: r for r in load_tsv(os.path.join(SUB, "Table_S5a_sigCpG_5CpG_pos.tsv"))}
s5b = load_tsv(os.path.join(SUB, "Table_S5b_MR_CpG_x_IgAN.tsv"))
s5c = {r["cpg"]: r for r in load_tsv(os.path.join(SUB, "Table_S5c_coloc_mQTL_x_IgAN.tsv"))}
s5d = load_tsv(os.path.join(SUB, "Table_S5d_coloc_mQTL_x_eQTL.tsv"))
t5 = [  # cpg, pos, kb, mrqtlP, beta, se, mr_beta, mr_p, pph4_igan, pph4_naive
 ("cg19603390","7222222","+25.7","3.1e-10",0.20,0.025,0.064,0.282,0.006,0.006),
 ("cg19473623","7224869","+28.3","2.2e-21",0.14,0.012,0.093,0.282,0.010,0.020),
 ("cg17994788","7261594","+65.0","2.4e-14",0.14,0.013,0.212,0.152,0.014,0.792),
 ("cg04827551","7268805","+72.2","1.3e-29",-0.34,0.023,-0.057,0.064,0.007,0.006),
 ("cg16101574","7291514","+94.9","<1e-300",0.42,0.010,0.023,0.279,0.007,0.823),
]
TSS = 7196565
for cpg,pos,kb,mp,b,se,mrb,mrp,p4i,p4n in t5:
    ra = s5a.get(cpg)
    chk(f"T5 {cpg} pos37", pos, ra["pos37"] if ra else "MISSING")
    if ra:
        kbu = (int(ra["pos37"]) - TSS)/1000
        chk(f"T5 {cpg} kbTSS", kb, f"+{kbu:.1f}")
    rc = s5c.get(cpg)
    chk(f"T5 {cpg} PPH4xIgAN", p4i, round(float(rc["PPH4"]),3), 0.0006)
    print(f"   S5c {cpg} PPH1={float(rc['PPH1']):.4f}")
    dn = [r for r in s5d if r["cpg"]==cpg and r["context"]=="OneK1K_B_naive"]
    chk(f"T5 {cpg} PPH4xNaive", p4n, round(float(dn[0]["PPH4"]),3), 0.0006)
    # MR: smallest p among Wald rs13226913 / Wald topcis / clumped IVW
    rows = [r for r in s5b if r["cpg"]==cpg]
    for r in rows:
        print(f"   S5b {cpg} {r['method']} {r['instrument']}: b={r['mr_beta']} p={r['mr_p']}")
    cand = [r for r in rows if r["mr_p"] and (r["instrument"]=="rs13226913" or "topcis" in r["method"] or "IVW" in r["method"])]
    best = min(cand, key=lambda r: float(r["mr_p"])) if cand else None
    if best:
        chk(f"T5 {cpg} MR beta(minP set)", mrb, round(float(best["mr_beta"]),3), 0.0006)
        chk(f"T5 {cpg} MR P(minP set)", mrp, round(float(best["mr_p"]),3), 0.0006)
    else:
        print("   (no cand rows)")
# mQTL beta/se/P from raw
raw = load_tsv(os.path.join(DATA, r"M4_mediation\lead_x_mQTL_sigCpG_raw.tsv"))
for cpg,pos,kb,mp,b,se,mrb,mrp,p4i,p4n in t5:
    r1 = [r for r in raw if r["cpg"]==cpg and r["lead_rsid"]=="rs13226913"]
    r2 = [r for r in raw if r["cpg"]==cpg and r["lead_rsid"]=="rs10238682"]
    if r1:
        chk(f"T5 {cpg} mQTL beta", b, round(float(r1[0]["beta_a1"]),2), 0.006)
        chk(f"T5 {cpg} mQTL se", se, round(float(r1[0]["se_are"]),3), 0.0006)
    cands = [float(x["pval_mre"]) for x in (r1+r2)]
    if mp.startswith("<"):
        print("PASS?" , f"T5 {cpg} mQTL P {mp} vs min_mre {min(cands):.3e}")
    else:
        chk(f"T5 {cpg} mQTL P(min mre)", float(mp), min(cands), float(mp)*0.06)

print("\n### S11 heterogeneity ###")
for r in load_tsv(os.path.join(SUB, "Table_S11_Bcell_eQTL_heterogeneity.tsv")):
    print(r)

print("\n### M2b summary.json (allelic score) ###")
for fn in ["summary.json","summary_m2b.json"]:
    p = os.path.join(DATA, r"M2b_susie_coloc", fn)
    if os.path.exists(p):
        with open(p, encoding="utf-8") as f: print(fn, f.read()[:1200])

print("\n### regulon empirical P (regulon_axis_composite.tsv) ###")
rac = load_tsv(os.path.join(RES, r"sc_regulatory\regulon_axis_composite.tsv"))
print("cols:", list(rac[0].keys()))
for tf,bs in [("DNMT1","B_plasma"),("KLF2","B_naive"),("SP1","B_naive")]:
    rows = [r for r in rac if r.get("tf")==tf and r.get("b_subset")==bs]
    for r in rows: print(tf, bs, r)
# pass count
passcol = [c for c in rac[0] if "emp" in c.lower() or "pass" in c.lower()]
print("emp/pass cols:", passcol)
npass = sum(1 for r in rac if passcol and r.get(passcol[0]) in ("True","1","TRUE"))
print("rows:", len(rac), "n_pass:", npass)

print("\n### S4b locus lookup min P ###")
s4b = load_tsv(os.path.join(SUB, "Table_S4b_cisQTL_topvar_IgAN_locus_lookup.tsv"))
for r in s4b: print({k:v for k,v in r.items() if k in ("gene","context","min_p","igan_min_p","region_min_p") or "p" in k.lower()})

print("\n### S12/S13 checks ###")
s12 = load_tsv(os.path.join(SUB, "Table_S12_EAS_eQTL_crossref.tsv"))
c1 = [r for r in s12 if r["gene"]=="C1GALT1"]
eg = sum(1 for r in c1 if r["is_eGene_FDR0.05"]=="yes")
print("S12 C1GALT1 subsets:", len(c1), "eGene yes:", eg)
c1c1 = [r for r in s12 if r["gene"]=="C1GALT1C1"]
print("S12 C1GALT1C1 subsets:", len(c1c1), "eGene yes:", sum(1 for r in c1c1 if r["is_eGene_FDR0.05"]=="yes"))
s13 = load_tsv(os.path.join(SUB, "Table_S13_lead_C1GALT1_eQTL_EAS.tsv"))
for rs in ["rs10238682","rs13226913"]:
    rows = [r for r in s13 if r["lead_rsid"]==rs and r["is_independent_signal"]=="yes"]
    bl = [r for r in rows if r["B_lineage"]=="yes"]
    print(rs, "indep subsets:", len(rows), "B-lineage:", len(bl), [r["cell_subset"] for r in bl])

print("\n### 42,259 location in manuscript ###")
import glob
for fp in glob.glob(os.path.join(MS, "*.md")):
    with open(fp, encoding="utf-8", errors="replace") as f:
        for i, ln in enumerate(f, 1):
            if "42,259" in ln or "23,458" in ln or "16,398" in ln:
                print(os.path.basename(fp), i, ln[:300])

print("\n### FAILS:", len(fails))
for f_ in fails: print("FAIL:", f_)
