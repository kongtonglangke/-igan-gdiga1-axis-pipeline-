# -*- coding: utf-8 -*-
"""预检: 修复列序后在多语境加载轴基因 eQTL"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os, pandas as pd
EQTL_CACHE = _paths.at("00_rawdata/eQTL_Catalogue/region_cache")

def load(qtd, chr_, ensg, tss, flank=1_000_000):
    path = EQTL_CACHE(f"{qtd}.{chr_}.tsv")
    if not os.path.exists(path):
        return None
    rows = []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            fld = line.rstrip("\n").split("\t")
            if len(fld) != 19 or fld[0] != ensg:
                continue
            rsid = fld[18]
            if not rsid or rsid == "NA":
                continue
            try:
                pos = int(fld[2]); beta = float(fld[9]); se = float(fld[10]); p = float(fld[8])
            except Exception:
                continue
            if abs(pos - tss) > flank:
                continue
            rows.append((rsid, pos, beta, se, p, fld[3], fld[4]))
    return pd.DataFrame(rows, columns=["rsid","pos","beta","se","p","ref","alt"]) if rows else None

for sym, eg, ch, tss in [("GALNT12","ENSG00000119514","9",98_808_381),
                          ("C1GALT1C1","ENSG00000171155","X",120_624_955)]:
    for qtd, ctx in [("QTD000356","GTEx_blood"),("QTD000261","GTEx_kidney")]:
        df = load(qtd, ch, eg, tss)
        if df is None or len(df) == 0:
            print(f"{sym} {ctx}: NO DATA")
            continue
        print(f"{sym} {ctx}: n={len(df)} p范围[{df['p'].min():.2e},{df['p'].max():.2e}] "
              f"rsid={df['rsid'].iloc[0]} alt={df['alt'].unique()[:6]}")
