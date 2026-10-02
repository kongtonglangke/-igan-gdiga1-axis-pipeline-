#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""M5-1: 从 GCST90018866 (IgAN 跨祖先 meta, GRCh38) 提取 4 个机制轴 lead 的效应+EAF。
输出: 阶段3/M5_东亚/lead_x_IgAN_eaf.tsv（补全阶段2表缺失的 effect_allele_frequency）
"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import gzip, csv, sys, time

SRC = _paths.legacy("00_rawdata/GCST90018866/GCST90018866.h.tsv.gz")
OUT = _paths.legacy("阶段3/M5_东亚/lead_x_IgAN_eaf.tsv")
LEADS = {"rs13226913", "rs10238682", "rs7856182", "rs5910940"}

t0 = time.time()
hits = []
with gzip.open(SRC, "rt", encoding="utf-8", errors="replace") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for i, row in enumerate(reader):
        rs = row.get("rsid", "")
        if rs in LEADS:
            hits.append({
                "lead_rsid": rs,
                "chr": row.get("chromosome"),
                "pos38": row.get("base_pair_location"),
                "EA": row.get("effect_allele"),
                "OA": row.get("other_allele"),
                "beta": row.get("beta"),
                "se": row.get("standard_error"),
                "EAF": row.get("effect_allele_frequency"),
                "p": row.get("p_value"),
            })
        if i % 5000000 == 0 and i:
            print("  scanned %dM rows..." % (i // 1000000), flush=True)

print("scanned %d rows in %.1fs; %d lead hits" % (i + 1, time.time() - t0, len(hits)))
with open(OUT, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(hits[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(hits)
print("saved:", OUT)
for h in hits:
    print(h)
