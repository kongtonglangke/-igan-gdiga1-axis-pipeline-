#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""M4-1: GoDMC API 定向导出 → M4 中介分析输入层
数据源: GoDMC REST API (api.godmc.org.uk/v0.1), 坐标 GRCh37
步骤:
  1) 4 个机制轴 lead SNP → assoc_meta/rsid/ 全 mQTL 关联
  2) 显著 CpG (mre/are P<5e-8) → info/cpg/ 取位置注释
  3) 显著 CpG → assoc_meta/cpg/ 全量 mQTL (供 IV 选择 / coloc)
输出: 阶段3/M4_中介/
  lead_x_mQTL_sigCpG.tsv     lead 级显著 CpG 汇总
  sigCpG_pos.tsv             CpG 位置注释
  cpg_full_assoc_{cpg}.tsv   每显著 CpG 全量关联(SNP工具变量集)
"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import json, csv, os, time, urllib.request, urllib.error

API = "http://api.godmc.org.uk/v0.1"
OUT = _paths.at("阶段3/M4_中介")
os.makedirs(OUT, exist_ok=True)
LEADS = {"rs13226913": "C1GALT1", "rs10238682": "C1GALT1",
         "rs7856182": "GALNT12", "rs5910940": "C1GALT1C1"}
P_THRESH = 5e-8

def get(path, timeout=120):
    url = "%s/%s" % (API, path)
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 research"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.loads(r.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            print("  HTTP %d on %s, retry" % (e.code, url), flush=True)
        except Exception as e:
            print("  ERR %s (%s), retry" % (e, url), flush=True)
        time.sleep(3)
    return None

def pval(x):
    try: return float(x)
    except: return 1.0

# ---- 1) lead 级关联 ----
lead_rows, sig_cpgs = [], {}
for rs, gene in LEADS.items():
    d = get("assoc_meta/rsid/%s" % rs)
    rows = (d or {}).get("assoc_meta", [])
    print("== %s (%s): %d assoc rows" % (rs, gene, len(rows)), flush=True)
    for x in rows:
        pm, pa = pval(x.get("pval_mre")), pval(x.get("pval_are"))
        lead_rows.append({"lead_rsid": rs, "lead_gene": gene, "cpg": x.get("cpg"),
                          "snp": x.get("snp"), "rsid": x.get("rsid"),
                          "beta_a1": x.get("beta_a1"), "se_are": x.get("se_are"),
                          "pval_mre": pm, "pval_are": pa,
                          "direction": x.get("direction"), "num_studies": x.get("num_studies"),
                          "cistrans": x.get("cistrans")})
        if pm < P_THRESH or pa < P_THRESH:
            sig_cpgs.setdefault(x.get("cpg"), set()).add(rs)

with open(OUT("lead_x_mQTL_sigCpG_raw.tsv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(lead_rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(lead_rows)
print("saved lead_x_mQTL_sigCpG_raw.tsv (%d rows)" % len(lead_rows))

# ---- 2) 显著 CpG 位置注释 ----
pos_rows = []
for cpg in sig_cpgs:
    d = get("info/cpg/%s" % cpg)
    if isinstance(d, list) and d:
        x = d[0]
        pos_rows.append({"cpg": cpg, "chr": x.get("chr"), "pos37": x.get("pos"),
                         "lead_hits": ",".join(sorted(sig_cpgs[cpg])),
                         "info": str(x)[:300]})
        print("  %s @ chr%s:%s  (leads: %s)" % (cpg, x.get("chr"), x.get("pos"), ",".join(sorted(sig_cpgs[cpg]))), flush=True)
    else:
        print("  %s info not found" % cpg, flush=True)

with open(OUT("sigCpG_pos.tsv"), "w", newline="", encoding="utf-8") as f:
    if pos_rows:
        w = csv.DictWriter(f, fieldnames=list(pos_rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(pos_rows)
print("saved sigCpG_pos.tsv (%d CpGs)" % len(pos_rows))

# ---- 3) 显著 CpG 全量关联 (IV 集) ----
for cpg in sig_cpgs:
    d = get("assoc_meta/cpg/%s" % cpg)
    rows = (d or {}).get("assoc_meta", [])
    if not rows:
        print("  %s: no full assoc" % cpg, flush=True); continue
    fn = OUT("cpg_full_assoc_%s.tsv" % cpg)
    with open(fn, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    n_iv = sum(1 for x in rows if pval(x.get("pval_mre")) < P_THRESH or pval(x.get("pval_are")) < P_THRESH)
    print("  %s full assoc %d rows -> %s (p<5e-8 IVs: %d)" % (cpg, len(rows), os.path.basename(fn), n_iv), flush=True)
    time.sleep(1)

print("[DONE] M4-1 GoDMC API export finished")
