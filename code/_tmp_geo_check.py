# -*- coding: utf-8 -*-
"""GEO original-paper backfill evidence gathering (step 4.4)."""
import json, time, urllib.request, urllib.parse

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

def get(url, params, retries=3):
    q = urllib.parse.urlencode(params)
    full = url + "?" + q
    for i in range(retries):
        try:
            req = urllib.request.Request(full, headers={"User-Agent": "research-script/0.1 (contact: user)"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            print("  retry", i, e)
            time.sleep(3)
    return None

# 1) esummary for PMID 28646076 (GSE93798 linked publication)
print("=== esummary PMID 28646076 ===")
out = get(BASE + "esummary.fcgi", {"db": "pubmed", "id": "28646076", "retmode": "json"})
if out:
    d = json.loads(out)
    for uid, rec in d.get("result", {}).items():
        if uid == "uids":
            continue
        print(json.dumps({
            "title": rec.get("title"),
            "authors": [a.get("name") for a in rec.get("authors", [])][:12],
            "source": rec.get("source"),
            "pubdate": rec.get("pubdate"),
            "volume": rec.get("volume"),
            "issue": rec.get("issue"),
            "pages": rec.get("pages"),
            "doi": [x.get("value") for x in rec.get("articleids", []) if x.get("idtype") == "doi"],
        }, ensure_ascii=False, indent=1))

# 2) gds efetch XML for the two GDS records — look for PubMed/Reference fields
print("\n=== gds efetch 200154046 (GSE73953) / 200115857 (GSE115857) ===")
for gds in ["200154046", "200115857"]:
    out = get(BASE + "efetch.fcgi", {"db": "gds", "id": gds, "retmode": "xml"})
    if out:
        import re
        # grab relevant XML fragments
        for tag in ["PubMedIds", "PubMedID", "Reference", "Citation", "PDAT", "title", "summary", "PDAT"]:
            for m in re.finditer(r"<" + tag + r"[^>]*>(.*?)</" + tag + r">", out, re.S):
                t = re.sub(r"<[^>]+>", " ", m.group(1)).strip()
                if t:
                    print(f"  [{gds}] {tag}: {t[:400]}")
        if "PubMed" not in out:
            print(f"  [{gds}] NO PubMed/Reference element in gds XML")

# 3) pubmed esearch for papers that mention the accessions
print("\n=== pubmed esearch GSE73953 / GSE115857 (all fields) ===")
for acc in ["GSE73953", "GSE115857"]:
    out = get(BASE + "esearch.fcgi", {"db": "pubmed", "term": f'"{acc}"[All Fields]', "retmode": "json", "retmax": "10"})
    if out:
        d = json.loads(out)
        ids = d.get("esearchresult", {}).get("idlist", [])
        print(f"  {acc}: {len(ids)} pubmed hits -> {ids}")
