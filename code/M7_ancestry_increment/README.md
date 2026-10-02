# M7 — ancestry-matched increment (added 2026-09-26, manuscript v1.9)

This module adds the four **ancestry-matched** checks introduced in x3-docx **v1.9**
(see `稿件x3_修订标黄说明_v1.9_20260926.md`). No main figure is affected; all
outputs feed **Additional file 1: Tables S14–S15** and the East Asian / glycan-
specificity statements of Methods and Results.

## A — Chinese population allele frequencies (Table S14, Table 6 column 5)
- Source: Westlake Biobank for Chinese (WBBC), `https://wbbc.westlake.edu.cn/data/WBBC.chr{7,9}.GRCh38.vcf.gz`
  (V20210103; sites + AF; 4,480 Chinese individuals; **autosomes only — no chrX**).
- Script: `step_m7_1_wbbc_chinese_af.py` → Chinese AF for the four Gd-IgA1 leads
  (rs13226913-C 0.9345; rs10238682-G 0.5289; rs7856182-T 0.02623; rs5910940 not covered).
- **Gotcha**: WBBC `CHROM` is written `chr7`/`chr9`; strip the `chr` prefix before matching.

## C — Asian-only IgAN meta-analysis, ancestry-matched look-up (Table S14)
- Source: Kiryluk Lab, `https://www.columbiamedicine.org/divisions/kiryluk/gwas/IgA/IgAN_Combined_metaanalysis_Asian_only.txt`
  (376,485,006 B; 4,590 biopsy-diagnosed cases / 7,573 controls; **hg19**; SNP key `CHR:BP_hg19`;
  chromosomes 1–22, **no chrX**).
- Script: `step_m7_2_axis_lookup.py` (full-file scan; the file is not sorted by position, so a
  complete scan is required). Result: rs13226913 P = 0.0766, rs10238682 P = 0.0397,
  rs7856182 P = 0.394 — **all far above 5 × 10⁻⁸**, consistent with the disease-negative result;
  all three betas are positive with the Gd-IgA1-increasing allele as the effect allele.
  rs5910940 untestable (no chrX).

## B — kidney-region eQTL context — **NOT DONE** (resource unavailable)
Three independent checks on 2026-09-26: (i) eQTL Catalogue release 8 lists only
`GTEx_v10 kidney_cortex` for kidney (already used in the manuscript); (ii) `nephqtl.org` query API
returned HTTP 503 throughout; (iii) `nephqtl2.org` redirects to a `hugeampkpn.org` single-page app
with no data/API endpoint. No result was written into the manuscript.

## D — serum total IgA specificity control (Table S15)
- Source: Kiryluk Lab, `https://www.columbiamedicine.org/divisions/kiryluk/gwas/IgA/META.IGA.LEVELS.ALL.COMBINED.txt`
  (596,831,842 B; Liu et al. 2022, *Nat Commun* 13:6859 [29]; 41,263 individuals; **hg19**; no chrX).
- **Critical gotcha**: the file is sorted **lexicographically by the SNP string `CHR:BP`**, not by
  numeric position. Chromosome order is therefore `10,11,…,19,1,20,21,22,2,3,…,9`, and within a
  chromosome the order is by the *string* (`7:100000274` < `7:6696565`). Locate regions by
  **binary search on the SNP string**, then filter records by numeric BP.
- Scripts: `step_m7_3a_local_block_extract.py` (chr1/chr17 taken from locally downloaded 1 MB blocks),
  `step_m7_3b_serum_iga_range_fetch.py` (chr7/chr9 fetched by HTTP Range + SNP-string binary search),
  `step_m7_4_serum_iga_control.py` (analysis + Table S15).
- Result: rs13226913 P = 0.917, rs10238682 P = 0.216, rs7856182 P = 0.170; regional minima
  5.2 × 10⁻⁵ (C1GALT1), 1.25 × 10⁻⁴ (GALNT12), 1.25 × 10⁻⁴ (GALNT2), 1.96 × 10⁻⁴ (ST6GALNAC2);
  **zero** genome-wide significant variants in any of the four testable regions → the axis signal is
  specific to the glycan trait rather than to total antibody quantity.

## Utility scripts
- `util_dl_kiryluk_blocks.py <url> <out> [threads] [blocksize]` — 1 MB-block, resumable,
  parallel downloader for the Kiryluk file server (the server throttles roughly 10–40 KB/s per
  connection, so ≥ 16 threads are needed).
- `util_remote_range_locate.py` — builds a coarse chromosome index over the remote file
  (`{chr, bp}` at a fixed byte interval) so that regions can be located without downloading the
  whole file. Uses `urllib` with `ssl.CERT_NONE`; plain `requests` fails against this host.

## Delivered tables (single source of truth)
`step_m7_5_build_tables.py` regenerates **`Table_S14_ancestry_matched_lookups.tsv`**
byte-for-byte from the four inputs (Asian-only summary statistics, WBBC chr7/chr9 hits,
the archived gnomAD/Han-Gd-IgA1/IgAN-meta columns under `data/M5_eastasia/`, and the
`LEADS` coordinate table inside the script):

```
python step_m7_5_build_tables.py                       # -> out/Table_S14_ancestry_matched_lookups.tsv
```

Formatting rules (fixed, so the table is reproducible rather than hand-typed):
gnomAD EAS/NFE `%.3f`; WBBC Chinese and its sub-populations `%.5f` below 0.1 else `%.4f`;
β/SE verbatim from the summary-statistics file; Asian-only P `%.3g`; Han-Gd-IgA1 P `%.2g`
(exponent normalised, `e-09` → `e-9`); European IgAN P `%.3f` (NA when untestable).
`Table_S15_serumIgA_axis_regions.tsv` is written by `step_m7_4_serum_iga_control.py`.

## Data mirror
`reproducibility/data/M7_ancestry_increment/` carries the delivered tables and the
lead-level look-up: `Table_S14_ancestry_matched_lookups.tsv`,
`Table_S15_serumIgA_axis_regions.tsv`, `D_lookup.tsv`.

## Reproduction order
```
python step_m7_1_wbbc_chinese_af.py                                    # A
python util_dl_kiryluk_blocks.py <IgAN_Asian_only_url>  IgAN_Asian_only.txt 24
python step_m7_2_axis_lookup.py IgAN_Asian_only.txt <serum_iga_axis.tsv> out/  # C
python util_dl_kiryluk_blocks.py <META.IGA.LEVELS_url> META_IGA_LEVELS.txt 24
python step_m7_3a_local_block_extract.py                               # chr1/chr17
python step_m7_3b_serum_iga_range_fetch.py                             # chr7/chr9
python step_m7_4_serum_iga_control.py                                  # D -> Table S15
python step_m7_5_build_tables.py                                       # -> Table S14
```
