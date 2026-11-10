# M4 — rebuilt MR layer (added 2026-11-10)

The methylation→IgAN Mendelian randomisation layer was **rebuilt from scratch** for
the *European Journal of Medical Research* revision (referee comment: "the MR is not
good enough"). This directory supersedes the MR/coloc half of
`code/M4_mediation/step_m4_2_中介分析.py`; the CpG discovery half
(`step_m4_1_fetch_mqtl.py`, `Table S5a`) is unchanged.

It produces **Figure 4 (panels A–F)**, **Figure S-Appendix A1**, and
**Tables S5b / S5c / S5d**. Results, Methods text and the figure legends that the
package ships are byte-identical to the outputs below.

## What changed versus the pre-revision layer

| Item | Before | Now |
|---|---|---|
| Instruments | the two published Gd-IgA1 leads, or a per-CpG top sentinel | selected **from the mQTL data itself**: P < 5e-8, GoDMC-style clumping (r² < 0.001 / 10 Mb), F ≥ 10 |
| Instrument count | 1–2 per CpG, F not reported | 9 independent signals → 7 harmonisable → **4 CpGs estimable**; F 62.8–4050 |
| Diagnostics | none | Cochran's Q, I², MR-Egger intercept, leave-one-out, Steiger, sample-overlap audit, positive control (eQTL), type-I-error calibration, power curve |
| Reverse MR | not attempted | attempted and reported as **not estimable** (no variant in the locus reaches P < 5e-8 for IgAN) |
| `Table S5d` colocalisation | computed with GoDMC `se_mre` | computed with GoDMC **`se`** (see §"One substantive correction") |
| Wald ratio rows | emitted for every CpG | emitted **only when k = 1**; for k ≥ 2 the Wald expression is numerically identical to fixed-effect IVW, so a second row under that label was both redundant and misnamed |

## One substantive correction (Table S5d)

GoDMC releases two standard errors. `se` comes from the 36-cohort **fixed-effect**
meta-analysis on which the significance calls and the reported P values rest;
`se_mre` comes from the **multiplicative random-effects** sensitivity analysis
(Min et al., *Nat Genet* 2021). The pre-revision `Table S5d` was built with
`se_mre`; re-running the identical 25 comparisons with the two columns side by side
reproduces the old table from `se_mre` **exactly** (`max |ΔPPH4| = 0` over 25 rows),
so the disagreement between the two versions was purely a column choice.

`out/coloc_se_diagnosis.tsv` holds both columns for all 25 comparisons. The ratio
`se_mre/se` is 1.44–3.07 (median 1.86): using `se_mre` inflates every mQTL standard
error roughly twofold, which systematically pushes PPH4 towards 0. The shipped
`Table S5d` therefore uses `se`. Seven comparisons change verdict band, and the sets
of "shared causal variant" calls under the two columns are **disjoint**:

| Comparison | `se_mre` (retired) | `se` (shipped) |
|---|---|---|
| cg19473623 × kidney cortex | 0.000 | **0.992** |
| cg19473623 × naive B | 0.020 | **0.825** |
| cg17994788 × naive B | 0.792 | **0.870** |
| cg17994788 × kidney cortex | 0.966 | 0.605 |
| cg16101574 × naive B | 0.823 | 0.780 |
| cg17994788 × whole blood | 0.617 | < 1e-6 |
| cg17994788 × memory B | 0.553 | 0.144 |

The `se` call at cg19473623 × kidney cortex is independently supported: the mQTL
lead rs10257835 carries an eQTL of the same sign (z = +8.5 vs mQTL z = +18.1), and
PPH4 stays at 0.86–0.9995 across the whole prior grid
(W = 0.01–0.10 × p12 = 1e-6–1e-4).

## Layout and how to run

```
code/M4_mediation/mr_rebuild_2026-11/   the eight R scripts + run_all.sh
data/M4_mediation/mr_rebuild_2026-11/
    inputs/   two IgAN extracts (chr7 locus window; genome-wide significant set)
    out/      19 result tables, including the three Table S5b/S5c/S5d sources
results/M4_mediation/                   the six figure files
```

```bash
cd code/M4_mediation/mr_rebuild_2026-11
./run_all.sh            # needs Rscript on PATH; override with RSCRIPT=/path/to/Rscript
```

`run_all.sh` only exports paths — every script reads them from
`00_config.R` (`MR_ROOT`, `MR_DATA_DIR`, `MQTL_DIR`, `EQTL_CACHE`, `MR_FIG_DIR`),
and no absolute path is hard-coded anywhere.

### External inputs you must supply

| Variable | What it needs | Where to get it |
|---|---|---|
| `MQTL_DIR` | the five `cpg_full_assoc_<CpG>.tsv` GoDMC cis-mQTL tables | already in `data/M4_mediation/`; `step_m4_1_fetch_mqtl.py` documents the query |
| `EQTL_CACHE` | eQTL Catalogue r8 region cache for the four expression contexts (GTEx v10 whole blood / kidney cortex, OneK1K naive / memory / intermediate B) | eQTL Catalogue FTP, `QTD000356`, `QTD000261`, `QTD000606–608` |

`03_coloc.R` is the only script that touches `EQTL_CACHE`; the other seven run
without it.

## Reproducibility notes

- R 4.3.3; `MendelianRandomization`, `MRPRESSO`, `RadialMR`, `coloc`, `ggplot2`,
  `patchwork`, `data.table`.
- Scripts are ASCII-only on purpose: on Windows, `Rscript` decodes source files with
  the local code page, so a Chinese path inside a script aborts `source()`.
- `07_simulation.R` uses a fixed seed; the reported 44-row power curve is stable
  across runs. The two `cairo_pdf` outputs embed a creation timestamp, so their
  md5 changes between runs while the PNG/TIFF and every TSV stay byte-identical —
  check reproducibility on the TSV/PNG, not on the PDF.
