## =====================================================================
## 06_reverse.R  -  is a reverse (IgAN -> methylation) MR possible here?
## A reverse MR needs IgAN risk variants as instruments for IgAN; at the
## C1GALT1 locus no variant reaches the instrument threshold, so the reverse
## direction is not estimable. This script documents that with numbers.
## out: out/11_reverse_feasibility.tsv
## =====================================================================
.this_dir <- dirname(normalizePath(sub("^--file=", "",
  commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]),
  winslash = "/", mustWork = FALSE))
source(file.path(.this_dir, "00_config.R"))

igan <- fread(file.path(DIR_DATA, "igan_chr7_region.tsv"), sep = "\t", header = TRUE, data.table = FALSE)
names(igan) <- c("chr", "pos", "ea", "oa", "beta", "se", "eaf", "pval", "vid", "hm_conv", "hm_code", "rsid")
igan <- igan[!is.na(igan$beta) & !is.na(igan$se) & igan$se > 0, ]

d <- fread(file.path(DIR_OUT, "02_harmonised.tsv"), sep = "\t", data.table = FALSE)
d <- d[d$ok, ]

## locus window in outcome (GRCh38) coordinates
WIN_LO <- 7100000; WIN_HI <- 7400000
w <- igan[igan$chr == 7 & igan$pos >= WIN_LO & igan$pos <= WIN_HI, ]

## IgAN-genome-wide-significant variants inside the window
w$z <- w$beta / w$se
w$F_equiv <- w$z^2

feas <- data.frame(
  item = c("IgAN GWAS variants in locus window (chr7:7.10-7.40 Mb, GRCh38)",
           "variants reaching P<5e-8 (instrument threshold)",
           "variants reaching P<1e-5 (relaxed threshold)",
           "smallest P in window",
           "F-statistic of the best in-window variant",
           "F-statistic of the best in-window variant at P<5e-8",
           "mQTL instruments with IgAN P<0.05"),
  value = c(nrow(w),
            sum(w$pval < 5e-8),
            sum(w$pval < 1e-5),
            signif(min(w$pval), 4),
            round(max(w$F_equiv), 2),
            ifelse(sum(w$pval < 5e-8) > 0,
                   round(max(w$F_equiv[w$pval < 5e-8]), 2), NA),
            sum(d$p_outcome < 0.05)),
  stringsAsFactors = FALSE)

cat("[06] reverse-MR feasibility:\n"); print(feas)

cat("\n[06] association of the 7 forward instruments with IgAN:\n")
print(d[, c("cpg", "rsid", "beta_outcome", "se_outcome", "p_outcome", "Fstat")], digits = 4)

## Steiger / direction evidence is reported in 04_directionality.R
verdict <- ifelse(sum(w$pval < 5e-8) == 0,
  "Reverse MR not estimable: no IgAN risk variant at the locus reaches P<5e-8, so no valid instrument for IgAN -> methylation exists. Direction is instead assessed with the Steiger test.",
  "Reverse MR estimable: IgAN risk variants are available at the locus.")
cat("\n[06]", verdict, "\n")

write.table(feas, file.path(DIR_OUT, "11_reverse_feasibility.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
