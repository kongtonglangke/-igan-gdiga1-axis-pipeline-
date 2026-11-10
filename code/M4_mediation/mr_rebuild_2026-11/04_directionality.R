## =====================================================================
## 04_directionality.R  -  Steiger directionality + sample-overlap diagnostics
## out: out/07_steiger.tsv   out/08_sample_overlap.tsv
## =====================================================================
.this_dir <- dirname(normalizePath(sub("^--file=", "",
  commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]),
  winslash = "/", mustWork = FALSE))
source(file.path(.this_dir, "00_config.R"))

d  <- fread(file.path(DIR_OUT, "02_harmonised.tsv"), sep = "\t", data.table = FALSE)
d  <- d[d$ok, ]
iv <- fread(file.path(DIR_OUT, "01_instrument_candidates.tsv"), sep = "\t", data.table = FALSE)

## ---------- helper: variance explained ----------
## r2 = beta^2 / (beta^2 + se^2 * 2 * N * MAF * (1 - MAF))  (TwoSampleMR get_r_from_bsen equivalent)
r2_from <- function(b, se, n, maf) {
  maf <- pmin(maf, 1 - maf)
  b^2 / (b^2 + se^2 * 2 * n * maf * (1 - maf))
}

steig <- list()
for (cpg in unique(d$cpg)) {
  x <- d[d$cpg == cpg, ]
  ne <- median(x$n_exposure, na.rm = TRUE)
  r2e <- sum(r2_from(x$beta_exposure, x$se_exposure, ne, x$eaf_exposure))
  r2o <- sum(r2_from(x$beta_outcome,  x$se_outcome,  N_IGAN_EFF, x$eaf_outcome))
  r2e <- min(r2e, 0.999); r2o <- min(r2o, 0.999)
  re <- sqrt(r2e); ro <- sqrt(r2o)
  ## Steiger test on Fisher-z scale
  zstat <- (atanh(re) - atanh(ro)) / sqrt(1 / (ne - 3) + 1 / (N_IGAN_EFF - 3))
  p_up  <- pnorm(zstat, lower.tail = FALSE)   # H1: exposure -> outcome
  steig[[cpg]] <- data.frame(cpg = cpg, k = nrow(x), n_exposure = ne, n_outcome = N_IGAN_EFF,
    R2_exposure = r2e, R2_outcome = r2o, ratio = r2e / r2o,
    steiger_z = zstat, p_direction = p_up,
    verdict = ifelse(p_up < 0.05, "exposure precedes outcome", "direction not supported"))
}
## CpGs without instruments
for (cpg in setdiff(CPG_INFO$cpg, unique(d$cpg)))
  steig[[cpg]] <- data.frame(cpg = cpg, k = 0, n_exposure = NA, n_outcome = N_IGAN_EFF,
    R2_exposure = NA, R2_outcome = NA, ratio = NA, steiger_z = NA, p_direction = NA,
    verdict = "not estimable (no instrument)")
st <- do.call(rbind, steig); rownames(st) <- NULL
cat("[04] Steiger directionality:\n"); print(st, digits = 4)

## ---------- sample overlap ----------
## (1) documented independence, (2) empirical residual correlation,
## (3) analytic statement about what overlap can and cannot do.
ov <- do.call(rbind, lapply(unique(d$cpg), function(cpg) {
  x <- d[d$cpg == cpg, ]
  bx <- x$beta_exposure; by <- x$beta_outcome
  b_ivw <- sum(bx * by / x$se_outcome^2) / sum(bx^2 / x$se_outcome^2)
  resid <- by - b_ivw * bx
  r_raw  <- if (nrow(x) > 2) cor(bx, by) else NA
  r_res  <- if (nrow(x) > 2) cor(bx, resid) else NA
  data.frame(cpg = cpg, k = nrow(x),
             cor_bx_by = r_raw, cor_bx_resid = r_res,
             p_cor_resid = if (nrow(x) > 2) cor.test(bx, resid)$p.value else NA,
             stringsAsFactors = FALSE)
}))
cat("\n[04] empirical exposure-outcome residual correlation (overlap screen):\n")
print(ov, digits = 4)

## ---------- variance explained per instrument (for STROBE-MR item 6) ----------
ve <- merge(d[, c("cpg", "rsid", "eaf_exposure", "beta_exposure", "se_exposure",
                  "n_exposure", "Fstat", "beta_outcome", "se_outcome")],
            iv[, c("cpg", "rsid", "maf")], by = c("cpg", "rsid"), all.x = TRUE)
ve$R2 <- r2_from(ve$beta_exposure, ve$se_exposure, ve$n_exposure, ve$maf)
ve$Fstat_alt <- (ve$beta_exposure / ve$se_exposure)^2

write.table(st, file.path(DIR_OUT, "07_steiger.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ov, file.path(DIR_OUT, "08_sample_overlap.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ve, file.path(DIR_OUT, "09_variance_explained.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("\n[04] N exposure (GoDMC) vs N outcome (IgAN): %s vs %s\n",
            paste(unique(round(ve$n_exposure)), collapse = "/"),
            format(N_IGAN, big.mark = ",")))
cat("[04] Design: GoDMC (36 European cohorts, methylation) and GCST90018866 (IgAN meta-analysis)\n")
cat("      share no cohort by design; overlap cannot be excluded for the European component.\n")
cat("      Under a true null, sample overlap inflates type I error but cannot create a null;\n")
cat("      a null estimate therefore cannot be an artefact of sample overlap.\n")
