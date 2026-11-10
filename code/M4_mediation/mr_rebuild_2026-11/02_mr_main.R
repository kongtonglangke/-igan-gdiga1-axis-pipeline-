## =====================================================================
## 02_mr_main.R  -  univariable MR, all sensitivity estimators, diagnostics
## out: out/03_mr_estimates.tsv        main estimates by method
##      out/04_mr_diagnostics.tsv      F, Q, I2, Egger intercept, LOO, PRESSO/Radial
##      out/05_single_snp.tsv          per-SNP Wald ratios
## =====================================================================
.this_dir <- dirname(normalizePath(sub("^--file=", "",
  commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]),
  winslash = "/", mustWork = FALSE))
source(file.path(.this_dir, "00_config.R"))

harm <- fread(file.path(DIR_OUT, "02_harmonised.tsv"), sep = "\t", data.table = FALSE)
d    <- harm[harm$ok, ]
d    <- d[order(d$cpg, d$pos), ]
cat(sprintf("[02] analysable instruments: %d SNPs / %d CpGs\n", nrow(d), length(unique(d$cpg))))
print(table(d$cpg))


## ---------- single-SNP Wald ratios ----------
ss <- do.call(rbind, lapply(split(d, d$cpg), function(x) {
  data.frame(cpg = x$cpg, rsid = x$rsid, pos = x$pos,
             bx = x$beta_exposure, bxse = x$se_exposure,
             by = x$beta_outcome,  byse = x$se_outcome,
             wald_beta = x$beta_outcome / x$beta_exposure,
             wald_se   = x$se_outcome / abs(x$beta_exposure),
             Fstat = x$Fstat, stringsAsFactors = FALSE)
}))
ss$wald_p <- 2 * pnorm(-abs(ss$wald_beta / ss$wald_se))
ss$wald_or <- exp(ss$wald_beta)
cat("\n[02] single-SNP Wald ratios:\n")
print(ss[, c("cpg", "rsid", "bx", "by", "wald_beta", "wald_se", "wald_p", "wald_or")], digits = 4)

## ---------- multi-method estimation ----------
run_cpg <- function(cpg) {
  x <- d[d$cpg == cpg, ]
  k <- nrow(x)
  bx <- x$beta_exposure; bxse <- x$se_exposure
  by <- x$beta_outcome;  byse <- x$se_outcome
  mrobj <- mr_input(bx = bx, bxse = bxse, by = by, byse = byse,
                    exposure = cpg, outcome = "IgAN", snps = x$rsid)
  out <- list(); note <- character(0)

  ## Wald ratio: defined only for ONE instrument.  With k >= 2 the inverse-variance
  ## expression below reproduces IVW (fixed) exactly, so emitting it under the label
  ## "Wald ratio" would duplicate the IVW row under a name that does not apply.
  if (k == 1) {
    b_w <- by[1] / bx[1]
    se_w <- byse[1] / abs(bx[1])
    p_w <- 2 * pnorm(-abs(b_w / se_w))
    out[[length(out) + 1]] <- data.frame(cpg = cpg, method = "Wald ratio", k = k,
      b = b_w, se = se_w, p = p_w, or = exp(b_w),
      lo = exp(b_w - 1.96 * se_w), hi = exp(b_w + 1.96 * se_w))
  }

  if (k >= 2) {
    for (md in c("fixed", "random")) {
      r <- mr_ivw(mrobj, model = md)
      out[[length(out) + 1]] <- data.frame(cpg = cpg, method = paste0("IVW (", md, ")"), k = k,
        b = r@Estimate, se = r@StdError, p = r@Pvalue, or = exp(r@Estimate),
        lo = exp(r@Estimate - 1.96 * r@StdError), hi = exp(r@Estimate + 1.96 * r@StdError))
    }
  }
  if (k >= 3) {
    eg <- mr_egger(mrobj, robust = FALSE)
    out[[length(out) + 1]] <- data.frame(cpg = cpg, method = "MR-Egger (slope)", k = k,
      b = eg@Estimate, se = eg@StdError.Est, p = eg@Pvalue.Est, or = exp(eg@Estimate),
      lo = exp(eg@Estimate - 1.96 * eg@StdError.Est), hi = exp(eg@Estimate + 1.96 * eg@StdError.Est))
    wm <- mr_median(mrobj, weighting = "weighted")
    out[[length(out) + 1]] <- data.frame(cpg = cpg, method = "Weighted median", k = k,
      b = wm@Estimate, se = wm@StdError, p = wm@Pvalue, or = exp(wm@Estimate),
      lo = exp(wm@Estimate - 1.96 * wm@StdError), hi = exp(wm@Estimate + 1.96 * wm@StdError))
    pm <- mr_median(mrobj, weighting = "penalized")
    out[[length(out) + 1]] <- data.frame(cpg = cpg, method = "Penalised median", k = k,
      b = pm@Estimate, se = pm@StdError, p = pm@Pvalue, or = exp(pm@Estimate),
      lo = exp(pm@Estimate - 1.96 * pm@StdError), hi = exp(pm@Estimate + 1.96 * pm@StdError))
    mb <- mr_mbe(mrobj, weighting = "weighted")
    out[[length(out) + 1]] <- data.frame(cpg = cpg, method = "Weighted mode", k = k,
      b = mb@Estimate, se = mb@StdError, p = mb@Pvalue, or = exp(mb@Estimate),
      lo = exp(mb@Estimate - 1.96 * mb@StdError), hi = exp(mb@Estimate + 1.96 * mb@StdError))
    ml <- mr_maxlik(mrobj)
    out[[length(out) + 1]] <- data.frame(cpg = cpg, method = "Maximum likelihood", k = k,
      b = ml@Estimate, se = ml@StdError, p = ml@Pvalue, or = exp(ml@Estimate),
      lo = exp(ml@Estimate - 1.96 * ml@StdError), hi = exp(ml@Estimate + 1.96 * ml@StdError))
  } else {
    note <- c(note, sprintf("k=%d: Egger / weighted median / mode / maxlik not estimable (need >=3 IVs)", k))
  }
  list(est = do.call(rbind, out), note = note, x = x, k = k)
}

res   <- lapply(unique(d$cpg), run_cpg)
names(res) <- unique(d$cpg)
est   <- do.call(rbind, lapply(res, `[[`, "est"))

## CpGs with zero instruments
zero <- setdiff(CPG_INFO$cpg, unique(d$cpg))
if (length(zero)) {
  est <- rbind(est, data.frame(cpg = zero, method = "not estimable", k = 0,
    b = NA, se = NA, p = NA, or = NA, lo = NA, hi = NA))
}

## interpretable scaling: OR per 0.1 unit of methylation M-value
est$or_per0.1 <- exp(est$b * 0.1)
est$lo_per0.1 <- exp((est$b - 1.96 * est$se) * 0.1)
est$hi_per0.1 <- exp((est$b + 1.96 * est$se) * 0.1)

cat("\n[02] MR estimates:\n")
print(est[order(est$cpg, est$method), c("cpg", "method", "k", "b", "se", "p", "or", "or_per0.1")], digits = 4)

## ---------- diagnostics ----------
diag_rows <- list()
for (cpg in names(res)) {
  r <- res[[cpg]]; x <- r$x; k <- r$k
  bx <- x$beta_exposure; bxse <- x$se_exposure
  by <- x$beta_outcome;  byse <- x$se_outcome
  ## instrument strength
  Fmean <- mean(bx^2 / bxse^2); Fmin <- min(bx^2 / bxse^2)
  ## variance explained (R2) on the exposure
  n_exp <- median(x$n_exposure, na.rm = TRUE)
  maf   <- x$eaf_exposure; maf <- pmin(maf, 1 - maf)
  r2    <- (bx^2) / (bx^2 + 2 * n_exp * maf * (1 - maf) * bxse^2)  # per-SNP R2 (TwoSampleMR formula)
  R2    <- sum(r2)
  ## IVW heterogeneity
  Q <- Qdf <- Qp <- I2 <- NA
  if (k >= 2) {
    ## Cochran's Q on the ratio scale (standard IVW heterogeneity statistic)
    b_ivw <- sum(bx * by / byse^2) / sum(bx^2 / byse^2)
    Q     <- sum((by - b_ivw * bx)^2 / byse^2)
    Qdf   <- k - 1
    Qp    <- pchisq(Q, Qdf, lower.tail = FALSE)
    I2    <- ifelse(Q > 0, max(0, (Q - Qdf) / Q), 0)
  }
  ## Egger intercept
  eint <- eintse <- eintp <- NA
  if (k >= 3) {
    eg <- mr_egger(mr_input(bx, bxse, by, byse), robust = FALSE)
    eint <- eg@Intercept; eintse <- eg@StdError.Int; eintp <- eg@Pvalue.Int
  }
  ## LOO (manual IVW refit; avoids package API drift)
  loo_min <- loo_max <- loo_minp <- NA
  if (k >= 3) {
    loo_b <- loo_p <- numeric(0)
    for (j in seq_len(k)) {
      bb <- sum(bx[-j] * by[-j] / byse[-j]^2) / sum(bx[-j]^2 / byse[-j]^2)
      se_loo <- sqrt(1 / sum(bx[-j]^2 / byse[-j]^2))
      loo_b <- c(loo_b, bb); loo_p <- c(loo_p, 2 * pnorm(-abs(bb / se_loo)))
    }
    loo_min <- min(loo_b); loo_max <- max(loo_b); loo_minp <- min(loo_p)
  }
  ## MR-PRESSO (needs >= 4 IVs)
  presso_p <- presso_nout <- NA
  if (k >= 4) {
    pr <- try(mr_presso(BetaOutcome = "by", BetaExposure = "bx", SdOutcome = "byse",
                        SdExposure = "bxse", OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                        data = data.frame(bx = bx, bxse = bxse, by = by, byse = byse),
                        NbDistribution = 2000, SignifThreshold = 0.05), silent = TRUE)
    if (!inherits(pr, "try-error")) {
      presso_p <- pr$`Main MR results`$`P-value`[1]
      presso_nout <- sum(pr$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices` > 0, na.rm = TRUE)
    }
  }
  ## Radial MR
  rad_p <- rad_nout <- NA
  if (k >= 3) {
    rr <- try(suppressWarnings(capture.output(
      obj <- ivw_radial(mr_input(bx, bxse, by, byse), alpha = 0.05, weights = 3))), silent = TRUE)
    if (exists("obj") && !inherits(rr, "try-error")) {
      rad_p    <- pchisq(obj$qstatistic, obj$df, lower.tail = FALSE)
      rad_nout <- sum(as.logical(obj$outliers), na.rm = TRUE)
      rm(obj)
    }
  }
  diag_rows[[cpg]] <- data.frame(cpg = cpg, k = k, F_mean = Fmean, F_min = Fmin,
    R2 = R2, n_exposure = n_exp, Q = Q, Q_df = Qdf, Q_p = Qp, I2 = I2,
    egger_intercept = eint, egger_intercept_se = eintse, egger_intercept_p = eintp,
    loo_beta_min = loo_min, loo_beta_max = loo_max, loo_p_min = loo_minp,
    presso_global_p = presso_p, presso_n_outlier = presso_nout,
    radial_Q_p = rad_p, radial_n_outlier = rad_nout)
}
diag <- do.call(rbind, diag_rows)
diag <- rbind(diag, data.frame(cpg = zero, k = 0, F_mean = NA, F_min = NA, R2 = NA,
  n_exposure = NA, Q = NA, Q_df = NA, Q_p = NA, I2 = NA, egger_intercept = NA,
  egger_intercept_se = NA, egger_intercept_p = NA, loo_beta_min = NA, loo_beta_max = NA,
  loo_p_min = NA, presso_global_p = NA, presso_n_outlier = NA, radial_Q_p = NA,
  radial_n_outlier = NA))
diag <- diag[match(CPG_INFO$cpg, diag$cpg), ]

## ---------- detectable effect at 80% power ----------
## se of the IVW estimate; minimum detectable |beta| = (z_0.975 + z_0.80) * se
for (i in seq_len(nrow(diag))) {
  cpg <- diag$cpg[i]
  e <- est[est$cpg == cpg & est$method %in% c("Wald ratio", "IVW (fixed)"), ]
  if (nrow(e) == 0) { diag$se_ivw[i] <- NA; next }
  se <- e$se[1]
  diag$se_ivw[i] <- se
  diag$md_beta[i] <- (qnorm(0.975) + qnorm(0.80)) * se
  diag$md_or[i]   <- exp(diag$md_beta[i])
  diag$md_or_per0.1[i] <- exp(diag$md_beta[i] * 0.1)
}
## empirical power at the observed point estimate
diag$power_obs <- ifelse(is.na(diag$se_ivw), NA, pnorm(abs(est$b[match(diag$cpg, est$cpg)]) / diag$se_ivw - qnorm(0.975)))

cat("\n[02] diagnostics:\n")
print(diag[, c("cpg", "k", "F_mean", "F_min", "R2", "Q", "Q_p", "I2", "egger_intercept_p",
               "presso_global_p", "md_or", "md_or_per0.1")], digits = 4)

## ---------- save ----------
ss <- as.data.frame(ss); rownames(ss) <- NULL
write.table(ss, file.path(DIR_OUT, "05_single_snp.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
est <- as.data.frame(est); rownames(est) <- NULL
write.table(est, file.path(DIR_OUT, "03_mr_estimates.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
diag$md_beta[diag$k < 1] <- NA; diag$md_or[diag$k < 1] <- NA; diag$md_or_per0.1[diag$k < 1] <- NA
diag <- as.data.frame(diag); rownames(diag) <- NULL
write.table(diag, file.path(DIR_OUT, "04_mr_diagnostics.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n[02] notes:\n")
for (cpg in names(res)) if (length(res[[cpg]]$note)) cat(" -", cpg, ":", res[[cpg]]$note, "\n")
cat(sprintf("\n[02] CpGs without any analysable instrument: %s\n",
            ifelse(length(zero), paste(zero, collapse = ", "), "none")))
