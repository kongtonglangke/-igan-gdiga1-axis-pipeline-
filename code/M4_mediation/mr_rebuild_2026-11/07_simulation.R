## =====================================================================
## 07_simulation.R  -  estimator calibration (type I error and power)
## Uses the OBSERVED instrument structure (bx, bxse) and regenerates the
## outcome under a known truth, so the operating characteristics of exactly
## this design are quantified rather than asserted.
## out: out/12_simulation.tsv
## =====================================================================
.this_dir <- dirname(normalizePath(sub("^--file=", "",
  commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]),
  winslash = "/", mustWork = FALSE))
source(file.path(.this_dir, "00_config.R"))
set.seed(20261110)

d <- fread(file.path(DIR_OUT, "02_harmonised.tsv"), sep = "\t", data.table = FALSE)
d <- d[d$ok, ]

NSIM <- 10000
sim_one <- function(x, true_b) {
  bx <- x$beta_exposure; bxse <- x$se_exposure
  byse <- x$se_outcome
  by <- true_b * bx + rnorm(length(bx), 0, byse)
  b  <- sum(bx * by / byse^2) / sum(bx^2 / byse^2)
  s  <- sqrt(1 / sum(bx^2 / byse^2))
  2 * pnorm(-abs(b / s))
}

res <- list()
for (cpg in unique(d$cpg)) {
  x <- d[d$cpg == cpg, ]
  ## type I error under the null
  p0 <- replicate(NSIM, sim_one(x, 0))
  ## power at several true effect sizes
  pw <- sapply(c(log(1.10), log(1.20), log(1.30), log(1.50)), function(tb)
    mean(replicate(NSIM / 2, sim_one(x, tb)) < 0.05))
  res[[cpg]] <- data.frame(cpg = cpg, k = nrow(x),
    type1_error = mean(p0 < 0.05),
    power_OR1.10 = pw[1], power_OR1.20 = pw[2],
    power_OR1.30 = pw[3], power_OR1.50 = pw[4], stringsAsFactors = FALSE)
}
sim <- do.call(rbind, res); rownames(sim) <- NULL
cat("[07] calibrated operating characteristics (10,000 / 5,000 replicates per cell):\n")
print(sim, digits = 3)

write.table(sim, file.path(DIR_OUT, "12_simulation.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## ---- power curve over a grid of true ORs (for plotting) ----
grid <- seq(1.0, 1.5, by = 0.05)
curve <- do.call(rbind, lapply(unique(d$cpg), function(cpg) {
  x <- d[d$cpg == cpg, ]
  data.frame(cpg = cpg, k = nrow(x), OR = grid,
             power = sapply(grid, function(o) mean(replicate(1000, sim_one(x, log(o))) < 0.05)))
}))
write.table(curve, file.path(DIR_OUT, "12b_power_curve.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("\n[07] power curve written (", nrow(curve), "rows )\n")
