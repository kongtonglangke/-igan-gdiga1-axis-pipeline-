## =====================================================================
## 01_instruments.R  -  instrument selection + allele harmonisation
## out: out/01_instrument_candidates.tsv  out/02_harmonised.tsv
## =====================================================================
.this_dir <- dirname(normalizePath(sub("^--file=", "",
  commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]),
  winslash = "/", mustWork = FALSE))
source(file.path(.this_dir, "00_config.R"))

## ---------- 1. read the 5 cis-mQTL tables ----------
fread_mqtl <- function(cpg) {
  f <- file.path(MQLT_D, sprintf("cpg_full_assoc_%s.tsv", cpg))
  d <- fread(f, sep = "\t", header = TRUE, data.table = FALSE,
             colClasses = list(character = c("rsid", "snp", "a1", "a2", "direction", "cistrans")))
  d$chr     <- as.integer(sub("^chr([0-9XY]+):.*$", "\\1", d$snp))
  d$pos     <- as.integer(sub("^chr[0-9XY]+:([0-9]+):.*$", "\\1", d$snp))
  d$vtype   <- sub("^chr[0-9XY]+:[0-9]+:", "", d$snp)
  d$cpg     <- cpg
  d$clumped <- as.logical(d$clumped)
  d$beta    <- as.numeric(d$beta_a1)   # GoDMC primary meta-analysis effect (methylation M-value)
  d$se      <- as.numeric(d$se)
  d$pval    <- as.numeric(d$pval)
  d$eaf     <- as.numeric(d$freq_a1)
  d$maf     <- pmin(d$eaf, 1 - d$eaf)
  d$Fstat   <- (d$beta / d$se)^2
  d$beta_mre <- as.numeric(d$beta_a1);     d$se_mre <- as.numeric(d$se_mre);     d$p_mre <- as.numeric(d$pval_mre)
  d$beta_are <- as.numeric(d$beta_are_a1); d$se_are <- as.numeric(d$se_are);     d$p_are <- as.numeric(d$pval_are)
  d$rsid[is.na(d$rsid) | d$rsid == "" | d$rsid == "NA"] <- NA_character_
  d
}

mqtl <- do.call(rbind, lapply(CPG_INFO$cpg, fread_mqtl))
mqtl <- merge(mqtl, CPG_INFO[, c("cpg", "gene", "relation")], by = "cpg", all.x = TRUE)
cat(sprintf("[01] mQTL records: %d rows / %d CpGs\n", nrow(mqtl), length(unique(mqtl$cpg))))

## ---------- 2. instrument filter: P<5e-8 AND independent (clumped) ----------
mqtl$is_sig     <- !is.na(mqtl$pval) & mqtl$pval < P_INSTR
mqtl$is_clumped <- !is.na(mqtl$clumped) & mqtl$clumped
mqtl$is_cis     <- mqtl$cistrans %in% c("True", TRUE)
mqtl$pass_F     <- !is.na(mqtl$Fstat) & mqtl$Fstat >= F_MIN

mqtl$stage <- with(mqtl, ifelse(!is_cis, "not_cis",
                         ifelse(!is_sig, "P_ge_5e-8",
                         ifelse(!is_clumped, "not_independent",
                         ifelse(!pass_F, "weak_F", "INCLUDED")))))

iv <- mqtl[mqtl$stage == "INCLUDED", ]
cat(sprintf("[01] qualifying instruments: %d SNPs / %d CpGs\n", nrow(iv), length(unique(iv$cpg))))
print(table(iv$cpg))
cat("\n[01] per-CpG filter funnel:\n"); print(table(mqtl$cpg, mqtl$stage))

cat(sprintf("\n[01] 5 CpGs span chr%s:%s-%s (%.1f kb)\n", unique(iv$chr)[1],
            format(min(CPG_INFO$pos37), big.mark = ","),
            format(max(CPG_INFO$pos37), big.mark = ","),
            (max(CPG_INFO$pos37) - min(CPG_INFO$pos37)) / 1000))
cat("[01] instruments:\n")
print(iv[order(iv$pos), c("cpg", "rsid", "snp", "pos", "beta", "se", "pval", "maf", "Fstat")])

## ---------- 3. outcome side: IgAN GWAS, chr7 region ----------
igan <- fread(file.path(DIR_DATA, "igan_chr7_region.tsv"), sep = "\t", header = TRUE, data.table = FALSE)
names(igan) <- c("chr", "pos", "ea", "oa", "beta", "se", "eaf", "pval", "vid", "hm_conv", "hm_code", "rsid")
igan <- igan[!is.na(igan$beta) & !is.na(igan$se) & igan$se > 0, ]

match_outcome <- function(r) {
  hit <- NULL
  if (!is.na(r$rsid)) hit <- igan[igan$rsid %in% r$rsid, , drop = FALSE]
  if (is.null(hit) || nrow(hit) == 0) hit <- igan[igan$chr == r$chr & abs(igan$pos - r$pos) <= 1, , drop = FALSE]
  if (nrow(hit) == 0) return(NULL)
  hit[1, , drop = FALSE]
}

## ---------- 4. harmonisation ----------
## Palindromic (A/T, C/G) variants are resolved by allele frequency, which is the
## standard approach (equivalent to TwoSampleMR action = 2). A variant is only
## retained when the frequency difference is unambiguous; otherwise it is dropped.
COMP  <- c(A = "T", T = "A", C = "G", G = "C")
FREQ_TOL <- 0.15

harmonise_one <- function(bx, bx_ea, bx_oa, fx, by, by_ea, by_oa, fy) {
  ea <- toupper(bx_ea); oa <- toupper(bx_oa)
  ea2 <- toupper(by_ea); oa2 <- toupper(by_oa)
  palin <- (ea == COMP[[oa]] || ea == oa)
  if (!palin) {
    if (ea == ea2 && oa == oa2)                 return(list(ok = TRUE, why = "aligned",               by = by))
    if (ea == oa2 && oa == ea2)                 return(list(ok = TRUE, why = "flipped",               by = -by))
    if (COMP[[ea]] == ea2 && COMP[[oa]] == oa2) return(list(ok = TRUE, why = "chain_complement",      by = by))
    if (COMP[[ea]] == oa2 && COMP[[oa]] == ea2) return(list(ok = TRUE, why = "chain_complement_flip", by = -by))
    return(list(ok = FALSE, why = "allele_mismatch", by = NA_real_))
  }
  ## palindromic: infer strand from allele frequencies
  d_same <- abs(fx - fy); d_flip <- abs(fx - (1 - fy))
  if (d_same < FREQ_TOL && d_same <= d_flip) return(list(ok = TRUE, why = "palindromic_same_strand", by = by))
  if (d_flip < FREQ_TOL && d_flip <  d_same) return(list(ok = TRUE, why = "palindromic_flipped",     by = -by))
  list(ok = FALSE, why = "palindromic_unresolved", by = NA_real_)
}

rows <- list()
for (i in seq_len(nrow(iv))) {
  r <- iv[i, ]; o <- match_outcome(r)
  if (is.null(o)) {
    rows[[length(rows) + 1]] <- data.frame(
      cpg = r$cpg, rsid = r$rsid, snp = r$snp, chr = r$chr, pos = r$pos,
      effect_allele = r$a1, other_allele = r$a2, eaf_exposure = r$eaf,
      beta_exposure = r$beta, se_exposure = r$se, p_exposure = r$pval,
      Fstat = r$Fstat, n_exposure = r$samplesize,
      beta_outcome = NA, se_outcome = NA, p_outcome = NA, eaf_outcome = NA,
      harmonise = "absent_in_outcome", stringsAsFactors = FALSE)
    next
  }
  h <- harmonise_one(r$beta, r$a1, r$a2, r$eaf, o$beta, o$ea, o$oa, o$eaf)
  rows[[length(rows) + 1]] <- data.frame(
    cpg = r$cpg, rsid = ifelse(is.na(r$rsid), o$rsid, r$rsid), snp = r$snp, chr = r$chr, pos = r$pos,
    effect_allele = r$a1, other_allele = r$a2, eaf_exposure = r$eaf,
    beta_exposure = r$beta, se_exposure = r$se, p_exposure = r$pval,
    Fstat = r$Fstat, n_exposure = r$samplesize,
    beta_outcome = if (h$ok) h$by else NA, se_outcome = if (h$ok) o$se else NA,
    p_outcome = if (h$ok) o$pval else NA, eaf_outcome = o$eaf,
    harmonise = h$why, stringsAsFactors = FALSE)
}
harm <- do.call(rbind, rows)
harm$ok <- !is.na(harm$beta_outcome)
cat(sprintf("\n[01] harmonisation: %d instruments, %d aligned, %d failed\n",
            nrow(harm), sum(harm$ok), sum(!harm$ok)))
print(table(harm$cpg, harm$harmonise))
cat("\n[01] harmonised detail:\n")
print(harm[order(harm$cpg, harm$pos),
           c("cpg", "rsid", "pos", "effect_allele", "other_allele", "beta_exposure",
             "se_exposure", "beta_outcome", "se_outcome", "Fstat", "harmonise")], digits = 4)

## ---------- 5. save ----------
fwrite(mqtl[, c("cpg", "gene", "relation", "chr", "pos", "vtype", "rsid", "a1", "a2",
                "beta", "se", "pval", "eaf", "maf", "Fstat",
                "beta_mre", "se_mre", "p_mre", "beta_are", "se_are", "p_are",
                "samplesize", "clumped", "is_sig", "is_clumped", "pass_F", "stage")],
       file.path(DIR_OUT, "01_instrument_candidates.tsv"), sep = "\t")
fwrite(harm, file.path(DIR_OUT, "02_harmonised.tsv"), sep = "\t")

## ---------- 6. can reverse MR be instrumented at this locus? ----------
reg <- igan[igan$chr == 7 & igan$pos >= 7150000 & igan$pos <= 7350000, ]
cat(sprintf("\n[01] in-locus IgAN association: n=%d, min P=%.4g, n(P<5e-8)=%d\n",
            nrow(reg), min(reg$pval, na.rm = TRUE), sum(reg$pval < 5e-8, na.rm = TRUE)))
gws <- fread(file.path(DIR_DATA, "igan_gws_p5e8.tsv"), sep = "\t", header = TRUE, data.table = FALSE)
cat(sprintf("[01] genome-wide significant IgAN variants in total: %d\n", nrow(gws)))

writeLines(c(sprintf("mQTL rows: %d", nrow(mqtl)),
             sprintf("instruments: %d", nrow(iv)),
             sprintf("harmonised ok: %d", sum(harm$ok)),
             sprintf("in-locus IgAN min P: %.6g", min(reg$pval, na.rm = TRUE)),
             sprintf("in-locus IgAN GWS variants: %d", sum(reg$pval < 5e-8, na.rm = TRUE)),
             sprintf("genome-wide IgAN GWS variants: %d", nrow(gws))),
           file.path(DIR_LOG, "01_instruments.log"))
