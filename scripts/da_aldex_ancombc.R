#!/usr/bin/env Rscript
# Differential abundance with TWO tools (section 11 binding requirement): ALDEx2 (Welch t on CLR) and
# ANCOMBC2 (ANCOM-BC2), for the three pairwise crop contrasts.  Consensus = significant in BOTH (BH<0.05).
# ASCII-only source.
suppressPackageStartupMessages({library(ALDEx2); library(ANCOMBC); library(phyloseq)})
root <- "D:/CodexWork/mic-re"; outv <- file.path(root, "out_analysis")
say <- function(...) cat(sprintf(...), "\n")
fm <- read.csv(file.path(outv, "feature_matrix.csv"), check.names = FALSE, stringsAsFactors = FALSE)
meta_cols <- c("project","base","crop","sample","all_reads","bacteria_reads","organelle_reads",
               "genus_reads","family_reads","unassigned_reads","n_asv")
feats <- setdiff(colnames(fm), meta_cols)
fm <- fm[fm$bacteria_reads >= 1000, ]
X <- t(as.matrix(fm[, feats])); colnames(X) <- paste(fm$project, fm$sample, sep = "|")
crop <- fm$crop
# drop features absent from the whole dataset or present in < 10 samples (ANCOM-BC2 needs estimable taxa)
present <- rowSums(X > 0)
keep <- present >= 10
say("features %d -> %d kept (present in >=10 samples)", length(present), sum(keep))
X <- X[keep, , drop = FALSE]
pairs_ <- list(c("maize","rice"), c("maize","wheat"), c("rice","wheat"))
say("ALDEx2 start: %s", format(Sys.time()))
ald_rds <- file.path(outv, "da_aldex2.rds")
if (file.exists(ald_rds)) {
  ald <- readRDS(ald_rds); say("reuse cached ALDEx2 results")
} else {
  ald <- list()
  for (pr in pairs_) {
    sub <- crop %in% pr
    cond <- as.character(crop[sub])   # ALDEx2 wants a character vector, not a factor
    cl <- aldex.clr(X[, sub, drop = FALSE], cond, mc.samples = 128, denom = "all", verbose = FALSE)
    tt <- aldex.ttest(cl, paired.test = FALSE, verbose = FALSE)
    ef <- aldex.effect(cl, verbose = FALSE)
    key <- paste(pr, collapse = "_")
    ald[[key]] <- data.frame(feature = rownames(tt), we.eBH = tt$we.eBH, wi.eBH = tt$wi.eBH, effect = ef$effect)
    say("  ALDEx2 %s: %d features with we.eBH<0.05", key, sum(ald[[key]]$we.eBH < 0.05, na.rm = TRUE))
  }
  saveRDS(ald, ald_rds)            # save BEFORE the next tool so a later failure cannot lose it
}
say("ANCOMBC2 start: %s", format(Sys.time()))
anc_rds <- file.path(outv, "da_ancombc2.rds")
anc <- if (file.exists(anc_rds)) readRDS(anc_rds) else list()
for (ref in c("maize", "rice", "wheat")) {
  if (!is.null(anc[[paste0("ref_", ref)]])) { say("  ANCOMBC2 ref=%s: cached", ref); next }
  mdl <- data.frame(crop = factor(crop, levels = c(ref, setdiff(c("maize","rice","wheat"), ref))))
  rownames(mdl) <- colnames(X)     # phyloseq requires sample names to match the OTU table
  ps <- phyloseq(otu_table(X, taxa_are_rows = TRUE), sample_data(mdl))
  key <- paste0("ref_", ref)
  res <- try(ancombc2(data = ps, fix_formula = "crop", p_adj_method = "BH", group = "crop",
                      struc_zero = FALSE, verbose = FALSE), silent = TRUE)
  if (inherits(res, "try-error")) { say("  ANCOMBC2 failed (ref %s): %s", ref, substr(as.character(res), 1, 140)); next }
  r <- res$res
  anc[[key]] <- r
  saveRDS(anc, anc_rds)
  say("  ANCOMBC2 ref=%s: columns %s", ref, paste(head(grep("^q_crop|^diff_crop", colnames(r), value = TRUE), 6), collapse = ", "))
}
# consensus per contrast
out <- c("# Differential abundance: ALDEx2 + ANCOMBC2 consensus (section 11)", "",
         sprintf("Features analysed: %d (present in >=10 samples; depth >=1000 bacteria reads).", nrow(X)),
         "ALDEx2 = Welch t on CLR (128 MC samples, denom=all); ANCOMBC2 = ANCOM-BC2 with BH adjustment.",
         "Consensus = significant in BOTH tools at BH < 0.05.", "",
         "| contrast | ALDEx2 we.eBH<0.05 | ANCOMBC2 q<0.05 | consensus |", "| --- | --- | --- | --- |")
for (pr in pairs_) {
  key <- paste(pr, collapse = "_")
  a <- ald[[key]]; a_sig <- a$feature[a$we.eBH < 0.05]
  k2 <- paste0("ref_", pr[1]); b_sig <- character(0); n_anc <- 0
  if (!is.null(anc[[k2]])) {
    r <- anc[[k2]]; lvl <- pr[2]
    qcol <- grep(paste0("^q_crop", lvl, "$"), colnames(r), value = TRUE)
    if (length(qcol)) { n_anc <- sum(r[[qcol]] < 0.05, na.rm = TRUE)
      b_sig <- r[[1]][which(r[[qcol]] < 0.05)] }
  }
  cons <- intersect(a_sig, b_sig)
  out <- c(out, sprintf("| %s vs %s | %d | %d | **%d** |", pr[1], pr[2], length(a_sig), n_anc, length(cons)))
  if (length(cons)) {
    aa <- a[a$feature %in% cons, ]
    write.csv(aa[order(-abs(aa$effect)), ], file.path(outv, sprintf("da_consensus_%s_%s.csv", pr[1], pr[2])), row.names = FALSE)
  }
}
if (length(ald)) saveRDS(ald, file.path(outv, "da_aldex2.rds"))
if (length(anc)) saveRDS(anc, file.path(outv, "da_ancombc2.rds"))
out <- c(out, "", "> Per-contrast consensus feature tables are written as da_consensus_<a>_<b>.csv.",
         "> Effect sizes: ALDEx2 effect (difference in CLR units, median-based) is reported; ANCOM-BC2 LFC is kept in the RDS.")
writeLines(out, file.path(outv, "da_report.md"))
say("wrote out_analysis/da_report.md")
cat(paste(out, collapse = "\n"), "\n")

