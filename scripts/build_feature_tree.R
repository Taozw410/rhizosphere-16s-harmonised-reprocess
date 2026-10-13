#!/usr/bin/env Rscript
# Align feature representatives and build an ML tree (DECIPHER) for UniFrac. ASCII-only source.
# usage: Rscript build_feature_tree.R
suppressPackageStartupMessages(library(DECIPHER)); suppressPackageStartupMessages(library(Biostrings))
outv <- "D:/CodexWork/mic-re/out_analysis"
say <- function(...) cat(sprintf(...), "\n")
al_f <- file.path(outv, "feature_aln.fa")
if (file.exists(al_f)) {
  al <- readDNAStringSet(al_f); say("reuse existing alignment: %d x %d", length(al), width(al)[1])
} else {
  seqs <- readDNAStringSet(file.path(outv, "feature_reps.fa"))
  say("feature representatives: %d ; length %d-%d", length(seqs), min(width(seqs)), max(width(seqs)))
  t0 <- Sys.time()
  al <- AlignSeqs(seqs, verbose = FALSE, processors = 2)
  say("alignment done in %.1f min ; width %d", as.numeric(difftime(Sys.time(), t0, units = "mins")), width(al)[1])
  writeXStringSet(al, al_f)
}
t0 <- Sys.time()
tr <- try(Treeline(myXStringSet = al, method = "ML", processors = 2, verbose = FALSE), silent = TRUE)
if (inherits(tr, "try-error")) {
  say("Treeline on the alignment failed (%s) -> using a distance-matrix route", substr(as.character(attr(tr, "condition")$message), 1, 90))
  d <- DistanceMatrix(al, type = "dist", verbose = FALSE)
  tr <- Treeline(myDist = d, method = "ML", processors = 2, verbose = FALSE)
}
say("tree built in %.1f min ; tips %d", as.numeric(difftime(Sys.time(), t0, units = "mins")), length(tr$tip.label))
ape::write.tree(tr, file.path(outv, "feature_tree.nwk"))
say("wrote out_analysis/feature_tree.nwk")
