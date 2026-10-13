#!/usr/bin/env Rscript
# Dispersion check (betadisper) re-run with a robust p extraction. ASCII-only source.
suppressPackageStartupMessages(library(vegan))
outv <- "D:/CodexWork/mic-re/out_analysis"
fm <- read.csv(file.path(outv, "feature_matrix.csv"), check.names = FALSE, stringsAsFactors = FALSE)
meta_cols <- c("project","base","crop","sample","all_reads","bacteria_reads","organelle_reads",
               "genus_reads","family_reads","unassigned_reads","n_asv")
feats <- setdiff(colnames(fm), meta_cols)
fm <- fm[fm$bacteria_reads >= 1000, ]
X <- as.matrix(fm[, feats]); rownames(X) <- paste(fm$project, fm$sample, sep = "|")
X <- X[, colSums(X) > 0, drop = FALSE]
P <- sweep(X, 1, rowSums(X), "/")
zf <- function(Pm) { q <- Pm
  for (i in seq_len(nrow(q))) { z <- q[i, ] == 0
    if (any(z)) { mn <- min(q[i, q[i, ] > 0]); q[i, z] <- 0.65 * mn; q[i, ] <- q[i, ] / sum(q[i, ]) } }
  log(q) - rowMeans(log(q)) }
crop <- factor(fm$crop, levels = c("maize","rice","wheat"))
L <- c("# Dispersion check (betadisper) - section 11 requirement", "",
       "Homogeneous dispersion is required for a PERMANOVA R2 to be interpreted as a location (composition) effect.",
       "`permutest` tests the crop factor for differences in within-group dispersion.", "",
       "| distance | F | p (permutest, 999) | mean distance to centroid (maize / rice / wheat) |", "| --- | --- | --- | --- |")
for (nm in c("Bray-Curtis","Aitchison")) {
  D <- if (nm == "Bray-Curtis") vegdist(P, method = "bray") else dist(zf(P))
  bd <- betadisper(D, crop)
  pt <- permutest(bd, permutations = 999)
  pv <- pt$tab$`Pr(>F)`[which(rownames(pt$tab) == "Groups")]
  fv <- pt$tab$F[which(rownames(pt$tab) == "Groups")]
  if (length(pv) == 0) { pv <- pt$tab[2, ncol(pt$tab)]; fv <- pt$tab[2, "F"] }
  L <- c(L, sprintf("| %s | %.3f | %.3f | %.4f / %.4f / %.4f |", nm, fv, pv,
                    bd$group.distances[1], bd$group.distances[2], bd$group.distances[3]))
}
writeLines(L, file.path(outv, "dispersion_check.md"))
cat(paste(L, collapse = "\n"), "\n")
