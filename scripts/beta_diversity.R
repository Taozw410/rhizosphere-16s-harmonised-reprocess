#!/usr/bin/env Rscript
# Beta diversity, final 45 projects (pre-registration section 11). ASCII-only source (GBK issue).
# Bray-Curtis + Aitchison (CLR+Euclidean); adonis2 PERMANOVA at sample and project level; betadisper;
# leave-one-project-out sensitivity (R2 range and share of runs with p<0.05).
suppressPackageStartupMessages(library(vegan))
root <- "D:/CodexWork/mic-re"; outv <- file.path(root, "out_analysis")
fm <- read.csv(file.path(outv, "feature_matrix.csv"), check.names = FALSE, stringsAsFactors = FALSE)
meta_cols <- c("project","base","crop","sample","all_reads","bacteria_reads","organelle_reads",
               "genus_reads","family_reads","unassigned_reads","n_asv")
feats <- setdiff(colnames(fm), meta_cols)
MINDEPTH <- 1000
n0 <- nrow(fm); fm <- fm[fm$bacteria_reads >= MINDEPTH, ]
X <- as.matrix(fm[, feats]); rownames(X) <- paste(fm$project, fm$sample, sep = "|")
X <- X[, colSums(X) > 0, drop = FALSE]
P <- sweep(X, 1, rowSums(X), "/")
cat(sprintf("samples %d -> %d kept (>= %d reads; dropped %d); features %d\n", n0, nrow(fm), MINDEPTH, n0-nrow(fm), ncol(P)))

zf <- function(Pm) { q <- Pm
  for (i in seq_len(nrow(q))) { z <- q[i, ] == 0
    if (any(z)) { mn <- min(q[i, q[i, ] > 0]); if (!is.finite(mn)) mn <- 1e-6
      q[i, z] <- 0.65 * mn; q[i, ] <- q[i, ] / sum(q[i, ]) } }
  log(q) - rowMeans(log(q)) }
D_bc <- vegdist(P, method = "bray"); D_ai <- dist(zf(P))
crop <- factor(fm$crop, levels = c("maize","rice","wheat")); proj <- fm$project

L <- c("# Beta diversity (final 45 projects; pre-registration section 11)", "",
       sprintf("Samples: %d (depth >= %d bacteria reads; %d dropped). Features: %d.", nrow(fm), MINDEPTH, n0-nrow(fm), ncol(P)),
       "Distances: Bray-Curtis and Aitchison (CLR after multiplicative zero replacement, delta = 0.65 * row-wise min positive).",
       "PERMANOVA = vegan::adonis2, 999 permutations.", "",
       "## 1. PERMANOVA", "", "| distance | level | R2 | F | p |", "| --- | --- | --- | --- | --- |")

agg <- aggregate(P, by = list(project = fm$project, crop = as.character(fm$crop)), FUN = mean)
Pa <- agg[, intersect(colnames(agg), colnames(P)), drop = FALSE]; Pa <- sweep(Pa, 1, rowSums(Pa), "/")
Da_bc <- vegdist(Pa, method = "bray"); Da_ai <- dist(zf(Pa))

for (nm in c("Bray-Curtis","Aitchison")) {
  D <- if (nm == "Bray-Curtis") D_bc else D_ai
  a <- adonis2(D ~ crop, data = fm, permutations = 999)
  L <- c(L, sprintf("| %s | sample | %.4f | %.2f | %s |", nm, a$R2[1], a$F[1], ifelse(a$`Pr(>F)`[1] <= 0.001, "<0.001", sprintf("%.3f", a$`Pr(>F)`[1]))))
  Da <- if (nm == "Bray-Curtis") Da_bc else Da_ai
  a2 <- adonis2(Da ~ crop, data = data.frame(crop = agg$crop), permutations = 999)
  L <- c(L, sprintf("| %s | project (1 point each) | %.4f | %.2f | %s |", nm, a2$R2[1], a2$F[1],
                    ifelse(a2$`Pr(>F)`[1] <= 0.001, "<0.001", sprintf("%.3f", a2$`Pr(>F)`[1]))))
}
L <- c(L, "", "## 2. Dispersion check (betadisper; required by section 11)", "",
       "| distance | permutest p | mean distance to centroid (maize / rice / wheat) |", "| --- | --- | --- |")
for (nm in c("Bray-Curtis","Aitchison")) {
  D <- if (nm == "Bray-Curtis") D_bc else D_ai
  bd <- betadisper(D, crop); bt <- permutest(bd, permutations = 999)
  L <- c(L, sprintf("| %s | %.3f | %.4f / %.4f / %.4f |", nm, bt$tab$`Pr(>F)`[2],
                    bd$group.distances[1], bd$group.distances[2], bd$group.distances[3]))
}
L <- c(L, "", "## 3. Leave-one-project-out sensitivity (sample level)", "",
       "| distance | R2 range | runs with p<0.05 | runs |", "| --- | --- | --- | --- |")
for (nm in c("Bray-Curtis","Aitchison")) {
  D <- if (nm == "Bray-Curtis") D_bc else D_ai
  M <- as.matrix(D); r2s <- c(); hits <- 0; k <- 0
  for (p in unique(proj)) {
    idx <- proj != p
    if (length(unique(crop[idx])) < 3) next
    aa <- adonis2(as.dist(M[idx, idx]) ~ crop[idx], permutations = 199)
    r2s <- c(r2s, aa$R2[1]); hits <- hits + (aa$`Pr(>F)`[1] <= 0.05); k <- k + 1
  }
  L <- c(L, sprintf("| %s | %.4f - %.4f | %d / %d | %d |", nm, min(r2s), max(r2s), hits, k, k))
}
writeLines(L, file.path(outv, "beta_diversity_report.md"))
cat(paste(L, collapse = "\n"), "\n")
