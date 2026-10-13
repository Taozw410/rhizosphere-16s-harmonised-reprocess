#!/usr/bin/env Rscript
# Variance decomposition: how much of the community variation is crop vs study, and how much of the
# crop signal is inseparable from study-level technical covariates (amplicon region / depth / ASV count /
# sample count / taxonomic resolution).  ASCII-only source.
suppressPackageStartupMessages(library(vegan))
root <- "D:/CodexWork/mic-re"; outv <- file.path(root, "out_analysis")
fm <- read.csv(file.path(outv, "feature_matrix.csv"), check.names = FALSE, stringsAsFactors = FALSE)
meta_cols <- c("project","base","crop","sample","all_reads","bacteria_reads","organelle_reads",
               "genus_reads","family_reads","unassigned_reads","n_asv")
feats <- setdiff(colnames(fm), meta_cols)
fm <- fm[fm$bacteria_reads >= 1000, ]
P <- as.matrix(fm[, feats]); P <- P[, colSums(P) > 0, drop = FALSE]
Rr <- sweep(P, 1, rowSums(P), "/")
zf <- function(Pm) { q <- Pm
  for (i in seq_len(nrow(q))) { z <- q[i, ] == 0
    if (any(z)) { mn <- min(q[i, q[i, ] > 0]); q[i, z] <- 0.65 * mn; q[i, ] <- q[i, ] / sum(q[i, ]) } }
  log(q) - rowMeans(log(q)) }
D_s_bc <- vegdist(Rr, method = "bray"); D_s_ai <- dist(zf(Rr))
say <- function(...) cat(sprintf(...), "\n")
L <- c("# Variance decomposition: crop vs study vs study-level technical covariates", "",
       sprintf("Sample level: %d samples from 45 projects. Study level: 45 project means.", nrow(fm)), "",
       "## 1. Sample level (context: how much variation sits between studies)", "",
       "| model | R2 | F | p |", "| --- | --- | --- | --- |")
for (nm in c("Bray-Curtis","Aitchison")) {
  D <- if (nm == "Bray-Curtis") D_s_bc else D_s_ai
  a <- adonis2(D ~ fm$project, permutations = 999)
  L <- c(L, sprintf("| %s ~ project | %.4f | %.2f | %s |", nm, a$R2[1], a$F[1],
                    ifelse(a$`Pr(>F)`[1] <= 0.001, "<0.001", sprintf("%.3f", a$`Pr(>F)`[1]))))
  b <- adonis2(D ~ fm$crop, permutations = 999)
  L <- c(L, sprintf("| %s ~ crop | %.4f | %.2f | %s |", nm, b$R2[1], b$F[1],
                    ifelse(b$`Pr(>F)`[1] <= 0.001, "<0.001", sprintf("%.3f", b$`Pr(>F)`[1]))))
}
# ---- study level
agg <- aggregate(Rr, by = list(project = fm$project, crop = fm$crop), FUN = mean)
Pa <- agg[, intersect(colnames(agg), colnames(Rr)), drop = FALSE]
Pa <- sweep(Pa, 1, rowSums(Pa), "/")
cov <- read.csv(file.path(outv, "study_covariates.csv"), stringsAsFactors = FALSE)
cov <- cov[match(agg$project, cov$project), ]
D_a_bc <- vegdist(Pa, method = "bray"); D_a_ai <- dist(zf(Pa))
X_crop <- model.matrix(~ crop, data = cov)[, -1, drop = FALSE]
X_tech <- model.matrix(~ log(depth_median) + asv_len_mode + log1p(n_asv_median) + n_samples +
                         unassigned_pct + genus_pct, data = cov)[, -1, drop = FALSE]
X_tech2 <- model.matrix(~ log(depth_median) + asv_len_mode + log1p(n_asv_median) + n_samples,
                        data = cov)[, -1, drop = FALSE]
L <- c(L, "", "## 2. Study level: crop alone vs technical covariates alone", "",
       "| distance | model | R2 | F | p |", "| --- | --- | --- | --- | --- |")
for (nm in c("Bray-Curtis","Aitchison")) {
  D <- if (nm == "Bray-Curtis") D_a_bc else D_a_ai
  for (tag in c("crop","technical (all)","technical (no resolution)")) {
    X <- if (tag == "crop") X_crop else if (tag == "technical (all)") X_tech else X_tech2
    a <- adonis2(D ~ X, permutations = 999)
    L <- c(L, sprintf("| %s | %s | %.4f | %.2f | %s |", nm, tag, a$R2[1], a$F[1],
                      ifelse(a$`Pr(>F)`[1] <= 0.001, "<0.001", sprintf("%.3f", a$`Pr(>F)`[1]))))
  }
}
L <- c(L, "", "## 3. Variance partitioning (db-RDA): crop [a] vs technical [c], shared [b]", "",
       "| distance | X2 = technical (all) | crop only [a] | shared [b] | technical only [c] | residual [d] |",
       "| --- | --- | --- | --- | --- | --- |")
for (nm in c("Bray-Curtis","Aitchison")) {
  D <- if (nm == "Bray-Curtis") D_a_bc else D_a_ai
  vp <- varpart(D, X_crop, X_tech)
  fr <- vp$part$indfract
  L <- c(L, sprintf("| %s | yes | %.3f | %.3f | %.3f | %.3f |", nm,
                    fr$Adj.R.square[1], fr$Adj.R.square[2], fr$Adj.R.square[3], fr$Adj.R.square[4]))
}
for (nm in c("Bray-Curtis","Aitchison")) {
  D <- if (nm == "Bray-Curtis") D_a_bc else D_a_ai
  vp <- varpart(D, X_crop, X_tech2)
  fr <- vp$part$indfract
  L <- c(L, sprintf("| %s | no resolution | %.3f | %.3f | %.3f | %.3f |", nm,
                    fr$Adj.R.square[1], fr$Adj.R.square[2], fr$Adj.R.square[3], fr$Adj.R.square[4]))
}
L <- c(L, "", "> Read-out: [a] = crop variance NOT shared with the technical covariates;",
       "> [b] = variance crop and technique explain jointly (inseparable in this design);",
       "> [c] = technical-only.  If [b] dominates, the crop effect cannot be separated from method/resolution.",
       "> The sample-level row `~ project` shows how much variation is between studies at all.")
writeLines(L, file.path(outv, "variance_decomposition.md"))
cat(paste(L, collapse = "\n"), "\n")
