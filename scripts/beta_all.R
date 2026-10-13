#!/usr/bin/env Rscript
# Unified beta-diversity analysis over FOUR distances (section 11): Bray-Curtis, Aitchison,
# unweighted UniFrac, weighted UniFrac.  For each: adonis2 at sample and project level (999 perms),
# betadisper + permutest (dispersion check), leave-one-project-out sensitivity (99 perms).
# UniFrac matrices come from unifrac_features.py; samples are aligned by name.
# ASCII-only source.
suppressPackageStartupMessages(library(vegan))
outv <- "D:/CodexWork/mic-re/out_analysis"
fm <- read.csv(file.path(outv, "feature_matrix.csv"), check.names = FALSE, stringsAsFactors = FALSE)
meta_cols <- c("project","base","crop","sample","all_reads","bacteria_reads","organelle_reads",
               "genus_reads","family_reads","unassigned_reads","n_asv")
feats <- setdiff(colnames(fm), meta_cols)
fm <- fm[fm$bacteria_reads >= 1000, ]
nm <- paste(fm$project, fm$sample, sep = "|")
X <- as.matrix(fm[, feats]); rownames(X) <- nm; X <- X[, colSums(X) > 0, drop = FALSE]
P <- sweep(X, 1, rowSums(X), "/")
zf <- function(Pm) { q <- Pm
  for (i in seq_len(nrow(q))) { z <- q[i, ] == 0
    if (any(z)) { mn <- min(q[i, q[i, ] > 0]); q[i, z] <- 0.65 * mn; q[i, ] <- q[i, ] / sum(q[i, ]) } }
  log(q) - rowMeans(log(q)) }
D <- list()
D[["Bray-Curtis"]] <- as.matrix(vegdist(P, method = "bray"))
D[["Aitchison"]] <- as.matrix(dist(zf(P)))
for (nmf in c("unweighted UniFrac", "weighted UniFrac")) {
  f <- file.path(outv, if (nmf == "unweighted UniFrac") "dist_unifrac_uw.csv" else "dist_unifrac_w.csv")
  if (file.exists(f)) {
    m <- as.matrix(read.csv(f, row.names = 1, check.names = FALSE))
    m <- m[rownames(m) %in% nm, colnames(m) %in% nm]
    D[[nmf]] <- m[nm, nm]
  } else cat("missing", f, "\n")
}
crop <- factor(fm$crop, levels = c("maize","rice","wheat")); proj <- fm$project
agg <- aggregate(P, by = list(project = fm$project, crop = as.character(fm$crop)), FUN = mean)
Pa <- sweep(agg[, intersect(colnames(agg), colnames(P)), drop = FALSE], 1,
            rowSums(agg[, intersect(colnames(agg), colnames(P)), drop = FALSE]), "/")
Da <- list()
Da[["Bray-Curtis"]] <- as.matrix(vegdist(Pa, method = "bray"))
Da[["Aitchison"]] <- as.matrix(dist(zf(Pa)))
for (nmf in c("unweighted UniFrac", "weighted UniFrac")) {
  if (is.null(D[[nmf]])) next
  # project-level UniFrac: mean composition per project -> recompute from representative distance?
  # not valid; instead average the sample-level distances per project pair
  M <- D[[nmf]]; labs <- proj
  up <- unique(proj); Pm <- matrix(0, length(up), length(up), dimnames = list(up, up))
  for (i in seq_along(up)) for (j in seq_along(up)) {
    ii <- labs == up[i]; jj <- labs == up[j]
    Pm[i, j] <- mean(M[ii, jj])
  }
  Da[[nmf]] <- Pm
}
say <- function(...) cat(sprintf(...), "\n")
L <- c("# Beta diversity over four distances (final 45 projects; pre-registration section 11)", "",
       sprintf("Samples: %d (depth >= 1000 bacteria reads). Features: %d. PERMANOVA = adonis2.", nrow(fm), ncol(P)),
       "UniFrac computed in the section-17 feature space on an NJ tree of feature representatives",
       "(midpoint-rooted in Python; 43 negative NJ branch lengths clamped to 0).", "",
       "## 1. PERMANOVA", "", "| distance | level | R2 | F | p |", "| --- | --- | --- | --- | --- |")
for (k in names(D)) {
  a <- adonis2(as.dist(D[[k]]) ~ crop, data = fm, permutations = 999)
  L <- c(L, sprintf("| %s | sample | %.4f | %.2f | %s |", k, a$R2[1], a$F[1],
                    ifelse(a$`Pr(>F)`[1] <= 0.001, "<0.001", sprintf("%.3f", a$`Pr(>F)`[1]))))
  a2 <- adonis2(as.dist(Da[[k]]) ~ crop, data = data.frame(crop = agg$crop), permutations = 999)
  L <- c(L, sprintf("| %s | project | %.4f | %.2f | %s |", k, a2$R2[1], a2$F[1],
                    ifelse(a2$`Pr(>F)`[1] <= 0.001, "<0.001", sprintf("%.3f", a2$`Pr(>F)`[1]))))
}
L <- c(L, "", "## 2. Dispersion check (betadisper, permutest 999)", "",
       "| distance | F | p | mean distance to centroid (maize/rice/wheat) |", "| --- | --- | --- | --- |")
for (k in names(D)) {
  bd <- betadisper(as.dist(D[[k]]), crop); pt <- permutest(bd, permutations = 999)
  pv <- pt$tab$`Pr(>F)`[which(rownames(pt$tab) == "Groups")]
  fv <- pt$tab$F[which(rownames(pt$tab) == "Groups")]
  L <- c(L, sprintf("| %s | %.3f | %.3f | %.4f / %.4f / %.4f |", k, fv, pv,
                    bd$group.distances[1], bd$group.distances[2], bd$group.distances[3]))
}
L <- c(L, "", "## 3. Leave-one-project-out (sample level, 99 perms)", "",
       "| distance | R2 range | runs p<0.05 | runs |", "| --- | --- | --- | --- |")
for (k in names(D)) {
  M <- D[[k]]; r2s <- c(); hits <- 0; n <- 0
  for (p in unique(proj)) {
    idx <- proj != p
    if (length(unique(crop[idx])) < 3) next
    aa <- adonis2(as.dist(M[idx, idx]) ~ crop[idx], permutations = 99)
    r2s <- c(r2s, aa$R2[1]); hits <- hits + (aa$`Pr(>F)`[1] <= 0.05); n <- n + 1
  }
  L <- c(L, sprintf("| %s | %.4f - %.4f | %d / %d | %d |", k, min(r2s), max(r2s), hits, n, n))
}
writeLines(L, file.path(outv, "beta_all_report.md"))
cat(paste(L, collapse = "\n"), "\n")
