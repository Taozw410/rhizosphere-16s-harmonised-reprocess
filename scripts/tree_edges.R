#!/usr/bin/env Rscript
# Export the feature tree as an edge table (+tip labels) so the vectorised UniFrac in Python can use it.
# Roots the tree at the midpoint (the feature tree has no natural outgroup). ASCII-only source.
suppressPackageStartupMessages(library(ape))
outv <- "D:/CodexWork/mic-re/out_analysis"
tr <- read.tree(file.path(outv, "feature_tree.nwk"))
if (!is.rooted(tr)) { tr <- midpoint(tr); cat("midpoint-rooted\n") }
cat(sprintf("tips %d ; nodes %d ; rooted %s\n", length(tr$tip.label), tr$Nnode, is.rooted(tr)))
df <- data.frame(edge_id = seq_len(nrow(tr$edge)), parent = tr$edge[, 1], child = tr$edge[, 2],
                 length = tr$edge.length,
                 tip_label = ifelse(tr$edge[, 2] <= length(tr$tip.label), tr$tip.label[tr$edge[, 2]], ""),
                 stringsAsFactors = FALSE)
write.csv(df, file.path(outv, "tree_edges.csv"), row.names = FALSE)
cat("wrote out_analysis/tree_edges.csv\n")
