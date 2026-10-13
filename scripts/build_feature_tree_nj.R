#!/usr/bin/env Rscript
# Fast neighbour-joining tree on the DECIPHER distance matrix (standard practice for UniFrac where an
# ML search over ~2000 tips is impractical). Also exports the edge table for the vectorised UniFrac.
# ASCII-only source.
suppressPackageStartupMessages(library(DECIPHER)); suppressPackageStartupMessages(library(ape)); suppressPackageStartupMessages(library(Biostrings))
outv <- "D:/CodexWork/mic-re/out_analysis"
say <- function(...) cat(sprintf(...), "\n")
al <- readDNAStringSet(file.path(outv, "feature_aln.fa"))
d <- DistanceMatrix(al, type = "dist", verbose = FALSE)
say("distance matrix done: %d x %d", attr(d, "Size"), attr(d, "Size"))
d[is.na(d)] <- max(d, na.rm = TRUE)   # overlap method leaves NAs; fill with the max observed distance
tr <- nj(d)
# rooting is done on the edge table in Python (ape::midpoint belongs to phangorn, not installed)
say("NJ tree: tips %d ; rooted %s", length(tr$tip.label), is.rooted(tr))
if (any(tr$edge.length < 0, na.rm = TRUE)) {
  cat(sprintf("clamping %d negative branch lengths to 0 (NJ after filling missing distances)\n", sum(tr$edge.length < 0, na.rm = TRUE)))
  tr$edge.length[tr$edge.length < 0] <- 0
}
tr$edge.length[is.na(tr$edge.length)] <- 0
write.tree(tr, file.path(outv, "feature_tree_nj.nwk"))
df <- data.frame(edge_id = seq_len(nrow(tr$edge)), parent = tr$edge[, 1], child = tr$edge[, 2],
                 length = tr$edge.length,
                 tip_label = ifelse(tr$edge[, 2] <= length(tr$tip.label), tr$tip.label[tr$edge[, 2]], ""),
                 stringsAsFactors = FALSE)
write.csv(df, file.path(outv, "tree_edges.csv"), row.names = FALSE)
say("wrote feature_tree_nj.nwk and tree_edges.csv")



