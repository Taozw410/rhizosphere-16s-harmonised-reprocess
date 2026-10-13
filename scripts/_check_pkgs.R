for (p in c("vegan","phyloseq","DECIPHER","randomForest","dplyr","ggplot2","iNEXT","ANCOMBC","ALDEx2","edgeR","DESeq2")) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("  %-14s %s\n", p, if (ok) "OK" else "MISSING"))
}
