z <- "D:/CodexWork/mic-re/out_analysis/pkgzips/ANCOMBC_2.12.1.zip"
install.packages(z, repos = NULL, type = "binary")
cat("--- load test ---\n")
r <- try(library(ANCOMBC), silent = TRUE)
if (inherits(r, "try-error")) cat("LOAD ERROR:", substr(as.character(attr(r, "condition")$message), 1, 200), "\n")
cat("ANCOMBC loads:", "ANCOMBC" %in% loadedNamespaces(), " version:", as.character(packageVersion("ANCOMBC")), "\n")
cat("--- final verification of the DA toolchain ---\n")
for (p in c("ALDEx2","ANCOMBC","zCompositions","iNEXT","vegan","phyloseq","randomForest"))
  cat(sprintf("%-14s %s %s\n", p, requireNamespace(p, quietly=TRUE),
              if (requireNamespace(p, quietly=TRUE)) as.character(packageVersion(p)) else ""))
