zips <- list.files("D:/CodexWork/mic-re/out_analysis/pkgzips", pattern = "\\.zip$", full.names = TRUE)
want <- c("microbiome","compositions","dplyr","Rtsne","tidyr","tensorA","robustbase","bayesm","purrr",
          "rhdf5","DEoptimR","Rhdf5lib","rhdf5filters")
sel <- zips[vapply(basename(zips), function(b) any(startsWith(b, paste0(want, "_"))), logical(1))]
cat("zips:", length(sel), "\n")
for (pass in 1:4) {
  before <- rownames(installed.packages())
  invisible(try(install.packages(sel, repos = NULL, type = "binary"), silent = TRUE))
  new <- setdiff(rownames(installed.packages()), before)
  cat(sprintf("pass %d: +%d (%s)\n", pass, length(new), paste(new, collapse = ", ")))
  if (!length(new)) break
}
cat("--- verify ---\n")
for (p in c("microbiome","compositions","rhdf5","Rhdf5lib","rhdf5filters","dplyr","tidyr","tensorA",
            "purrr","Rtsne","robustbase","bayesm","DEoptimR","ANCOMBC","ALDEx2"))
  cat(sprintf("%-14s %s %s\n", p, requireNamespace(p, quietly=TRUE),
              if (requireNamespace(p, quietly=TRUE)) as.character(packageVersion(p)) else ""))
