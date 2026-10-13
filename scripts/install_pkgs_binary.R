options(repos = c(CRAN = "https://cloud.r-project.org",
                  BioCsoft = "https://bioconductor.org/packages/3.21/bioc",
                  BioCann  = "https://bioconductor.org/packages/3.21/data/annotation",
                  BioCexp  = "https://bioconductor.org/packages/3.21/data/experiment"),
        download.file.method = "libcurl", timeout = 3600, install.packages.check.source = "no")
say <- function(...) cat(sprintf(...), "\n")
want <- c("zCompositions", "ALDEx2", "ANCOMBC", "iNEXT")
say("binary retry start: %s", format(Sys.time()))
for (attempt in 1:10) {
  left <- want[!vapply(want, function(p) requireNamespace(p, quietly = TRUE), TRUE)]
  if (!length(left)) break
  say("--- attempt %d (binary): %s", attempt, paste(left, collapse = ", "))
  r <- try(install.packages(left, repos = getOption("repos"), type = "binary", Ncpus = 2), silent = TRUE)
  if (inherits(r, "try-error")) say("install.packages raised: %s", substr(as.character(r), 1, 120))
  Sys.sleep(20)
}
say("--- verify ---")
for (p in want) say("%-14s %s %s", p, requireNamespace(p, quietly = TRUE),
                    if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "")
say("binary retry end: %s", format(Sys.time()))
