options(repos = c(CRAN = "https://cloud.r-project.org"), download.file.method = "libcurl", timeout = 1800)
say <- function(...) cat(sprintf(...), "\n")
say("install start: %s", format(Sys.time()))
# bypass BiocManager's config.yaml fetch (it is flaky through the local proxy): probe BioC versions directly
bio <- NULL
for (v in c("3.22", "3.21", "3.20")) {
  u <- sprintf("https://bioconductor.org/packages/%s/bioc/src/contrib/PACKAGES.gz", v)
  ok <- !inherits(try(download.file(u, tempfile(), quiet = TRUE), silent = TRUE), "try-error")
  say("bioc %s reachable: %s", v, ok)
  if (ok) { bio <- v; break }
}
if (is.null(bio)) { say("no Bioconductor release reachable -> abort"); quit(status = 1) }
repos <- c(CRAN = "https://cloud.r-project.org",
           BioCsoft = sprintf("https://bioconductor.org/packages/%s/bioc", bio),
           BioCann  = sprintf("https://bioconductor.org/packages/%s/data/annotation", bio),
           BioCexp  = sprintf("https://bioconductor.org/packages/%s/data/experiment", bio))
pkgs <- c("zCompositions", "ALDEx2", "ANCOMBC", "iNEXT")
for (attempt in 1:6) {
  say("--- attempt %d: %s", attempt, paste(pkgs, collapse = ", "))
  install.packages(pkgs, repos = repos, Ncpus = 2)
  left <- pkgs[!vapply(pkgs, function(p) requireNamespace(p, quietly = TRUE), TRUE)]
  if (!length(left)) break
  say("remaining after attempt %d: %s", attempt, paste(left, collapse = ", "))
  pkgs <- left
  Sys.sleep(15)
}
say("--- verify ---")
for (p in c("zCompositions", "ALDEx2", "ANCOMBC", "iNEXT"))
  say("%-14s %s %s", p, requireNamespace(p, quietly = TRUE),
      if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "")
say("install end: %s", format(Sys.time()))
