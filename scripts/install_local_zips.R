zdir <- "D:/CodexWork/mic-re/out_analysis/pkgzips"
say <- function(...) cat(sprintf(...), "\n")
zips <- list.files(zdir, pattern = "\\.zip$", full.names = TRUE)
say("local zips: %d", length(zips))
for (pass in 1:5) {
  before <- rownames(installed.packages())
  r <- try(install.packages(zips, repos = NULL, type = "binary"), silent = TRUE)
  if (inherits(r, "try-error")) say("pass %d raised: %s", pass, substr(as.character(r), 1, 100))
  after <- rownames(installed.packages())
  new <- setdiff(after, before)
  say("pass %d installed %d new packages", pass, length(new))
  if (!length(new)) break
}
say("--- verify ---")
for (p in c("zCompositions", "ALDEx2", "ANCOMBC", "iNEXT"))
  say("%-14s %s %s", p, requireNamespace(p, quietly = TRUE),
      if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "")
