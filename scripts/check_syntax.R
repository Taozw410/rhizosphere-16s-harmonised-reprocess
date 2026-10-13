p <- "D:/CodexWork/mic-re/scripts/run_dada2_parallel.R"
ok <- tryCatch({ parse(p); "PARSE OK" }, error = function(e) paste("PARSE ERROR:", conditionMessage(e)))
cat(ok, "\n")
