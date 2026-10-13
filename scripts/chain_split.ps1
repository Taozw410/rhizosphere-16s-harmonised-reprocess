$root = 'D:\CodexWork\mic-re'
$log  = "$root\out50k\stage2_split.log"
"split chain start $(Get-Date -Format 'HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8
# A) 短读批次（raw50k = 批次1 剩余15 + 批次2 共17 = 32 个项目；已完成项目会 SKIP）
Add-Content $log "--- A short-read batch start $(Get-Date -Format 'HH:mm:ss') ---" -Encoding utf8
conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript "$root\scripts\run_dada2_parallel.R" ALL "$root\data\raw50k" "$root\out_stage2" 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz' 3 1e7 6 "$root\meta\primer_table_stage2.csv" *>> $log
Add-Content $log "--- A done $(Get-Date -Format 'HH:mm:ss') ---" -Encoding utf8
# B) 长读单端批次
Add-Content $log "--- B long-read batch start $(Get-Date -Format 'HH:mm:ss') ---" -Encoding utf8
conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript "$root\scripts\run_dada2_parallel.R" ALL "$root\data\raw50k_longread" "$root\out_stage2" 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz' 3 1e7 6 "$root\meta\primer_table_stage2.csv" *>> $log
Add-Content $log "--- B done $(Get-Date -Format 'HH:mm:ss') ---" -Encoding utf8
# C) 报告
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\summarize_w1.py" "$root\out_stage2" *>> $log
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\classify_w1.py" "$root\out_stage2" *>> $log
Add-Content $log "all done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
