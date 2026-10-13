$root = 'D:\CodexWork\mic-re'
$log  = "$root\out50k\batch2_chain.log"
"chain start $(Get-Date -Format 'HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8
# 1) 等下载
while (@(Get-Process curl -ErrorAction SilentlyContinue).Count -gt 0) { Start-Sleep -Seconds 60 }
Add-Content $log "download done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\verify_downloads.py" *>> $log
# 2) 抽样
$projs = (Get-Content "$root\meta\batch2_selected_projects.txt") -join ','
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\subsample_fastq.py" "$root\data\raw" "$root\data\raw50k" 50000 20261006 $projs *>> $log
Add-Content $log "subsample done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
# 3) 等批次1管线结束（无 Rscript）
while (@(Get-Process Rscript -ErrorAction SilentlyContinue).Count -gt 0) { Start-Sleep -Seconds 60 }
Add-Content $log "batch1 pipeline done; starting batch2 pipeline $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript "$root\scripts\run_dada2_parallel.R" ALL "$root\data\raw50k" "$root\out_stage2" 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz' 3 1e7 6 "$root\meta\primer_table_batch2.csv" *>> $log
Add-Content $log "pipeline done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\summarize_w1.py" "$root\out_stage2" *>> $log
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\classify_w1.py" "$root\out_stage2" *>> $log
Add-Content $log "all done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
