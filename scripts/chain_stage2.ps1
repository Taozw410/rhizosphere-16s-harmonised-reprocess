$root = 'D:\CodexWork\mic-re'
$log  = "$root\out50k\stage2_chain.log"
"chain start $(Get-Date -Format 'HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8

# 1) 等下载结束（无 curl 进程）
while (@(Get-Process curl -ErrorAction SilentlyContinue).Count -gt 0) { Start-Sleep -Seconds 60 }
Add-Content $log "download finished $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

# 2) 完整性校验
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\verify_downloads.py" *>> $log
Add-Content $log "verify done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

# 3) 50k 抽样（仅批次1的 22 个项目）
$projs = (Get-Content "$root\meta\stage2_batch1_projects.txt") -join ','
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\subsample_fastq.py" "$root\data\raw" "$root\data\raw50k" 50000 20261006 $projs *>> $log
Add-Content $log "subsample done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

# 4) 跑管线（独立输出目录 out_stage2）
New-Item -ItemType Directory -Force -Path "$root\out_stage2" | Out-Null
Copy-Item "$root\out50k\_shared\errF.rds" "$root\out_stage2\errF.rds" -Force -ErrorAction SilentlyContinue
New-Item -ItemType Directory -Force -Path "$root\out_stage2\_shared" | Out-Null
Copy-Item "$root\out50k\_shared\errF.rds" "$root\out_stage2\_shared\errF.rds" -Force -ErrorAction SilentlyContinue
Copy-Item "$root\out50k\_shared\errR.rds" "$root\out_stage2\_shared\errR.rds" -Force -ErrorAction SilentlyContinue
Add-Content $log "pipeline start $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript "$root\scripts\run_dada2_parallel.R" ALL "$root\data\raw50k" "$root\out_stage2" 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz' 3 1e7 6 "$root\meta\primer_table_stage2.csv" *>> $log
Add-Content $log "pipeline done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

# 5) 报告
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\summarize_w1.py" "$root\out_stage2" *>> $log
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\classify_w1.py" "$root\out_stage2" *>> $log
Add-Content $log "all done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
