# 等第一批（7 个项目）跑完，自动启动第二批（3 个补充项目）
$log = 'D:\CodexWork\mic-re\out50k\batch2.log'
"watcher start $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8
while (@(Get-Process Rscript -ErrorAction SilentlyContinue).Count -gt 0) { Start-Sleep -Seconds 60 }
Add-Content -LiteralPath $log -Value "batch1 done $(Get-Date -Format 'HH:mm:ss') -> launching batch2" -Encoding utf8
conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript 'D:\CodexWork\mic-re\scripts\run_dada2_parallel.R' 'PRJNA430104,PRJNA677863,PRJNA1479906' 'D:\CodexWork\mic-re\data\raw50k' 'D:\CodexWork\mic-re\out50k' 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz' 3 1e7 6 'D:/CodexWork/mic-re/meta/primer_table.csv' *>> $log
Add-Content -LiteralPath $log -Value "batch2 done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
