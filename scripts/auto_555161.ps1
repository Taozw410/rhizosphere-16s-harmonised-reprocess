$log = 'D:\CodexWork\mic-re\out50k\rerun_555161.log'
"chain start $(Get-Date -Format 'HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8
while (@(Get-Process Rscript -ErrorAction SilentlyContinue).Count -gt 0) { Start-Sleep -Seconds 60 }
Add-Content -LiteralPath $log -Value "previous run done $(Get-Date -Format 'HH:mm:ss'); archiving old 555161" -Encoding utf8
$o = 'D:\CodexWork\mic-re\out50k'
$arch = Join-Path $o '_archive_555161_trunclen250'
if (Test-Path (Join-Path $o 'PRJNA555161')) {
  New-Item -ItemType Directory -Force -Path $arch | Out-Null
  Move-Item -LiteralPath (Join-Path $o 'PRJNA555161') -Destination (Join-Path $arch 'PRJNA555161') -Force
}
Add-Content -LiteralPath $log -Value "relaunching PRJNA555161 truncLen 280/250 $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript 'D:\CodexWork\mic-re\scripts\run_dada2_parallel.R' 'PRJNA555161' 'D:\CodexWork\mic-re\data\raw50k' $o 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz' 3 1e7 6 'D:/CodexWork/mic-re/meta/primer_table.csv' *>> $log
Add-Content -LiteralPath $log -Value "555161 done $(Get-Date -Format 'HH:mm:ss'); regen summary+tiers" -Encoding utf8
& 'D:\bioenv\envs\riz16s\python.exe' 'D:\CodexWork\mic-re\scripts\summarize_w1.py' *>> 'D:\CodexWork\mic-re\out50k\summary.log'
& 'D:\bioenv\envs\riz16s\python.exe' 'D:\CodexWork\mic-re\scripts\classify_w1.py' *>> 'D:\CodexWork\mic-re\out50k\summary.log'
Add-Content -LiteralPath $log -Value "all done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
