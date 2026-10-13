# Phase 5c：PRJNA1476760 的单端抢救（R2 归档损坏、R1 完好）
$root = 'D:\CodexWork\mic-re'
$log  = "$root\out50k\phase5c_chain.log"
$py   = 'D:\bioenv\envs\riz16s\python.exe'
$pwsh = 'C:\Users\hthome\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\powershell\pwsh.exe'
$ref  = 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz'
"phase5c start $(Get-Date -Format 'HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8
while ($true) {
  Start-Sleep -Seconds 300
  $done = (Test-Path "$root\out50k\lr_resume_chain.log") -and ((Get-Content "$root\out50k\lr_resume_chain.log" -Raw -ErrorAction SilentlyContinue) -match 'LR RESUME ALL DONE')
  if ($done) { break }
}
Add-Content -LiteralPath $log -Value "lr resume finished $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
$before = "$root\out50k\_archive_single_end_before\PRJNA1476760"
New-Item -ItemType Directory -Force -Path (Split-Path $before) | Out-Null
if (Test-Path "$root\out_stage2\PRJNA1476760") { Move-Item -LiteralPath "$root\out_stage2\PRJNA1476760" -Destination $before -Force }
$p = Start-Process -FilePath $pwsh -WindowStyle Hidden -PassThru -ArgumentList @('-NoProfile','-Command',
  "conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\run_dada2_parallel.R' 'PRJNA1476760' '$root\data\raw50k_se' '$root\out_stage2' '$ref' 2 1e7 6 '$root\meta\primer_table_se_1476760.csv' *>> '$root\out50k\se_stream.log'")
Add-Content -LiteralPath $log -Value "single-end stream launched $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
$p.WaitForExit()
Add-Content -LiteralPath $log -Value "stream done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
& $py "$root\scripts\master_scan.py" *>> $log
& $py "$root\scripts\marker_gene_qc.py" *>> $log
& $py "$root\scripts\summarize_w1.py" "$root\out_stage2" *>> $log
& $py "$root\scripts\classify_w1.py" "$root\out_stage2" *>> $log
& $py "$root\scripts\merged_report.py" *>> $log
& $py "$root\scripts\taxonomy_diagnosis.py" *>> $log
& $py "$root\scripts\select_final45.py" --out "$root\meta\final45_selection.csv" --report "$root\out50k\选集报告.md" *>> $log
Add-Content -LiteralPath $log -Value "PHASE5C ALL DONE $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
