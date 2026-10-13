$root = 'D:\CodexWork\mic-re'
$log  = "$root\out50k\batch8_chain.log"
"batch8 chain start $(Get-Date -Format 'HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8
$pwsh = 'C:\Users\hthome\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\powershell\pwsh.exe'

while ($true) {
  Start-Sleep -Seconds 180
  $a = (Test-Path "$root\out50k\batch7_chain.log") -and ((Get-Content "$root\out50k\batch7_chain.log" -Raw) -match 'ALL DONE')
  $b = (Test-Path "$root\out50k\batch8_fetch.log") -and ((Get-Content "$root\out50k\batch8_fetch.log" -Raw) -match 'BATCH8 FETCH DONE')
  if ($a -and $b) { break }
}
Add-Content $log "batch7 finished + batch8 data ready $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

$projs = Get-Content "$root\meta\batch8_selected.txt" -Encoding UTF8 | Where-Object { $_.Trim() } | ForEach-Object { $_.Trim() }
$info = foreach ($p in $projs) {
  $f = @(Get-ChildItem "$root\data\raw50k\$p" -File -ErrorAction SilentlyContinue)
  [pscustomobject]@{ project = $p; n = $f.Count }
}
$A = @(); $B = @(); $sA = 0; $sB = 0
foreach ($x in ($info | Sort-Object -Property n -Descending)) {
  if ($sA -le $sB) { $A += $x.project; $sA += $x.n } else { $B += $x.project; $sB += $x.n }
}
Add-Content $log "b8 stream A=$($A -join ',') ($sA)" -Encoding utf8
Add-Content $log "b8 stream B=$($B -join ',') ($sB)" -Encoding utf8

$pA = Start-Process -FilePath $pwsh -WindowStyle Hidden -PassThru -ArgumentList @('-NoProfile','-Command',
  "conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\run_dada2_parallel.R' '$($A -join ',')' '$root\data\raw50k' '$root\out_stage2' 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz' 2 1e7 6 '$root\meta\primer_table_batch8.csv' *>> '$root\out50k\batch8_streamA.log'")
Start-Sleep -Seconds 15
$pB = Start-Process -FilePath $pwsh -WindowStyle Hidden -PassThru -ArgumentList @('-NoProfile','-Command',
  "conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\run_dada2_parallel.R' '$($B -join ',')' '$root\data\raw50k' '$root\out_stage2' 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz' 2 1e7 6 '$root\meta\primer_table_batch8.csv' *>> '$root\out50k\batch8_streamB.log'")
$pA.WaitForExit(); $pB.WaitForExit()
Add-Content $log "b8 streams done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\master_scan.py" *>> $log
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\marker_gene_qc.py" *>> $log
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\summarize_w1.py" "$root\out_stage2" *>> $log
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\classify_w1.py" "$root\out_stage2" *>> $log
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\merged_report.py" *>> $log
Add-Content $log "ALL DONE $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
