# batch7 chain v2: start ONE stream immediately (pending long-read stream A is already finished, so the
# "max two streams" rule still holds), then start the second stream once the pending chain reports ALL DONE.
$root = 'D:\CodexWork\mic-re'
$log  = "$root\out50k\batch7_chain.log"
$py   = 'D:\bioenv\envs\riz16s\python.exe'
$pwsh = 'C:\Users\hthome\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\powershell\pwsh.exe'
$ref  = 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz'
"batch7 chain v2 start $(Get-Date -Format 'HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8

# primer table (already prepared with per-project truncLen + primers)
Add-Content -LiteralPath $log -Value "primer table ready $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

$projs = @(Get-Content "$root\meta\batch7_selected.txt" -Encoding UTF8 | Where-Object { $_.Trim() } | ForEach-Object { $_.Trim() })
# order by sample count descending, split into two halves
$info = foreach ($p in $projs) {
  $f = @(Get-ChildItem "$root\data\raw50k\$p" -File -ErrorAction SilentlyContinue)
  [pscustomobject]@{ project = $p; n = $f.Count }
}
$sorted = @($info | Sort-Object -Property n -Descending)
$half = [int][math]::Ceiling($sorted.Count / 2)
$A = @($sorted[0..($half-1)] | Select-Object -ExpandProperty project)
$B = @($sorted[$half..($sorted.Count-1)] | Select-Object -ExpandProperty project)
Add-Content -LiteralPath $log -Value "stream A (start now) = $($A -join ',')" -Encoding utf8
Add-Content -LiteralPath $log -Value "stream B (after pending ALL DONE) = $($B -join ',')" -Encoding utf8

$pA = Start-Process -FilePath $pwsh -WindowStyle Hidden -PassThru -ArgumentList @('-NoProfile','-Command',
  "conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\run_dada2_parallel.R' '$($A -join ',')' '$root\data\raw50k' '$root\out_stage2' '$ref' 2 1e7 6 '$root\meta\primer_table_batch7.csv' *>> '$root\out50k\batch7_streamA.log'")
Add-Content -LiteralPath $log -Value "stream A launched $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

while ($true) {
  Start-Sleep -Seconds 180
  $done = (Test-Path "$root\out50k\pend_chain.log") -and ((Get-Content "$root\out50k\pend_chain.log" -Raw -ErrorAction SilentlyContinue) -match 'ALL DONE')
  if ($done) { break }
}
Add-Content -LiteralPath $log -Value "pending ALL DONE $(Get-Date -Format 'HH:mm:ss'); launching stream B" -Encoding utf8
$pB = Start-Process -FilePath $pwsh -WindowStyle Hidden -PassThru -ArgumentList @('-NoProfile','-Command',
  "conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\run_dada2_parallel.R' '$($B -join ',')' '$root\data\raw50k' '$root\out_stage2' '$ref' 2 1e7 6 '$root\meta\primer_table_batch7.csv' *>> '$root\out50k\batch7_streamB.log'")
$pA.WaitForExit(); $pB.WaitForExit()
Add-Content -LiteralPath $log -Value "b7 streams done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

& $py "$root\scripts\master_scan.py" *>> $log
& $py "$root\scripts\marker_gene_qc.py" *>> $log
& $py "$root\scripts\summarize_w1.py" "$root\out_stage2" *>> $log
& $py "$root\scripts\classify_w1.py" "$root\out_stage2" *>> $log
& $py "$root\scripts\merged_report.py" *>> $log
Add-Content -LiteralPath $log -Value "ALL DONE $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
