# Strategy-A rerun chain (11 merge-limited projects): wait for batch8, archive old results,
# set the planned truncLen, rerun in two streams, refresh gates/reports, write the A/B table.
$root = 'D:\CodexWork\mic-re'
$log  = "$root\out50k\trunclen2_chain.log"
$py   = 'D:\bioenv\envs\riz16s\python.exe'
$pwsh = 'C:\Users\hthome\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\powershell\pwsh.exe'
$ref  = 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz'
$projs = @('PRJNA996548','PRJNA596848','PRJNA665860','PRJNA753114','PRJNA1334676','PRJNA816337',
           'PRJNA1189047','PRJNA548437','PRJNA565000','PRJNA781096','PRJNA887236')
$stamp = '20261008'

"strategy-A rerun chain start $(Get-Date -Format 'HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8

while ($true) {
  Start-Sleep -Seconds 300
  $done = (Test-Path "$root\out50k\batch8_chain.log") -and ((Get-Content "$root\out50k\batch8_chain.log" -Raw -ErrorAction SilentlyContinue) -match 'ALL DONE')
  if (-not $done) { continue }
  $r = @(Get-CimInstance Win32_Process -Filter "Name='Rscript.exe'" -ErrorAction SilentlyContinue |
         Where-Object { $_.CommandLine -like '*run_dada2_parallel.R*' })
  if ($r.Count -eq 0) { break }
  Add-Content -LiteralPath $log -Value "waiting: $($r.Count) R process(es) $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
}
Add-Content -LiteralPath $log -Value "batch8 finished and R idle $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

# archive superseded results (never delete)
$beforeRoot = "$root\out50k\_archive_trunclen_before"
New-Item -ItemType Directory -Force -Path $beforeRoot | Out-Null
foreach ($p in $projs) {
  $d = "$root\out_stage2\$p"
  if (Test-Path $d) {
    $arch = Join-Path $beforeRoot $p
    if (Test-Path $arch) { $arch = "${arch}_$(Get-Date -Format 'HHmmss')" }
    Move-Item -LiteralPath $d -Destination $arch -Force
    Add-Content -LiteralPath $log -Value "archived $p" -Encoding utf8
  }
}

# update primer table from the plan (also clears primers so the planned truncLen stays valid)
& $py "$root\scripts\apply_trunclen_plan.py" "$root\meta\trunclen_plan_rerun11.csv" "$root\meta\primer_table_stage2.csv" *>> $log

# two balanced streams by sample count
$info = foreach ($p in $projs) {
  $f = @(Get-ChildItem "$root\data\raw50k\$p" -File -ErrorAction SilentlyContinue)
  [pscustomobject]@{ project = $p; n = $f.Count }
}
$A = @(); $B = @(); $sA = 0; $sB = 0
foreach ($x in ($info | Sort-Object -Property n -Descending)) {
  if ($sA -le $sB) { $A += $x.project; $sA += $x.n } else { $B += $x.project; $sB += $x.n }
}
Add-Content -LiteralPath $log -Value "streams A=$($A -join ',') ($sA)  B=$($B -join ',') ($sB)" -Encoding utf8

$pA = Start-Process -FilePath $pwsh -WindowStyle Hidden -PassThru -ArgumentList @('-NoProfile','-Command',
  "conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\run_dada2_parallel.R' '$($A -join ',')' '$root\data\raw50k' '$root\out_stage2' '$ref' 2 1e7 6 '$root\meta\primer_table_stage2.csv' *>> '$root\out50k\trunclen2_streamA.log'")
Start-Sleep -Seconds 15
$pB = Start-Process -FilePath $pwsh -WindowStyle Hidden -PassThru -ArgumentList @('-NoProfile','-Command',
  "conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\run_dada2_parallel.R' '$($B -join ',')' '$root\data\raw50k' '$root\out_stage2' '$ref' 2 1e7 6 '$root\meta\primer_table_stage2.csv' *>> '$root\out50k\trunclen2_streamB.log'")
Add-Content -LiteralPath $log -Value "streams launched $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
$pA.WaitForExit(); $pB.WaitForExit()
Add-Content -LiteralPath $log -Value "streams done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

& $py "$root\scripts\master_scan.py" *>> $log
& $py "$root\scripts\marker_gene_qc.py" *>> $log
& $py "$root\scripts\summarize_w1.py" "$root\out_stage2" *>> $log
& $py "$root\scripts\classify_w1.py" "$root\out_stage2" *>> $log
& $py "$root\scripts\merged_report.py" *>> $log
& $py "$root\scripts\trunclen_ab_report.py" "$beforeRoot" "$root\out_stage2" ($projs -join ',') --out "$root\out50k\trunclen_ab_report.md" *>> $log
Add-Content -LiteralPath $log -Value "ALL DONE $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
