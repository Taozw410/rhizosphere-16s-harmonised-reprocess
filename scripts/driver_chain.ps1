# 单一驱动链：batch7 -> batch8 -> A 重跑，每阶段两条流（各 2 核），阶段间跑闸门。
# 长读组（1261944/607209/1157226/509515）本轮暂停：其探针特征与已排除项目同型、单项目 10 小时量级，
# 数据与中间缓存保留在磁盘，可随时用同一条命令恢复。
$root = 'D:\CodexWork\mic-re'
$log  = "$root\out50k\driver_chain.log"
$py   = 'D:\bioenv\envs\riz16s\python.exe'
$pwsh = 'C:\Users\hthome\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\powershell\pwsh.exe'
$ref  = 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz'

function Launch-Stream {
  param([string]$projs, [string]$table, [string]$streamlog)
  return Start-Process -FilePath $pwsh -WindowStyle Hidden -PassThru -ArgumentList @('-NoProfile','-Command',
    "conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\run_dada2_parallel.R' '$projs' '$root\data\raw50k' '$root\out_stage2' '$ref' 2 1e7 6 '$table' *>> '$streamlog'")
}

"driver chain start $(Get-Date -Format 'HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8

# ---------- Phase 1: batch7 (20 projects) ----------
$b7 = @(Get-Content "$root\meta\batch7_selected.txt" -Encoding UTF8 | Where-Object { $_.Trim() } | ForEach-Object { $_.Trim() })
$info = foreach ($p in $b7) { [pscustomobject]@{ project = $p; n = @(Get-ChildItem "$root\data\raw50k\$p" -File -ErrorAction SilentlyContinue).Count } }
$sorted = @($info | Sort-Object -Property n -Descending)
$half = [int][math]::Ceiling($sorted.Count / 2)
$A = @($sorted[0..($half-1)] | Select-Object -ExpandProperty project)
$B = @($sorted[$half..($sorted.Count-1)] | Select-Object -ExpandProperty project)
Add-Content -LiteralPath $log -Value "PHASE1 batch7  A=$($A -join '|')  B=$($B -join '|')" -Encoding utf8
$pA = Launch-Stream ($A -join ',') "$root\meta\primer_table_batch7.csv" "$root\out50k\b7_streamA.log"
Start-Sleep -Seconds 15
$pB = Launch-Stream ($B -join ',') "$root\meta\primer_table_batch7.csv" "$root\out50k\b7_streamB.log"
$pA.WaitForExit(); $pB.WaitForExit()
Add-Content -LiteralPath $log -Value "PHASE1 streams done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
& $py "$root\scripts\master_scan.py" *>> $log
& $py "$root\scripts\marker_gene_qc.py" *>> $log
& $py "$root\scripts\summarize_w1.py" "$root\out_stage2" *>> $log
& $py "$root\scripts\classify_w1.py" "$root\out_stage2" *>> $log
& $py "$root\scripts\merged_report.py" *>> $log
Add-Content -LiteralPath $log -Value "PHASE1 gates done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

# ---------- Phase 2: batch8 (6 projects) ----------
$b8 = @(Get-Content "$root\meta\batch8_selected.txt" -Encoding UTF8 | Where-Object { $_.Trim() } | ForEach-Object { $_.Trim() })
$info = foreach ($p in $b8) { [pscustomobject]@{ project = $p; n = @(Get-ChildItem "$root\data\raw50k\$p" -File -ErrorAction SilentlyContinue).Count } }
$sorted = @($info | Sort-Object -Property n -Descending)
$half = [int][math]::Ceiling($sorted.Count / 2)
$C = @($sorted[0..($half-1)] | Select-Object -ExpandProperty project)
$D = @($sorted[$half..($sorted.Count-1)] | Select-Object -ExpandProperty project)
Add-Content -LiteralPath $log -Value "PHASE2 batch8  A=$($C -join '|')  B=$($D -join '|')" -Encoding utf8
$pC = Launch-Stream ($C -join ',') "$root\meta\primer_table_batch8.csv" "$root\out50k\b8_streamA.log"
Start-Sleep -Seconds 15
$pD = Launch-Stream ($D -join ',') "$root\meta\primer_table_batch8.csv" "$root\out50k\b8_streamB.log"
$pC.WaitForExit(); $pD.WaitForExit()
Add-Content -LiteralPath $log -Value "PHASE2 streams done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
& $py "$root\scripts\master_scan.py" *>> $log
& $py "$root\scripts\marker_gene_qc.py" *>> $log
& $py "$root\scripts\summarize_w1.py" "$root\out_stage2" *>> $log
& $py "$root\scripts\classify_w1.py" "$root\out_stage2" *>> $log
& $py "$root\scripts\merged_report.py" *>> $log
Add-Content -LiteralPath $log -Value "PHASE2 gates done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

# ---------- Phase 3: strategy-A reruns (11 projects) ----------
$projs = @('PRJNA996548','PRJNA596848','PRJNA665860','PRJNA753114','PRJNA1334676','PRJNA816337','PRJNA1189047','PRJNA548437','PRJNA565000','PRJNA781096','PRJNA887236')
$beforeRoot = "$root\out50k\_archive_trunclen_before"
New-Item -ItemType Directory -Force -Path $beforeRoot | Out-Null
foreach ($p in $projs) {
  $d = "$root\out_stage2\$p"
  if (Test-Path $d) {
    $arch = Join-Path $beforeRoot $p
    if (Test-Path $arch) { $arch = "${arch}_$(Get-Date -Format 'HHmmss')" }
    Move-Item -LiteralPath $d -Destination $arch -Force
  }
}
& $py "$root\scripts\apply_trunclen_plan.py" "$root\meta\trunclen_plan_rerun11.csv" "$root\meta\primer_table_stage2.csv" "$root\meta\primer_pairs_20261008.csv" *>> $log
$info = foreach ($p in $projs) { [pscustomobject]@{ project = $p; n = @(Get-ChildItem "$root\data\raw50k\$p" -File -ErrorAction SilentlyContinue).Count } }
$sorted = @($info | Sort-Object -Property n -Descending)
$half = [int][math]::Ceiling($sorted.Count / 2)
$E = @($sorted[0..($half-1)] | Select-Object -ExpandProperty project)
$F = @($sorted[$half..($sorted.Count-1)] | Select-Object -ExpandProperty project)
Add-Content -LiteralPath $log -Value "PHASE3 reruns  A=$($E -join '|')  B=$($F -join '|')" -Encoding utf8
$pE = Launch-Stream ($E -join ',') "$root\meta\primer_table_stage2.csv" "$root\out50k\rerun_streamA.log"
Start-Sleep -Seconds 15
$pF = Launch-Stream ($F -join ',') "$root\meta\primer_table_stage2.csv" "$root\out50k\rerun_streamB.log"
$pE.WaitForExit(); $pF.WaitForExit()
Add-Content -LiteralPath $log -Value "PHASE3 streams done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
& $py "$root\scripts\master_scan.py" *>> $log
& $py "$root\scripts\marker_gene_qc.py" *>> $log
& $py "$root\scripts\summarize_w1.py" "$root\out_stage2" *>> $log
& $py "$root\scripts\classify_w1.py" "$root\out_stage2" *>> $log
& $py "$root\scripts\merged_report.py" *>> $log
& $py "$root\scripts\trunclen_ab_report.py" "$beforeRoot" "$root\out_stage2" ($projs -join ',') --out "$root\out50k\trunclen_ab_report.md" *>> $log
& $py "$root\scripts\select_final45.py" --out "$root\meta\final45_selection.csv" --report "$root\out50k\选集报告.md" *>> $log
Add-Content -LiteralPath $log -Value "PIPELINE ALL DONE $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
