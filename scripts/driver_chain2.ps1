# 续跑驱动链：接管 batch7 之后的全部阶段。Phase 1 的两条流已在运行，本链只等它们结束，不重复启动。
$root = 'D:\CodexWork\mic-re'
$log  = "$root\out50k\driver2_chain.log"
$py   = 'D:\bioenv\envs\riz16s\python.exe'
$pwsh = 'C:\Users\hthome\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\powershell\pwsh.exe'
$ref  = 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz'

function Launch-Stream {
  param([string]$projs, [string]$table, [string]$streamlog)
  Start-Process -FilePath $pwsh -WindowStyle Hidden -PassThru -ArgumentList @('-NoProfile','-Command',
    "conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\run_dada2_parallel.R' '$projs' '$root\data\raw50k' '$root\out_stage2' '$ref' 2 1e7 6 '$table' *>> '$streamlog'")
}
function Run-Gates {
  & $py "$root\scripts\master_scan.py" *>> $log
  & $py "$root\scripts\marker_gene_qc.py" *>> $log
  & $py "$root\scripts\summarize_w1.py" "$root\out_stage2" *>> $log
  & $py "$root\scripts\classify_w1.py" "$root\out_stage2" *>> $log
  & $py "$root\scripts\merged_report.py" *>> $log
}
function Split-Half {
  param([string[]]$Projects)
  $info = foreach ($p in $Projects) { [pscustomobject]@{ project = $p; n = @(Get-ChildItem "$root\data\raw50k\$p" -File -ErrorAction SilentlyContinue).Count } }
  $sorted = @($info | Sort-Object -Property n -Descending)
  $half = [int][math]::Ceiling($sorted.Count / 2)
  return @{ A = @($sorted[0..($half-1)] | Select-Object -ExpandProperty project); B = @($sorted[$half..($sorted.Count-1)] | Select-Object -ExpandProperty project) }
}

"driver2 start $(Get-Date -Format 'HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8

# ---- 0) 等 Phase 1 的两条流结束（它们已在运行）----
while ($true) {
  Start-Sleep -Seconds 120
  $a = (Test-Path "$root\out50k\b7_streamA.log") -and ((Get-Content "$root\out50k\b7_streamA.log" -Raw -ErrorAction SilentlyContinue) -match 'all done')
  $b = (Test-Path "$root\out50k\b7_streamB.log") -and ((Get-Content "$root\out50k\b7_streamB.log" -Raw -ErrorAction SilentlyContinue) -match 'all done')
  if ($a -and $b) { break }
}
Add-Content -LiteralPath $log -Value "PHASE1 streams finished $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
Run-Gates
Add-Content -LiteralPath $log -Value "PHASE1 gates done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

# ---- 1) 换入支持"单端去引物"的管线（此刻没有 R 在跑，文件未锁）----
Copy-Item -LiteralPath "$root\scripts\_bak\run_dada2_parallel.patched4.R" -Destination "$root\scripts\run_dada2_parallel.R" -Force
conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript "$root\scripts\_bak\parse_check4.R" *>> $log
Add-Content -LiteralPath $log -Value "pipeline swapped (single-end primer trimming enabled) $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

# ---- 2) Phase 2: batch8 ----
$b8 = @(Get-Content "$root\meta\batch8_selected.txt" -Encoding UTF8 | Where-Object { $_.Trim() } | ForEach-Object { $_.Trim() })
$sp = Split-Half $b8
Add-Content -LiteralPath $log -Value "PHASE2 batch8 A=$($sp.A -join '|')  B=$($sp.B -join '|')" -Encoding utf8
$pA = Launch-Stream ($sp.A -join ',') "$root\meta\primer_table_batch8.csv" "$root\out50k\b8_streamA.log"
Start-Sleep -Seconds 15
$pB = Launch-Stream ($sp.B -join ',') "$root\meta\primer_table_batch8.csv" "$root\out50k\b8_streamB.log"
$pA.WaitForExit(); $pB.WaitForExit()
Run-Gates
Add-Content -LiteralPath $log -Value "PHASE2 gates done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

# ---- 3) Phase 3: strategy-A reruns (short-overlap rule) ----
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
$sp = Split-Half $projs
Add-Content -LiteralPath $log -Value "PHASE3 reruns A=$($sp.A -join '|')  B=$($sp.B -join '|')" -Encoding utf8
$pC = Launch-Stream ($sp.A -join ',') "$root\meta\primer_table_stage2.csv" "$root\out50k\rerun_streamA.log"
Start-Sleep -Seconds 15
$pD = Launch-Stream ($sp.B -join ',') "$root\meta\primer_table_stage2.csv" "$root\out50k\rerun_streamB.log"
$pC.WaitForExit(); $pD.WaitForExit()
Run-Gates
& $py "$root\scripts\trunclen_ab_report.py" "$beforeRoot" "$root\out_stage2" ($projs -join ',') --out "$root\out50k\trunclen_ab_report.md" *>> $log
Add-Content -LiteralPath $log -Value "PHASE3 gates done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

# ---- 4) Phase 4: short-overlap A/B on 4 projects (item 3, queued after Phase 3) ----
$ab = @('PRJDB8720','PRJEB45341','PRJEB45605','PRJNA1060287')
$abBefore = "$root\out50k\_archive_shortov_ab"
New-Item -ItemType Directory -Force -Path $abBefore | Out-Null
foreach ($p in $ab) {
  $d = "$root\out_stage2\$p"
  if (Test-Path $d) {
    $arch = Join-Path $abBefore $p
    if (Test-Path $arch) { $arch = "${arch}_$(Get-Date -Format 'HHmmss')" }
    Move-Item -LiteralPath $d -Destination $arch -Force
  }
}
$sp = Split-Half $ab
Add-Content -LiteralPath $log -Value "PHASE4 short-overlap A/B A=$($sp.A -join '|')  B=$($sp.B -join '|')" -Encoding utf8
$pE = Launch-Stream ($sp.A -join ',') "$root\meta\primer_table_ab_short.csv" "$root\out50k\ab_streamA.log"
Start-Sleep -Seconds 15
$pF = Launch-Stream ($sp.B -join ',') "$root\meta\primer_table_ab_short.csv" "$root\out50k\ab_streamB.log"
$pE.WaitForExit(); $pF.WaitForExit()
Run-Gates
& $py "$root\scripts\trunclen_ab_report.py" "$abBefore" "$root\out_stage2" ($ab -join ',') --out "$root\out50k\shortoverlap_ab_report.md" *>> $log
Add-Content -LiteralPath $log -Value "PHASE4 gates done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

# ---- 5) Phase 5: 预注册的"救回"复测（4 个项目） ----
$p5 = @('PRJNA413594','PRJNA555161','PRJNA1476760','PRJEB34506')
$p5Before = "$root\out50k\_archive_salvage_before"
New-Item -ItemType Directory -Force -Path $p5Before | Out-Null
foreach ($p in $p5) {
  $d = "$root\out_stage2\$p"
  if (Test-Path $d) {
    $arch = Join-Path $p5Before $p
    if (Test-Path $arch) { $arch = "${arch}_$(Get-Date -Format 'HHmmss')" }
    Move-Item -LiteralPath $d -Destination $arch -Force
  }
}
$sp = Split-Half $p5
Add-Content -LiteralPath $log -Value "PHASE5 salvage A=$($sp.A -join '|')  B=$($sp.B -join '|')" -Encoding utf8
$pG = Launch-Stream ($sp.A -join ',') "$root\meta\primer_table_phase5.csv" "$root\out50k\salvage_streamA.log"
Start-Sleep -Seconds 15
$pH = Launch-Stream ($sp.B -join ',') "$root\meta\primer_table_phase5.csv" "$root\out50k\salvage_streamB.log"
$pG.WaitForExit(); $pH.WaitForExit()
Run-Gates
& $py "$root\scripts\apply_salvage_rule.py" --projects ($p5 -join ',') *>> $log
& $py "$root\scripts\trunclen_ab_report.py" "$p5Before" "$root\out_stage2" ($p5 -join ',') --out "$root\out50k\salvage_ab_report.md" *>> $log
Add-Content -LiteralPath $log -Value "PHASE5 gates done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

& $py "$root\scripts\select_final45.py" --out "$root\meta\final45_selection.csv" --report "$root\out50k\选集报告.md" *>> $log
Add-Content -LiteralPath $log -Value "PIPELINE ALL DONE $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
