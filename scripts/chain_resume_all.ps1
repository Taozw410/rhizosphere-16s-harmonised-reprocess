# 续跑链：接上被中断的剩余队列（Phase5b 收尾 → 长读补跑 → Phase5c → Phase5d → 最终选集）
$root = 'D:\CodexWork\mic-re'
$log  = "$root\out50k\resume_all.log"
$py   = 'D:\bioenv\envs\riz16s\python.exe'
$pwsh = 'C:\Users\hthome\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\powershell\pwsh.exe'
$ref  = 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz'

function Launch {
  param([string]$script, [string]$projs, [string]$rawroot, [string]$table, [string]$streamlog)
  Start-Process -FilePath $pwsh -WindowStyle Hidden -PassThru -ArgumentList @('-NoProfile','-Command',
    "conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\$script' '$projs' '$rawroot' '$root\out_stage2' '$ref' 2 1e7 6 '$table' *>> '$streamlog'")
}
function Gates {
  & $py "$root\scripts\master_scan.py" *>> $log
  & $py "$root\scripts\marker_gene_qc.py" *>> $log
  & $py "$root\scripts\summarize_w1.py" "$root\out_stage2" *>> $log
  & $py "$root\scripts\classify_w1.py" "$root\out_stage2" *>> $log
  & $py "$root\scripts\merged_report.py" *>> $log
}

"resume-all start $(Get-Date -Format 'HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8

# ---- A) Phase 5b 收尾：596848 + 665860 的 keep-all 变体 ----
$pA = Launch 'run_dada2_parallel_keepall.R' 'PRJNA596848,PRJNA665860' "$root\data\raw50k" "$root\meta\primer_table_phase5b.csv" "$root\out50k\keepall_stream.log"
$pA.WaitForExit()
Add-Content -LiteralPath $log -Value "PHASE5B streams done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
Gates
& $py "$root\scripts\pick_phase5b_winner.py" *>> $log
& $py "$root\scripts\apply_salvage_rule.py" --projects 'PRJNA1476760,PRJNA1310568' *>> $log
& $py "$root\scripts\trunclen_ab_report.py" "$root\out50k\_archive_phase5b_before" "$root\out_stage2" 'PRJNA1310568,PRJNA596848,PRJNA665860' --out "$root\out50k\phase5b_ab_report.md" *>> $log
Add-Content -LiteralPath $log -Value "PHASE5B ALL DONE $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

# ---- B) 长读补跑：3 个水稻单端项目 ----
$b1 = Launch 'run_dada2_parallel.R' 'PRJNA1261944' "$root\data\raw50k_longread" "$root\meta\primer_table_lr_resume.csv" "$root\out50k\lr_streamA.log"
Start-Sleep -Seconds 15
$b2 = Launch 'run_dada2_parallel.R' 'PRJNA1157226,PRJNA509515' "$root\data\raw50k_longread" "$root\meta\primer_table_lr_resume.csv" "$root\out50k\lr_streamB.log"
$b1.WaitForExit(); $b2.WaitForExit()
Add-Content -LiteralPath $log -Value "LR streams done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
Gates
& $py "$root\scripts\taxonomy_diagnosis.py" *>> $log
Add-Content -LiteralPath $log -Value "LR RESUME ALL DONE $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

# ---- C) Phase 5c：PRJNA1476760 单端抢救 ----
$beforeC = "$root\out50k\_archive_single_end_before"
New-Item -ItemType Directory -Force -Path $beforeC | Out-Null
if (Test-Path "$root\out_stage2\PRJNA1476760") { Move-Item -LiteralPath "$root\out_stage2\PRJNA1476760" -Destination (Join-Path $beforeC 'PRJNA1476760') -Force }
$c1 = Launch 'run_dada2_parallel.R' 'PRJNA1476760' "$root\data\raw50k_se" "$root\meta\primer_table_se_1476760.csv" "$root\out50k\se_stream.log"
$c1.WaitForExit()
Add-Content -LiteralPath $log -Value "PHASE5C stream done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
Gates
& $py "$root\scripts\taxonomy_diagnosis.py" *>> $log
Add-Content -LiteralPath $log -Value "PHASE5C ALL DONE $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

# ---- D) Phase 5d：596848 / 665860 的 244,244 变体（短读长样本组修正） ----
$beforeD = "$root\out50k\_archive_phase5d_before"
New-Item -ItemType Directory -Force -Path $beforeD | Out-Null
foreach ($p in @('PRJNA596848','PRJNA665860')) {
  if (Test-Path "$root\out_stage2\$p") {
    $arch = Join-Path $beforeD $p
    if (Test-Path $arch) { $arch = "${arch}_$(Get-Date -Format 'HHmmss')" }
    Move-Item -LiteralPath "$root\out_stage2\$p" -Destination $arch -Force
  }
}
$d1 = Launch 'run_dada2_parallel_keepall.R' 'PRJNA596848' "$root\data\raw50k" "$root\meta\primer_table_phase5d.csv" "$root\out50k\f5d_streamA.log"
Start-Sleep -Seconds 15
$d2 = Launch 'run_dada2_parallel_keepall.R' 'PRJNA665860' "$root\data\raw50k" "$root\meta\primer_table_phase5d.csv" "$root\out50k\f5d_streamB.log"
$d1.WaitForExit(); $d2.WaitForExit()
Add-Content -LiteralPath $log -Value "PHASE5D streams done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
Gates
& $py "$root\scripts\pick_phase5b_winner.py" *>> $log
& $py "$root\scripts\trunclen_ab_report.py" "$root\out50k\_archive_phase5b_before" "$root\out_stage2" 'PRJNA596848,PRJNA665860' --out "$root\out50k\phase5d_ab_report.md" *>> $log
Add-Content -LiteralPath $log -Value "PHASE5D ALL DONE $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

# ---- E) 最终选集 ----
& $py "$root\scripts\select_final45.py" --out "$root\meta\final45_selection.csv" --report "$root\out50k\选集报告.md" *>> $log
Add-Content -LiteralPath $log -Value "RESUME ALL DONE $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
