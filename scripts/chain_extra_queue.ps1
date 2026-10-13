# 追加队列（用户 2026-10-10 决定）：等主链 RESUME ALL DONE 后跑三个补充项目，再刷新选集。
#   PRJNA607209  小麦 单端长读  数据齐但漏跑             -> run_dada2_parallel.R
#   PRJNA554261  水稻 单端      混合文库，需去 515F 引物 -> run_dada2_parallel_keepall.R
#   PRJNA578725  小麦 双端      去 27F/534R + truncLen 271,271 -> run_dada2_parallel.R
# 资源纪律：全程最多两条流各 2 核；主链未结束前不动手。
$root = 'D:\CodexWork\mic-re'
$log  = "$root\out50k\extra_queue.log"
$py   = 'D:\bioenv\envs\riz16s\python.exe'
$pwsh = 'C:\Users\hthome\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\powershell\pwsh.exe'
$ref  = 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz'
$tab  = "$root\meta\primer_table_extra.csv"

function Say { param([string]$m) Add-Content -LiteralPath $log -Value "$((Get-Date -Format 'HH:mm:ss'))  $m" -Encoding utf8 }
function Launch {
  param([string]$script, [string]$projs, [string]$rawroot, [string]$streamlog)
  Start-Process -FilePath $pwsh -WindowStyle Hidden -PassThru -ArgumentList @('-NoProfile','-Command',
    "conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\$script' '$projs' '$rawroot' '$root\out_stage2' '$ref' 2 1e7 6 '$tab' *>> '$streamlog'")
}

"extra queue start $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8

# ---- 1) 等主链收尾（标记唯一，只在真正结束时出现）----
while ($true) {
  Start-Sleep -Seconds 300
  $done = (Test-Path "$root\out50k\resume_all.log") -and ((Get-Content "$root\out50k\resume_all.log" -Raw -ErrorAction SilentlyContinue) -match 'RESUME ALL DONE')
  if ($done) { break }
  $alive = [bool](Get-Process -Id 3736 -ErrorAction SilentlyContinue)
  $still = [bool](Get-Process -Name Rscript -ErrorAction SilentlyContinue)
  if ((-not $alive) -and (-not $still)) { Say "MAIN CHAIN GONE and no R running, and resume_all.log has no RESUME ALL DONE -> ABORT extra queue (needs manual resume first)"; exit 1 }
}
Say "main chain finished, extra queue starting"

# ---- 2) 归档旧结果（§12.3：重跑前必须移走旧目录，否则管线会因 project_qc+asv_taxonomy 同时存在而整项目跳过）----
$before = "$root\out50k\_archive_extra_before"
New-Item -ItemType Directory -Force -Path $before | Out-Null
foreach ($p in @('PRJNA554261','PRJNA578725')) {
  if (Test-Path "$root\out_stage2\$p") {
    $dst = Join-Path $before $p
    if (Test-Path $dst) { $dst = "${dst}_$(Get-Date -Format 'HHmmss')" }
    Move-Item -LiteralPath "$root\out_stage2\$p" -Destination $dst -Force
    Say "archived old result: $p"
  }
}

# ---- 3) 两条流：单端长读 607209 ｜ 单端混合文库 554261 ----
$a = Launch 'run_dada2_parallel.R'        'PRJNA607209' "$root\data\raw50k_longread" "$root\out50k\ex_streamA.log"
Start-Sleep -Seconds 15
$b = Launch 'run_dada2_parallel_keepall.R' 'PRJNA554261' "$root\data\raw50k"           "$root\out50k\ex_streamB.log"
$a.WaitForExit(); $b.WaitForExit()
Say "streams 607209/554261 exited"

# ---- 4) 单流：双端 578725（去引物 + truncLen 271,271）----
$c = Launch 'run_dada2_parallel.R' 'PRJNA578725' "$root\data\raw50k" "$root\out50k\ex_streamC.log"
$c.WaitForExit()
Say "stream 578725 exited"

# ---- 5) 刷新闸门与报告 ----
& $py "$root\scripts\master_scan.py"        *>> $log
& $py "$root\scripts\marker_gene_qc.py"     *>> $log
& $py "$root\scripts\summarize_w1.py" "$root\out_stage2" *>> $log
& $py "$root\scripts\classify_w1.py" "$root\out_stage2"  *>> $log
& $py "$root\scripts\merged_report.py"      *>> $log
& $py "$root\scripts\select_final45.py" --out "$root\meta\final45_selection.csv" --report "$root\out50k\选集报告.md" *>> $log
& $py "$root\scripts\tax_bins.py" --selection "$root\meta\final45_selection.csv" --out "$root\out50k\tax_bin_report.md" *>> $log
Add-Content -LiteralPath $log -Value "$(Get-Date -Format 'HH:mm:ss')  EXTRA QUEUE ALL DONE" -Encoding utf8
