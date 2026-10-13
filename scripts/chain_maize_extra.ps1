# 玉米扩池：等追加队列（EXTRA QUEUE ALL DONE）后，抽样并处理 3 个已筛通过的玉米候选，再刷新闸门与选集。
# 资源纪律：全程最多两条流各 2 核；与追加队列串行，不并行。
$root = 'D:\CodexWork\mic-re'
$log  = "$root\out50k\maize_extra.log"
$py   = 'D:\bioenv\envs\riz16s\python.exe'
$pwsh = 'C:\Users\hthome\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\powershell\pwsh.exe'
$ref  = 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz'
$projs = 'PRJNA633917,PRJNA662201,PRJNA967263'
$tab  = "$root\meta\primer_table_maize.csv"

function Say { param([string]$m) Add-Content -LiteralPath $log -Value "$((Get-Date -Format 'HH:mm:ss'))  $m" -Encoding utf8 }
function Launch { param([string]$script,[string]$pjs,[string]$raw,[string]$slog)
  Start-Process -FilePath $pwsh -WindowStyle Hidden -PassThru -ArgumentList @('-NoProfile','-Command',
    "conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\$script' '$pjs' '$raw' '$root\out_stage2' '$ref' 2 1e7 6 '$tab' *>> '$slog'") }

"maize extra start $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8

# ---- 1) 等追加队列收尾（标记唯一）----
while ($true) {
  Start-Sleep -Seconds 300
  $done = (Test-Path "$root\out50k\extra_queue.log") -and ((Get-Content "$root\out50k\extra_queue.log" -Raw -ErrorAction SilentlyContinue) -match 'EXTRA QUEUE ALL DONE')
  if ($done) { break }
  $alive = [bool](Get-Process -Id 7568 -ErrorAction SilentlyContinue)
  $r = [bool](Get-Process -Name Rscript -ErrorAction SilentlyContinue)
  if ((-not $alive) -and (-not $r)) { Say 'EXTRA QUEUE GONE without ALL DONE -> ABORT (manual resume needed)'; exit 1 }
}
Say 'extra queue finished; maize expansion starting'

# ---- 2) gzip 完整性校验（§12.7：体积不够不算成功）----
& $py "$root\scripts\verify_gzip_raw.py" "$root\data\raw" $projs --workers 2 *>> $log
if ($LASTEXITCODE -ne 0) { Say "gzip check FAILED -> ABORT"; exit 2 }
Say 'gzip check passed'

# ---- 3) 50k 抽样 ----
& $py "$root\scripts\subsample_fastq.py" "$root\data\raw" "$root\data\raw50k" 50000 20261007 $projs *>> $log
Say 'subsampled 50k'

# ---- 4) 引物检测 -> truncLen 规划 -> 管线参数表 ----
& $py "$root\scripts\prep_maize_primers.py" *>> $log
Copy-Item -LiteralPath "$root\meta\primer_table_maize_raw.csv" -Destination $tab -Force
& $py "$root\scripts\plan_trunclen.py" "$root\data\raw50k" $projs --asv-root '' --primer-table $tab --tag maize --out "$root\meta\trunclen_plan_maize.csv" *>> $log
& $py "$root\scripts\apply_trunclen_plan.py" "$root\meta\trunclen_plan_maize.csv" $tab "$root\meta\primer_table_maize_raw.csv" *>> $log
Say 'primer table ready'
Get-Content -LiteralPath $tab -Encoding UTF8 | ForEach-Object { Say ("  tab: " + $_) }

# ---- 5) 两条流 ----
$mzA = Launch 'run_dada2_parallel.R' 'PRJNA633917' "$root\data\raw50k" "$root\out50k\mz_streamA.log"
Start-Sleep -Seconds 15
$mzB = Launch 'run_dada2_parallel.R' 'PRJNA662201,PRJNA967263' "$root\data\raw50k" "$root\out50k\mz_streamB.log"
$mzA.WaitForExit(); $mzB.WaitForExit()
Say 'maize streams exited'

# ---- 6) 闸门 + 选集 + 回退规模报告 ----
& $py "$root\scripts\master_scan.py" *>> $log
& $py "$root\scripts\marker_gene_qc.py" *>> $log
& $py "$root\scripts\summarize_w1.py" "$root\out_stage2" *>> $log
& $py "$root\scripts\classify_w1.py" "$root\out_stage2" *>> $log
& $py "$root\scripts\merged_report.py" *>> $log
& $py "$root\scripts\select_final45.py" --out "$root\meta\final45_selection.csv" --report "$root\out50k\选集报告.md" *>> $log
& $py "$root\scripts\tax_bins.py" --selection "$root\meta\final45_selection.csv" --out "$root\out50k\tax_bin_report.md" *>> $log
Add-Content -LiteralPath $log -Value "$(Get-Date -Format 'HH:mm:ss')  MAIZE EXTRA ALL DONE" -Encoding utf8
