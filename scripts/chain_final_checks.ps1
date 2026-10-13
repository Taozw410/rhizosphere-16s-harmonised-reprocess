# 收尾三件（用户 2026-10-11 决定）：A' 正式裁决 + IDTAXA 交叉验证 + 最终 45 项目的功效模拟。
# 串行占用资源：A' 阶段两条流各 2 核；随后 IDTAXA 单核；全程不与其他链并行。
$root = 'D:\CodexWork\mic-re'
$log  = "$root\out50k\final_checks.log"
$py   = 'D:\bioenv\envs\riz16s\python.exe'
$pwsh = 'C:\Users\hthome\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\powershell\pwsh.exe'
$ref  = 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz'
$tab  = "$root\meta\primer_table_abprime.csv"
$abroot = "$root\out_abprime"

function Say { param([string]$m) Add-Content -LiteralPath $log -Value "$((Get-Date -Format 'HH:mm:ss'))  $m" -Encoding utf8 }
function LaunchAB { param([string]$pjs,[string]$slog)
  Start-Process -FilePath $pwsh -WindowStyle Hidden -PassThru -ArgumentList @('-NoProfile','-Command',
    "conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\run_dada2_parallel_ab.R' '$pjs' '$root\data\raw50k' '$abroot' '$ref' 2 1e7 6 '$tab' *>> '$slog'") }
function Run { param([string]$exe,[string[]]$argv,[string]$label)
  & $exe @argv *>> $log
  $rc = $LASTEXITCODE
  Say ("{0} exit={1}" -f $label, $rc)
  return $rc }

"final checks start $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8
New-Item -ItemType Directory -Force -Path $abroot | Out-Null

# ---- 1) A' 正式裁决：4 个项目跑 A' 臂到独立输出根（不覆盖 A 臂）----
Say 'Aprime arm starting (2 streams)'
$s1 = LaunchAB 'PRJNA715843,PRJNA593708,PRJNA427813' "$root\out50k\ab_stream1.log"
Start-Sleep -Seconds 15
$s2 = LaunchAB 'PRJNA596848' "$root\out50k\ab_stream2.log"
$s1.WaitForExit(); $s2.WaitForExit()
Say 'Aprime streams exited'
[void](Run $py @("$root\scripts\ab_prime_report.py", "--abprime", $abroot, "--out", "$root\out50k\abprime_report.md") 'ab_prime_report')

# ---- 2) 最终 45 项目的功效模拟（§8.3 要求）----
[void](Run $py @("$root\scripts\power_simulation.py") 'power_simulation')

# ---- 3) IDTAXA 交叉验证（§3；DECIPHER 官方预训练模型不可达 -> 用本机 SILVA 138.1 自建）----
$train = "$root\ref_region\idtaxa_train.fa"
if (-not (Test-Path $train)) {
  [void](Run $py @("$root\scripts\make_idtaxa_train.py", "--cap", "6000", "--out", $train) 'make_idtaxa_train')
}
Say 'IDTAXA starting (single core)'
$idlog = "$root\out50k\idtaxa.log"
Start-Process -FilePath $pwsh -WindowStyle Hidden -PassThru -ArgumentList @('-NoProfile','-Command',
  "conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\idtaxa_crossval.R' '$train' '$root\out50k\idtaxa' 300 '$root\meta\final45_selection.csv' *>> '$idlog'") | Out-Null
Wait-Process -Name Rscript -Timeout 21600 -ErrorAction SilentlyContinue
Say "IDTAXA step finished (see out50k/idtaxa.log)"

Add-Content -LiteralPath $log -Value "$(Get-Date -Format 'HH:mm:ss')  FINAL CHECKS ALL DONE" -Encoding utf8
