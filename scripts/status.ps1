# 一键状态：各链进度、项目数、磁盘
$root = 'D:\CodexWork\mic-re'
"===== 时间 $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') ====="
"--- 运行中的链 ---"
Get-CimInstance Win32_Process | Where-Object { $_.Name -eq 'pwsh.exe' -and $_.CommandLine -match '\-File' } | ForEach-Object {
  $f = [regex]::Match($_.CommandLine, '([A-Za-z0-9_]+\.ps1)').Groups[1].Value
  if ($f) { "  $f  (pid $($_.ProcessId))" }
}
"--- R 进程结构（master / PSOCK worker） ---"
$rps = @(Get-CimInstance Win32_Process -Filter "Name='Rscript.exe'")
$masters = @($rps | Where-Object { $_.CommandLine -and $_.CommandLine -notmatch 'workRSOCK|slaveRSOCK' })
$workers = @($rps | Where-Object { $_.CommandLine -match 'workRSOCK|slaveRSOCK' })
"  master = $($masters.Count)   worker = $($workers.Count)   curl = $(@(Get-Process curl -ErrorAction SilentlyContinue).Count)"
foreach ($m in $masters) {
  $g = Get-Process -Id $m.ProcessId -ErrorAction SilentlyContinue
  $rss = if ($g) { [math]::Round($g.WorkingSet64/1MB,0) } else { 0 }
  if ($rss -lt 20) { continue }
  $pl = [regex]::Match($m.CommandLine, '(PRJ[A-Z]*\d+(?:,PRJ[A-Z]*\d+)*)').Groups[1].Value
  if ($pl.Length -gt 70) { $pl = $pl.Substring(0,70) + '...' }
  "    master pid=$($m.ProcessId) RSS=${rss}MB  正在处理=$pl"
}
foreach ($w in $workers) {
  $g = Get-Process -Id $w.ProcessId -ErrorAction SilentlyContinue
  $rss = if ($g) { [math]::Round($g.WorkingSet64/1MB,0) } else { 0 }
  $port = [regex]::Match($w.CommandLine,'PORT=(\d+)').Groups[1].Value
  "    worker pid=$($w.ProcessId) RSS=${rss}MB  port=$port"
}
"  (worker 只在 filterAndTrim/dada/mergePairs 干活；learnErrors、去嵌合、分类学只在 master 跑，"
"   所以会出现 worker 0% CPU 却占着内存的正常现象)"
"--- 链日志末行 ---"
foreach ($n in 'driver_chain.log', 'pend_chain.log', 'batch7_chain.log', 'batch8_chain.log', 'trunclen2_chain.log') {
  $p = "$root\out50k\$n"
  if (Test-Path $p) {
    $last = (Get-Content $p -Encoding UTF8 | Where-Object { $_.Trim() } | Select-Object -Last 1)
    "  {0,-22} {1}" -f $n, $last
  }
}
"--- 处理结果目录 ---"
"  out50k   = $(@(Get-ChildItem "$root\out50k" -Directory | Where-Object { $_.Name -match '^PRJ' }).Count) 个已完成"
"  out_stage2 = $(@(Get-ChildItem "$root\out_stage2" -Directory | Where-Object { $_.Name -match '^PRJ' }).Count) 个已完成"
"--- 数据 ---"
"  data/raw    = $(@(Get-ChildItem "$root\data\raw" -Directory).Count) 项目"
"  data/raw50k = $(@(Get-ChildItem "$root\data\raw50k" -Directory).Count) 项目"
"--- 磁盘/内存 ---"
Get-PSDrive C, D | ForEach-Object { "  $($_.Name): 剩余 $([math]::Round($_.Free/1GB,1)) GB" }
"  可用内存 $([math]::Round((Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory/1MB,1)) GB"
"--- 最近一次闸门结果 ---"
if (Test-Path "$root\out50k\marker_gene_qc.csv") {
  Import-Csv "$root\out50k\marker_gene_qc.csv" | Group-Object verdict | ForEach-Object { "  $($_.Name) = $($_.Count)" }
}
