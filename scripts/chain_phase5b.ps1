# Phase 5b：等主链结束后跑 3 个变体测试
#   PRJNA1310568 → A'（maxEE_R=5、不截断）
#   PRJNA596848 / PRJNA665860 → 保留全部样本变体（去掉 --discard-untrimmed）
$root = 'D:\CodexWork\mic-re'
$log  = "$root\out50k\phase5b_chain.log"
$py   = 'D:\bioenv\envs\riz16s\python.exe'
$pwsh = 'C:\Users\hthome\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\powershell\pwsh.exe'
$ref  = 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz'
"phase5b start $(Get-Date -Format 'HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8
while ($true) {
  Start-Sleep -Seconds 300
  $done = (Test-Path "$root\out50k\driver2_chain.log") -and ((Get-Content "$root\out50k\driver2_chain.log" -Raw -ErrorAction SilentlyContinue) -match 'PIPELINE ALL DONE')
  if ($done) { break }
}
Add-Content -LiteralPath $log -Value "main pipeline finished $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
$beforeRoot = "$root\out50k\_archive_phase5b_before"
New-Item -ItemType Directory -Force -Path $beforeRoot | Out-Null
foreach ($p in @('PRJNA1310568','PRJNA596848','PRJNA665860')) {
  $d = "$root\out_stage2\$p"
  if (Test-Path $d) {
    $arch = Join-Path $beforeRoot $p
    if (Test-Path $arch) { $arch = "${arch}_$(Get-Date -Format 'HHmmss')" }
    Move-Item -LiteralPath $d -Destination $arch -Force
    Add-Content -LiteralPath $log -Value "archived $p" -Encoding utf8
  }
}
$pA = Start-Process -FilePath $pwsh -WindowStyle Hidden -PassThru -ArgumentList @('-NoProfile','-Command',
  "conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\run_dada2_parallel_ab.R' 'PRJNA1310568' '$root\data\raw50k' '$root\out_stage2' '$ref' 2 1e7 6 '$root\meta\primer_table_phase5b.csv' *>> '$root\out50k\ab_stream.log'")
Start-Sleep -Seconds 15
$pB = Start-Process -FilePath $pwsh -WindowStyle Hidden -PassThru -ArgumentList @('-NoProfile','-Command',
  "conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\run_dada2_parallel_keepall.R' 'PRJNA596848,PRJNA665860' '$root\data\raw50k' '$root\out_stage2' '$ref' 2 1e7 6 '$root\meta\primer_table_phase5b.csv' *>> '$root\out50k\keepall_stream.log'")
Add-Content -LiteralPath $log -Value "streams launched $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
$pA.WaitForExit(); $pB.WaitForExit()
Add-Content -LiteralPath $log -Value "streams done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
& $py "$root\scripts\master_scan.py" *>> $log
& $py "$root\scripts\marker_gene_qc.py" *>> $log
& $py "$root\scripts\summarize_w1.py" "$root\out_stage2" *>> $log
& $py "$root\scripts\classify_w1.py" "$root\out_stage2" *>> $log
& $py "$root\scripts\merged_report.py" *>> $log
& $py "$root\scripts\pick_phase5b_winner.py" *>> $log
& $py "$root\scripts\apply_salvage_rule.py" --projects 'PRJNA413594,PRJNA555161,PRJNA1476760,PRJEB34506,PRJNA1310568' *>> $log
& $py "$root\scripts\trunclen_ab_report.py" "$beforeRoot" "$root\out_stage2" 'PRJNA1310568,PRJNA596848,PRJNA665860' --out "$root\out50k\phase5b_ab_report.md" *>> $log
& $py "$root\scripts\select_final45.py" --out "$root\meta\final45_selection.csv" --report "$root\out50k\选集报告.md" *>> $log
Add-Content -LiteralPath $log -Value "PHASE5B ALL DONE $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
