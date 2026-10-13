# 长读水稻补跑：等 Phase 5b 结束后，处理 3 个此前暂停的长读项目（单端，raw 在 data/raw50k_longread）
$root = 'D:\CodexWork\mic-re'
$log  = "$root\out50k\lr_resume_chain.log"
$py   = 'D:\bioenv\envs\riz16s\python.exe'
$pwsh = 'C:\Users\hthome\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\powershell\pwsh.exe'
$ref  = 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz'
"lr resume chain start $(Get-Date -Format 'HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8
while ($true) {
  Start-Sleep -Seconds 300
  $done = (Test-Path "$root\out50k\phase5b_chain.log") -and ((Get-Content "$root\out50k\phase5b_chain.log" -Raw -ErrorAction SilentlyContinue) -match 'PHASE5B ALL DONE')
  if ($done) { break }
}
Add-Content -LiteralPath $log -Value "phase5b finished $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
$projs = @('PRJNA1157226','PRJNA1261944','PRJNA509515')
# 两个流：1261944(最大值) 单独一流，另两个一流
$pA = Start-Process -FilePath $pwsh -WindowStyle Hidden -PassThru -ArgumentList @('-NoProfile','-Command',
  "conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\run_dada2_parallel.R' 'PRJNA1261944' '$root\data\raw50k_longread' '$root\out_stage2' '$ref' 2 1e7 6 '$root\meta\primer_table_lr_resume.csv' *>> '$root\out50k\lr_streamA.log'")
Start-Sleep -Seconds 15
$pB = Start-Process -FilePath $pwsh -WindowStyle Hidden -PassThru -ArgumentList @('-NoProfile','-Command',
  "conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\run_dada2_parallel.R' 'PRJNA1157226,PRJNA509515' '$root\data\raw50k_longread' '$root\out_stage2' '$ref' 2 1e7 6 '$root\meta\primer_table_lr_resume.csv' *>> '$root\out50k\lr_streamB.log'")
Add-Content -LiteralPath $log -Value "streams launched $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
$pA.WaitForExit(); $pB.WaitForExit()
Add-Content -LiteralPath $log -Value "streams done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
& $py "$root\scripts\master_scan.py" *>> $log
& $py "$root\scripts\marker_gene_qc.py" *>> $log
& $py "$root\scripts\summarize_w1.py" "$root\out_stage2" *>> $log
& $py "$root\scripts\classify_w1.py" "$root\out_stage2" *>> $log
& $py "$root\scripts\merged_report.py" *>> $log
& $py "$root\scripts\taxonomy_diagnosis.py" *>> $log
& $py "$root\scripts\select_final45.py" --out "$root\meta\final45_selection.csv" --report "$root\out50k\选集报告.md" *>> $log
Add-Content -LiteralPath $log -Value "LR RESUME ALL DONE $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
