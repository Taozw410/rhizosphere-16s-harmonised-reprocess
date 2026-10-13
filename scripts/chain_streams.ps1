$root = 'D:\CodexWork\mic-re'
$log  = "$root\out50k\stream_chain.log"
"stream chain start $(Get-Date -Format 'HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8
$A = Get-Content "$root\meta\stream_A.txt" -Raw
$B = Get-Content "$root\meta\stream_B.txt" -Raw
$common = "$root\data\raw50k", "$root\out_stage2", 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz', '2', '1e7', '6', "$root\meta\primer_table_stage2.csv"

Add-Content $log "launch A ($(($A -split ',').Count) projects, ncore=2) $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
$pA = Start-Process -FilePath 'C:\Users\hthome\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\powershell\pwsh.exe' -WindowStyle Hidden -PassThru -ArgumentList @('-NoProfile','-Command',"conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\run_dada2_parallel.R' '$A' '$root\data\raw50k' '$root\out_stage2' 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz' 2 1e7 6 '$root\meta\primer_table_stage2.csv' *>> '$root\out50k\streamA.log'")

Start-Sleep -Seconds 10
Add-Content $log "launch B ($(($B -split ',').Count) projects, ncore=2) $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
$pB = Start-Process -FilePath 'C:\Users\hthome\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\powershell\pwsh.exe' -WindowStyle Hidden -PassThru -ArgumentList @('-NoProfile','-Command',"conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\run_dada2_parallel.R' '$B' '$root\data\raw50k' '$root\out_stage2' 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz' 2 1e7 6 '$root\meta\primer_table_stage2.csv' *>> '$root\out50k\streamB.log'")

Add-Content $log "waiting for both streams ... pidA=$($pA.Id) pidB=$($pB.Id)" -Encoding utf8
$pA.WaitForExit(); Add-Content $log "stream A done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
$pB.WaitForExit(); Add-Content $log "stream B done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

Add-Content $log "--- long-read batch start $(Get-Date -Format 'HH:mm:ss') ---" -Encoding utf8
conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript "$root\scripts\run_dada2_parallel.R" ALL "$root\data\raw50k_longread" "$root\out_stage2" 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz' 3 1e7 6 "$root\meta\primer_table_stage2.csv" *>> $log
Add-Content $log "--- long-read done $(Get-Date -Format 'HH:mm:ss') ---" -Encoding utf8
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\summarize_w1.py" "$root\out_stage2" *>> $log
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\classify_w1.py" "$root\out_stage2" *>> $log
Add-Content $log "ALL DONE $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
