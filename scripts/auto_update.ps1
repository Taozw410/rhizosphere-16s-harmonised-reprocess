$log = 'D:\CodexWork\mic-re\out50k\auto_update.log'
"watcher start $(Get-Date -Format 'HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8
# 等所有 R 进程结束（1479906 收尾 + PRJNA996221）
while (@(Get-Process Rscript -ErrorAction SilentlyContinue).Count -gt 0) { Start-Sleep -Seconds 60 }
Add-Content -LiteralPath $log -Value "R finished $(Get-Date -Format 'HH:mm:ss'); regenerating" -Encoding utf8
& 'D:\bioenv\envs\riz16s\python.exe' 'D:\CodexWork\mic-re\scripts\summarize_w1.py' *>> $log
& 'D:\bioenv\envs\riz16s\python.exe' 'D:\CodexWork\mic-re\scripts\classify_w1.py' *>> $log
Add-Content -LiteralPath $log -Value "done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
