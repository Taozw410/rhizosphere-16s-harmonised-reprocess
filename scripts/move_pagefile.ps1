$log = 'D:\CodexWork\mic-re\out50k\pagefile_move.log'
"pagefile move start $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8
try {
    # 1) 关闭自动管理
    $cs = Get-CimInstance Win32_ComputerSystem
    if ($cs.AutomaticManagedPagefile) {
        Set-CimInstance -InputObject $cs -Property @{AutomaticManagedPagefile = $false} -ErrorAction Stop
        Add-Content $log "automatic management disabled" -Encoding utf8
    } else {
        Add-Content $log "automatic management already off" -Encoding utf8
    }
    # 2) 删除现有页面文件设置（即 C 盘那个）
    Get-CimInstance Win32_PageFileSetting -ErrorAction SilentlyContinue | ForEach-Object {
        Add-Content $log "removing setting: $($_.Name)" -Encoding utf8
        Remove-CimInstance -InputObject $_
    }
    # 3) 在 D 盘新建：初始 8 GB / 最大 24 GB
    $pf = New-CimInstance -ClassName Win32_PageFileSetting -Property @{
        Name = 'D:\pagefile.sys'; InitialSize = 8192; MaximumSize = 24576
    } -ErrorAction Stop
    Add-Content $log "created: $($pf.Name) init=$($pf.InitialSize)MB max=$($pf.MaximumSize)MB" -Encoding utf8
    # 4) 复核
    Add-Content $log "--- verify ---" -Encoding utf8
    Get-CimInstance Win32_PageFileSetting | ForEach-Object { Add-Content $log "  $($_.Name) init=$($_.InitialSize) max=$($_.MaximumSize)" -Encoding utf8 }
    $c2 = Get-CimInstance Win32_ComputerSystem
    Add-Content $log "AutomaticManagedPagefile = $($c2.AutomaticManagedPagefile)" -Encoding utf8
    Add-Content $log "SUCCESS" -Encoding utf8
} catch {
    Add-Content $log "ERROR: $($_.Exception.Message)" -Encoding utf8
}
