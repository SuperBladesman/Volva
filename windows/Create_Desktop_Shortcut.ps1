# Creates a desktop (and Start menu) shortcut to Volva.exe in this folder.
$exe = Join-Path $PSScriptRoot 'Volva.exe'
if (-not (Test-Path $exe)) { Write-Host "Volva.exe not found next to this script."; exit 1 }
$shell = New-Object -ComObject WScript.Shell
foreach ($dir in @([Environment]::GetFolderPath('Desktop'), (Join-Path ([Environment]::GetFolderPath('Programs')) ''))) {
    $lnk = $shell.CreateShortcut((Join-Path $dir 'Volva 4.1.lnk'))
    $lnk.TargetPath = $exe
    $lnk.WorkingDirectory = [Environment]::GetFolderPath('MyDocuments')
    $lnk.IconLocation = "$exe,0"
    $lnk.Description = 'Volva 4.1 - XRD crystallite size and strain analysis'
    $lnk.Save()
    Write-Host "Shortcut created in $dir"
}
