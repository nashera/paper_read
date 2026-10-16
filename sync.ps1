$time=(Get-Date -UFormat "%A %m/%d/%Y %R %Z")
$currentPath=(Get-Location).path
Set-Location -Path $PSScriptRoot
git add --all;
git commit -m "$time"
git pull --rebase
expect.exe auto.lua
Set-Location $currentPath
