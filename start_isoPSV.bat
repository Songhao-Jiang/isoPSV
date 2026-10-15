@echo off
setlocal
chcp 65001 >nul
set LC_CTYPE=
set LC_ALL=
set LANG=
cd /d "%~dp0"

echo ============================================
echo  isoPSV multiengine - local web shell launcher
echo  http://127.0.0.1:8780
echo  (closing this window stops the web shell)
echo ============================================

rem bundled portable runtimes take precedence
set "PY="
set "PYARG="
if exist "%~dp0runtime\python\python.exe" set "PY=%~dp0runtime\python\python.exe"
if not defined PY (where python >nul 2>nul && set "PY=python")
if not defined PY (
  where py >nul 2>nul && set "PY=py" && set "PYARG=-3"
)
if not defined PY (where python3 >nul 2>nul && set "PY=python3")

if exist "%~dp0runtime\R\bin\x64\Rscript.exe" (
  set "PATH=%~dp0runtime\R\bin\x64;%PATH%"
  set "RSCRIPT=%~dp0runtime\R\bin\x64\Rscript.exe"
  echo using bundled portable R
)

if not defined PY (
  echo.
  echo [isoPSV] Python not found. The web shell needs Python 3.8+.
  echo Keep runtime\python in the package folder, or install Python
  echo from www.python.org and put it on PATH.
  pause
  exit /b 1
)

rem v1.0 kill the stale server still holding port 8780
rem after an upgrade - otherwise the new server fails to bind and the browser
rem keeps talking to the OLD process ("page looks new but buttons dead").
powershell -NoProfile -Command "Get-NetTCPConnection -LocalPort 8780 -State Listen -ErrorAction SilentlyContinue | ForEach-Object { Write-Host ('[isoPSV] stopping stale server PID ' + $_.OwningProcess); Stop-Process -Id $_.OwningProcess -Force -ErrorAction SilentlyContinue }"

rem open the browser 2 seconds after the server starts
start "" /b powershell -NoProfile -WindowStyle Hidden -Command "Start-Sleep -Seconds 2; Start-Process 'http://127.0.0.1:8780'"

"%PY%" %PYARG% "%~dp0web\server.py"
pause
