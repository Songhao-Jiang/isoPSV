@echo off
rem tests/run_all.bat - Windows one-click test runner (P3-3, 2026-09-25)
rem Runs every test suite via tests/run_all.R and writes tests/last_run.log
set LC_CTYPE=
set LC_ALL=
set LANG=
cd /d "%~dp0"
Rscript run_all.R
echo.
echo exit code: %errorlevel%
pause
