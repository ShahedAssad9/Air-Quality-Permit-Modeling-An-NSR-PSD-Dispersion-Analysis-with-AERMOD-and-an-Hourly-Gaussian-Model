@echo off
title AERMOD project - step C: AERMOD runs
cd /d "%~dp0"
echo Step C: running AERMOD (19 runs, several at a time, below-normal priority).
echo Progress is written to logs\C_aermod.log
echo.
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0tools\run_aermod.ps1" %*
echo.
echo Step C finished. You can close this window.
pause >nul
