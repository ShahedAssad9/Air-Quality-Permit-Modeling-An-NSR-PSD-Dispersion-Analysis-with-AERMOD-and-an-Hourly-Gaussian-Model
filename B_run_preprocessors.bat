@echo off
title AERMOD project - step B: prepare the data and run AERSURFACE, AERMINUTE, AERMET, AERMAP, BPIPPRM
cd /d "%~dp0"
echo Step B: preparing the data and running the EPA preprocessors.
echo Progress is written to logs\B_preprocessors.log
echo.
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0tools\run_preprocessors.ps1" %*
echo.
echo Step B finished. You can close this window.
pause >nul
