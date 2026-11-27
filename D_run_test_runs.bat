@echo off
title AERMOD project - step D: AERMOD test runs (downwash and terrain switched off)
cd /d "%~dp0"
echo Step D (optional): 10 extra AERMOD runs that explain why AERMOD and Part A differ.
echo   SO2 1-hr, NO2 1-hr and PM2.5 24-hr, each with building downwash off, terrain off, and both off,
echo   plus SO2 1-hr with the stack raised to its 75 m GEP height. About 15 minutes.
echo Progress is written to logs\C_aermod_diagnostic.log
echo Do not start any other step while this window is open.
echo.
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0tools\run_aermod.ps1" diagnostic
echo.
echo Step D finished. You can close this window.
pause >nul
