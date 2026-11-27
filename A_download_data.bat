@echo off
title AERMOD project - step A: download the EPA programs and the data
cd /d "%~dp0"
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0tools\stop_downloads.ps1"
echo Downloading the EPA programs and the data (about 600 MB), 24 files at a time; interrupted files are resumed.
echo Progress is written to logs\A_download_data.log
echo.
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0tools\download_data.ps1"
echo.
echo Step A finished. You can close this window.
pause >nul
