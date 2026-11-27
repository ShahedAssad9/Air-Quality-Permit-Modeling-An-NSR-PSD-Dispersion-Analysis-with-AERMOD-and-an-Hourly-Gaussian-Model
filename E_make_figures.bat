@echo off
title AERMOD project - step E: result tables and figures (Python)
cd /d "%~dp008_Python_Analysis"
if not exist "%~dp0logs" mkdir "%~dp0logs"
set "LOG=%~dp0logs\E_make_figures.log"
set PYTHONIOENCODING=utf-8
set PY=
where python >nul 2>nul && set PY=python
if not defined PY where py >nul 2>nul && set PY=py
if not defined PY goto nopython
%PY% -c "import numpy, matplotlib" >nul 2>nul
if errorlevel 1 goto nopackages

echo Step E started %date% %time% > "%LOG%"
echo Step E: making the result tables and figures. Progress is written to logs\E_make_figures.log
echo.
echo 1/4  concept illustrations
%PY% make_illustrations.py >> "%LOG%" 2>&1 || goto failed
echo 2/4  Part A with the synthetic demo weather - about 1 to 2 minutes
%PY% run_part_a.py >> "%LOG%" 2>&1 || goto failed
echo 3/4  AERMOD design values and regulatory tests
%PY% aermod_results.py >> "%LOG%" 2>&1 || goto failed
echo 4/4  Part A against AERMOD - the first time, this also runs Part A on the AERMET data
%PY% compare_with_aermod.py >> "%LOG%" 2>&1 || goto failed
echo Step E finished OK %date% %time% >> "%LOG%"
echo.
echo Step E finished. The figures are in the images folder. You can close this window.
pause >nul
exit /b 0

:failed
echo FAILED >> "%LOG%"
echo.
echo Something went wrong. Read logs\E_make_figures.log
pause >nul
exit /b 1

:nopython
echo Python was not found. Install Python 3.9 or newer from python.org, then run this file again.
pause >nul
exit /b 1

:nopackages
echo Python needs numpy and matplotlib. Open a command window in this folder and run:
echo     %PY% -m pip install -r requirements.txt
pause >nul
exit /b 1
