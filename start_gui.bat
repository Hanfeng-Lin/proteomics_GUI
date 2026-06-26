@echo off
setlocal
cd /d "%~dp0"
title DIA-NN Volcano Explorer

REM ==================================================================
REM  1) Locate conda and the "proteomics" environment (clear errors,
REM     never silently close -- every error path pauses).
REM ==================================================================
set "CONDA_ROOT="
set "CONDA_BASE="
for %%R in (
  "%USERPROFILE%\anaconda3"
  "%USERPROFILE%\miniconda3"
  "%USERPROFILE%\Anaconda3"
  "%USERPROFILE%\Miniconda3"
  "%LOCALAPPDATA%\anaconda3"
  "%LOCALAPPDATA%\miniconda3"
  "C:\ProgramData\anaconda3"
  "C:\ProgramData\miniconda3"
  "C:\ProgramData\Anaconda3"
  "C:\ProgramData\Miniconda3"
) do (
  if not defined CONDA_BASE if exist "%%~R\Scripts\activate.bat" set "CONDA_BASE=%%~R"
  if not defined CONDA_ROOT if exist "%%~R\envs\proteomics\python.exe" set "CONDA_ROOT=%%~R"
)

REM Fallback: ask a conda on PATH for its base install (no output if conda absent).
if not defined CONDA_ROOT (
  for /f "delims=" %%i in ('conda info --base 2^>nul') do (
    if not defined CONDA_BASE set "CONDA_BASE=%%i"
    if exist "%%i\envs\proteomics\python.exe" set "CONDA_ROOT=%%i"
  )
)

if not defined CONDA_ROOT (
  echo ============================================================
  if not defined CONDA_BASE (
    echo [ERROR] Anaconda/Miniconda is not installed, or could not be found.
    echo.
    echo Install Anaconda from https://www.anaconda.com/download
    echo ^(or Miniconda^), then run this script again.
  ) else (
    echo [ERROR] Conda was found but the "proteomics" environment is missing.
    echo   conda location: %CONDA_BASE%
    echo.
    echo Create the environment once with:
    echo     conda env create -f "%~dp0environment.yml"
  )
  echo ============================================================
  echo.
  pause
  exit /b 1
)
echo Found "proteomics" environment under: %CONDA_ROOT%

REM ==================================================================
REM  2) Update check via the VERSION file (best-effort; never blocks launch).
REM     Compares the local VERSION to the one on GitHub (works for git and ZIP).
REM ==================================================================
set "LOCAL_VER="
for /f "usebackq delims=" %%V in (`powershell -NoProfile -Command "if(Test-Path 'VERSION'){(Get-Content -Raw 'VERSION').Trim()}"`) do set "LOCAL_VER=%%V"
set "REMOTE_VER="
for /f "usebackq delims=" %%V in (`powershell -NoProfile -Command "try{$ProgressPreference='SilentlyContinue';[Net.ServicePointManager]::SecurityProtocol=[Net.SecurityProtocolType]::Tls12;(Invoke-RestMethod -UseBasicParsing -TimeoutSec 6 'https://raw.githubusercontent.com/Hanfeng-Lin/proteomics_GUI/main/VERSION').Trim()}catch{}"`) do set "REMOTE_VER=%%V"

if not defined REMOTE_VER goto :launch
if not defined LOCAL_VER set "LOCAL_VER=(unknown)"
if /i "%LOCAL_VER%"=="%REMOTE_VER%" goto :launch

echo.
echo ------------------------------------------------------------
echo A newer version of proteomics_GUI is available on GitHub.
echo     installed: %LOCAL_VER%
echo     latest   : %REMOTE_VER%
echo ------------------------------------------------------------
set "DOUPD="
set /p DOUPD="Update now before launching? [Y/N] "
if /i "%DOUPD%"=="Y" (
  call "%~dp0update.bat"
  echo.
  echo Update finished. Please run start_gui.bat again to launch the updated app.
  pause
  exit /b 0
)

:launch
REM ==================================================================
REM  3) Activate the environment and launch the GUI.
REM ==================================================================
echo Activating...
call "%CONDA_ROOT%\Scripts\activate.bat" proteomics
if errorlevel 1 (
  echo [ERROR] Failed to activate the "proteomics" environment.
  pause
  exit /b 1
)

cd /d "%~dp0"
echo Starting the Volcano Explorer GUI...
python gui.py
if errorlevel 1 (
  echo.
  echo [ERROR] The GUI exited with an error. See the messages above.
  pause
  exit /b 1
)
