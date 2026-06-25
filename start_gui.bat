@echo off
setlocal
cd /d "%~dp0"
title DIA-NN Volcano Explorer

REM ------------------------------------------------------------------
REM  Find a conda installation that contains the "proteomics" env.
REM ------------------------------------------------------------------
set "CONDA_ROOT="
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
  if exist "%%~R\envs\proteomics\python.exe" (
    set "CONDA_ROOT=%%~R"
    goto :found
  )
)

REM Fallback: ask conda on PATH for its base install, then check there.
where conda >nul 2>nul
if %errorlevel%==0 (
  for /f "delims=" %%i in ('conda info --base 2^>nul') do (
    if exist "%%i\envs\proteomics\python.exe" set "CONDA_ROOT=%%i"
  )
)

:found
if not defined CONDA_ROOT (
  echo [ERROR] Could not find a conda environment named "proteomics".
  echo Looked in the common Anaconda/Miniconda locations and on PATH.
  echo If your conda is installed elsewhere, add its path to this script.
  echo.
  pause
  exit /b 1
)

echo Found "proteomics" environment under: %CONDA_ROOT%
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
