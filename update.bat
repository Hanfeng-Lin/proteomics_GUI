@echo off
setlocal
cd /d "%~dp0"
title proteomics_GUI - Update

set "REPO=https://github.com/Hanfeng-Lin/proteomics_GUI"
set "ZIPURL=%REPO%/archive/refs/heads/main.zip"

REM If this is a git clone and git is available, prefer a clean git pull.
where git >nul 2>nul
if %errorlevel%==0 if exist ".git" (
  echo Updating via git pull...
  git pull
  if errorlevel 1 (
    echo.
    echo [ERROR] git pull failed. If you have local edits, run one and retry:
    echo   git stash            keep your changes
    echo   git checkout -- .    discard them
    pause
    exit /b 1
  )
  echo.
  echo Up to date. Launch with start_gui.bat.
  pause
  exit /b 0
)

REM Otherwise download and apply the latest ZIP.
set "TMP=%TEMP%\proteomics_GUI_update"
set "ZIP=%TMP%\latest.zip"
if exist "%TMP%" rmdir /s /q "%TMP%"
mkdir "%TMP%"

echo Downloading the latest version from GitHub...
powershell -NoProfile -Command "$ProgressPreference='SilentlyContinue'; [Net.ServicePointManager]::SecurityProtocol=[Net.SecurityProtocolType]::Tls12; Invoke-WebRequest -UseBasicParsing -Uri '%ZIPURL%' -OutFile '%ZIP%'"
if not exist "%ZIP%" (
  echo [ERROR] Download failed. Check your internet connection, or download manually:
  echo   %ZIPURL%
  pause
  exit /b 1
)

echo Extracting...
powershell -NoProfile -Command "Expand-Archive -Force -LiteralPath '%ZIP%' -DestinationPath '%TMP%'"
set "SRC=%TMP%\proteomics_GUI-main"
if not exist "%SRC%" (
  echo [ERROR] Unexpected archive layout - update aborted.
  pause
  exit /b 1
)

echo Updating files ^(your data files and proteomics_GUI_output are kept^)...
REM Copy everything except this running script ^(replacing it mid-run can corrupt it^).
robocopy "%SRC%" "%~dp0." /E /XF "%~nx0" /NFL /NDL /NJH /NJS /NP >nul

rmdir /s /q "%TMP%"
echo.
echo Updated to the latest version. Launch with start_gui.bat.
echo ^(If update.bat itself changed upstream, re-download the ZIP once.^)
pause
