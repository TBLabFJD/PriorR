@echo off
rem Starts PriorR (Windows). Needs R installed: https://cran.r-project.org/bin/windows/base/
setlocal
set "RSCRIPT="
where Rscript >nul 2>nul && set "RSCRIPT=Rscript"
rem R installed for all users or only for this user (the last folder found is the newest version)
if not defined RSCRIPT (
  for /d %%D in ("%ProgramFiles%\R\R-*") do set "RSCRIPT=%%D\bin\Rscript.exe"
)
if not defined RSCRIPT (
  for /d %%D in ("%LOCALAPPDATA%\Programs\R\R-*") do set "RSCRIPT=%%D\bin\Rscript.exe"
)
if not defined RSCRIPT (
  echo R is not installed. Install it from https://cran.r-project.org/bin/windows/base/
  pause
  exit /b 1
)
"%RSCRIPT%" "%~dp0launch.R"
if errorlevel 1 pause
