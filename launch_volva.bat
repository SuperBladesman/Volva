@echo off
REM Launch Völva 4.1 from its conda environment.
REM Edit CONDA_ROOT if Anaconda is installed elsewhere.
set CONDA_ROOT=%USERPROFILE%\anaconda3
call "%CONDA_ROOT%\Scripts\activate.bat" volva41
cd /d "%~dp0"
python -m volva
if errorlevel 1 pause
