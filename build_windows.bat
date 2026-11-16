@echo off
setlocal
REM =====================================================================
REM  Build the stand-alone Windows version of Volva 4.1.
REM  Run from Anaconda Prompt, inside the volva41 environment:
REM      conda activate volva41
REM      build_windows.bat
REM  Output goes to %USERPROFILE%\Volva_build (outside Google Drive):
REM      dist\Volva\Volva.exe               the application folder
REM      dist\Volva\Create desktop shortcut.bat
REM      Volva_4.1.3_Setup.exe              installer (only if Inno Setup 6 is installed)
REM =====================================================================
cd /d "%~dp0"
set OUT=%USERPROFILE%\Volva_build

echo [1/3] Running tests...
python -m pytest -q || (echo Tests failed - not building. & exit /b 1)

echo [2/3] Building application folder (takes a few minutes)...
REM Build with a clean PATH so no DLLs from other Python/Qt installs (e.g. base Anaconda) are bundled
set "PATH=%CONDA_PREFIX%;%CONDA_PREFIX%\Library\bin;%CONDA_PREFIX%\Scripts;%SystemRoot%\System32;%SystemRoot%"
pyinstaller --noconfirm --clean --distpath "%OUT%\dist" --workpath "%OUT%\work" volva.spec || (echo PyInstaller failed. & exit /b 1)
copy /y "windows\Create_Desktop_Shortcut.ps1" "%OUT%\dist\Volva\" >nul
copy /y "windows\Create desktop shortcut.bat" "%OUT%\dist\Volva\" >nul
REM The work folder holds an intermediate Volva.exe that cannot run on its own; remove it to avoid confusion
rmdir /s /q "%OUT%\work"

echo [3/3] Installer...
set ISCC=
if exist "%ProgramFiles(x86)%\Inno Setup 6\ISCC.exe" set ISCC=%ProgramFiles(x86)%\Inno Setup 6\ISCC.exe
if exist "%LOCALAPPDATA%\Programs\Inno Setup 6\ISCC.exe" set ISCC=%LOCALAPPDATA%\Programs\Inno Setup 6\ISCC.exe
if defined ISCC (
    "%ISCC%" /Q /DDistDir="%OUT%\dist\Volva" /O"%OUT%" "windows\Volva_installer.iss" && echo Installer: %OUT%\Volva_4.1.3_Setup.exe
) else (
    echo Inno Setup 6 not found - skipping the installer. The app folder is ready to use.
)

echo.
echo Done. Test it by double-clicking %OUT%\dist\Volva\Volva.exe
explorer "%OUT%\dist\Volva"
endlocal
