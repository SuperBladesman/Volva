@echo off
REM Double-click to put a Volva shortcut on the desktop and in the Start menu.
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0Create_Desktop_Shortcut.ps1"
pause
