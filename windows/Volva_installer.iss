; Inno Setup script for Volva 4.1 (https://jrsoftware.org/isinfo.php, free).
; Produces a single Volva_4.1.2_Setup.exe that installs per user (no admin
; rights needed), adds desktop and Start menu shortcuts and an uninstaller.
; build_windows.bat compiles this automatically if Inno Setup 6 is installed.
#define AppVer "4.1.2"
#ifndef DistDir
  #define DistDir "..\dist\Volva"
#endif

[Setup]
AppId={{6F1C2B9E-4A7D-4C1E-9B8A-56A0E0C4D411}
AppName=Volva
AppVersion={#AppVer}
AppVerName=Volva {#AppVer}
AppPublisher=University of Sheffield
DefaultDirName={localappdata}\Programs\Volva
DefaultGroupName=Volva
PrivilegesRequired=lowest
DisableProgramGroupPage=yes
OutputBaseFilename=Volva_{#AppVer}_Setup
OutputDir=.
SetupIconFile=..\volva\gui\volva.ico
UninstallDisplayIcon={app}\Volva.exe
Compression=lzma2
SolidCompression=yes
WizardStyle=modern

[Tasks]
Name: "desktopicon"; Description: "Create a desktop shortcut"; GroupDescription: "Shortcuts:"

[Files]
Source: "{#DistDir}\*"; DestDir: "{app}"; Flags: recursesubdirs ignoreversion

[Icons]
Name: "{group}\Volva {#AppVer}"; Filename: "{app}\Volva.exe"; WorkingDir: "{userdocs}"
Name: "{userdesktop}\Volva {#AppVer}"; Filename: "{app}\Volva.exe"; WorkingDir: "{userdocs}"; Tasks: desktopicon

[Run]
Filename: "{app}\Volva.exe"; Description: "Start Volva now"; Flags: nowait postinstall skipifsilent
