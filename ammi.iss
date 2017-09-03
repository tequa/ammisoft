; -- Example1.iss --
; Demonstrates copying 3 files and creating an icon.

; SEE THE DOCUMENTATION FOR DETAILS ON CREATING .ISS SCRIPT FILES!


[Setup]
AppName=AMMISOFT
AppVersion=1.0
DefaultDirName={pf}\AMMISOFT
DefaultGroupName=AMMISOFT
UninstallDisplayIcon={app}\uninstall.exe
OutputBaseFilename=AMMISOFT-installer

[Files]
Source: "C:\Users\DMoran\Desktop\ammi\ammimain\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs
Source: "C:\Users\DMoran\Desktop\ammi\ammisoftextra\*"; DestDir: "{userdocs}\AMMISOFT"
Source: "C:\Users\DMoran\Desktop\ammi\ammiappdata\*"; DestDir: "{localappdata}\AMMISOFT"

[Icons]
Name: "{group}\AMMISOFT"; Filename: "{app}\WinPython-64bit-2.7.13.1Zero\python-2.7.13.amd64\pythonw.exe"; WorkingDir: "{app}"; Parameters: """{app}\WinPython-64bit-2.7.13.1Zero\python-2.7.13.amd64\Scripts\ammisoft.pyw"""; IconFilename: "{app}\WinPython-64bit-2.7.13.1Zero\settings\wheat.ico"
Name: "{commondesktop}\AMMISOFT"; Filename: "{app}\WinPython-64bit-2.7.13.1Zero\python-2.7.13.amd64\pythonw.exe"; WorkingDir: "{app}"; Parameters: """{app}\WinPython-64bit-2.7.13.1Zero\python-2.7.13.amd64\Scripts\ammisoft.pyw"""; IconFilename: "{app}\WinPython-64bit-2.7.13.1Zero\settings\wheat.ico"

