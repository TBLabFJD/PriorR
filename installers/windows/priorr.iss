; Installer of PriorR for Windows (Inno Setup 6: https://jrsoftware.org/isinfo.php).
; Built by .github/workflows/installers.yml:  iscc installers\windows\priorr.iss
; Installs for the current user (no administrator needed), so PriorR can install its R packages
; next to it. Asks for the shared folder of the lab and writes priorr.cfg. R is not installed:
; when it is missing the installer says so and offers to open the download page.

#define AppVersion "4.2"
#define Root "..\.."

[Setup]
AppId={{6F2C1E5A-4B7D-4F1E-9C3A-8E2D5B7A1C40}
AppName=PriorR
AppVersion={#AppVersion}
AppVerName=PriorR v{#AppVersion}
AppPublisher=Bioinformatics Unit - Genetics and Genomics Department, IIS-Fundación Jiménez Díaz
AppPublisherURL=https://github.com/TBLabFJD/PriorR
DefaultDirName={localappdata}\Programs\PriorR
DisableProgramGroupPage=yes
PrivilegesRequired=lowest
OutputDir=Output
OutputBaseFilename=PriorR-{#AppVersion}-setup
SetupIconFile=priorr.ico
UninstallDisplayIcon={app}\priorr.ico
LicenseFile=..\LICENSE.txt
WizardStyle=modern
Compression=lzma2
SolidCompression=yes

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"
Name: "spanish"; MessagesFile: "compiler:Languages\Spanish.isl"

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"

[Files]
Source: "{#Root}\launch.R"; DestDir: "{app}"; Flags: ignoreversion
Source: "{#Root}\PriorR.bat"; DestDir: "{app}"; Flags: ignoreversion
Source: "{#Root}\priorr.cfg.example"; DestDir: "{app}"; Flags: ignoreversion
Source: "{#Root}\README.md"; DestDir: "{app}"; Flags: ignoreversion
Source: "priorr.ico"; DestDir: "{app}"; Flags: ignoreversion
Source: "{#Root}\app\*"; DestDir: "{app}\app"; Flags: ignoreversion recursesubdirs createallsubdirs; Excludes: "cache\*,shiny_bookmarks\*"

[Icons]
Name: "{userprograms}\PriorR"; Filename: "{app}\PriorR.bat"; WorkingDir: "{app}"; IconFilename: "{app}\priorr.ico"
Name: "{userdesktop}\PriorR"; Filename: "{app}\PriorR.bat"; WorkingDir: "{app}"; IconFilename: "{app}\priorr.ico"; Tasks: desktopicon

[Run]
Filename: "{app}\PriorR.bat"; Description: "Start PriorR now"; Flags: postinstall nowait skipifsilent shellexec; Check: RInstalled

[Code]
var
  SharedPage: TInputOptionWizardPage;
  FolderPage: TInputDirWizardPage;

function RInstalled(): Boolean;
var
  Found: TFindRec;
begin
  Result := RegKeyExists(HKLM, 'SOFTWARE\R-core\R') or RegKeyExists(HKCU, 'SOFTWARE\R-core\R') or
            RegKeyExists(HKLM, 'SOFTWARE\R-core\R64') or RegKeyExists(HKCU, 'SOFTWARE\R-core\R64');
  if not Result then
    Result := FindFirst(ExpandConstant('{commonpf}\R\R-*'), Found) or
              FindFirst(ExpandConstant('{localappdata}\Programs\R\R-*'), Found);
end;

procedure InitializeWizard();
begin
  SharedPage := CreateInputOptionPage(wpSelectDir, 'Data of the lab',
    'Where does PriorR keep the users, comments, classifications and profiles?',
    'Every computer of the lab must use the same shared folder to see the same users and comments.',
    True, False);
  SharedPage.Add('Only this computer (a folder "shared" inside PriorR)');
  SharedPage.Add('A shared folder of the lab (network folder)');
  SharedPage.SelectedValueIndex := 0;

  FolderPage := CreateInputDirPage(SharedPage.ID, 'Shared folder of the lab',
    'Choose the shared folder of the lab', 'For example \\server\genetics\priorr_data or Z:\priorr_data.',
    False, '');
  FolderPage.Add('');
end;

function ShouldSkipPage(PageID: Integer): Boolean;
begin
  Result := (PageID = FolderPage.ID) and (SharedPage.SelectedValueIndex = 0);
end;

procedure CurStepChanged(CurStep: TSetupStep);
var
  Cfg, Folder: String;
  ErrorCode: Integer;
begin
  if CurStep = ssPostInstall then
  begin
    Cfg := ExpandConstant('{app}\priorr.cfg');
    // an update keeps the configuration of the lab
    if not FileExists(Cfg) then
    begin
      if SharedPage.SelectedValueIndex = 0 then
        Folder := 'shared'
      else
        Folder := FolderPage.Values[0];
      SaveStringToFile(Cfg, '# PriorR configuration of the lab (see priorr.cfg.example)' + #13#10 +
                            'shared_folder: ' + Folder + #13#10, False);
    end;
    if not RInstalled() then
      if MsgBox('PriorR needs R, which is not installed on this computer.' + #13#10#13#10 +
                'Install R from CRAN and then start PriorR. Open the download page now?',
                mbConfirmation, MB_YESNO) = IDYES then
        ShellExec('open', 'https://cran.r-project.org/bin/windows/base/', '', '', SW_SHOWNORMAL, ewNoWait, ErrorCode);
  end;
end;
