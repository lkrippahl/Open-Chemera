{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 12.3.2011
Purpose:
  Configuration files and folders. Multiplatform support
Requirements:
Revisions:
To do:
*******************************************************************************}

unit chemeraconfig;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Forms, FileUtil;

type
  TChemConfig=record
    ConfigFile, DataFolder, LigandsFolder:string
  end;

var
  Config:TChemConfig;

implementation

procedure DefaultConfig;

begin
  with Config do
    begin
    {$ifdef win32}
    ConfigFile := ExtractFilePath(Application.EXEName) + 'chemera.ini';
    DataFolder:=ExtractFilePath(Application.EXEName)+'data';
    {$endif}
    {$ifdef Unix}
    ConfigFile := GetAppConfigFile(False);
    DataFolder:=GetAppConfigDir(False);
    {$endif}
    LigandsFolder:=AppendPathDelim(DataFolder)+'ligands'
    end;
end;

initialization
DefaultConfig;

end.

