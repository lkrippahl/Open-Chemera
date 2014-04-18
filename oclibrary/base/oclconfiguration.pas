{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 12.3.2011
Purpose:
  Configuration files and folders. Also handles basic molecule data.
  Data stored in the configuration folders (for data used by several
  applications; data for specific applications should be handled by
  the application and kept in appfolder)
Requirements:
  Tested on Windows7 and Kubuntu
  Files:
    atomdata.txt:
      tab separated file with element symbols and atomic numbers
      ** more efficient if more common elements are first **

Revisions:
To do:
  perhaps function AAOneLetterCode(TLC:string):string; should not be here
    (in sequence?)
*******************************************************************************}

unit oclconfiguration;

{$mode objfpc}{$H+}

interface

uses
  Interfaces, Classes, SysUtils, Forms, FileUtil, basetypes, stringutils,LCLProc;

type
  TOCLConfig=record
    OCLPath, AppPath, AppConfig, OCLConfig, MonomersPath:string;
    DefaultAtomicRadius:TFloat;
  end;
  TOCLAtomData=record
    //elements sorted by atomic number (= to index +1)
    Symbol,Name:string;
    CovalentRadius,VdWRadius,UnitedRadius:TFloat;
    CPKColor:TCoord;
  end;
  TOCLAAData=record
    Name,TLCode:string;
    OLCode,Category:Char;
    Hydropathy:Single;
  end;

var
  Config:TOCLConfig;
  AtomData:array of TOCLAtomData;
  AAData:array of TOCLAAData;
const
  AtomDataFile='atomdata.txt';
  AADataFile='aminoaciddata.txt';

procedure LoadAtomData;
procedure LoadAAData;
function AtomicNumber(Symbol:string):Integer;
function AAIndex(Code:string):Integer;
function AAOneLetterCode(TLC:string):string;

implementation

procedure RemoveComments(var Sl:TStringList);

var
  f,ix:Integer;
  s:string;

begin
  f:=0;
  while f<Sl.Count do
    begin
    ix:=Pos('#',Sl.Strings[f]);
    if ix<1 then Inc(f)
    else
    if ix=1 then Sl.Delete(f)
    else
      begin
      s:=Sl.Strings[f];
      Delete(s,ix,Length(s));
      Sl.Strings[f]:=s;
      Inc(f);
      end;
    end;
end;

procedure DefaultConfig;

var exename:string;

begin
  DecimalSeparator:='.';
    { TODO : Improve this }

  exename:=ChangeFileExt(ExtractFileName(Application.Exename),'');
  with Config do
    begin
    OCLPath:=GetAppConfigDir(False);
    if Pos(exename,OCLPath)=Length(OCLPath)-Length(exename) then
      Delete(OCLPath,Length(OCLPath)-Length(exename),Length(exename)+1);

    OCLPath:=AppendPathDelim(OCLPath+'oclibrary');
    OCLConfig:=OCLPath+'oclibrary.ini';

    AppPath:=AppendPathDelim(OCLPath+exename);
    AppConfig := AppPath+ChangeFileExt(exename,'.ini');
    MonomersPath:=AppendPathDelim(OCLPath+'monomers');

    //TODO: this should go into a file
    DefaultAtomicRadius:=1.5;

    end;
end;

procedure LoadAtomData;

var
  sl:TStringList;
  f:Integer;
  tmp:TSimpleStrings;

begin
  try
  AtomData:=nil;
  sl:=TStringList.Create;
  sl.LoadFromFile(Config.OCLPath+AtomDataFile);
  RemoveComments(sl);
  SetLength(AtomData,sl.Count);
  for f:=0 to sl.Count-1 do
    begin
    tmp:=SplitString(sl.Strings[f],#9);
    with AtomData[f] do
      begin
      Symbol:=tmp[0];
      CovalentRadius:=StringToFloat(tmp[1]);
      VdWRadius:=StringToFloat(tmp[2]);
      UnitedRadius:=StringToFloat(tmp[3]);
      Name:=tmp[4];
      CPKColor:=Coord(StrToInt(tmp[5])/255, StrToInt(tmp[6])/255, StrToInt(tmp[7])/255);
      end
    end;
  sl.Free;
  except
  sl.Free;
  raise Exception.Create('Error loading '+Config.OCLPath+AtomDataFile)
  end;
end;

procedure LoadAAData;
var
  sl:TStringList;
  f:Integer;
  tmp:TSimpleStrings;

begin
  try
  AAData:=nil;
  sl:=TStringList.Create;
  sl.LoadFromFile(Config.OCLPath+AADataFile);
  RemoveComments(sl);
  SetLength(AAData,sl.Count);
  for f:=0 to sl.Count-1 do
    begin
    tmp:=SplitString(sl.Strings[f],#9);
    with AAData[f] do
      begin
      Name:=tmp[0];
      TLCode:=tmp[1];
      OLCode:=tmp[2,1];
      Category:=tmp[3,1];
      Hydropathy:=StringToFloat(tmp[4]);
      end
    end;
  sl.Free;
  except
  sl.Free;
  raise Exception.Create('Error loading '+Config.OCLPath+AADataFile)
  end;
end;

function AtomicNumber(Symbol: string): Integer;

begin
  Result:=High(AtomData);
  while (Result>=0) and (AtomData[Result].Symbol<>Symbol) do
    Dec(Result);
  if Result>=0 then Inc(Result);
end;

function AAIndex(Code: string): Integer;
begin
  Result:=High(AAData);
  if Length(Code)=3 then
    while (Result>=0) and (AAData[Result].TLCode<>Code) do
      Dec(Result)
  else
    while (Result>=0) and (AAData[Result].OLCode<>Code) do
      Dec(Result);
end;

function AAOneLetterCode(TLC: string): string;

var ix:Integer;

begin
  ix:=AAIndex(TLC);
  if ix>=0 then
    Result:=AAData[ix].OLCode
  else Result:='';
end;

initialization
DefaultConfig;

end.

