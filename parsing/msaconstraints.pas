{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 23.6.2011
Purpose:
  creates constraint files for MSA adjustments with CaSPER
Requirements:
Revisions:
To do:
*******************************************************************************}

unit msaconstraints;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, alignment, clustalparser, basetypes, ocstringutils, debugutils;

type

  { TCPMSA }

  TCPMSA=class
    private
      FSubMatrix:TSubMatrix;
      FMSA:TMSA;
      procedure GenerateConstraints(FileName:string;Threshold:TOCFLoat);
      procedure GenerateDomains(FileName:string;Window:Integer);
      procedure GenerateQueries(FileName:string;PDBMap:TOCIntegers);
    public
      MinShiftable:Integer; //minimum number of shiftable lines in block
                            //default is 20% of MSA, set on AssignMSA
      procedure ReadSubMatrix(FileName:string);
      procedure ReadClustalMSA(FileName:string);
      procedure AssignMSA(AMSA:TMSA);
      procedure GenerateCP(FileName:string;PDBMap:TOCIntegers); //creates all files
      procedure GenerateBlocks(FileName:string);
  end;


implementation


{ TCPMSA }

procedure TCPMSA.ReadSubMatrix(FileName: string);
begin
  FSubMatrix:=ReadBLASTMatrix(FileName);
end;

procedure TCPMSA.ReadClustalMSA(FileName: string);
begin
  FMSA:=ReadClustal(FileName);
end;

procedure TCPMSA.AssignMSA(AMSA: TMSA);
begin
  FMSA:=AMSA;
  MinShiftable:=Round(Length(FMSA.Alignment)/5);
end;

procedure TCPMSA.GenerateCP(FileName: string;PDBMap:TOCIntegers);
begin
  GenerateConstraints(FileName+'.const',0);
  //GenerateDomains(FileName+'.doms',Window);
  if PDBMap<>nil then GenerateQueries(FileName+'.query',PDBMap);
  SaveMSAToFile(FMSA,FileName+'.msa');
  GenerateBlocks(FileName+'.blocks');
end;

procedure TCPMSA.GenerateConstraints(FileName: string; Threshold: TOCFLoat);

var
  buf:TStringList;
  f,g:Integer;

begin
  buf:=TStringList.Create;
  buf.Add('# Forbidden combinations in any column');
  for f:=0 to 18 do //only 20 regular AAs
    for g:=f+1 to 19 do
      if FSubMatrix.Matrix[f,g]<=Threshold then
            buf.Add(FSubMatrix.MonomerIndex[f+1]+FSubMatrix.MonomerIndex[g+1]);
  buf.SaveToFile(FileName);
  buf.Free;
end;

procedure TCPMSA.GenerateDomains(FileName: string;Window:Integer);
//TO DO: assumes query sequence is the first sequence in the MSA. Revise...


var
  buf:TStringList;
  domains:TOCStrings;
  maxlength:Integer;

procedure SetDomains(Ix:Integer);

function Combinations(Gaps,AAs:string):TOCStrings;

var
  s:string;
  tmp:TOCStrings;
  f:Integer;

begin
  Result:=nil;
  if Gaps<>'' then
    begin
    s:=Gaps[1];
    tmp:=Combinations(Copy(Gaps,2,Length(Gaps)),AAs);
    if tmp<>nil then
      for f:=0 to High(tmp) do
        AddUniqueToArray(s+tmp[f],Result);
    end;
  if AAs<>'' then
    begin
    s:=AAs[1];
    tmp:=Combinations(Gaps,Copy(AAs,2,Length(AAs)));
    for f:=0 to High(tmp) do
      AddUniqueToArray(s+tmp[f],Result);
    end;
  //check if no recursive calls returned results
  //In that case there is only one non empty (gaps or aas), at most
  if (Result=nil) and ((Gaps+AAs)<>'') then
    AddToArray(Gaps+AAs,Result)
end;

var
  f,g:Integer;
  s,gaps,aas:string;
  tmpvals:TOCStrings;

begin
  maxlength:=0;

  for f:=0 to High(FMSA.Alignment) do
    begin
    //split residues and gaps to generate domains recursively
    s:=Copy(FMSA.Alignment[f],Ix-Window,2*Window+1);
    aas:='';
    gaps:='';
    for g:=1 to Length(s) do
      if s[g]=FMSA.GapMarker then gaps:=gaps+s[g]
      else aas:=aas+s[g];
    tmpvals:=Combinations(gaps,aas);
    domains[f]:=FlattenStrings(tmpvals,' ');
    if Length(domains[f])>maxlength then
      maxlength:=Length(domains[f]);
    end;
end;

var
  f,g:Integer;
  s:string;

begin
  buf:=TStringList.Create;
  s:=FMSA.Alignment[0];
  SetLength(domains,Length(FMSA.Alignment));
  for f:=Window+1 to Length(s)-Window do
    if s[f]<>FMSA.GapMarker then
      begin
      SetDomains(f);
        buf.Add('Center Position: '+IntToStr(f));
        for g:=0 to High(domains) do
          buf.Add(domains[g]);
      end;
  buf.SaveToFile(FileName);
  buf.Free;
end;

procedure TCPMSA.GenerateQueries(FileName: string;PDBMap:TOCIntegers);

//assumes query is first on MSA

var
  buf:TStringList;
  f,ix:Integer;
  s:string;

begin
  buf:=TStringList.Create;
  buf.Add('Positions to test (first position is number 1):');
  ix:=-1;
  s:=FMSA.Alignment[0];
  for f:=1 to Length(s) do
    if s[f]<>FMSA.GapMarker then
      begin
      Inc(ix);
      if PDBMap[ix]>=0 then
        buf.Add(IntToStr(f));
      end;
  buf.SaveToFile(FileName);
  buf.Free;
end;

procedure TCPMSA.GenerateBlocks(FileName: string);

const

  Window=5;
var
  blankcounts:TOCIntegers;
  sl:TStringList;

function IsAnchor(Start,Len:Integer):Boolean;

//checks if this block cannot be changed

var
  f,cchar,countshift:Integer;
  s:string;

begin
  Result:=True;
  countshift:=0;
  for f:=1 to High(FMSA.Alignment) do //assumes 0 is fixed query
      begin
      s:=Copy(FMSA.Alignment[f],Start,Len);
      cchar:=CountInString(s,FMSA.GapMarker);
      if (cchar>0) and (cchar<Length(s)) then
        begin
        Inc(countshift);
        if countshift>=MinShiftable then
          Exit(False);
        end;
      end;
end;

procedure AddToList(Start,Count:Integer);

var
  g,f:Integer;
  pref,suf:string;

begin
  if Start<1 then
    begin
      pref:=FMSA.GapMarker;
      Start:=1;
    end
    else pref:='';
  if Count+Start>Length(FMSA.Alignment[0]) then
    begin
      suf:=FMSA.GapMarker;
      Count:=Length(FMSA.Alignment[0])-Start-1;
    end
    else suf:='';
  if Count>0 then
    begin
    sl.Add('***Block***');
    sl.Add('Start:'+IntToStr(Start));
    sl.Add('Count:'+IntToStr(Count));
    for f:=0 to High(FMSA.Alignment) do
      sl.Add(pref+Copy(FMSA.Alignment[f],Start,Count)+suf);
    end;
end;

procedure BuildList;

//TO DO: messy code...

var start,finish:Integer;

begin
  start:=1;
  while start<=Length(FMSA.Alignment[0]) do
    begin
    //seek to first block that is not anchor
    while (start<=Length(FMSA.Alignment[0])-Window) and IsAnchor(start,Window) do
      Inc(start);
    if start>1 then start:=start+Window-1; //at least one block before failing
    finish:=start+Window;
    while (finish<=Length(FMSA.Alignment[0])-Window) and
      not IsAnchor(finish,Window) do
      Inc(finish);
    AddToList(start-1,finish+2-start);
    start:=finish+1;
    end;
end;

begin
  if FMSA.Alignment=nil then Exit;
  sl:=TStringList.Create;
  BuildList;
  sl.SaveToFile(FileName);
  sl.Free;
end;

end.

