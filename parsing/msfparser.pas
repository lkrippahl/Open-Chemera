{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 23.6.2011
Purpose:
  Parser for MSF multiple sequence alignment files.
Requirements:
Revisions:
To do:
  This is a quick fix for Balibase. Redo the parsing according to the correct
  format (see, eg, http://tcoffee.vital-it.ch/Doc/doc3.html)
*******************************************************************************}
unit msfparser;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, alignment, ocstringutils;


function ReadMSF(fromfile:string):TMSA;
procedure SaveMSF(const ToFile:string; const MSA:TMSA);

implementation

function ReadMSF(fromfile: string): TMSA;

var
  sl:TStringList;

  function ParseHeader:Integer;
  //returns the index of the first line after the header
  //or 0 if there is no header
  //TO DO: this is only skipping the header, if it exists

  var f:Integer;

  begin
    Result:=0;
    for f:=0 to sl.Count-1 do
      if Pos('//',sl.Strings[f])>0 then
        begin
        Result:=f+1;
        Break;
        end;
  end;

var
  f,Ix,g:Integer;
  nam,seq,s,accesscode:string;


begin
  sl:=TStringList.Create;
  sl.LoadFromFile(fromfile);
  Result.GapMarker:='-';
  for f:=ParseHeader to sl.Count-1 do   // Process header
    if sl.Strings[f]<>''  then          // skip blank lines
    begin
    s:=sl.Strings[f];
    nam:=SnipString(s,' ');  // Assume a space
                             // between description and sequence
                             // TO DO: header should have the sequence
                             // identifiers
    if nam<>'' then
      begin
      seq:=Deblank(s);
      for g:=1 to Length(seq) do
        //TO DO: improve this to convert all gaps used in mds
        if seq[g]='.' then seq[g]:=Result.GapMarker;

      Ix:=LastIndexOf(nam,Result.SequenceIDs);
      if Ix<0 then
        begin
        AddToArray(nam,Result.SequenceIDs);
        AddToArray(seq,Result.Alignment);
        end
      else
        Result.Alignment[Ix]:=Result.Alignment[Ix]+seq;
      end;
    end;

  sl.Free;
end;


procedure SaveMSF(const ToFile:string; const MSA:TMSA);

procedure WriteHeader(Sl:TStringList);

var
  f:Integer;
  s,t:string;

begin
  t:=IntToStr(Length(MSA.Alignment[0]));
  Sl.Add('MSF: '+t+'Type: P  Check: 0 ..');
  Sl.Add('');
  for f:=0 to High(MSA.SequenceIDs) do
    begin
    Sl.Add('Name: '+MSA.SequenceIDs[f]+' Length: '+t+' Check: 0 Weight: 1.00');
    end;
  Sl.Add('//');
end;

function MaxIdLength:Integer;

var f:Integer;

begin
  Result:=6;
  for f:=0 to High(MSA.SequenceIds) do
    if Length(MSA.SequenceIds[f])>Result then
      Result:=Length(MSA.SequenceIds[f]);
end;

var
  sl:TStringList;
  f,g,cix:Integer;
  idlen,seqlen:Integer;
  s:string;

begin
  if MSA.Alignment=nil then
    Exit;
  sl:=TStringList.Create;
  WriteHeader(sl);
  idlen:=MaxIdLength+1; //length of string for id, including spaces at the end
  seqlen:=Length(MSA.Alignment[0]);
  cix:=1;              //current column in alignment
  while cix<=seqlen do
    begin
    for f:=0 to High(MSA.Alignment) do
      begin
      s:=AppendToLength(MSA.SequenceIDs[f],' ',idlen);
      for g:=0 to 4 do // 5 blocks of 10
        s:=s+' '+Copy(MSA.Alignment[f],cix+10*g,10);
      sl.Add(s);
      end;
    sl.Add('');
    cix:=cix+50;       // 50 per line
    end;
  sl.SaveToFile(ToFile);
  sl.Free;
end;

end.

