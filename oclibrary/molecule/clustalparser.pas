{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 9.1.2011
Purpose:
  Parser for clustall MSA files.
Requirements:
Revisions:
To do:
*******************************************************************************}

unit clustalparser;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, alignment, stringutils;


function ReadClustal(fromfile:string):TMSA;
procedure WriteClustal(const ToFile:string; const MSA:TMSA);

implementation

function ReadClustal(fromfile: string): TMSA;

var
  sl:TStringList;
  f,Ix:Integer;
  nam,seq,s,accesscode:string;

begin
  sl:=TStringList.Create;
  sl.LoadFromFile(fromfile);
  Result.GapMarker:='-';
  for f:=1 to sl.Count-1 do   // first line is header, "CLUSTAL"
    begin
    s:=sl.Strings[f];
    nam:=SnipString(s,' ');   // Clustal files should have a space
                              // between description and sequence
    if nam<>'' then
      begin
      seq:=TrimmedBlanks(s);  //check if there is a number at the end of the line
      if Pos(' ',seq)>0 then    //separated by a space (sequence number)
        seq:=SnipString(seq,' ');

      if Pos('/',nam)>0 then   // EBI clustal has sequence info in description
        nam:=SnipString(nam, '/');

      while Pos('|',nam)>0 do  // EBI clustal has accession code etc before ide
        accesscode:=SnipString(nam,'|');
      if accesscode<>'' then   // compose name with accession code to prevent duplicates
        nam:=accesscode+'|'+nam;
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

procedure WriteClustal(const ToFile:string; const MSA:TMSA);

var
  sl:TStringList;

begin

end;

end.

