{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 7.4.2011
Purpose:
  Parser for sequence query data from the European Bioinformatics Institute
  in fasta format. Only reads sequences, not alignments
    http://www.ebi.ac.uk/
Requirements:
Revisions:
To do:
*******************************************************************************}

unit ebifastaparser;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, sequence, alignment, stringutils;

function ReadEBIFasta(Filename:string):TOCSequences;

implementation

function ReadEBIFasta(Filename:string):TOCSequences;

var
  sl:TStringList;
  f:Integer;
  s:string;
  tmp:TOCSequence;

procedure ParseDescription;

//assumes the following example:
//>sp|P08877|PTHP_BACSU Phosphocarrier protein HPr OS=Bacillus subtilis GN=ptsH PE=1 SV=3
//    *code* * id    ** *Description              *   *** organism ***  ** gene *`*
begin
  with tmp do
    begin
    SnipString(s,'|');
    Code:=SnipString(s,'|');
    ID:=SnipString(s,' ');
    Description:=SnipString(S,'OS=');
    Organism:=TrimmedBlanks(SnipString(S,'GN='));
    Gene:=TrimmedBlanks(SnipString(S,'PE='));
    Evidence:=TrimmedBlanks(SnipString(S,'SV='));
    end;
end;

begin
  Result:=nil;
  tmp:=EmptySequence;
  sl:=TStringList.Create;
  sl.LoadFromFile(FileName);
  for f:=0 to sl.Count-1 do
    begin
    s:=sl.Strings[f];
    if Pos('>',s)=1 then  //new sequence found
      begin
      if tmp.Sequence<>'' then
        AddSequence(tmp,Result);
      tmp:=EmptySequence;
      ParseDescription;
      end
    else tmp.Sequence:=tmp.Sequence+TrimmedBlanks(s);
    end;
  if tmp.Sequence<>'' then
    AddSequence(tmp,Result);
  sl.Free;
end;

end.

