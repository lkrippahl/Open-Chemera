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
  Classes, SysUtils, sequence, alignment, ocstringutils,fasta, lclproc;

function ReadEBIFasta(Filename:string):TOCSequences;

implementation

function ReadEBIFasta(Filename:string):TOCSequences;

var
  f:Integer;
  s:string;
  tmp:TOCSequence;
  freader:TFastaReader;

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
    if Pos('GN=',S)>0 then              //in some cases there is no gene info...
      begin
      Organism:=TrimmedBlanks(SnipString(S,'GN='));
      Gene:=TrimmedBlanks(SnipString(S,'PE='));
      end
    else
      begin
      Organism:=TrimmedBlanks(SnipString(S,'PE='));
      Gene:='';
      end;
    Evidence:=TrimmedBlanks(SnipString(S,'SV='));
    end;
end;

begin
  freader:=TFastaReader.Create(FileName);
  Result:=nil;
  tmp:=EmptySequence;
  for f:=0 to High(freader.Seqs) do
    begin
    s:=freader.Seqs[f].Id;
    ParseDescription;
    tmp.Sequence:=freader.Seqs[f].Sequence;
    AddSequence(tmp,Result);
    end;
  freader.Free;
end;

end.

