{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 9.1.2011
Purpose:
  Basic sequence data structures and utils
  Based on UniProtKB and EBI results
Requirements:
Revisions:
To do:
*******************************************************************************}

unit sequence;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes;


type
  //these fields are based on the EBI results
  TOCSequence=record
    Sequence:string;
    Organism:string;
    Database:string;    // e.g. uniprot
    ID:string;          //protein identifier in the database. E.g. Q6AW43_9ANNE
    Code:string;         //protein accession code, e.g. Q6AW43

    Evidence:string;    //UniProtKB, evidence level, can be a number 1-5
                        //see http://www.uniprot.org/manual/protein_existence
                        //TODO: should it always be a number?
    Gene:string;
    Description:string;

    //these are used when the set of sequences results from database queries

     Score, Bits, Expectation, Identity, Positives:TFloat;
    end;
  TOCSequences=array of TOCSequence;
  TOCSequencesSet=array of TOCSequences;

  TOCCompressedSequenceRec=record
    Symbol:Char;
    First,Last:Integer;
  end;
  TOCCompressedSequence=array of TOCCompressedSequenceRec;
  TOCCompressedSequences=array of TOCCompressedSequence;



procedure AddSequence(se:TOCSequence; var ss:TOCSequences);
function EmptySequence:TOCSequence;
function ListOrganisms(const Sequences:TOCSequences):TSimpleStrings;
  //Lists all different organisms, by order of occurrence in the sequence array
function FilterByCommonOrganisms(SequencesSet:TOCSequencesSet):TOCSequencesSet;
  //the resulting set contains the same set of sequence arrays but with each array
  //containing only thoses sequences from organisms present in all sequence arrays,
  //and only the first of each.
function FirstIxByOrganism(const Sequences:TOCSequences;Organism:string):Integer;

  //Compressed sequences store repeated letters as first..last segments
  //These are mostly useful for storing alignment columns, as these often have repeats
function CompressSequence(Sequence:string):TOCCompressedSequence;
function CompressSequences(Sequences:TOCSequences):TOCCompressedSequences;
function CompressSequences(Sequences:TSimpleStrings):TOCCompressedSequences;


implementation

procedure AddSequence(se: TOCSequence; var ss: TOCSequences);

var l:Integer;

begin
  l:=Length(ss);
  SetLength(ss,l+1);
  ss[l]:=se;
end;

function EmptySequence:TOCSequence;

begin
  with Result do
    begin
    Sequence:='';
    Organism:='';
    Database:='';
    ID:='';
    Evidence:='';
    Code:='';
    Gene:='';
    Description:='';
    Score:=-1;
    Bits:=-1;
    Expectation:=-1;
    Identity:=-1;
    Positives:=-1;
    end;
end;

function ListOrganisms(const Sequences:TOCSequences):TSimpleStrings;

var
  f:Integer;
begin
  Result:=nil;
  for f:=0 to High(Sequences) do
    if not IsInArray(Sequences[f].Organism,Result) then
      AddToArray(Sequences[f].Organism,Result);
end;

function FilterByCommonOrganisms(SequencesSet:TOCSequencesSet):TOCSequencesSet;

var
  organisms:TSimpleStrings;
  inall:array of Boolean;
  indices:array of TIntegers;
  f,g,count:Integer;

begin
  organisms:=ListOrganisms(SequencesSet[0]);
  SetLength(inall,Length(organisms));
  SetLength(indices,Length(SequencesSet),Length(organisms));
  for f:=0 to High(inall) do
    inall[f]:=True;
  for f:=0 to High(organisms) do
    begin
    g:=0;
    while inall[f] and (g<=High(SequencesSet)) do
      begin
      indices[g,f]:=FirstIxByOrganism(SequencesSet[g],organisms[f]);
      inall[f]:=indices[g,f]>=0;
      Inc(g);
      end;
    end;
  count:=0;
  for f:=0 to High(inall) do if inall[f] then Inc(count);
  SetLength(Result,Length(SequencesSet),count);
  count:=0;
  for f:=0 to High(inall) do
    if inall[f] then
      begin
      for g:=0 to High(Result) do
        Result[g,count]:=SequencesSet[g,indices[g,f]];
      Inc(count);
      end;
end;

function FirstIxByOrganism(const Sequences:TOCSequences;Organism:string):Integer;

begin
  Result:=0;
  while (Result<=High(Sequences)) and (Organism<>Sequences[Result].Organism) do
    Inc(Result);
  if Result>High(Sequences) then Result:=-1;
end;

function CompressSequence(Sequence: string): TOCCompressedSequence;

var
  f:Integer;
  lastix:Integer;
  lastsymbol:string;

begin
  lastix:=-1;
  lastsymbol:='';
  SetLength(Result,Length(Sequence));
  for f:=1 to Length(Sequence) do
    if Sequence[f]<>lastsymbol then
      begin
      lastsymbol:=Sequence[f];
      if lastix>=0 then
        Result[Lastix].Last:=f-1;
      Inc(lastix);
      Result[lastix].First:=f;
      Result[lastix].Symbol:=Sequence[f];
      end;
  if lastix>=0 then
     Result[lastix].Last:=Length(Sequence);
  SetLength(Result,lastix+1);
end;

function CompressSequences(Sequences: TOCSequences): TOCCompressedSequences;

var f:Integer;

begin
  SetLength(Result,Length(Sequences));
  for f:=0 to High(Result) do
    Result[f]:=CompressSequence(Sequences[f].Sequence);
end;

function CompressSequences(Sequences: TSimpleStrings): TOCCompressedSequences;

var f:Integer;

begin
  SetLength(Result,Length(Sequences));
  for f:=0 to High(Result) do
    Result[f]:=CompressSequence(Sequences[f]);
end;

end.

