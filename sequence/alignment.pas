{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 9.1.2011
Purpose:
  Functions for alignment of sequences
Requirements:
Revisions:
To do:
  SinglesToMSA is ignoring the TrimToQuery parameter. Rewrite
*******************************************************************************}

unit alignment;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils,basetypes, ocstringutils, debugutils, LCLProc;

const

// Default marker for sequence gaps
  DefaultGapMarker='-';

type

{
  Sequences stored in these data structures are meant to be
  the alignment results, not the original sequences.
  The SequenceID fields should be used to identify the original sequences
}

  TMSA=record
    // this record holds multiple alignment matrixes
    // it may also hold a set of single alignments to one query sequence
    SequenceIDs:TOCStrings;
    Alignment:TOCStrings;
    GapMarker:char;
  end;

  TMSAs=array of TMSA;

  TSingleAlignment=record
    QueryId,MatchID:string;
    Score,Bits,Expectation,Probability:TOCFloat;
    Identity,Positives:Integer;
    QueryStart,QueryEnd,MatchStart,MatchEnd:Integer;
    AlignedQuery,AlignedMatch:string;
  end;

  TSingleAlignments=array of TSingleAlignment;

  TSubMatrix=record
    //Substitution matrices
    MonomerIndex:string;      //one letter code for monomers, same length as matrix
    Matrix:TOCMatrix;         //substitution values
    Comments:TOCStrings;
  end;

function TrimMSA(msa:TMSA; seqid:string):TMSA;
  // removes all gaps from the identified sequence and trims all others
  // helps to find the residue variations from one sequence to all others
  // if the sequence is not found returns an empty MSA

function SinglesToMSA(query:string;sas:TSingleAlignments;
  TrimToQuery:Boolean=False; gapmarker:char=DefaultGapMarker):TMSA;
  // Builds a MSA with the query sequence as the first sequence and
  // all matches in the array as aligned sequences.
  // Gaps in the query sequence at each alignment are removed from the matches,
  // if TrimToQuery is true.
  // The gapmarker must match the marker used in the source sequences.
  // Warning: this function doesn't check if all alignments correspond to the
  // same query sequence.

procedure AddAlignment(al:TSingleAlignment; var als:TSingleAlignments);

function ReadBLASTMatrix(FileName:string):TSubMatrix;
  //Reads a substitution matrix in the BLAST format
  //see ftp://ftp.ncbi.nih.gov/blast/matrices for examples

procedure SaveMSAToFile(MSA:TMSA; FileName:string);

implementation

function TrimMSA(msa: TMSA; seqid: string):TMSA;

var
  f,six:Integer;
  positionlist:TOCIntegers;

procedure GetList(s:string; gm:char);

var f:Integer;

begin
  positionlist:=nil;
  for f:=1 to Length(s) do
    if s[f]<>gm then AddToArray(f,positionlist);
end;

function FromList(s:string):string;

var f:Integer;

begin
  Result:='';
  for f:=0 to High(positionlist) do
    Result:=Result+s[positionlist[f]];
end;


begin
  Result.SequenceIDs:=nil;
  Result.Alignment:=nil;
  Result.GapMarker:=msa.GapMarker;
  six:=LastIndexOf(seqid,msa.SequenceIDs);
  if six>=0 then
    begin
    GetList(msa.Alignment[six],msa.GapMarker);
    SetLength(Result.Alignment,Length(msa.Alignment));
    SetLength(Result.SequenceIDs,Length(msa.SequenceIDs));
    for f:=0 to High(msa.Alignment) do
      begin
      Result.SequenceIDs[f]:=msa.SequenceIDs[f];
      Result.Alignment[f]:=FromList(msa.Alignment[f]);
      if f<10 then DebugLn(Result.Alignment[f]);
      end;
    end;
end;

function SinglesToMSA(query: string; sas: TSingleAlignments;
  TrimToQuery:Boolean=False; gapmarker:char=DefaultGapMarker): TMSA;

//TO DO: implement TrimToQuery. Currently being ignored, and trimming everything

var
  f,l:Integer;

function BuildMatch(ix:Integer):string;

var
  f,qst:Integer;
  sq,sm:string;

begin
  Result:='';
  qst:=sas[ix].QueryStart;
  sq:=sas[ix].AlignedQuery;
  sm:=sas[ix].AlignedMatch;
  while Length(Result)<qst-1 do
    Result:=Result+gapmarker;
  for f:=1 to Length(sq) do
    if sq[f]<>gapmarker then
      Result:=Result+sm[f];
  while Length(Result)<Length(query) do
    Result:=Result+GapMarker;
end;

begin
  Result.SequenceIDs:=nil;
  Result.Alignment:=nil;
  l:=Length(sas);
  if l>1 then
    begin
    SetLength(Result.SequenceIDs,l);
    SetLength(Result.Alignment,l);       
    for f:=0 to High(sas) do
      begin
      Result.Alignment[f]:=BuildMatch(f);
      Result.SequenceIDs[f]:=sas[f].MatchID;
      end;
    end;
end;

procedure AddAlignment(al: TSingleAlignment; var als: TSingleAlignments);

var l:Integer;

begin
  l:=Length(als);
  SetLength(als,l+1);
  als[l]:=al;
end;

function ReadBLASTMatrix(FileName:string):TSubMatrix;

var
  buf:TStringList;
  f,ix:Integer;
  s:string;

begin
  Result.Comments:=nil;
  Result.Matrix:=nil;
  Result.MonomerIndex:='';
  buf:=TStringList.Create;
  buf.LoadFromFile(FileName);
  for f:=0 to buf.Count-1 do
    begin
    s:=buf.Strings[f];

    // snip comments
    ix:=Pos('#',s);
    if ix>0 then
      begin
      AddToArray(Copy(s,ix+1,Length(s)),Result.Comments);
      Delete(s,ix,Length(s));
      end;
    if s<>'' then
      begin
      if Result.MonomerIndex='' then    //first data line must be headings
        begin
        for ix:=1 to Length(s) do
          if s[ix]>' ' then
            Result.MonomerIndex:=Result.MonomerIndex+s[ix];
            SetLength(Result.Matrix,Length(Result.MonomerIndex));
        end
      else                               //other lines are value lines
        begin
        ix:=Pos(s[1],Result.MonomerIndex);
        if ix>0 then
          begin
          Dec(ix);
          Delete(s,1,1);                 //remove monomer code
          Result.Matrix[ix]:=StringToFloats(s);
          end;
        end;

      end;

    end;
  buf.Free;
end;

procedure SaveMSAToFile(MSA:TMSA; FileName:string);

var
  f:Integer;
  buf:TStringList;

begin
  buf:=TStringList.Create;
  buf.Add('**IDENTIFIERS**');
  for f:=0 to High(MSA.SequenceIDs) do
    buf.Add(MSA.SequenceIds[f]);
  buf.Add('**ALIGNMENT**');
  for f:=0 to High(MSA.Alignment) do
    buf.Add(MSA.Alignment[f]);
  buf.SaveToFile(FileName);
  buf.Free;
end;

end.

