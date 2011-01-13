unit alignment;

{
  General alignment utilities and basetypes
}

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils,basetypes, ocstringutils, debugutils;

const

// Default marker for sequence gaps
  DefaultGapMarker='-';

type

{
  Sequences stored in these data structures are meant to be
  the alignment results, not the original sequences.
  The SequenceID fields shoule be used to identify the original sequences
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

function TrimMSA(msa:TMSA; seqid:string):TMSA;
  // removes all gaps from the identified sequence and trims all others
  // helps to find the residue variations from one sequence to all others
  // if the sequence is not found returns an empty MSA

function SinglesToMSA(query:string;sas:TSingleAlignments; gapmarker:char=DefaultGapMarker):TMSA;
  // Builds a MSA with the query sequence as the first sequence and
  // all matches in the array as aligned sequences.
  // Gaps in the query sequence at each alignment are removed from the matches.
  // The gapmarker must match the marker used in the source sequences.
  // Warning: this function doesn't check if all alignments correspond to the
  // same query sequence.

procedure AddAlignment(al:TSingleAlignment; var als:TSingleAlignments);

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
    if s[f]<>gm then AddInteger(f,positionlist);
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

  if six>0 then
    begin
    GetList(msa.Alignment[six],msa.GapMarker);
    SetLength(Result.Alignment,Length(msa.Alignment));
    SetLength(Result.SequenceIDs,Length(msa.SequenceIDs));
    for f:=0 to High(msa.Alignment) do
      begin
      Result.SequenceIDs[f]:=msa.SequenceIDs[f];
      Result.Alignment[f]:=FromList(msa.Alignment[f]);
      end;
    end;
end;

function SinglesToMSA(query: string; sas: TSingleAlignments; gapmarker:char=DefaultGapMarker): TMSA;

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


end.

