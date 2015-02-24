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
  SumOfPairScores is not calculating gap scores
  NeedlemanWunschAlign raises exception on invalid residues (like CA HOH, etc)
  May best be merged with sequence unit...
*******************************************************************************}

unit alignment;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils,basetypes, stringutils, progress, debugutils, fasta, sequence;

const

// Default marker for sequence gaps
  DefaultGapMarker='-';
// Default marker used in substitution matrixes. This is the index used in
  DefaultSubstitutionGapMarker='*';

type

  //Record for returning pairwise alignment results (mapping and score)
  TPairwiseAlignment=record
    Map:TIntegers;
      //indicates for each residue of molecule 1 which residue of 2 is a match.
      //-1 if not matched
    Score:TFloat;
  end;


{
  Sequences stored in these data structures are meant to be
  the alignment results, not the original sequences.
  The SequenceID fields should be used to identify the original sequences
}

  TMSA=record
    // this record holds multiple alignment matrixes
    // it may also hold a set of single alignments to one query sequence
    SequenceIDs:TSimpleStrings;
    Alignment:TSimpleStrings;
    GapMarker:char;
  end;

  TMSAs=array of TMSA;

  TSingleAlignment=record
    QueryId,MatchID:string;
    Score,Bits,Expectation,Probability:TFloat;
    Identity,Positives:Integer;
    QueryStart,QueryEnd,MatchStart,MatchEnd:Integer;
    AlignedQuery,AlignedMatch:string;
  end;

  TSingleAlignments=array of TSingleAlignment;

  TSubMatrix=record
    //Substitution matrices
    MonomerIndex:string;    //one letter code for monomers, same length as matrix
    Matrix:TMatrix;         //substitution values
    //NOTE: MonomerIndex is 1 based. Must subtract 1 to obtain index in Matrix
    Comments:TSimpleStrings;
    GapPenalty,GapExtension:TFloat;
  end;

function TrimMSA(msa:TMSA; seqid:string):TMSA;
  // removes all gaps from the identified sequence and trims all others
  // helps to find the residue variations from one sequence to all others
  // if the sequence is not found returns an empty MSA

function TransposeMSA(const MSA:TMSA):TSimpleStrings;
  // Returns strings with the MSA columns

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
function SumOfPairScores(const MSA:TMSA;const SubMat:TSubMatrix):TFLoats;
function GetSubScore(const SubMat:TSubMatrix;A1,A2:Char;SubGapMarker:Char=DefaultSubstitutionGapMarker):TFloat;
  //if A1 or A2 are not found, uses SubGapMarker.

function ScoreSequencePair(const Seq1,Seq2:string; const SubMat:TSubMatrix;
  SubGapMarker:Char=DefaultSubstitutionGapMarker):TFloat;overload;
function ScoreSequencePair(const Seq1,Seq2:TOCCompressedSequence;
  const SubMat:TSubMatrix;SubGapMarker:Char=DefaultSubstitutionGapMarker):TFloat;overload;




function GapPenalty(const SubMat:TSubMatrix; GapLen: Integer): TFloat;
function NeedlemanWunschAlign(Seq1Ix, Seq2Ix: TIntegers; const SubMat:TSubMatrix):TPairwiseAlignment;overload;
  //Sequences are the indexes in the SubMatrix scores.
  // Result is the mapping from Sequence 1 to 2

function NeedlemanWunschAlign(Seq1, Seq2: string; const SubMat:TSubMatrix):TPairwiseAlignment;overload;
  //Converts sequences to indeces, calls NeedlemanWunschAlign(Seq1Ix, Seq2Ix...

function IndexSequence(Seq:string;SubMat:TSubMatrix):TIntegers;
  //Returns 0 based array with indexes for SubMat. -1 for not found

function AlignmentToSequence(const Align,Seq:string; GapMarker:Char):TIntegers;
  //Indexes alignment to sequence, starting at index 1 (not index 0), for indexing strings

function SequenceFromIndex(Seq:string; Filter:TIntegers):string;

function ReadMSA(FileName:string):TMSA;
  //Fasta only. TODO:detect file types

function MapSequenceToMSA(Seq, MSASeq: string; GapMarker: char): TIntegers;
function FindMatch(Seq: string; MSA: TMSA; out Map:TIntegers): Integer;

implementation

function TrimMSA(msa: TMSA; seqid: string):TMSA;

var
  f,six:Integer;
  positionlist:TIntegers;

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
      end;
    end;
end;

function TransposeMSA(const MSA:TMSA):TSimpleStrings;

var
  f,g:Integer;

begin
  if MSA.Alignment=nil then
    Result:=nil
  else
    begin
    SetLength(Result,Length(MSA.Alignment[0]));
    for f:=1 to Length(MSA.Alignment[0]) do
      begin
      SetLength(Result[f-1],Length(MSA.Alignment));
      for g:=0 to High(MSA.Alignment) do
        Result[f-1,g+1]:=MSA.Alignment[g,f];
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
  //TO DO: gap penalties should come from the file...
  Result.GapPenalty:=-10;
  Result.GapExtension:=-1;

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

function SumOfPairScores(const MSA:TMSA;const SubMat:TSubMatrix):TFLoats;

// TO DO: Calculate gap scores


var
  col,l1,l2:Integer;
begin
  Result:=nil;
  if MSA.Alignment<>nil then
      begin

      SetLength(Result,Length(MSA.Alignment[0]));
      for col:=0 to High(Result) do
        begin
        Result[col]:=0;
        for l1:=1 to High(MSA.Alignment)-1 do
          for l2:=l1+1 to High(MSA.Alignment) do
              Result[col]:=Result[col]+GetSubScore(SubMat,
                MSA.Alignment[l1,col+1],MSA.Alignment[l2,col+1]);
        end;
      end;

end;

function GetSubScore(const SubMat:TSubMatrix;A1,A2:Char;SubGapMarker:Char=DefaultSubstitutionGapMarker):TFloat;

var i1,i2:Integer;

begin
  i1:=Pos(A1,SubMat.MonomerIndex)-1;
  if i1<0 then i1:=Pos(SubGapMarker,SubMat.MonomerIndex)-1;
  i2:=Pos(A2,SubMat.MonomerIndex)-1;
  if i2<0 then i1:=Pos(SubGapMarker,SubMat.MonomerIndex)-1;
  if (i1>=0) and (i2>=0) then
    Result:=SubMat.Matrix[i1,i2]
  else Result:=0;
end;

function ScoreSequencePair(const Seq1, Seq2: string; const SubMat: TSubMatrix;
  SubGapMarker:Char=DefaultSubstitutionGapMarker): TFloat;

var f:Integer;

begin
  Result:=0;
  Assert(Length(Seq1)=Length(Seq2),'Error in scoring sequence pair: sequences must have equal lengths');
  for f:=1 to Length(Seq1) do
    Result:=Result+GetSubScore(SubMat,Seq1[f],Seq2[f],SubGapMarker);
end;

function ScoreSequencePair(const Seq1, Seq2: TOCCompressedSequence;
  const SubMat: TSubMatrix;SubGapMarker:Char=DefaultSubstitutionGapMarker): TFloat;

var
  ix1,ix2:Integer;
  st,en:Integer;
  tmp:TFloat;

begin
  Assert(Seq1[High(Seq1)].Last=Seq2[High(Seq2)].Last,'Error in scoring sequence pair: sequences must have equal lengths');
  Result:=0;
  ix1:=0;
  ix2:=0;
  while (ix1<Length(Seq1)) do
    begin
    st:=Max(Seq1[ix1].First,Seq2[ix2].First);
    en:=Min(Seq1[ix1].Last,Seq2[ix2].Last)+1;
    tmp:=GetSubScore(SubMat,Seq1[ix1].Symbol,Seq2[ix2].Symbol,SubGapMarker)*(en-st);
    Result:=Result+tmp;
    if en>Seq1[ix1].Last then Inc(ix1);
    if en>Seq2[ix2].Last then Inc(ix2);
    end;
end;

function GapPenalty(const SubMat:TSubMatrix; GapLen: Integer): TFloat;
begin
  if GapLen<=0 then Result:=0
  else Result:=SubMat.GapPenalty+GapLen*SubMat.GapExtension;
end;

function NeedlemanWunschAlign(Seq1, Seq2: string; const SubMat: TSubMatrix
  ): TPairwiseAlignment;

var
   seq1ix,seq2ix:TIntegers;

begin
  seq1ix:=IndexSequence(Seq1,SubMat);
  seq2ix:=IndexSequence(Seq2,SubMat);
  Result:=NeedlemanWunschAlign(seq1ix, seq2ix, SubMat);
end;

function IndexSequence(Seq:string;SubMat:TSubMatrix):TIntegers;

var f:Integer;

begin
  SetLength(Result,Length(Seq));
  for f:=0 to High(Result) do
    Result[f]:=Pos(Seq[f+1],SubMat.MonomerIndex)-1;
end;

function NeedlemanWunschAlign(Seq1Ix, Seq2Ix: TIntegers; const SubMat:TSubMatrix):TPairwiseAlignment;


var
  scoremat:TMatrix;
  total:Integer;
  task:TRunningTask;   //for progress report
  taskstep:Single;

function TotalVal(S1,E1,S2,E2:Integer):TFloat;

var
  f:Integer;
  foundgap:Boolean;

begin
  task.Step(taskstep);
  Result:=scoremat[E1,E2];
  if (S1<E1) and (S1>=0) then
    Result:=Result+GapPenalty(SubMat,E1-S1-1);
  if (S2<E2) and (S2>=0) then
    Result:=Result+GapPenalty(SubMat,E2-S2-1);
end;

procedure EvalMat(I1,I2:Integer);
// indexes are for the matrix, starting at 0

var
  x,y:Integer;
  t,tot,max:TFLoat;

begin
  if (Seq1Ix[I1]>=0) and (Seq2Ix[I2]>=0) then
    tot:=SubMat.Matrix[Seq1Ix[I1],Seq2Ix[I2]]
  //if one residue is not present in the substitution matrix
  //then the score starts at 0.
  else tot:=0;
  if (I1<High(Seq1Ix)) and (I2<High(Seq2Ix)) then
    begin
    max:=-10000;
    for x:=I1+1 to High(Seq1Ix) do
      begin
      t:=TotalVal(I1,x,I2+1,I2+1);
      if t>max then max:=t;
      end;
    for y:=I2+1 to High(Seq2Ix) do
      begin
      t:=TotalVal(I1+1,I1+1,I2,y);
      if t>max then max:=t;
      end;
    tot:=tot+max;
    end;
  ScoreMat[I1,I2]:=tot;
end;

procedure CalcLines(I1,I2:Integer);

var f:Integer;

begin
  for f:=0 to I1 do EvalMat(f,I2);
  for f:=0 to I2-1 do EvalMat(I1,f);
end;

procedure BuildMat;

var i1,i2:Integer;

begin
  SetLength(ScoreMat,Length(Seq1Ix),Length(Seq2Ix));
  i1:=High(ScoreMat);
  i2:=High(ScoreMat[0]);
  repeat
    CalcLines(i1,i2);
    Dec(i1);
    Dec(i2);
  until(i1<0) or (i2<0);
end;

procedure Align;

var
  I1,I2,NI1,NI2,f:Integer;

procedure NextLine(var I1,I2:Integer);

var
  x,y,s1,s2:Integer;
  t,max:TFLoat;

begin
  max:=0;
  s1:=I1+1;
  s2:=I2+1;
  I1:=Length(ScoreMat);
  I2:=Length(ScoreMat[0]);
  try
  for x:=s1 to High(ScoreMat) do
    begin
    t:=TotalVal(s1-1,x,s2,s2);
    if max<t then
      begin
      max:=t;
      I1:=x;
      I2:=s2;
      end;
    end;
  for y:=s2 to High(ScoreMat[0]) do
    begin
    t:=TotalVal(s1,s1,s2-1,y);
    if max<t then
      begin
      max:=t;
      I1:=s1;
      I2:=y;
      end;
    end;
  except
    writeln(s1,' ',s2,' ',length(scoremat),' ',length(scoremat[0]));
  end;
  if max>Result.Score then
    Result.Score:=max;
end;


begin
  NI1:=-1;
  NI2:=-1;
  repeat
    NextLine(NI1,NI2);
    if (NI1<=High(Seq1Ix)) and (NI2<=High(Seq2Ix)) then
      Result.Map[NI1]:=NI2;
  until (NI1>=High(Seq1Ix)) or (NI2>=High(Seq2Ix));
end;

var f:Integer;

begin
  SetLength(Result.Map,Length(Seq1Ix));
  for f:=0 to High(Result.Map) do Result.Map[f]:=-1;
  if (Seq1Ix<>nil) and (Seq2Ix<>nil) then
    begin
    try
      taskstep:=2/Length(Seq1Ix)/Length(Seq2Ix)/(Length(Seq1Ix)+Length(Seq2Ix));
      task:=NewTask(False,'Aligning');
      BuildMat;
      Result.Score:=ScoreMat[0,0];
      Align;
    finally
      FreeTask(task);
    end;
    ScoreMat:=nil;
    end;
end;

function AlignmentToSequence(const Align,Seq:string; GapMarker:Char):TIntegers;

var c,f:Integer;
    tmp:string; //debug
begin
  c:=0;
  tmp:='';
  SetLength(Result,Length(Align));
  for f:=1 to Length(Align) do
    if Align[f]<>GapMarker then
      begin
      Inc(c);
      Result[f-1]:=c;
      tmp:=tmp+Align[f]
      end
    else
      Result[f-1]:=-1;
  DebugReport('AlignToSequence');
  DebugReport(Seq);
  DebugReport(tmp);
end;

function SequenceFromIndex(Seq:string; Filter:TIntegers):string;

var f:Integer;

begin
  Result:='';
  for f:=0 to High(Filter) do
    if Filter[f]>0 then
      Result:=Result+Seq[Filter[f]];
end;

function ReadMSA(FileName: string): TMSA;

var
  reader:TFastaReader;
  f:Integer;

begin
  Result.GapMarker:=DefaultGapMarker;
  reader:=TFastaReader.Create(FileName);
  SetLength(Result.Alignment,Length(Reader.Seqs));
  SetLength(Result.SequenceIds,Length(Reader.Seqs));
  for f:=0 to High(Reader.Seqs) do
    begin
    Result.Alignment[f]:=Reader.Seqs[f].Sequence;
    Result.SequenceIds[f]:=Reader.Seqs[f].ID;
    end;
end;

function MapSequenceToMSA(Seq, MSASeq: string; GapMarker: char): TIntegers;

var
  s,m,ls,lm:Integer;


begin
  lm:=Length(MSASeq);
  ls:=Length(Seq);
  SetLength(Result,ls);
  s:=1;
  m:=1;
  while (s<=ls) and (m<=lm) do
    begin
    if Seq[s]=MSASeq[m] then
      begin
      Result[s-1]:=m;
      Inc(s);
      Inc(m);
      end
    else if MSASeq[m]=GapMarker then Inc(m)
    else
      begin
      Result:=nil;
      Break;
      end;
    end;
end;

function FindMatch(Seq: string; MSA: TMSA; out Map:TIntegers): Integer;

var f:Integer;

begin
  Result:=-1;
  for f:=0 to High(MSA.Alignment) do
    begin
    Map:=MapSequenceToMsa(Seq,MSA.Alignment[f],MSA.GapMarker);
    if Map<>nil then
      begin
      Result:=f;
      Break
      end;
    end;
end;

end.


