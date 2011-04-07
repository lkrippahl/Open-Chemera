{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 9.1.2011
Purpose:
  sequence evolution correlations and structure comparisons
Requirements:
Revisions:
To do:
 1.Single function to extract sequences by matching organisms instead
 of having each cross* function do its own
 2. Then use OrganismUnknown to check if organism was found
 3. Then change ImportEBIBlastOrgIDs to remove the creation of unique
 organism IDS. Also ImportEBIFastaOrgIDs
 4. there is some code duplication on storing alignments. Use TMSA?
*******************************************************************************}
unit seqstruccomp;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, fasta, pdbparser, clustalparser, alignment,
  ocstringutils, ebixmlparser, ebifastaparser, progress, chartutils, graphics, threedcalc, sequence;

const
  //default organism name
  OrganismUnknown='Organism unknown';

  //Result types
  //function RTAsLabel returns a text label for each mode
  RSLTSelfDistance=0;
  RSLTSelfMI=1;
  RSLTCrossDistance=2;
  RSLTCrossMI=3;
  RSLTSelfInteract=4;  //Not implemented
  RSLTCrossInteract=5; // implemented with DifferentInteraction

type
  //Information on query sequence and structure
  TQueryData=record
    ID:string;
    Sequence:string;
    ResCoords:TOCCoords;
    MultiAlign:TOCStrings;
    //multialign is trimmed or expanded so the matrix matches the
    //length of the query sequence.
    SeqIDs:TOCStrings;
    OrganismIDs:TOCStrings;
    end;

  TComparisonResult=record
    ID1,ID2:string;
    ResultType:Integer;
    Data:TOCMatrix;
    end;
  TResults = array of TComparisonResult;

  TQueries = array of TQueryData;

  { TComparer }

  TComparer=class
  protected
    FQueries:TQueries;
    FResults:TResults;
    FScaleMI:Boolean; // if true scales MutualInformation by the size of symbol strings

    //Data for interaction propensities
    FInteractionPropensities:TOCMatrix;
    FInteractionPropensitiesAAs:string;


    function RTAsLabel(m:Integer):string;
    function Transpose(sl:TOCStrings):TOCStrings;
      //transposes an array of strings.
      //all strings are assumed to be the same length
    procedure AddQuery(nq:TQueryData);
    function QueryIxById(id:string):Integer;
    function MutualInformation(s1,s2:string):TOCFloat; //s1 and s2 must be of equal length and aligned

    function Interaction(s1,s2:string):TOCFloat;       //s1 and s2 must be of equal length and aligned
      //average of interaction scores for all pairings

    function DifferentInteraction(s1,s2:string):TOCFloat;       //s1 and s2 must be of equal length and aligned
      //average interaction scores for the different pairings
      //with 0 if only one pairing
      //****DOESN'T WORK WELL. TO REMOVE*****

    procedure AddResult(id1,id2:string;rt:Integer;d:TOCMatrix);
    function ListSelfResults(id:string):TOCIntegers;
    function ListCrossResults(id1,id2:string):TOCIntegers;
    function FindFirstMatrix(id1,id2:string;dt:Integer):TOCMatrix;
      //finds the first data matrix with these ids and type

    function NeighbourCount(QIx:Integer;Cutoff:Single):TOCIntegers;

  public
    property ScaleMI:Boolean read FScaleMI write FScaleMI;
    constructor Create(DataFolder:string=''); //data folder without terminating slash
    procedure LoadInteractionPropensities(FileName:string);
    procedure ClearData;
    procedure LoadQuery(id,seqfile,pdbfile:string);
    procedure LoadClustalMSA(queryid,alnfile:string;TrimToQuery:Boolean = False);
      // queryid is both the query sequence id and the id of the corresponding
      // sequence in the .aln file
    procedure LoadEBIBlast(queryid,filename:string; TrimToQuery:Boolean = False);
    procedure ImportEBIBlastOrgIDs(queryid,filename:string);
      // imports only the organism ids from the EBI blast xml file
      // into a previously loaded msa
    procedure ImportEBIFastaOrgIDs(QueryId,Filename:string);
      // imports only the organism ids from the EBI fasta file
      // into a previously loaded msa

    procedure SelfMutualInf(seqid:string);
    procedure SelfDistance(seqid:string);
    procedure SelfPlot(seqid:string; sl:TStringList);
    procedure CrossPlot(seq1,seq2:string; sl:TStringList);
    procedure CrossDistance(sid1,sid2:string);
    procedure CrossMutualInf(sid1,sid2:string;Cutoff:Single; MaxNeibs:Integer);
    procedure CrossInteract(sid1,sid2:string;Cutoff:Single;MaxNeibs:Integer);
    procedure ClearResults;
    procedure FilteredHistogram(sid1, sid2: string; filter,
              data,bins: Integer; filtermin, filtermax: Real; sl: TStringList;
              BinColumn:Boolean);
    procedure LowestAndCount(sid1,sid2:string;low,count:Integer;
      countmin,countmax:Real;sl:TStringList);
    procedure PlotData(sid1,sid2:string;datax,datay:Integer;bmp:TBitmap);
    function GetMSA(Id:string):TMSA;
    function MergeAsMSA(Id1,ID2:string):TMSA;
      //merges by organism id, concatenating the sequences

  end;

implementation

{ TComparer }

procedure TComparer.LoadInteractionPropensities(FileName: string);

var
  sl:TStringList;
  f,g:Integer;
  s:TOCStrings;

begin
  sl:=TStringList.Create;
  sl.LoadFromFile(FileName);
  FInteractionPropensitiesAAs:=sl.Strings[0];
  SetLength(FInteractionPropensities,20,20);
  for f:=1 to 20 do
    begin
    s:=SplitString(sl.Strings[f],#9);
    for g:=0 to 19 do
      FInteractionPropensities[f-1,g]:=StrToFloat(s[g]);
    end;
  sl.Free;
end;

function TComparer.RTAsLabel(m: Integer): string;
begin
  case m of
    RSLTSelfDistance:Result:='Distance';
    RSLTSelfMI:Result:='Mut.Inf.';
    RSLTCrossDistance:Result:='X Distance';
    RSLTCrossMI:Result:='X MI';
  else Result:='Unknown';
  end;
end;

function TComparer.Transpose(sl: TOCStrings): TOCStrings;

var
  f,g:Integer;
  s:string;
begin
  Result:=nil;
  if sl<>nil then
    begin
    SetLength(Result,Length(sl[0]));
    for f:=1 to Length(sl[0]) do
      begin
      s:='';
      for g:=0 to High(sl) do
        s:=s+sl[g,f];
      Result[f-1]:=s;
      end;
    end;
end;

procedure TComparer.AddQuery(nq: TQueryData);
begin
  SetLength(FQueries,Length(FQueries)+1);
  FQueries[High(FQueries)]:=nq;
end;

function TComparer.QueryIxById(id: string): Integer;
begin
  Result:=High(FQueries);
  while (Result>=0) and (FQueries[Result].ID<>id) do
    Dec(Result);
end;

function TComparer.MutualInformation(s1, s2: string): TOCFloat;

const
  Ln2=0.69314718056; //ln(2) to convert to log2

var
  f,g:Integer;
  symb1,symb2:string;
  step,min,mult:TOCFloat;
  marg1,marg2:array of TOCFloat;
  cross:array of array of TOCFLoat;

procedure BuildSymbols;

var f,g:Integer;

begin
  symb1:='';
  symb2:='';
  for f:=1 to Length(s1) do
    if Pos(s1[f],symb1)<1 then symb1:=symb1+s1[f];
  for f:=1 to Length(s2) do
    if Pos(s2[f],symb2)<1 then symb2:=symb2+s2[f];
  SetLength(marg1,Length(symb1));
  SetLength(marg2,Length(symb2));
  SetLength(cross,Length(symb1),Length(symb2));
  for f:=0 to Length(symb1)-1 do
    marg1[f]:=0;
  for f:=0 to Length(symb2)-1 do
    marg2[f]:=0;
  for f:=0 to High(cross) do
    for g:=0 to High(cross[0]) do
      cross[f,g]:=0;
end;

procedure Count;

var f,i1,i2:Integer;

begin
  step:=1/Length(s1);
  for f:=1 to Length(s1) do
    begin
    i1:=Pos(s1[f],symb1)-1;
    i2:=Pos(s2[f],symb2)-1;
    marg1[i1]:=marg1[i1]+step;
    marg2[i2]:=marg2[i2]+step;
    cross[i1,i2]:=cross[i1,i2]+step;
    end;
end;

begin
  Result:=0;
  BuildSymbols;
  Count;
  min:=step/2;
  for f:=0 to High(cross) do
    for g:=0 to High(cross[0]) do
      begin
      mult:=marg1[f]*marg2[g];
      if cross[f,g]>min then
        Result:=Result+cross[f,g]*Ln(cross[f,g]/mult);
      end;
  Result:=Result/Ln2;
  if FScaleMI then Result:=Result*Ln(Sqrt(Length(symb1)*Length(symb2)));
end;

function TComparer.Interaction(s1, s2: string): TOCFloat;

var
  f,g,p1,p2:Integer;

begin
  Result:=0;
  for f:=1 to Length(s1) do
    begin
    p1:=Pos(s1[f],FInteractionPropensitiesAAs)-1;
    p2:=Pos(s2[f],FInteractionPropensitiesAAs)-1;
    if (p1>=0) and (p2>=0) then
      Result:=Result+FInteractionPropensities[p1,p2];
    end;
  Result:=Result/Length(s1);
end;

function TComparer.DifferentInteraction(s1, s2: string): TOCFloat;
var
  f,g,c,p1,p2:Integer;
  diffs:array of array of Integer;

begin
  Result:=0;
  SetLength(diffs,Length(FInteractionPropensitiesAAs),Length(FInteractionPropensitiesAAs));
  for f:=0 to High(diffs) do
    for g:=0 to High(diffs) do
      diffs[f,g]:=0;
  c:=0;
  for f:=1 to Length(s1) do
    begin
    p1:=Pos(s1[f],FInteractionPropensitiesAAs)-1;
    p2:=Pos(s2[f],FInteractionPropensitiesAAs)-1;
    if (p1>=0) and (p2>=0) and (diffs[p1,p2]=0) then
      begin
      Inc(c);
      diffs[p1,p2]:=1;
      Result:=Result+FInteractionPropensities[p1,p2];
      end;
    end;
  if c>1 then
    Result:=Result/sqrt(c)/4
  else Result:=0;

end;

procedure TComparer.AddResult(id1,id2:string; rt: Integer; d: TOCMatrix);

var ix:Integer;

begin
  ix:=Length(FResults);
  SetLength(FResults,ix+1);
  FResults[ix].ID1:=id1;
  FResults[ix].ID2:=id2;
  FResults[ix].Data:=d;
  FResults[ix].ResultType:=rt;
end;

function TComparer.ListSelfResults(id: string): TOCIntegers;

begin
  Result:=ListCrossResults(id,'');
end;

function TComparer.ListCrossResults(id1, id2: string): TOCIntegers;

var f:Integer;

begin
  Result:=nil;
  for f:=0 to High(FResults) do
    if (FResults[f].ID1=id1) and (FResults[f].ID2=id2) then
      AddToArray(f,Result);
end;

function TComparer.FindFirstMatrix(id1, id2: string; dt: Integer): TOCMatrix;

var f:Integer;

begin
  Result:=nil;
  for f:=0 to High(FResults) do
    if (FResults[f].ID1=id1) and
      (FResults[f].ID2=id2) and
      (FResults[f].ResultType=dt) then
      begin
      Result:=FResults[f].Data;
      Break;
    end;
end;

function TComparer.NeighbourCount(QIx: Integer; Cutoff: Single): TOCIntegers;

var
  f,g:Integer;
  d:TOCFloat;
begin
  Result:=FilledInts(Length(FQueries[QIx].ResCoords),0);
  for f:=0 to High(Result) do
    for g:=f+1 to High(Result) do
      begin
      d:=Distance(FQueries[QIx].ResCoords[f],FQueries[QIx].ResCoords[g]);
      if d<=Cutoff then
        begin
        Inc(Result[f]);
        Inc(Result[g]);
        end;
      end;
end;

constructor TComparer.Create(DataFolder:string='');
begin
  inherited Create;
  FScaleMI:=False;
  if DataFolder<>'' then
    begin
    LoadInteractionPropensities(DataFolder+PathDelim+'AAPropensities.txt');
    end;
end;

procedure TComparer.ClearData;
begin
  FQueries:=nil;
end;

procedure TComparer.LoadQuery(id, seqfile, pdbfile: string);

var
  nq:TQueryData;
  fr:TFastaReader;
  pr:TPDBReader;
  f,six,rix,ac:Integer;
  c:TOCCoord;

begin
  nq.ID:=id;
  c:=NullVector;
  fr:=TFastaReader.Create(seqfile);
  if fr.Seqs<>nil then
    nq.Sequence:=fr.Seqs[0].Sequence;
  fr.Free;
  pr:=TPDBReader.Create(pdbfile);
  SetLength(nq.ResCoords,Length(nq.Sequence));
  if pr.Atoms<>nil then
    begin
    six:=-1;
    ac:=0;
    rix:=-1;
    for f:=0 to High(pr.Atoms) do
      begin
      if rix<>pr.Atoms[f].ResSeq then
        begin
        if ac>0 then
          nq.ResCoords[six]:=ScaleVector(c,1/ac);
        ac:=0;
        Inc(six);
        rix:=pr.Atoms[f].ResSeq;
        c:=NullVector;
        end;
      inc(ac);
      c:=AddVectors(c,pr.Atoms[f].Coords);
      end;
    if ac>0 then
      nq.ResCoords[six]:=ScaleVector(c,1/ac);
    end;
  pr.Free;
  AddQuery(nq);
end;

procedure TComparer.LoadClustalMSA(queryid, alnfile: string; TrimToQuery:Boolean = False);

var
  msa:TMSA;
  ix:Integer;
begin
  ix:=QueryIxById(queryid);
  if ix>=0 then
    begin
    msa:=ReadClustal(alnfile);
    if TrimToQuery then msa:=TrimMSA(msa,queryid);
    FQueries[ix].MultiAlign:=msa.Alignment;
    FQueries[ix].SeqIDs:=msa.SequenceIDs;
    end;
end;

procedure TComparer.LoadEBIBlast(queryid, filename: string; TrimToQuery:Boolean=False);

var
  blast:TBlastPData;
  msa:TMSA;
  qseq:string;
  f,six:Integer;

begin
  six:=QueryIxById(queryid);
  if six>=0 then
    begin
    qseq:=FQueries[six].Sequence;
    blast:=ReadBlastPFile(filename);
    msa:=SinglesToMSA(qseq,blast.Alignments);
    FQueries[six].MultiAlign:=msa.Alignment;
    SetLength(FQueries[six].SeqIDs,Length(blast.Alignments));
    SetLength(FQueries[six].OrganismIDs,Length(blast.Alignments));
    for f:=0 to High(blast.Alignments) do
      begin
      FQueries[six].SeqIDs[f]:=blast.SequenceInfo[f].ID;
      FQueries[six].OrganismIDs[f]:=blast.SequenceInfo[f].Organism;
      end;
    end;
end;

procedure TComparer.ImportEBIBlastOrgIDs(queryid, filename: string);
var
  blast:TBlastPData;
  f,six,oix:Integer;

begin
  six:=QueryIxById(queryid);
  if six>=0 then
    begin
    blast:=ReadBlastPFile(filename);
    SetLength(FQueries[six].OrganismIDs,Length(FQueries[six].SeqIDs));

    //remove this part after fixing the unknown organism check.
    //this is a quick fix for forcing mismatches between unknown organisms
    for f:=0 to High(FQueries[six].OrganismIDs) do
      FQueries[six].OrganismIDs[f]:=OrganismUnknown+
        queryid+IntToStr(f);

    for f:=0 to High(blast.Alignments) do
      begin
      oix:=LastIndexOf(blast.SequenceInfo[f].ID,FQueries[six].SeqIDs);
      if oix>=0 then
        FQueries[six].OrganismIDs[oix]:=blast.SequenceInfo[f].Organism;
      end;
    end;
end;

procedure TComparer.ImportEBIFastaOrgIDs(QueryId, Filename: string);

var
  seqs:TOCSequences;
  f,six,oix:Integer;

begin
  six:=QueryIxById(QueryId);
  if six>=0 then
    begin
    seqs:=ReadEBIFasta(FileName);
    SetLength(FQueries[six].OrganismIDs,Length(FQueries[six].SeqIDs));

      //remove this part after fixing the unknown organism check.
      //this is a quick fix for forcing mismatches between unknown organisms
      for f:=0 to High(FQueries[six].OrganismIDs) do
        FQueries[six].OrganismIDs[f]:=OrganismUnknown+
          queryid+IntToStr(f);

    for f:=0 to High(seqs) do
      begin
      oix:=LastIndexOf(seqs[f].ID,FQueries[six].SeqIDs);
      if oix>=0 then
        FQueries[six].OrganismIDs[oix]:=seqs[f].Organism;
      end;
  end;
end;

procedure TComparer.SelfMutualInf(seqid: string);

var
  six:Integer;
  mis:TOCMatrix;

procedure CalcMis;

var
  f,g:Integer;
  seqs:TOCStrings;
begin
  seqs:=Transpose(FQueries[six].MultiAlign);
  for f:=0 to High(seqs) do
    for g:=f to High(seqs) do
      if g=f then mis[f,g]:=0
      else
        begin
        mis[f,g]:=MutualInformation(seqs[f],seqs[g]);
        mis[g,f]:=mis[f,g];
        end;
end;

begin
  six:=QueryIxById(seqid);
  if six>=0 then
    begin
    SetLength(mis,Length(FQueries[six].Sequence),Length(FQueries[six].Sequence));
    CalcMIs;
    AddResult(seqid,'',RSLTSelfMI,mis);
    end;
end;

procedure TComparer.SelfDistance(seqid: string);
var
  six:Integer;
  dists:TOCMatrix;

procedure CalcDists;

var
  f,g:Integer;

begin
  for f:=0 to High(dists) do
    begin
    dists[f,f]:=0;
    if f<High(dists) then
      for g:=f+1 to High(dists) do
        begin
        dists[f,g]:=Distance(FQueries[six].ResCoords[f],FQueries[six].ResCoords[g]);
        dists[g,f]:=dists[f,g];
        end;
    end;
end;

begin
  six:=QueryIxById(seqid);
  if six>=0 then
    begin
    SetLength(dists,Length(FQueries[six].Sequence),Length(FQueries[six].Sequence));
    CalcDists;
    AddResult(seqid,'',RSLTSelfDistance,dists);
    end;
end;

procedure TComparer.SelfPlot(seqid:string; sl: TStringList);

var
  six:Integer;
  f,g,h,hi:Integer;
  ixs:TOCIntegers;
  s:string;

begin
  ixs:=ListSelfResults(seqid);
  six:=QueryIxById(seqid);
  if ixs<>nil then
    begin
    s:='Res1'+#9+'Res2';
    for f:=0 to High(Ixs) do
      s:=s+#9+RTAsLabel(FResults[ixs[f]].ResultType);
    hi:=High(FResults[ixs[0]].Data);
    sl.Add(s);
    for g:=0 to hi-1 do
      for h:=g+1 to hi do
        begin
        if six>=0 then
          s:=FQueries[six].Sequence[g+1]+#9+
            FQueries[six].Sequence[h+1]
        else s:='?'+#9+'?';
        for f:=0 to High(ixs) do
          s:=s+#9+FloatToStr(FResults[ixs[f]].Data[g,h]);
        sl.Add(s);
        end;
    end;
end;

procedure TComparer.CrossPlot(seq1, seq2: string; sl: TStringList);
var
  six1,six2:Integer;
  f,g,h,h1,h2:Integer;
  ixs:TOCIntegers;
  s:string;

begin
  ixs:=ListCrossResults(seq1,seq2);
  six1:=QueryIxById(seq1);
  six2:=QueryIxById(seq2);
  if ixs<>nil then
    begin
    s:='Res1'+#9+'AA'+#9+'Res2'+#9+'AA';
    for f:=0 to High(Ixs) do
      s:=s+#9+RTAsLabel(FResults[ixs[f]].ResultType);
    h1:=High(FResults[ixs[0]].Data);
    h2:=High(FResults[ixs[0]].Data[0]);
    sl.Add(s);
    for g:=0 to h1 do
      for h:=0 to h2 do
        begin
        s:=IntToStr(g+1)+#9+FQueries[six1].Sequence[g+1]+#9+
            IntToStr(h+1)+#9+FQueries[six2].Sequence[h+1];
        for f:=0 to High(ixs) do
          s:=s+#9+FloatToStr(FResults[ixs[f]].Data[g,h]);
{        //FIX THIS!******
        if FResults[ixs[High(ixs)]].Data[g,h]>2 then
        //FIX THIS!******}
          sl.Add(s);
        end;
    end;
end;

procedure TComparer.CrossDistance(sid1, sid2: string);
var
  six1,six2:Integer;
  dists:TOCMatrix;
  task:TRunningTask;

procedure CalcDists;

var
  f,g:Integer;
  s:Real;

begin
  if dists<>nil then
    s:=1/Length(dists)
  else s:=0;
  for f:=0 to High(dists) do
    begin
    for g:=0 to High(dists[f]) do
      dists[f,g]:=Distance(FQueries[six1].ResCoords[f],FQueries[six2].ResCoords[g]);
    task.Step(s);
    end;
end;

begin
  six1:=QueryIxById(sid1);
  six2:=QueryIxById(sid2);
  if (six1>=0) and (six2>=0) then
    begin
    task:=NewTask(False,'Cross distances '+sid1+'/'+sid2);
    try
      SetLength(dists,Length(FQueries[six1].Sequence),Length(FQueries[six2].Sequence));
      CalcDists;
      AddResult(sid1,sid2,RSLTCrossDistance,dists);
    finally
      FreeTask(task);
    end;
    end;
end;

procedure TComparer.CrossMutualInf(sid1, sid2: string;Cutoff:Single; MaxNeibs:Integer);

var
  six1,six2:Integer;
  seqs1,seqs2:TOCStrings;
  mis:TOCMatrix;
  task:TRunningTask;
  ne1,ne2:TOCIntegers;

procedure ListSequences;

var
  f,ix:Integer;

begin

  seqs1:=nil;
  seqs2:=nil;

  //TO DO: This is duplicated in MergeAsMSA. Use that method instead

  for f:=0 to High(FQueries[six1].OrganismIDs) do
    begin
    ix:=LastIndexOf(FQueries[six1].OrganismIDs[f],FQueries[six2].OrganismIDs);
    if ix>=0 then
      begin
      AddToArray(FQueries[six1].MultiAlign[f],seqs1);
      AddToArray(FQueries[six2].MultiAlign[ix],seqs2);
      end;
    end;
  seqs1:=Transpose(seqs1);
  seqs2:=Transpose(seqs2);

end;

procedure CalcMis;

var
  f,g:Integer;
  s:Real;
begin
  if seqs1<>nil then s:=1/Length(seqs1)
  else s:=0;
  for f:=0 to High(seqs1) do
    begin
    for g:=0 to High(seqs2) do
      if (Cutoff<0) or ((ne1[f]<=MaxNeibs) and (ne2[g]<=MaxNeibs)) then
        mis[f,g]:=MutualInformation(seqs1[f],seqs2[g])
      else mis[f,g]:=0;
    task.Step(s);
    end;
end;


begin
  six1:=QueryIxById(sid1);
  six2:=QueryIxById(sid2);
  if Cutoff>0 then
    begin
    ne1:=NeighbourCount(six1,Cutoff);
    ne2:=NeighbourCount(six2,Cutoff);
    end;
  if (six1>=0) and (six2>=0) then
    begin
    task:=NewTask(False,'Mutual information '+sid1+'/'+sid2);
    try
      ListSequences;
      if seqs1<>nil then
        begin
        SetLength(mis,Length(FQueries[six1].Sequence),Length(FQueries[six2].Sequence));
        CalcMis;
        AddResult(sid1,sid2,RSLTCrossMI,mis);
        end;
      finally
        FreeTask(task);
      end;
    end;
end;

procedure TComparer.CrossInteract(sid1, sid2: string;Cutoff:Single; MaxNeibs:Integer);

var
  six1,six2:Integer;
  seqs1,seqs2:TOCStrings;
  mis:TOCMatrix;
  task:TRunningTask;
  ne1,ne2:TOCIntegers; //neighbour count

procedure ListSequences;

var
  f,ix:Integer;

begin
  seqs1:=nil;
  seqs2:=nil;
  for f:=0 to High(FQueries[six1].OrganismIDs) do
    begin
    ix:=LastIndexOf(FQueries[six1].OrganismIDs[f],FQueries[six2].OrganismIDs);
    if ix>=0 then
      begin
      AddToArray(FQueries[six1].MultiAlign[f],seqs1);
      AddToArray(FQueries[six2].MultiAlign[ix],seqs2);
      end;
    end;
  seqs1:=Transpose(seqs1);
  seqs2:=Transpose(seqs2);

end;

procedure CalcInters;

var
  f,g:Integer;
  s:Real;
begin
  if seqs1<>nil then s:=1/Length(seqs1)
  else s:=0;
  for f:=0 to High(seqs1) do
    begin
    for g:=0 to High(seqs2) do
      begin
      if (Cutoff<0) or ((ne1[f]<=MaxNeibs) and (ne2[g]<=MaxNeibs)) then
        mis[f,g]:=Interaction(seqs1[f],seqs2[g])
      else
        mis[f,g]:=0;
      end;
    task.Step(s);
    end;
end;


begin
  six1:=QueryIxById(sid1);
  six2:=QueryIxById(sid2);
  if Cutoff>0 then
    begin
    ne1:=NeighbourCount(six1,Cutoff);
    ne2:=NeighbourCount(six2,Cutoff);
    end;

  if (six1>=0) and (six2>=0) then
    begin
    task:=NewTask(False,'Cross Interaction '+sid1+'/'+sid2);
    try
      ListSequences;
      if seqs1<>nil then
        begin
        SetLength(mis,Length(FQueries[six1].Sequence),Length(FQueries[six2].Sequence));
        CalcInters;
        AddResult(sid1,sid2,RSLTCrossInteract,mis);
        end;
      finally
        FreeTask(task);
      end;
    end;
end;

procedure TComparer.ClearResults;
begin
  FResults:=nil;
end;

procedure TComparer.FilteredHistogram(sid1, sid2: string; filter,
  data,bins: Integer; filtermin, filtermax: Real; sl: TStringList; BinColumn:Boolean);

var
  f:Integer;
  lobin,upbin:TOCFloats;
  lo,up:TOCFloat;
  datam,filterm,filtered:TOCMatrix;
  hist:TOCIntegers;
  s:string;

procedure SetFiltered;

var
  f,g:Integer;
  r:TOCFloat;

begin
  Assert((length(datam)=length(filterm)) and (length(filterm[0])=Length(datam[0])),
    'Incorrect matrices');
  SetLength(filtered,Length(datam),Length(datam[0]));
  for f:=0 to High(datam) do
    for g:=0 to High(datam[0]) do
      begin
      r:=filterm[f,g];
      if (r<filtermin) or (r>filtermax) then
        filtered[f,g]:=lo-1 // to be outside the range of the lower bin
      else
        filtered[f,g]:=datam[f,g];
      end;
end;


begin
  lobin:=nil;
  upbin:=nil;
  datam:=FindFirstMatrix(sid1,sid2,data);
  filterm:=FindFirstMatrix(sid1,sid2,filter);
  if (filterm<>nil) and (datam<>nil) then
    begin
    lo:=Min(datam);
    up:=Max(datam);
    SetFiltered;
    CreateBins(lo,up,bins,lobin,upbin);
    hist:=Histogram(filtered,lobin,upbin);
    s:=FloatToStr(filtermin)+'-'+FloatToStr(filtermax);
    if BinColumn then s:='Bins'+#9+s;
    sl.Add(s);
    end;
    for f:=0 to High(hist) do
      begin
      s:=IntToStr(hist[f]);
      if BinColumn then
        s:=FloatToStr(lobin[f])+#9+s;
      sl.Add(s);
    end;
end;

procedure TComparer.LowestAndCount(sid1, sid2: string; low, count: Integer;
  countmin, countmax: Real; sl: TStringList);
var
  f:Integer;
  lom,countm:TOCMatrix;
  lowval:TOCFloat;
  countval:Integer;


begin
  lom:=FindFirstMatrix(sid1,sid2,low);
  countm:=FindFirstMatrix(sid1,sid2,count);
  if (lom<>nil) and (countm<>nil) then
    for f:=0 to High(lom) do
      begin
      lowval:=Min(lom[f]);
      countval:=CountBetween(countm[f],countmin,countmax);
      sl.Add(FloatToStr(lowval)+#9+IntToStr(countval));
    end;
end;

procedure TComparer.PlotData(sid1, sid2: string; datax, datay: Integer;
  bmp: TBitmap);

var
  f:Integer;
  xm,ym:TOCMatrix;
  minx,maxx,miny,maxy:TOCFloat;
  r:TRect;


begin
  xm:=FindFirstMatrix(sid1,sid2,datax);
  ym:=FindFirstMatrix(sid1,sid2,datay);
  if (ym<>nil) and (xm<>nil) then
    begin
    minx:=Min(xm);
    maxx:=Max(xm);
    miny:=Min(ym);
    maxy:=Max(ym);
    r:=Rect(0,0,bmp.Width-1,bmp.Height-1);
    for f:=0 to High(xm) do
      PlotXY(xm[f],ym[f],minx,maxx,miny,maxy,bmp,r,clBlue);
    end;
end;

function TComparer.GetMSA(Id: string): TMSA;
var
  queryix:Integer;

begin
  queryix:=QueryIxById(Id);
  if queryix>=0 then
    with FQueries[queryix] do
      begin
      Result.Alignment:=Copy(MultiAlign,0,High(MultiAlign));
      Result.SequenceIDs:=Copy(SeqIDs,0,High(SeqIDs));
      end;
end;

function TComparer.MergeAsMSA(Id1, Id2: string): TMSA;
var
  f,six,qix1,qix2:Integer;

begin
  qix1:=QueryIxById(Id1);
  qix2:=QueryIxById(Id2);
  Result.Alignment:=nil;
  Result.SequenceIds:=nil;
  if (qix1>=0) and (qix2>=0) then
  for f:=0 to High(FQueries[qix1].OrganismIDs) do
    begin
    six:=LastIndexOf(FQueries[qix1].OrganismIDs[f],FQueries[qix2].OrganismIDs);
    if six>=0 then
      begin
      AddToArray(FQueries[qix1].MultiAlign[f]+FQueries[qix2].MultiAlign[six],
        Result.Alignment);
      AddToArray(FQueries[qix1].SeqIds[f]+'+'+FQueries[qix2].SeqIds[six],
        Result.SequenceIds);
      end;
    end;
end;

end.

