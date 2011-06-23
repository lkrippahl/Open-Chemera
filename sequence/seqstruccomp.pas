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
 5. AlignToPDB should receive the SubMat, not the file name
 6. Indeces ids1 ids2 on TResult are only used for residues in PlotResults
    and crossdists Redo this...
*******************************************************************************}
unit seqstruccomp;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, fasta, pdbparser, clustalparser, alignment,
  ocstringutils, ebixmlparser, ebifastaparser, progress, chartutils, graphics,
  threedcalc, sequence, dictionaries, resinteraction, quicksort, debugutils, capsparser;

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
  RSLTPlot=6;

type
  //Information on query sequence and structure
  TQueryData=record
    ID:string; //assumed to be the ID of the query sequence in MSAs
    Sequence:string;
    ResCoords:TOCCoords;
    ResNames:TOCStrings;
    ResIDs:TOCIntegers;
    PDBSequence:string;
    PDBAlignment:TOCIntegers;
    MultiAlign:TOCStrings;
    //multialign is trimmed or expanded so the matrix matches the
    //length of the query sequence.  TODO: correct this.
    SeqIDs:TOCStrings;
    OrganismIDs:TOCStrings;
    GapMarker:char;
    end;

  TComparisonResult=record
    ID1,ID2:string;
    IDs1,IDs2:TOCIntegers; //identifiers for rows and columns, needed for plotresults
    ResultType:Integer;
    Data:TOCMatrix;
    Text:TOCStrings;
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
    FInterMat:TSubMatrix;


    function RTAsLabel(m:Integer):string;
    function Transpose(sl:TOCStrings):TOCStrings;
      //transposes an array of strings.
      //all strings are assumed to be the same length
    procedure AddQuery(nq:TQueryData);
    function QueryIxById(id:string):Integer;
    function DifferentInteraction(s1,s2:string):TOCFloat;       //s1 and s2 must be of equal length and aligned
      //average interaction scores for the different pairings
      //with 0 if only one pairing
      //****DOESN'T WORK WELL. TO REMOVE*****

    procedure AddResult(id1,id2:string;rt:Integer;d:TOCMatrix;Text:TOCStrings;IDs1,IDs2:TOCIntegers);
    function ListSelfResults(id:string):TOCIntegers;
    function ListCrossResults(id1,id2:string):TOCIntegers;
    function FindFirstMatrix(id1,id2:string;dt:Integer):TOCMatrix;
      //finds the first data matrix with these ids and type

    function NeighbourCount(QIx:Integer;Cutoff:Single):TOCIntegers;
    procedure SaveResult(Result:TComparisonResult; FileName:string);
    function AlignmentMap(Ix:Integer):TOCIntegers;
    function SequenceByColumn(ID:string;Filter:TOCIntegers):TOCStrings;
  public
    property ScaleMI:Boolean read FScaleMI write FScaleMI;
    constructor Create(DataFolder:string=''); //data folder without terminating slash
    procedure LoadInteractionPropensities(FileName:string);
    procedure ClearData;
    procedure LoadQuery(Id,SeqFile:string);
    procedure LoadPDB(Id, PdbFile: string);
    procedure AlignToPDB(Id,SubMatFile:string);
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
    procedure CrossMutualInf(sid1,sid2:string);
    procedure CrossInteract(sid1,sid2:string;IntMat:TSubMatrix);
    procedure ClearResults;
    procedure FilteredHistogram(sid1, sid2: string; filter,
              data,bins: Integer; filtermin, filtermax: Real; sl: TStringList;
              BinColumn:Boolean);
    procedure LowestAndCount(sid1,sid2:string;low,count:Integer;
      countmin,countmax:Real;sl:TStringList);
    procedure PlotResults(ResIx1,ResIx2:Integer;Min1,Max1,Min2,Max2:TOCFloat);
    procedure PlotMinDistAverageInt(DistIx,IntIx:Integer);
    procedure PlotData(sid1,sid2:string;datax,datay:Integer;bmp:TBitmap);
    function GetMSA(Id:string):TMSA;
    function GetFilteredMSA(ID:string; const Filter:TOCIntegers):TMSA;
    procedure MergeAsNewQuery(Id1,ID2,NewQuery:string);
      //merges by organism id, concatenating the sequences and the structure coords
    procedure MatchByOrganism(Id1,Id2:string);
      //filters both MSAs and sorts sequences to match organisms.
    procedure FilterSequencesByLength(ID:string; MinLen,MaxLen:Integer);
    procedure CleanGaps(ID:string);
    procedure CleanEmptySequences(ID:string);
    function GetPDBAlignment(ID:string):TOCIntegers;
    procedure BlankOutsideQuery(ID:string; Distance:Integer);
      //places gaps on all sequences in the alignment that are further than
      //distance from a non gap query position
     procedure BlankOutsidePDB(ID: string; Distance: Integer; PDBMap: TOCIntegers);
      //Same as BlankOutsideQuery, but using those that correspond to pdb
    procedure FlushResults(BaseFileName:string);
    procedure ParseCAPS(CapsFile,Pdb1,Pdb2, SubMatFile,ScriptFile:string); // TODO: quickfix; improve later

  end;

implementation

function EmptyQuery:TQueryData;

begin
  with Result do
    begin
    ID:='';
    Sequence:='';
    ResCoords:=nil;
    PDBSequence:='';
    PDBAlignment:=nil;
    MultiAlign:=nil;
    SeqIDs:=nil;
    OrganismIDs:=nil;
    GapMarker:='-';
    end;
end;

{ TComparer }

procedure TComparer.LoadInteractionPropensities(FileName: string);

begin
  FInterMat:=ReadBLASTMatrix(FileName);
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

function TComparer.DifferentInteraction(s1, s2: string): TOCFloat;
var
  f,g,c,p1,p2:Integer;
  diffs:array of array of Integer;

begin
  Result:=0;
  SetLength(diffs,Length(FInterMat.MonomerIndex),Length(FInterMat.MonomerIndex));
  for f:=0 to High(diffs) do
    for g:=0 to High(diffs) do
      diffs[f,g]:=0;
  c:=0;
  for f:=1 to Length(s1) do
    begin
    p1:=Pos(s1[f],FInterMat.MonomerIndex)-1;
    p2:=Pos(s2[f],FInterMat.MonomerIndex)-1;
    if (p1>=0) and (p2>=0) and (diffs[p1,p2]=0) then
      begin
      Inc(c);
      diffs[p1,p2]:=1;
      Result:=Result+FInterMat.Matrix[p1,p2];
      end;
    end;
  if c>1 then
    Result:=Result/sqrt(c)/4
  else Result:=0;

end;

procedure TComparer.AddResult(id1,id2:string; rt: Integer; d: TOCMatrix; Text:TOCStrings;IDs1,IDs2:TOCIntegers);

var ix:Integer;

begin
  ix:=Length(FResults);
  SetLength(FResults,ix+1);
  FResults[ix].ID1:=id1;
  FResults[ix].ID2:=id2;
  FResults[ix].IDs1:=IDs1;
  FResults[ix].IDs2:=IDs2;
  FResults[ix].Data:=d;
  FResults[ix].ResultType:=rt;
  FResults[ix].Text:=Text;
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
//TODO: deprecated, should not be used because of alignment between sequence and structure
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

procedure TComparer.SaveResult(Result: TComparisonResult; FileName: string);

var
  sl:TStringList;
  f,g:Integer;
  t,s:string;

begin
  sl:=TStringList.Create;
  with Result do
    begin
    sl.Add(ID1);
    sl.Add(ID2);
    sl.Add(IntToStr(ResultType));;
    for f:=0 to High(Data) do
      begin
      s:='';
      for g:=0 to High(Data[0]) do
        begin
        //TO DO: this is a quick fix
        t:=FloatToStr(Round(Data[f,g]*1000)/1000)+#9;
        if Pos('.',t)>0 then t[Pos('.',t)]:=',';
        s:=s+t;
        end;
      sl.Add(s);
      end;
    for f:=0 to High(Text) do
      sl.Add(Text[f]);
    end;
  sl.SaveToFile(FileName);
  sl.Free;
end;

function TComparer.AlignmentMap(Ix:Integer): TOCIntegers;

var
  seqix,f:Integer;
  tmp:TOCIntegers;
begin
  with FQueries[Ix] do
    begin
    seqix:=LastIndexOf(ID,SeqIDs);
    if seqix<0 then seqix:=0; //if not found, assume 0
    tmp:=AlignmentToSequence(MultiAlign[seqix],Sequence, GapMarker);
    Result:=nil;
    if PDBAlignment<>nil then
      begin
      for f:=0 to High(tmp) do
        if (tmp[f]>0) and (PDBAlignment[tmp[f]-1]>=0) then
          AddToArray(f+1,Result);
      end
    else
      for f:=0 to High(tmp) do
        if (tmp[f]>0) then
          AddToArray(f+1,Result);
    end;
end;

function TComparer.SequenceByColumn(ID: string; Filter: TOCIntegers): TOCStrings;

var
  tmpmsa:TMSA;

begin
  Result:=nil;
  tmpmsa:=GetFilteredMSA(ID,Filter);
  Result:=Transpose(tmpmsa.Alignment);
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

procedure TComparer.LoadQuery(Id, SeqFile: string);

var
  nq:TQueryData;
  fr:TFastaReader;

begin
  nq:=EmptyQuery;
  nq.ID:=id;
  nq.ResCoords:=nil;
  fr:=TFastaReader.Create(seqfile);
  if fr.Seqs<>nil then
    nq.Sequence:=fr.Seqs[0].Sequence;
  fr.Free;
  AddQuery(nq);
end;

procedure TComparer.LoadPDB(Id,PdbFile:string);

var
  pr:TPDBReader;
  ix,f,rix,ac:Integer;
  c:TOCCoord;
  resid:Integer;
  resname:string;

begin
  ix:=QueryIxById(Id);
  if ix>=0 then
  with FQueries[ix] do
    begin
    pr:=TPDBReader.Create(PdbFile);
    c:=NullVector;
    ResCoords:=nil;
    ResNames:=nil;
    ResIDs:=nil;
    if pr.Atoms<>nil then
      begin
      ac:=0;
      rix:=-1;
      for f:=0 to High(pr.Atoms) do
        begin
        if rix<>pr.Atoms[f].ResSeq then
          begin
          PDBSequence:=PDBSequence+OneLetterCode(pr.Atoms[f].ResName);
          if ac>0 then
            begin
            AddToArray(ScaleVector(c,1/ac),ResCoords);
            AddToArray(resname, ResNames);
            AddToArray(resid, ResIDs);
            end;
          ac:=0;
          rix:=pr.Atoms[f].ResSeq;
          c:=NullVector;
          resname:=pr.Atoms[f].ResName;
          resid:=pr.Atoms[f].ResSeq;

          end;
        inc(ac);
        c:=AddVectors(c,pr.Atoms[f].Coords);
        end;
      if ac>0 then
         begin
         AddToArray(ScaleVector(c,1/ac),ResCoords);
         AddToArray(resname, ResNames);
         AddToArray(resid, ResIDs);
         end;
      end;
    pr.Free;
    end;

end;

procedure TComparer.AlignToPDB(Id, SubMatFile: string);

var
  submat:TSubMatrix;
  queryseqixs,pdbseqixs:TOCIntegers;
  ix:Integer;

begin
  ix:=QueryIxById(Id);
  if ix>=0 then
  with FQueries[ix] do
    begin
    submat:=ReadBLASTMatrix(SubMatFile);
    queryseqixs:=IndexSequence(Sequence,SubMat);
    pdbseqixs:=IndexSequence(PDBSequence,SubMat);
    PDBAlignment:=NeedlemanWunschAlign(queryseqixs,pdbseqixs,submat);
    end;

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
    FQueries[ix].GapMarker:=msa.GapMarker;
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
      oix:=LastIndexOf(seqs[f].Code+'|'+seqs[f].ID,FQueries[six].SeqIDs);
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
    AddResult(seqid,'',RSLTSelfMI,mis, nil,nil,nil);
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
    AddResult(seqid,'',RSLTSelfDistance,dists, nil,nil,nil);
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
  cix1,cix2:TOCIntegers;
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
      dists[f,g]:=Distance(FQueries[six1].ResCoords[cix1[f]],
        FQueries[six2].ResCoords[cix2[g]]);
    task.Step(s);
    end;
end;

function IndexCoords(Align:TOCIntegers):TOCIntegers;

var f:Integer;

begin
  Result:=nil;
  for f:=0 to High(Align) do
    if Align[f]>=0 then AddToArray(Align[f],Result);
end;

var f:Integer;

begin
  six1:=QueryIxById(sid1);
  six2:=QueryIxById(sid2);
  if (six1>=0) and (six2>=0) then
    begin
    task:=NewTask(False,'Cross distances '+sid1+'/'+sid2);
    try
      cix1:=IndexCoords(FQueries[six1].PDBAlignment);
      cix2:=IndexCoords(FQueries[six2].PDBAlignment);
      SetLength(dists,Length(cix1),Length(cix2));
      CalcDists;
      for f:=0 to High(cix1) do cix1[f]:=FQueries[six1].ResIDs[cix1[f]];
      for f:=0 to High(cix2) do cix2[f]:=FQueries[six2].ResIDs[cix2[f]];
      AddResult(sid1,sid2,RSLTCrossDistance,dists, nil,cix1,cix2);
    finally
      FreeTask(task);
    end;
    end;
end;

procedure TComparer.CrossMutualInf(sid1, sid2: string);

var
  six1,six2:Integer;
  filter1,filter2:TOCIntegers;
  seqs1,seqs2:TOCStrings;
  mis:TOCMatrix;
  task:TRunningTask;

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
      mis[f,g]:=MutualInformation(seqs1[f],seqs2[g]);
    task.Step(s);
    end;
end;

begin
  six1:=QueryIxById(sid1);
  six2:=QueryIxById(sid2);
  if (six1>=0) and (six2>=0) then
    begin
    seqs1:=SequenceByColumn(sid1,AlignmentMap(six1));
    seqs2:=SequenceByColumn(sid2,AlignmentMap(six2));
    task:=NewTask(False,'Mutual information '+sid1+'/'+sid2);
    try
      if seqs1<>nil then
        begin
        SetLength(mis,Length(seqs1),Length(seqs2));
        CalcMis;
        seqs1:=Transpose(seqs1);
        seqs2:=Transpose(seqs2);
        AddToArray(sid2+':',seqs1);
        seqs1:=Concatenate(seqs1,seqs2);
        AddResult(sid1,sid2,RSLTCrossMI,mis,seqs1,nil,nil);
        end;
      finally
        FreeTask(task);
      end;
    end;
end;

procedure TComparer.CrossInteract(sid1, sid2: string;IntMat:TSubMatrix);

var
  six1,six2:Integer;
  seqs1,seqs2:TOCStrings;
  mis:TOCMatrix;
  task:TRunningTask;
  ne1,ne2:TOCIntegers; //neighbour count

function RandomSeq(s:string):string;

var
  f,ix:Integer;
begin
  Result:=s;
  for f:=1 to Length(s) do
    begin
    ix:=Random(Length(s)-f+1)+f;
    Result[f]:=s[ix];
    s[ix]:=s[f];
    end;
end;

procedure CalcInters;

var
  f,g,c:Integer;
  s,tmp:TOCFloat;
begin
  if seqs1<>nil then s:=1/Length(seqs1)
  else s:=0;
  for f:=0 to High(seqs1) do
    begin
    for g:=0 to High(seqs2) do
      begin
      mis[f,g]:=AverageInteraction(seqs1[f],seqs2[g],IntMat,FQueries[six1].GapMarker);
      //Quick fix for subtracting random baseline
      tmp:=0;
      for c:=1 to 10 do
        tmp:=tmp+AverageInteraction(seqs1[f],RandomSeq(seqs2[g]),IntMat,FQueries[six1].GapMarker);
      mis[f,g]:=mis[f,g]-tmp/10;
      //DEBUG
      if FResults[0].Data[f,g]<10 then
        begin
        DebugReport('Seq1:'+seqs1[f]);
        DebugReport('NotR:'+seqs2[g]);
        DebugReport('Rand:'+RandomSeq(seqs2[g]));
        DebugReport(tmp/10);
        DebugReport(mis[f,g]);
        end;
      end;
    task.Step(s);
    end;
end;


begin
  six1:=QueryIxById(sid1);
  six2:=QueryIxById(sid2);
  if (six1>=0) and (six2>=0) then
    begin
    seqs1:=SequenceByColumn(sid1,AlignmentMap(six1));
    seqs2:=SequenceByColumn(sid2,AlignmentMap(six2));
    task:=NewTask(False,'Cross Interaction '+sid1+'/'+sid2);
    try
      if seqs1<>nil then
        begin
        SetLength(mis,Length(seqs1),Length(seqs2));
        CalcInters;
        AddResult(sid1,sid2,RSLTCrossInteract,mis,nil,nil,nil);
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

procedure TComparer.PlotResults(ResIx1, ResIx2: Integer;Min1,Max1,Min2,Max2:TOCFloat);
//Result matrixes must have equal dimensions


var
  d1,d2,cross:TOCMatrix;
  f,g,c:Integer;
  x1,x2:TOCFloat;

begin
  d1:=FResults[ResIx1].Data;
  d2:=FResults[ResIx2].Data;

  c:=0;
  SetLength(cross,Length(d1)*Length(d1[0]),4);
  for f:=0 to High(d1) do
    for g:=0 to High(d1[0]) do
      begin
      x1:=d1[f,g];
      x2:=d2[f,g];
      if (x1>=Min1) and (x1<=Max1) and (x2>=Min2) and (x2<=Max2) then
        begin
        cross[c,0]:=f+1;
        cross[c,1]:=g+1;
        if FResults[ResIx1].IDs1<>nil then
          with FResults[ResIx1] do
            begin
            cross[c,0]:=Ids1[f];
            cross[c,1]:=IDs2[g];
            end;
        cross[c,2]:=x1;
        cross[c,3]:=x2;
        Inc(c);
        end;
      end;
  SetLength(cross,c,4);
  AddResult(FResults[ResIx1].ID1,FResults[ResIx1].ID1,RSLTPlot,cross,nil,nil,nil);
end;

procedure TComparer.PlotMinDistAverageInt(DistIx, IntIx: Integer);
//Result matrixes must have equal dimensions
//calculates minumum distance and average interaction score for each position

var
  d1,d2,cross:TOCMatrix;
  f,g,c:Integer;
  xd,xi:TOCFloat;
  vals:TOCFloats;
  sortedixs:TOCIntegers;

begin
  d1:=FResults[DistIx].Data;
  d2:=FResults[IntIx].Data;
  c:=0;
  SetLength(cross,Length(d1)+Length(d1[0]),4);
  //horizontal
  for f:=0 to High(d1) do
    begin
    xd:=d1[f,0];
    xi:=d2[f,0];
    vals:=nil;
    AddToArray(-xi,vals);
    for g:=1 to High(d1[0]) do
      begin
      if d1[f,g]<xd then xd:=d1[f,g];
      //xi:=xi+d2[f,g];
      AddToArray(-d2[f,g],vals);
      end;
    cross[c,0]:=0;
    cross[c,2]:=xd;
    //cross[c,3]:=xi/Length(d1[0]);
    sortedixs:=QSAscendingIndex(vals);
    cross[c,3]:=0;
    for g:=0 to 10 do
      cross[c,3]:=cross[c,3]-vals[sortedixs[g]];
    cross[c,1]:=FQueries[0].ResIds[FQueries[0].PDBAlignment[f]]; //FIX THIS!!!
    if FResults[DistIx].IDs1<>nil then
      cross[c,1]:=FResults[DistIx].IDs1[f];
    Inc(c);
    end;
    //vertical

  for g:=0 to High(d1[0]) do
    begin
    xd:=d1[0,g];
    xi:=d2[0,g];
    vals:=nil;
    AddToArray(-xi,vals);
    for f:=1 to High(d1) do
      begin
      if d1[f,g]<xd then xd:=d1[f,g];
      xi:=xi+d2[f,g];
      AddToArray(-d2[f,g],vals);
      end;
    cross[c,0]:=1;
    cross[c,2]:=xd;
    //cross[c,3]:=xi/Length(d1);
    sortedixs:=QSAscendingIndex(vals);
    cross[c,3]:=0;
    for f:=0 to 10 do
      cross[c,3]:=cross[c,3]-vals[sortedixs[f]];
    cross[c,1]:=FQueries[1].ResIds[FQueries[1].PDBAlignment[g]]; //FIX THIS!!!
    if FResults[DistIx].IDs2<>nil then
      cross[c,1]:=FResults[DistIx].IDs2[g];
    Inc(c);
    end;
  AddResult(FResults[DistIx].ID1,FResults[DistIx].ID2,RSLTPlot,cross,nil,nil,nil);
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
      Result.GapMarker:=GapMarker;
      Result.Alignment:=Copy(MultiAlign,0,Length(MultiAlign));
      Result.SequenceIDs:=Copy(SeqIDs,0,Length(SeqIDs));
      end;
end;

function TComparer.GetFilteredMSA(ID: string; const Filter:TOCIntegers): TMSA;
//filter is the indexes for each sequence, base 1

var
  f:Integer;

begin
  Result:=GetMSA(ID);
  for f:=0 to High(Result.Alignment) do
    Result.Alignment[f]:=SequenceFromIndex(Result.Alignment[f],Filter);
end;

procedure TComparer.MergeAsNewQuery(Id1,ID2,NewQuery:string);
var
  f,six,nqix,qix1,qix2:Integer;
  nq:TQueryData;

begin
  qix1:=QueryIxById(Id1);
  qix2:=QueryIxById(Id2);
  nqix:=QueryIxById(NewQuery);
  if (qix1>=0) and (qix2>=0) and (nqix<0) then
  begin
  nq.ID:=NewQuery;
  nq.ResCoords:=Concatenate(FQueries[qix1].ResCoords,FQueries[qix2].ResCoords);
  nq.GapMarker:=FQueries[qix1].GapMarker;
  nq.MultiAlign:=nil;
  nq.OrganismIDs:=nil;
  nq.SeqIDs:=nil;
  nq.Sequence:=FQueries[qix1].Sequence+FQueries[qix2].Sequence;
  for f:=0 to High(FQueries[qix1].OrganismIDs) do
    begin
    six:=LastIndexOf(FQueries[qix1].OrganismIDs[f],FQueries[qix2].OrganismIDs);
    if six>=0 then
      begin
      AddToArray(FQueries[qix1].OrganismIDs[f],nq.OrganismIDs);
      AddToArray(FQueries[qix1].MultiAlign[f]+FQueries[qix2].MultiAlign[six],
        nq.MultiAlign);
      AddToArray(FQueries[qix1].SeqIds[f]+'+'+FQueries[qix2].SeqIds[six],
        nq.SeqIds);
      end;
    end;
  AddQuery(nq);
  end;
end;

procedure TComparer.MatchByOrganism(Id1, Id2: string);

var
  qix1,qix2:Integer;
  ix1,ix2:TOCIntegers;

procedure BuildIndexes;

var
  qs1,qs2,f,g,six:Integer;

begin
  //place queryid in first position
  qs1:=LastIndexOf(FQueries[qix1].ID,FQueries[qix1].SeqIDs);
  qs2:=LastIndexOf(FQueries[qix2].ID,FQueries[qix2].SeqIDs);
  if (qs1>=0) and (qs1>=0) then
    begin
    AddToArray(qs1,ix1);
    AddToArray(qs2,ix2);
    end;
  for f:=0 to High(FQueries[qix1].OrganismIDs) do
    if f<>qs1 then
    begin
    six:=LastIndexOf(FQueries[qix1].OrganismIDs[f],FQueries[qix2].OrganismIDs);
    if (six>=0) and (six<>qs2) then
      begin
      AddToArray(f,ix1);
      AddToArray(six,ix2)
      end;
    end;
end;

procedure FilterAlignment(Qix:Integer;Index:TOCIntegers);

var
  f:Integer;
  orgs,aligns,ids:TOCStrings;

begin
  SetLength(orgs,Length(Index));
  SetLength(aligns,Length(Index));
  SetLength(ids,Length(Index));
  with FQueries[Qix] do
    begin
    for f:=0 to High(Index) do
      begin
      orgs[f]:=OrganismIds[Index[f]];
      aligns[f]:=MultiAlign[Index[f]];
      ids[f]:=SeqIDs[Index[f]];
      end;
    OrganismIds:=orgs;
    MultiAlign:=aligns;
    SeqIDs:=ids;
    end;
end;

begin
  qix1:=QueryIxById(Id1);
  qix2:=QueryIxById(Id2);
  if (qix1<0) or (qix2<0) then Exit;
  ix1:=nil;
  ix2:=nil;
  BuildIndexes;
  FilterAlignment(qix1,ix1);
  FilterAlignment(qix2,ix2);
end;

procedure TComparer.FilterSequencesByLength(ID: string; MinLen,
  MaxLen: Integer);

var
  qix,f,count:Integer;
  align,seqs,orgs:TOCStrings;

begin
  qix:=QueryIxById(ID);
  align:=nil;
  seqs:=nil;
  orgs:=nil;
  if qix>=0 then
    with FQueries[qix] do
      begin
      for f:=0 to High(MultiAlign) do
        begin
        count:=Length(MultiAlign[f])-CountInString(MultiAlign[f],GapMarker);
        if (count<=MaxLen) and (count>=MinLen) then
          begin
          AddToArray(MultiAlign[f],align);
          AddToArray(SeqIds[f],seqs);
          if OrganismIds<>nil then AddToArray(OrganismIds[f],orgs);
          end;
        end;
      MultiAlign:=align;
      SeqIds:=seqs;
      OrganismIds:=orgs;
      end;
end;

procedure TComparer.CleanGaps(ID: string);

var
  qix,f:Integer;
  colstokeep:TOCIntegers;

function IsToKeep(Ix:Integer):Boolean;

var g:Integer;

begin
  Result:=False;
  with FQueries[qix] do
    for g:=0 to High(MultiAlign) do
      if MultiAlign[g,Ix]<>GapMarker then
        Exit(True);
end;

function Selecteds(S:string):string;

var g:Integer;

begin
  SetLength(Result,Length(colstokeep));
  for g:=0 to High(colstokeep) do
    Result[g+1]:=S[colstokeep[g]];
end;

begin
  qix:=QueryIxById(ID);
  colstokeep:=nil;
  if (qix>=0) and (FQueries[qix].MultiAlign<>nil) then
    with FQueries[qix] do
      begin
      for f:=1 to Length(MultiAlign[0]) do
        if IsToKeep(f) then
          AddToArray(f,colstokeep);
      for f:=0 to High(MultiAlign) do
        MultiAlign[f]:=Selecteds(MultiAlign[f]);
      end;
end;

procedure TComparer.CleanEmptySequences(ID: string);

var
  newalign,neworgs,newids:TOCStrings;
  qix,f:Integer;


begin
  qix:=QueryIxById(ID);
    if qix>=0 then
      with FQueries[qix] do
        begin
        newalign:=nil;
        neworgs:=nil;
        newids:=nil;
        for f:=0 to High(MultiAlign) do
          if CountInString(MultiAlign[f],GapMarker)<>Length(MultiAlign[f]) then
            begin
            AddToArray(MultiAlign[f],newalign);
            AddToArray(SeqIDs[f],newids);
            AddToArray(OrganismIDs[f],neworgs);
            end;
        MultiAlign:=newalign;
        SeqIDs:=newids;
        OrganismIDs:=neworgs;
        end;
end;

function TComparer.GetPDBAlignment(ID: string): TOCIntegers;

var qix:Integer;

begin
  Result:=nil;
  qix:=QueryIxById(ID);
  if qix>=0 then
    Result:=FQueries[qix].PDBAlignment;
end;

procedure TComparer.BlankOutsideQuery(ID: string; Distance: Integer);

var
  queryix,ix:Integer;
  s:string;

function IsToBlank(Position:Integer):Boolean;

var f1,f2,f:Integer;

begin
  Result:=True;
  f1:=Position-Distance;
  if f1<1 then f1:=1;
  f2:=Position+Distance;
  if f2>Length(s) then f2:=Length(s);
  for f:=f1 to f2 do
    if s[f]<>FQueries[ix].GapMarker then
      Exit(False);
end;

var f,g:Integer;

begin
  ix:=QueryIxById(ID);
  if ix>=0 then
    with FQueries[ix] do
      begin
      queryix:=LastIndexOf(ID,SeqIDs);

      if queryix<0 then queryix:=0;
        //if ID doesn´t match, assume it's the first one

      s:=MultiAlign[queryix];
      for f:=1 to Length(s) do
        if IsToBlank(f) then
          for g:=0 to High(MultiAlign) do
            MultiAlign[g,f]:=GapMarker;
      end;

end;

procedure TComparer.BlankOutsidePDB(ID: string; Distance: Integer; PDBMap:TOCIntegers);

var
  queryix,ix:Integer;
  s:string;
  amap:TOCIntegers;

function IsToBlank(Position:Integer):Boolean;

var f1,f2,f:Integer;

begin
  Result:=True;
  f1:=Position-Distance;
  if f1<1 then f1:=1;
  f2:=Position+Distance;
  if f2>Length(s) then f2:=Length(s);
  for f:=f1 to f2 do
    if (amap[f-1]>=0) and (PDBMap[amap[f-1]]>=0) then
      Exit(False);
end;

procedure MapAlignment(Seq:string;GapMarker:Char);

var f,seqix:Integer;

begin
  SetLength(amap,Length(Seq));
  seqix:=0;
  for f:=1 to Length(Seq) do
    if Seq[f]=GapMarker then
      amap[f-1]:=-1
    else
      begin
      amap[f-1]:=seqix;
      Inc(seqix);
      end;
end;

var f,g:Integer;

begin
  ix:=QueryIxById(ID);
  if ix>=0 then
    with FQueries[ix] do
      begin
      queryix:=LastIndexOf(ID,SeqIDs);
      if queryix<0 then queryix:=0;
        //if ID doesn´t match, assume it's the first one

      s:=MultiAlign[queryix];
      MapAlignment(s,GapMarker);
      for f:=1 to Length(s) do
        if IsToBlank(f) then
          for g:=0 to High(MultiAlign) do
            MultiAlign[g,f]:=GapMarker;
      end;
end;

procedure TComparer.FlushResults(BaseFileName: string);

var f:Integer;

begin
  for f:=0 to High(FResults) do
    SaveResult(FResults[f],BaseFileName+IntToStr(f)+'.rslt');
end;

procedure TComparer.ParseCAPS(CapsFile, Pdb1, Pdb2, SubMatFile, ScriptFile: string);

var submat:TSubMatrix;

function AlignPDB(PdbFile:string;Seq:string):TOCIntegers;

var
  pr:TPDBReader;
  ix,f,rix,ac:Integer;
  resnames:TOCStrings;
  resid:Integer;
  resids:TOCIntegers;
  pdbsequence,resname:string;
  pdbseqixs,seqixs:TOCIntegers;

begin
  pr:=TPDBReader.Create(PdbFile);
  resnames:=nil;
  pdbsequence:='';
  resids:=nil;
  if pr.Atoms<>nil then
      begin
      ac:=0;
      rix:=-1;
      for f:=0 to High(pr.Atoms) do
        begin
        if rix<>pr.Atoms[f].ResSeq then
          begin
          pdbsequence:=pdbsequence+OneLetterCode(pr.Atoms[f].ResName);
          if ac>0 then
            begin
            AddToArray(resname, resnames);
            AddToArray(resid, resids);
            end;
          ac:=0;
          rix:=pr.Atoms[f].ResSeq;
          resname:=pr.Atoms[f].ResName;
          resid:=pr.Atoms[f].ResSeq;
          end;
        inc(ac);
        end;
      if ac>0 then
         begin
         AddToArray(resname, resnames);
         AddToArray(resid, resids);
         end;
      end;
    pr.Free;
    seqixs:=IndexSequence(Seq,SubMat);
    pdbseqixs:=IndexSequence(pdbsequence,SubMat);
    Result:=NeedlemanWunschAlign(seqixs,pdbseqixs,submat);
    //change mapping to residue IDs
    for f:=0 to High(Result) do
      if Result[f]>=0 then Result[f]:=resids[Result[f]];
end;


var
  sl:TStringList;
  q1,q2:TQueryData;
  ixs1,ixs2:TOCIntegers;
  s1,s2:string;
  map1,map2:TOCIntegers;
  f:Integer;

begin
  QuickCAPSParser(CapsFile,s1,s2,ixs1,ixs2);
  q1.Sequence:=s1;
  q2.Sequence:=s2;
  DebugReport(ixs1);
  DebugReport(ixs2);
  submat:=ReadBLASTMatrix(SubMatFile);
  map1:=AlignPdb(Pdb1,q1.Sequence);
  map2:=AlignPdb(Pdb2,q2.Sequence);
  DebugReport(map1);
  DebugReport(map2);

  sl:=TStringList.Create;
  //assumes layers are loaded in same order
  sl.Add('Desel');
  for f:=0 to high(ixs1) do
    if map1[ixs1[f]]>=0 then
    sl.Add('sel l:1 r:'+IntToStr(map1[ixs1[f]])+'  //'+
      IntToStr(ixs1[f]+1)+':'+q1.Sequence[ixs1[f]+1]);
  sl.Add('sel store:layer1');
  sl.Add('Desel');
  for f:=0 to high(ixs2) do
    if map2[ixs2[f]]>=0 then
    sl.Add('sel l:2 r:'+IntToStr(map2[ixs2[f]])+'  //'+
    IntToStr(ixs2[f]+1)+':'+q2.Sequence[ixs2[f]+1]);
  sl.Add('sel store:layer2');
  sl.Add('Desel');
  sl.SaveToFile(ScriptFile);
  sl.Free;
end;



end.

