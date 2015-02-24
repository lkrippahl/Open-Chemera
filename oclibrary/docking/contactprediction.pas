{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain. Enjoy.
********************************************************************************
Author: Ludwig Krippahl
Date: 5.9.2013
Purpose:
  Prediction of residue contacts in protein complexes.
  Currently generation of potential descriptors.
  Can also be used to encapsulate surface analysis

Requirements:
Revisions:
To do:
  Implement prediction after training SVM
*******************************************************************************}
unit contactprediction;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, INIFiles, basetypes, surface, pdbmolecules,
  molecules, geomhash, quicksort, pdbparser, oclconfiguration,
  alignment, sequence, stringutils, progress, pdbsurface;

const

  TargetIndex=0;
  ProbeIndex=1;

  //descriptor text
  MinSC='Minimum Sidechain ASA';
  MaxSC='Maximum Sidechain ASA';
  MinFull='Minimum  Full ASA';
  MaxFull='Maximum Ful ASA';
  CrossSCContact='Sidechain Contact ASA';
  CrossFullContact='Full Contact ASA';
  MaxGapFraction='Max Gap Fraction';
  MinGapFraction='Min Gap Fraction';
  MaxRelGapFraction='Max Relative Gap Fraction';
  MinRelGapFraction='Min Relative Gap Fraction';
  MaxSubFraction='Max Substitution Fraction';
  MinSubFraction='Min Substitution Fraction';
  MaxRelSubFraction='Max Rel Substitution Fraction';
  MinRelSubFraction='Min Rel Substitution Fraction';
  ContactFraction='Contact Score Fraction';
  RelContactFraction='Relative Contact Score Fraction';
  ScotchScore='SCOTCH Score'; //Madaoui, H. & Guerois, R. (2008), 'Coevolution at protein complex interfaces can be detected by the complementarity trace with important impact for predictive docking

type

  TResiduePairs=array[0..1] of TIntegers;

  TCDParameters=record
    TargetPDBFile,ProbePDBFile,TargetPDBID,ProbePDBID:string;

    TargetIDs,ProbeIDs:TSimpleStrings;
    TargetIndexes,ProbeIndexes:TIntegers;
      //note: indexes take precedence; ids are ignored if indexes are not nil

    ResidueMSAFiles:array of array[0..1] of TSimpleStrings;
      //each Must have same number of elements as IDs or Indexes
    MSADatabaseNames:TSimpleStrings;

    ContactMatrixFiles:TSimpleStrings;
    SubMatrixFile:string;

    SurfaceCutoff:TFloat;
      //Used to determine if residue is surface (minimum All surface)

    //Surface analysis options (ignored in descriptor generation)
    DoInnerContacts,DoCrossContacts:Boolean;

    //Report options
    WriteHeaders:Boolean;
    DescriptorFile,RawFile:string;

    //Distance computation, true to use distances instead of surface
    ComputeContactsByDistance:Boolean;
    ContactDistance:TFloat;


  { TODO : Add lists for excluded and useonly amino acids }
  end;

  TContactDescriptor=record
    Values:TMatrix;
    Id:string;
  end;
  TContactDescriptors=array of TContactDescriptor;

  TDescriptorResults=record
    Chains:array[0..1] of TMolecules;
    Residues:array[0..1] of TMolecules;
    SelectedIXs:array[0..1] of TIntegers;
      //Indexes of selected residues in original sequences.
    Selected:array[0..1] of TMolecules;
      //Selected residues
    Neighbours:array[0..1] of array of TIntegers;
      //These arrays index the selected residues only
    ResidueMSAs:array of array[0..1] of TSimpleStrings;
      //MSA columns for selected residues, for different databases
    CompressedMSAs:array of array [0..1] of TOCCompressedSequences;
    MSADBNames:TSimpleStrings;
      //Suffix to add to contact descriptor name for each database
    ContactMatrix,SubMatrix:TSubMatrix;
    Descriptors:TContactDescriptors;
  end;

  { TDescriptorGen }
  //Generates contact descriptors
  //also serves for surface analysis, neighbours and so forth
  //always assumes 2 structures at most
  //Sequences in MSA files must be correctly sorted and matched by organism
  //Also must correspond to the pdb sequences with X filling missing residue IDs
  //Residue IDs are assumed to be sequential when mapping to sequence
  { TODO : Extend to N structures, for N*(N-1) cross contacts. Should not be hard to do}
  TDescriptorGen=class
  protected
    FParameters:TCDParameters;
    FModelMan:TPDBModelMan;
    FPDBSurf:TPDBSurface;
    //All surface and contact values are indexed to the full set of used residues
    //(the residue lists in FSurfaceData)
    FSurfaceData:TSurfaceReports;
      //target and probe, lists only used residues in surface computations

    FCrossContacts:TCrossContactReport;
    //These are nil if no cross contacts requested

    // output values
    FDescriptors:TDescriptorResults;

    procedure CompressMSAs;
    procedure SetChainSurface(Index:Integer);
    procedure SelectSurfaceResidues(Index:Integer);
    procedure BuildNeighbourIndexes(Index:Integer);
    procedure SetInnerContacts(Index:Integer);
    procedure SetCrossContacts;
    procedure BuildMSAStrings(DBIndex, Index: Integer; Files: TSimpleStrings);
    procedure SetContactDescriptor(Ix:Integer;const Desc: TContactDescriptor);
    procedure AddContactDescriptor(const Desc:TContactDescriptor);
    procedure AddMinMaxSurface;
    procedure AddCrossContacts;

    // sequence scores. Index identifies the database for the sequences
    { TODO : Add weighted versions of each }
    procedure AddSequenceConservation(DBIndex:Integer);
    procedure AddGapFraction(DBIndex:Integer);
    procedure AddContactScores(DBIndex:Integer;CMLabel:string);
    procedure AddScotchScore(DBIndex:Integer);

    // patch (neighbour) scores
    procedure AddPatchScores;
      //adds patch version of all descriptors


  public
    property Descriptors:TDescriptorResults read FDescriptors;
    procedure ClearAll;
    procedure LoadModels(Parameters:TCDParameters);
      //clears everything, loads models, sets parameters;

    function GetSurfaceData(StructureIx:Integer):TSurfaceReport;
      //admitted residues are all used for surface calculations

    function FullInnerContacts(StructureIx:Integer):TMatrix;
    function SidechainInnerContacts(StructureIx:Integer):TMatrix;
    function FullCrossContacts:TMatrix;
    function SidechainCrossContacts:TMatrix;


    constructor Create(NumSurfacePoints:Integer);
    destructor Destroy;override;
    procedure GenerateAllDescriptors( const Parameters:TCDParameters; Verbose:Boolean=True);
      //Parameters MUST include valid files for both Target and Probe
    procedure ReportDescriptors(Report:TStrings;Headers:Boolean=false;OnlyContacts:Boolean=False);
    procedure SaveDescriptors(Headers:Boolean=false;OnlyContacts:Boolean=False);
      //Saves to output file specified in parameters
    procedure SaveRaw;
      //Saves raw data to output file specified in parameters

    procedure DistanceDescriptors(const Parameters: TCDParameters);
  end;

function ReadParameters(FileName:string):TCDParameters;
function EmptyDescriptorResults:TDescriptorResults;

implementation



function ReadParameters(FileName:string):TCDParameters;

const
  Separator=',';

var
  f:Integer;
  ini:TINIFile;
  pdbfolder:string;
  sequencefolders:TSimpleStrings;

  procedure SelectedChains(
      IDString,IndexString:string;
      out ChainIds:TSimpleStrings;
      out ChainIndexes:TIntegers);

  var
    ixs:TSimpleStrings;
    f:Integer;

  begin
    ChainIds:=SplitString(IDString,Separator);
    ixs:=SPlitString(IndexString,Separator);
    SetLength(ChainIndexes,Length(ixs));
    for f:=0 to High(ixs) do
      ChainIndexes[f]:=StrToInt(ixs[f])-1;
  end;

  procedure AddMSAs(Index:Integer;SL:TStringList;Folders:TSimpleStrings);

  var f,g,ix:Integer;
      names,vals, files:TSimpleStrings;
  begin
    ReadKeysValues(AsSimpleStrings(SL),names,vals);
    for f:=0 to High(names) do
      begin
      ix:=FirstIndexOf(names[f],Result.MSADatabaseNames);
      files:=SplitString(vals[f],Separator);
      for g:=0 to High(files) do
        files[g]:=IncludeTrailingPathDelimiter(Folders[ix])+files[g];
      Result.ResidueMSAFiles[ix,Index]:=files;
      end;
  end;



var s:string;
    tmp:TStringList;
    folders,names,vals:TSimpleStrings;

begin
  ini:=TINIFile.Create(FileName);

  pdbfolder:=ini.ReadString('Structure','PDBFolder','');
  if pdbfolder<>'' then pdbfolder:=IncludeTrailingPathDelimiter(pdbfolder);

  tmp:=TStringList.Create;

  with Result do
    begin
    SurfaceCutoff:=ini.ReadFloat('Structure','SurfaceCutoff',0);
    TargetPDBFile:=pdbfolder+ini.ReadString('Structure','TargetFile','');
    ProbePDBFile:=pdbfolder+ini.ReadString('Structure','ProbeFile','');
    TargetPDBID:=pdbfolder+ini.ReadString('Structure','TargetID','');
    ProbePDBID:=pdbfolder+ini.ReadString('Structure','ProbeID','');

    ComputeContactsByDistance:=
          ini.ReadString('Structure','ComputeContactsByDistance','0')<>'0';
    ContactDistance:=
          StrToFloat(ini.ReadString('Structure','ContactDistance','0'));


    SelectedChains(ini.ReadString('Structure','TargetChainIDs',''),
                   ini.ReadString('Structure','TargetChainIndexes',''),
                   Result.TargetIDs,Result.TargetIndexes);
    SelectedChains(ini.ReadString('Structure','ProbeChainIDs',''),
                   ini.ReadString('Structure','ProbeChainIndexes',''),
                   Result.ProbeIDs,Result.ProbeIndexes);

    tmp:=TStringList.Create;
    ini.ReadSectionValues('SequenceFolders',tmp);
    ReadKeysValues(AsSimpleStrings(tmp),MSADatabaseNames,sequencefolders);
    SetLength(ResidueMSAFiles,Length(sequencefolders));

    ini.ReadSectionValues('TargetSequenceFiles',tmp);
    AddMSAs(0,tmp,sequencefolders);
    ini.ReadSectionValues('ProbeSequenceFiles',tmp);
    AddMSAs(1,tmp,sequencefolders);

    ContactMatrixFiles:=SplitString(ini.ReadString('ScoreFiles','Contacts',''),Separator);
    SubMatrixFile:=ini.ReadString('ScoreFiles','Substitutions','');
    DescriptorFile:=ini.ReadString('Output','Descriptors','');
    RawFile:=ini.ReadString('Output','Stats','');
    end;
  ini.Free;
end;

function EmptyDescriptorResults:TDescriptorResults;

var f:Integer;

begin
  with Result do
    begin
    for f:=0 to 1 do
      begin
      Chains[f]:=nil;
      SelectedIXs[f]:=nil;
      Selected[f]:=nil;
      Neighbours[f]:=nil;
      end;
    Descriptors:=nil;
    end;
end;


{ TDescriptorGen }

procedure TDescriptorGen.CompressMSAs;

var f,g:Integer;


begin
  with FDescriptors do
    begin
    SetLength(CompressedMSAs,Length(ResidueMSAs));
    for f:=0 to High(CompressedMSAs) do
      for g:=0 to High(CompressedMSAs[f]) do
        CompressedMSAs[f,g]:=CompressSequences(ResidueMSAs[f,g]);
    end;
end;

procedure TDescriptorGen.SetChainSurface(Index: Integer);

begin
  FSurfaceData[Index]:=FPDBSurf.SurfaceStatistics(FDescriptors.Chains[Index]);
end;

procedure TDescriptorGen.LoadModels(Parameters:TCDParameters);

function GetChains(Layer:Integer;IXs:TIntegers;IDs:TSimpleStrings):TMolecules;

begin
  if IXs<>nil then
    Result:=FModelMan.GetChains(Layer,IXs)
  else Result:=FModelMan.GetChains(Layer,IDs);
end;

var f:Integer;

begin
  ClearAll;
  FParameters:=Parameters;
  FModelMan.LoadLayer(Parameters.TargetPDBFile);
  FDescriptors.Chains[0]:=GetChains(0,Parameters.TargetIndexes,Parameters.TargetIDs);
  if Parameters.ProbePDBFile<>'' then
    begin
    if Parameters.TargetPDBFile=Parameters.ProbePDBFile then
      FDescriptors.Chains[1]:=GetChains(0,Parameters.ProbeIndexes,Parameters.ProbeIDs)
    else
      begin
      FModelMan.LoadLayer(Parameters.ProbePDBFile);
      FDescriptors.Chains[1]:=GetChains(1,Parameters.ProbeIndexes,Parameters.ProbeIDs);
      end;
    end;
end;

function TDescriptorGen.GetSurfaceData(StructureIx: Integer): TSurfaceReport;
begin
  Result:=FSurfaceData[StructureIx];
end;

function TDescriptorGen.FullInnerContacts(StructureIx: Integer): TMatrix;
begin
  Result:=FSurfaceData[StructureIx].FullInnerContacts;
end;

function TDescriptorGen.SidechainInnerContacts(StructureIx: Integer): TMatrix;
begin
  Result:=FSurfaceData[StructureIx].SidechainInnerContacts;

end;

function TDescriptorGen.FullCrossContacts: TMatrix;
begin
  Result:=FCrossContacts.FullCrossContacts;
end;

function TDescriptorGen.SidechainCrossContacts: TMatrix;
begin
  Result:=FCrossContacts.SidechainCrossContacts;
end;

procedure TDescriptorGen.SelectSurfaceResidues(Index:Integer);

var f,curr:Integer;

begin
  with FDescriptors do
    begin
    curr:=0;
    SetLength(Selected[Index],Length(Residues[Index]));
    SetLength(SelectedIXs[Index],Length(Residues[Index]));
    for f:=0 to High(Residues[Index]) do
      if FSurfaceData[Index].SASurface[f]>FParameters.SurfaceCutoff then
        begin
        Selected[Index,curr]:=Residues[Index,f];
        SelectedIXs[Index,curr]:=f;
        Inc(curr);
        end;
    SetLength(Selected[Index],curr);
    SetLength(SelectedIXs[Index],curr);
    end;
end;

procedure TDescriptorGen.BuildNeighbourIndexes(Index: Integer);

var
  f,n,curr,r1,r2:Integer;
  tmpneighs:TIntegers;

begin
  with FDescriptors do
  with FSurfaceData[Index] do
    begin
    SetLength(Neighbours[Index],Length(SelectedIXs[Index]));
    SetLength(tmpneighs,Length(SelectedIXs[Index]));
    for f:=0 to High(SelectedIxs[Index]) do
      begin
      curr:=0;
      r1:=SelectedIXs[Index,f];
      for n:=0 to High(SelectedIXs[Index]) do
        begin
        r2:=SelectedIXs[Index,n];
        if (r1<>r2) and (FullInnerContacts[r1,r2]>0) then
          begin
          tmpneighs[curr]:=n;
          Inc(curr);
          end;
        end;
      if curr=0 then
        Neighbours[Index,f]:=nil
      else Neighbours[Index,f]:=Copy(tmpneighs,0,curr);
      end;
    end;
end;



procedure TDescriptorGen.SetInnerContacts(Index: Integer);
begin
  if FParameters.ComputeContactsByDistance then
    FPDBSurf.CalcInnerDistanceContacts(FSurfaceData[Index],FParameters.ContactDistance)
  else FPDBSurf.CalcInnerContacts(FSurfaceData[Index]);
end;

procedure TDescriptorGen.SetCrossContacts;
begin
  if FParameters.ComputeContactsByDistance then
    FCrossContacts:=FPDBSurf.CrossDistanceContactReport(
        FSurfaceData[TargetIndex],FSurfaceData[ProbeIndex],FParameters.ContactDistance)
  else FCrossContacts:=FPDBSurf.CrossContactReport(FSurfaceData[TargetIndex],FSurfaceData[ProbeIndex]);

  FCrossContacts.TargetIx:=TargetIndex;
  FCrossContacts.ProbeIx:=ProbeIndex;
end;

procedure TDescriptorGen.BuildMSAStrings(DBIndex, Index: Integer; Files: TSimpleStrings);

var
  f,rix:Integer;
  msa:TMSA;
  tmp:TSimpleStrings;

begin
  tmp:=nil;
  for f:=0 to High(Files) do
    begin
    msa:=ReadMSA(Files[f]);
    msa.GapMarker:='-';
    msa:=TrimMSA(msa, msa.SequenceIDs[0]);//remove all columns with gaps on query
    msa.GapMarker:='X';
    msa:=TrimMSA(msa, msa.SequenceIDs[0]);//remove all X, which are gaps on pdb
    tmp:=Concatenate(tmp,TransposeMSA(msa));
    end;
  with FDescriptors do
    begin
    SetLength(ResidueMSAs[DBIndex,Index],Length(Selected[Index]));
    for f:=0 to High(Selected[Index]) do
      begin
      rix:=SelectedIxs[Index,f];
      ResidueMSAs[DBIndex,Index,f]:=tmp[rix];
      if AAOneLetterCode(Residues[Index,rix].Name)<>tmp[rix,1] then
        raise Exception.Create('Residue mismatch on '+Residues[Index,rix].Parent.Name+' '+IntToStr(Residues[Index,rix].Id)+
          ' '+Residues[Index,rix].Name+' : '+tmp[rix,1]+' | '+IntToStr(rix));
      end;
    end;
end;

procedure TDescriptorGen.SetContactDescriptor(Ix:Integer;const Desc: TContactDescriptor);

var
  f:Integer;

begin
  with FDescriptors.Descriptors[Ix] do
    begin
    ID:=Desc.ID;
    SetLength(Values,Length(Desc.Values));
    for f:=0 to High(Desc.Values) do
      Values[f]:=Copy(Desc.Values[f],0,Length(Desc.Values[f]));
    end;
end;

procedure TDescriptorGen.AddContactDescriptor(const Desc: TContactDescriptor);

begin
  with FDescriptors do
    begin
    SetLength(Descriptors,Length(Descriptors)+1);
    SetContactDescriptor(High(Descriptors),Desc);
    end;
end;

procedure TDescriptorGen.AddMinMaxSurface;

var
  mins,maxs,minf,maxf:TContactDescriptor;
  f,g,r1,r2:Integer;

begin
  with FDescriptors do
    begin
    SetLength(mins.Values,Length(Selected[0]),Length(Selected[1]));
    SetLength(maxs.Values,Length(Selected[0]),Length(Selected[1]));
    SetLength(minf.Values,Length(Selected[0]),Length(Selected[1]));
    SetLength(maxf.Values,Length(Selected[0]),Length(Selected[1]));
    mins.Id:=MinSC;
    maxs.Id:=MaxSC;
    minf.Id:=MinFull;
    maxf.Id:=MaxFull;

    for f:=0 to High(SelectedIxs[0]) do
      for g:=0 to High(SelectedIxs[1]) do
        begin
        r1:=SelectedIxs[0,f];
        r2:=SelectedIxs[1,g];
        minf.Values[f,g]:=Min(FSurfaceData[0].SASurface[r1],FSurfaceData[1].SASurface[r2]);
        maxf.Values[f,g]:=Max(FSurfaceData[0].SASurface[r1],FSurfaceData[1].SASurface[r2]);
        mins.Values[f,g]:=Min(FSurfaceData[0].SASidechainSurface[r1],FSurfaceData[1].SASidechainSurface[r2]);
        maxs.Values[f,g]:=Max(FSurfaceData[0].SASidechainSurface[r1],FSurfaceData[1].SASidechainSurface[r2]);
        end;
    AddContactDescriptor(mins);
    AddContactDescriptor(maxs);
    AddContactDescriptor(minf);
    AddContactDescriptor(maxf);
    end;
end;


procedure TDescriptorGen.AddCrossContacts;


var
  cs,full:TContactDescriptor;
  f,g,r1,r2:Integer;

begin
  with FDescriptors do
    begin
    SetLength(cs.Values,Length(Selected[0]),Length(Selected[1]));
    SetLength(full.Values,Length(Selected[0]),Length(Selected[1]));
    cs.Id:=CrossSCContact;
    full.Id:=CrossFullContact;
    for f:=0 to High(SelectedIxs[0]) do
      for g:=0 to High(SelectedIxs[1]) do
        begin
        r1:=SelectedIxs[0,f];
        r2:=SelectedIxs[1,g];
        cs.Values[f,g]:=FCrossContacts.SidechainCrossContacts[r1,r2];
        full.Values[f,g]:=FCrossContacts.FullCrossContacts[r1,r2];
        end;
    AddContactDescriptor(cs);
    AddContactDescriptor(full);
    end;
end;



procedure TDescriptorGen.AddSequenceConservation(DBIndex:Integer);


procedure ScoreSubstitutions(Index:Integer; out Scores:TFloats; out Average:TFloat);

var
  f,g:Integer;
  len:Integer;
  seq:string;

begin
  with FDescriptors do
    begin
    SetLength(Scores,Length(Selected[Index]));
    Average:=0;
    len:=Length(ResidueMSAs[DBIndex,Index,0]);
    for f:=0 to High(Selected[Index]) do
      begin
      seq:=ResidueMSAs[DBIndex,Index,f];
      Scores[f]:=0;
      for g:=2 to len do
        Scores[f]:=Scores[f]+GetSubScore(SubMatrix,seq[1],seq[g]);
      Scores[f]:=Scores[f]/len;
      Average:=Average+Scores[f];
      end;
    Average:=Average/Length(ResidueMSAs[DBIndex,Index]);
    end;
end;

var
  mins,maxs,minrs,maxrs:TContactDescriptor;
  scores:array [0..1] of TFloats;
  averages:array[0..1] of TFloat;
  f,g,r1,r2:Integer;

begin
  with FDescriptors do
    begin
    SetLength(mins.Values,Length(Selected[0]),Length(Selected[1]));
    SetLength(maxs.Values,Length(Selected[0]),Length(Selected[1]));
    SetLength(minrs.Values,Length(Selected[0]),Length(Selected[1]));
    SetLength(maxrs.Values,Length(Selected[0]),Length(Selected[1]));


    ScoreSubstitutions(0,scores[0],averages[0]);
    ScoreSubstitutions(1,scores[1],averages[1]);


    maxs.Id:=MaxSubFraction+'('+MSADBNames[DBIndex]+')';
    mins.Id:=MinSubFraction+'('+MSADBNames[DBIndex]+')';
    maxrs.Id:=MaxRelSubFraction+'('+MSADBNames[DBIndex]+')';
    minrs.Id:=MinRelSubFraction+'('+MSADBNames[DBIndex]+')';

    for f:=0 to High(SelectedIxs[0]) do
      for g:=0 to High(SelectedIxs[1]) do
        begin
        mins.Values[f,g]:=Min(scores[0,f],scores[1,g]);
        maxs.Values[f,g]:=Max(scores[0,f],scores[1,g]);
        minrs.Values[f,g]:=Min(scores[0,f]-averages[0],scores[1,g]-averages[1]);
        maxrs.Values[f,g]:=Max(scores[0,f]-averages[0],scores[1,g]-averages[1]);
        end;
    AddContactDescriptor(mins);
    AddContactDescriptor(maxs);
    AddContactDescriptor(minrs);
    AddContactDescriptor(maxrs);
    end;

end;

procedure TDescriptorGen.AddGapFraction(DBIndex:Integer);


procedure CountGaps(Index:Integer; out Counts:TFloats; out Average:TFloat);

var
  f:Integer;
  len:Integer;

begin
  with FDescriptors do
    begin
    SetLength(Counts,Length(Selected[Index]));
    Average:=0;
    len:=Length(ResidueMSAs[DBIndex,Index,0]);
    for f:=0 to High(Selected[Index]) do
      begin
      Counts[f]:=CountInString(ResidueMSAs[DBIndex,Index,f],'-')/len;
      Average:=Average+Counts[f];
      end;
    Average:=Average/Length(ResidueMSAs[DBIndex,Index])+1;
    end;
end;

var
  ming,maxg,minrg,maxrg:TContactDescriptor;
  gapfractions:array [0..1] of TFloats;
  averages:array[0..1] of TFloat;
  f,g,r1,r2:Integer;

begin
  with FDescriptors do
    begin
    SetLength(ming.Values,Length(Selected[0]),Length(Selected[1]));
    SetLength(maxg.Values,Length(Selected[0]),Length(Selected[1]));
    SetLength(minrg.Values,Length(Selected[0]),Length(Selected[1]));
    SetLength(maxrg.Values,Length(Selected[0]),Length(Selected[1]));

    CountGaps(0,gapfractions[0],averages[0]);
    CountGaps(1,gapfractions[1],averages[1]);

    maxg.Id:=MaxGapFraction+'('+MSADBNames[DBIndex]+')';
    ming.Id:=MinGapFraction+'('+MSADBNames[DBIndex]+')';
    maxrg.Id:=MaxRelGapFraction+'('+MSADBNames[DBIndex]+')';
    minrg.Id:=MinRelGapFraction+'('+MSADBNames[DBIndex]+')';

    for f:=0 to High(SelectedIxs[0]) do
      for g:=0 to High(SelectedIxs[1]) do
        begin
        ming.Values[f,g]:=Min(gapfractions[0,f],gapfractions[1,g]);
        maxg.Values[f,g]:=Max(gapfractions[0,f],gapfractions[1,g]);
        minrg.Values[f,g]:=Min(gapfractions[0,f]/averages[0],gapfractions[1,g]/averages[1]);
        maxrg.Values[f,g]:=Max(gapfractions[0,f]/averages[0],gapfractions[1,g]/averages[1]);
        end;
    AddContactDescriptor(ming);
    AddContactDescriptor(maxg);
    AddContactDescriptor(minrg);
    AddContactDescriptor(maxrg);
    end;

end;

procedure TDescriptorGen.AddContactScores(DBIndex:Integer;CMLabel:string);

procedure ScoreContact(TRix,PRix:Integer;out Score,Expected:TFloat);

var
  f,g:Integer;
  len:Integer;
  tseq,pseq:TOCCompressedSequence;
  count,tpart,part:Integer;

begin
  with FDescriptors do
    begin
    tseq:=CompressedMSAs[DBIndex,0,TRix];
    pseq:=CompressedMSAs[DBIndex,1,PRix];
    Score:=ScoreSequencePair(tseq,pseq,ContactMatrix);
    Score:=Score/tseq[High(tseq)].Last;
    Expected:=0;
    count:=0;
    for f:=0 to High(tseq) do
      begin
      tpart:=tseq[f].Last-tseq[f].First+1;
      for g:=0 to High(pseq) do
        begin
        part:=tpart*(pseq[g].Last-pseq[g].First+1);
        Expected:=Expected+
          GetSubScore(ContactMatrix,tseq[f].Symbol,pseq[g].Symbol)*part;
        count:=count+part;
        end;
      end;
    Expected:=Expected/count;
    end;
end;


var
  cs,rcs:TContactDescriptor;
  score,expected:TFloat;
  f,g,r1,r2:Integer;

begin
  with FDescriptors do
    begin
    SetLength(cs.Values,Length(Selected[0]),Length(Selected[1]));
    SetLength(rcs.Values,Length(Selected[0]),Length(Selected[1]));
    cs.Id:=ContactFraction+'('+MSADBNames[DBIndex]+')'+CMLabel;
    rcs.Id:=RelContactFraction+'('+MSADBNames[DBIndex]+')'+CMLabel;
    WriteLn(Length(SelectedIxs[0]));
    for f:=0 to High(SelectedIxs[0]) do
      begin
      for g:=0 to High(SelectedIxs[1]) do
        begin
        ScoreContact(f,g,score,expected);
        cs.Values[f,g]:=score;
        rcs.Values[f,g]:=score-expected;
        end;
      if f mod 100 = 0 then
        Writeln(f);
      end;
    AddContactDescriptor(cs);
    AddContactDescriptor(rcs);
    end;
end;

procedure TDescriptorGen.AddScotchScore(DBIndex: Integer);

const
  Neutral='GAVLIMCFPWY';
  Polar='STNQ';
  Positive='KRH';
  Negative='DE';

function ResType(Res:string):Integer;

begin
  Result:=-1;
  if Pos(Res,Neutral)>0 then
    Result:=0
  else if Pos(Res,Polar)>0 then
    Result:=1
  else if Pos(Res,Positive)>0 then
    Result:=2
  else if Pos(Res,Negative)>0 then
    Result:=3;
end;

function ScoreColumns(Col1,Col2:string):TFloat;

var
  f,t1,t2:Integer;

begin
  Result:=0;
  for f:=1 to Length(Col1) do
    begin
    t1:=ResType(Col1[f]);
    t2:=ResType(Col2[f]);
    if ((t1=0) and (t2=0)) or
       ((t1=1) and (t2=1)) or
       ((t1=2) and (t2=3)) or
       ((t1=3) and (t2=2)) then
         Result:=Result+1
    end;
  Result:=Result/Length(Col1);
end;

var
  cs:TContactDescriptor;
  score:TFloat;
  f,g,r1,r2:Integer;

begin
  with FDescriptors do
    begin
    SetLength(cs.Values,Length(Selected[0]),Length(Selected[1]));
    cs.Id:=ScotchScore+'('+MSADBNames[DBIndex]+')';
    for f:=0 to High(SelectedIxs[0]) do
      for g:=0 to High(SelectedIxs[1]) do
        begin
        score:=ScoreColumns(ResidueMSAs[DBIndex,0,f],ResidueMSAs[DBIndex,1,g]);
        cs.Values[f,g]:=score;
        end;
    AddContactDescriptor(cs);
    end;
end;

procedure TDescriptorGen.AddPatchScores;

var
  allone,allall:TContactDescriptor;
  descount:Integer;
  total:TFloat;
  f,g,r1,r2,n,n1,n2,d:Integer;

begin
  with FDescriptors do
    begin
    SetLength(allone.Values,Length(Selected[0]),Length(Selected[1]));
    SetLength(allall.Values,Length(Selected[0]),Length(Selected[1]));
    descount:=High(Descriptors);
    SetLength(Descriptors,Length(Descriptors)*3);

    for d:=0 to descount do
      begin
      allone.ID:=Descriptors[d].ID+' (all to one)';
      allall.ID:=Descriptors[d].ID+' (all to all)';
      for r1:=0 to High(Selected[0]) do
        for r2:=0 to High(Selected[1]) do
          begin
          total:=Descriptors[d].Values[r1,r2];
          for n:=0 to High(Neighbours[0,r1]) do
            total:=total+Descriptors[d].Values[Neighbours[0,r1,n],r2];
          for n:=0 to High(Neighbours[1,r2]) do
            total:=total+Descriptors[d].Values[r1,Neighbours[1,r2,n]];
          allone.Values[r1,r2]:=total;
          for n1:=0 to High(Neighbours[0,r1]) do
            for n2:=0 to High(Neighbours[1,r2]) do
              total:=total+Descriptors[d].Values[Neighbours[0,r1,n1],Neighbours[1,r2,n2]];
          allall.Values[r1,r2]:=total;
          end;
      SetContactDescriptor(d+descount+1,allone);
      SetContactDescriptor(d+2*descount+2,allall);
      end;
    end;
end;

procedure TDescriptorGen.ClearAll;
begin
  FDescriptors:=EmptyDescriptorResults;
  SetLength(FSurfaceData,2);
  FSurfaceData[TargetIndex]:=EmptySurfaceReport;
  FSurfaceData[ProbeIndex]:=EmptySurfaceReport;
  FCrossContacts:=EmptyCrossContactReport;
end;


constructor TDescriptorGen.Create(NumSurfacePoints:Integer);
begin
  inherited Create;
  FPDBSurf:=TPDBSurface.Create(NumSurfacePoints);
  { TODO : Set parameters }
  FPDBSurf.ProbeRadius:=1.4;//H atoms are missing but using united atom radii
  FModelMan:=TPDBModelMan.Create(Config.MonomersPath);
end;

destructor TDescriptorGen.Destroy;
begin
  FPDBSurf.Free;
  FModelMan.Free;
  inherited Destroy;
end;

procedure TDescriptorGen.GenerateAllDescriptors(const Parameters: TCDParameters; Verbose:Boolean=True);

var f,g:Integer;

procedure Report(s:string);
begin
  if Verbose then
    WriteLn(s);
end;

begin
  ClearAll;
  Report('Loading');
  LoadModels(Parameters);
  FDescriptors.SubMatrix:=ReadBLASTMatrix(Parameters.SubMatrixFile);

  Report('Computing surfaces');
  SetChainSurface(TargetIndex);
  SetChainSurface(ProbeIndex);

  Report('Computing contacts');
  SetInnerContacts(TargetIndex);
  SetInnerContacts(ProbeIndex);

  //selection and descriptors use FDescriptor data
  //Copying to FDescriptors makes it easier to report results (?)
  Report('Selecting residues');
  FDescriptors.Residues[TargetIndex]:=FSurfaceData[TargetIndex].Residues;
  FDescriptors.Residues[ProbeIndex]:=FSurfaceData[ProbeIndex].Residues;
  SelectSurfaceResidues(TargetIndex);
  SelectSurfaceResidues(ProbeIndex);
  BuildNeighbourIndexes(TargetIndex);
  BuildNeighbourIndexes(ProbeIndex);


  Report('Building sequence data');
  FDescriptors.MSADBNames:=Parameters.MSADatabaseNames;
  SetLength(FDescriptors.ResidueMSAs,Length(Parameters.MSADatabaseNames));

  for f:=0 to High(Parameters.MSADatabaseNames) do
    begin
    Report(Parameters.MSADatabaseNames[f]);
    BuildMSAStrings(f,TargetIndex,Parameters.ResidueMSAFiles[f,TargetIndex]);
    BuildMSAStrings(f,ProbeIndex,Parameters.ResidueMSAFiles[f,ProbeIndex]);
    end;
  CompressMSAs;

  Report('Surface Scores');
  AddMinMaxSurface;

  Report('Sequence scores');
  for f:=0 to High(Parameters.MSADatabaseNames) do
    begin
    Report('  '+Parameters.MSADatabaseNames[f]);
    AddSequenceConservation(f);
    writeln('add');
    for g:=0 to High(Parameters.ContactMatrixFiles) do
      begin
      writeln(g);
      Report('    '+Parameters.ContactMatrixFiles[g]);
      FDescriptors.ContactMatrix:=ReadBLASTMatrix(Parameters.ContactMatrixFiles[g]);
      AddContactScores(f,ExtractFileName(Parameters.ContactMatrixFiles[g]));
      end;
    AddGapFraction(f);
    AddScotchScore(f);
    end;

  Report('Contacts');
  if Parameters.TargetPDBFile=Parameters.ProbePdbFile then
    begin
    SetCrossContacts;
    AddCrossContacts;
    end;

  Report('Patches');
  //after all scores to compute as patches (all to one, all to all)
  AddPatchScores;

end;

procedure TDescriptorGen.DistanceDescriptors(const Parameters: TCDParameters);

begin
  ClearAll;
  LoadModels(Parameters);

  SetChainSurface(TargetIndex);
  SetChainSurface(ProbeIndex);
  SetInnerContacts(TargetIndex);
  SetInnerContacts(ProbeIndex);

  FDescriptors.Residues[TargetIndex]:=FSurfaceData[TargetIndex].Residues;
  FDescriptors.Residues[ProbeIndex]:=FSurfaceData[ProbeIndex].Residues;
  SelectSurfaceResidues(TargetIndex);
  SelectSurfaceResidues(ProbeIndex);

  BuildNeighbourIndexes(TargetIndex);
  BuildNeighbourIndexes(ProbeIndex);

  AddMinMaxSurface;

  if Parameters.TargetPDBFile=Parameters.ProbePdbFile then
    begin
    SetCrossContacts;
    AddCrossContacts;
    end;
  AddPatchScores;

end;


procedure TDescriptorGen.ReportDescriptors(Report: TStrings; Headers: Boolean; OnlyContacts:Boolean=False);

var
  f,r1,r2,d,contactindex:Integer;
  s:string;

begin
  with FDescriptors do
  if Descriptors<>nil then
    begin
      if OnlyContacts then
        begin
        contactindex:=-1;
        for f:=0 to High(Descriptors) do
          if Descriptors[f].ID=CrossSCContact then
            begin
            contactindex:=f;
            Break;
            end;
        end;
      if Headers then
        begin
        s:='Target ID-Chain-Res Name-Res ID:'+
        'Probe ID-Chain-Res Name-Res ID';
        for f:=0 to High(Descriptors) do s:=s+#9+Descriptors[f].ID;
        Report.Add(s);
        end;
      for r1:=0 to High(Descriptors[0].Values) do
        for r2:=0 to High(Descriptors[0].Values[0]) do
          begin
          with Selected[0,r1] do
            s:=FParameters.TargetPDBID+'-'+Parent.Name+'-'+
               Name+'-'+IntToStr(ID);
          with Selected[1,r2] do
            s:=s+':'+FParameters.ProbePDBID+'-'+Parent.Name+'-'+
               Name+'-'+IntToStr(ID);
          for d:=0 to High(Descriptors) do
            s:=s+#9+FloatToStrF(Descriptors[d].Values[r1,r2],ffFixed,6,4);
          if OnlyContacts then
            begin
            if Descriptors[contactindex].Values[r1,r2]>0 then
               Report.Add(s);
            end
          else Report.Add(s);
          end;
    end;
end;

procedure TDescriptorGen.SaveDescriptors(Headers:Boolean=false;OnlyContacts:Boolean=False);

var sl:TStringList;

begin
  if FParameters.DescriptorFile<>'' then
    begin
    sl:=TStringList.Create;
    ReportDescriptors(sl,Headers,OnlyContacts);
    sl.SaveToFile(FParameters.DescriptorFile);
    sl.Free;
    end;
end;

procedure TDescriptorGen.SaveRaw;

var
  sl:TStringList;
  f,msadb,r,g,contactindex:Integer;
  s,pdbid:string;


begin
  if FParameters.RawFile<>'' then
    begin
    sl:=TStringList.Create;
    with FDescriptors do
      if Descriptors<>nil then
        begin
        for f:=0 to 1 do
          begin
          if f=0 then
            pdbid:=FParameters.TargetPDBID
          else
            pdbid:=FParameters.ProbePDBID;
          for msadb:=0 to High(FParameters.MSADatabaseNames) do
            begin
            sl.Add(FParameters.MSADatabaseNames[msadb]);
            for r:=0 to High(Selected[f]) do
              begin
              with Selected[f,r] do
                s:=pdbid+'-'+Parent.Name+'-'+Name+'-'+IntToStr(ID);
              sl.Add(s+#9+ResidueMSAs[msadb,f,r]);
              end;
            end;
          for r:=0 to High(Neighbours[f]) do
            begin
            with Selected[f,r] do
              s:=pdbid+'-'+Parent.Name+'-'+Name+'-'+IntToStr(ID)+':';
            for g:=0 to High(Neighbours[f,r]) do
              with Selected[f,Neighbours[f,r,g]] do
                s:=s+' '+Parent.Name+'-'+Name+'-'+IntToStr(ID);
            sl.Add(s);
            end;
          end;
        end;
    sl.SaveToFile(FParameters.RawFile);
    sl.Free;
    end;
end;

end.



