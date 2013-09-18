{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain. Enjoy.
********************************************************************************
Author: Ludwig Krippahl
Date: 5.9.2013
Purpose:
  Prediction of residue contacts in protein complexes.
  Currently generation of potential descriptors.
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
  alignment, stringutils, progress, pdbsurface;

const
  NumSurfacePoints=50;
    //number of points used per atom to compute surface
    //this is a constant because once the classifier is trained
    //this parameter should not change

type

  TResiduePairs=array[0..1] of TIntegers;

  TCDParameters=record
    TargetPDB,ProbePDB:string;

    TargetIDs,ProbeIDs:TSimpleStrings;
    TargetIndexes,ProbeIndexes:TIntegers;
      //note: indexes take precedence; ids are ignored if indexes are not nil

    TargetMSAFiles,ProbeMSAFiles:TSimpleStrings;
    TargetOrganismFiles,ProbeOrganismFiles:TSimpleStrings;
      //Must have same number of elements as IDs or Indexes

  { TODO : Add lists for excluded and useonly amino acids }
  end;

  TContactDescriptor=record
    Values:TMatrix;
    Id:string;
  end;
  TContactDescriptors=array of TContactDescriptor;

  TDescriptorResults=record
    TargetResidues,ProbeResidues:TMolecules;
    Descriptors:TContactDescriptors;
  end;

  { TDescriptorGen }

  TDescriptorGen=class
  protected
    FParameters:TCDParameters;
    FModelMan:TPDBModelMan;
    FPDBSurf:TPDBSurface;
    FChains:array[0..1] of TMolecules;
      //target and probe chain sets

    //All surface and contact values are indexed to the full set of used residues
    //(the residue lists in FSurfaceData)
    FSurfaceData:array[0..1] of TSurfaceReport;
      //target and probe, lists only used residues in surface computations
    FFullInnerContacts:array[0..1] of TMatrix;
    FSidechainInnerContacts:array[0..1] of TMatrix;

    FFullCrossContacts:TMatrix;
    FSideChainCrossContacts:TMatrix;
      //These are nil if no cross contacts requested

    //The remaining data indexes the SelectedResidues
    FSelectedResidues:array[0..1] of TMolecules;
      //Exposed residues at surface
    FNeighbourLists:array[0..1] of array of TIntegers;
      //TIntegers corresponds to one selected residue.


    // output values
    FDescriptors:TDescriptorResults;

    procedure SetChainSurface(Index:Integer);
    procedure FindSurfaceResidues;
    function GetChains(Layer:Integer;Indexes:TIntegers;IDs:TSimpleStrings):TMolecules;

  public
    property Descriptors:TDescriptorResults read FDescriptors;
    procedure ClearAll;
    procedure LoadModels(Parameters:TCDParameters);
      //clears everything, loads models, sets parameters;

    function GetSurfaceData(StructureIx:Integer):TSurfaceReport;
      //admitted residues are all used for surface calculations

    constructor Create;
    destructor Destroy;override;
    procedure GenerateAllDescriptors( const Parameters:TCDParameters);
  end;

function ReadParameters(FileName:string):TCDParameters;

implementation

const
  TargetIndex=0;
  ProbeIndex=1;


function ReadParameters(FileName:string):TCDParameters;

var
  f:Integer;
  ini:TINIFile;
  pdbfolder,sequencefolder:string;

  procedure SelectedChains(
      IDString,IndexString:string;
      out ChainIds:TSimpleStrings;
      out ChainIndexes:TIntegers);

  var
    ixs:TSimpleStrings;
    f:Integer;

  begin
    ChainIds:=SplitString(IDString,',');
    ixs:=SPlitString(IndexString,',');
    SetLength(ChainIndexes,Length(ixs));
    for f:=0 to High(ixs) do
      ChainIndexes[f]:=StrToInt(ixs[f])-1;
  end;

begin
  ini:=TINIFile.Create(FileName);

  pdbfolder:=ini.ReadString('Structure','PDBFolder','');
  if pdbfolder<>'' then pdbfolder:=IncludeTrailingPathDelimiter(pdbfolder);

  sequencefolder:=ini.ReadString('Structure','SequenceFolder','');
  if sequencefolder<>'' then pdbfolder:=IncludeTrailingPathDelimiter(sequencefolder);

  with Result do
    begin
    TargetPDB:=pdbfolder+ini.ReadString('Structure','Target','');
    ProbePDB:=pdbfolder+ini.ReadString('Structure','Target','');
    SelectedChains(ini.ReadString('Structure','TargetChainIDs',''),
                   ini.ReadString('Structure','TargetChainIndexes',''),
                   Result.TargetIDs,Result.TargetIndexes);
   SelectedChains(ini.ReadString('Structure','ProbeChainIDs',''),
                   ini.ReadString('Structure','ProbeChainIndexes',''),
                   Result.ProbeIDs,Result.ProbeIndexes);
    end;
  ini.Free;
end;

{ TDescriptorGen }

procedure TDescriptorGen.SetChainSurface(Index: Integer);

begin
  FSurfaceData[Index]:=FPDBSurf.SurfaceStatistics(FChains[Index]);
end;

procedure TDescriptorGen.LoadModels(Parameters:TCDParameters);


var f:Integer;

begin
  ClearAll;
  FParameters:=Parameters;
  FModelMan.LoadLayer(Parameters.TargetPDB);
  FChains[0]:=GetChains(0,Parameters.TargetIndexes,Parameters.TargetIDs);
  if Parameters.TargetPDB=Parameters.ProbePDB then
    FChains[1]:=GetChains(0,Parameters.ProbeIndexes,Parameters.ProbeIDs)
  else
    begin
    FModelMan.LoadLayer(Parameters.ProbePDB);
    FChains[1]:=GetChains(1,Parameters.ProbeIndexes,Parameters.ProbeIDs);
    end;
end;

function TDescriptorGen.GetSurfaceData(StructureIx: Integer): TSurfaceReport;
begin
  Result:=FSurfaceData[StructureIx];
end;

procedure TDescriptorGen.FindSurfaceResidues;
begin

end;

function TDescriptorGen.GetChains(Layer: Integer; Indexes: TIntegers;
  IDs: TSimpleStrings): TMolecules;

var f:Integer;

begin
  if Indexes<>nil then
    begin
    SetLength(Result,Length(Indexes));
    for f:=0 to High(Result) do
      Result[f]:=FModelMan.LayerByIx(0).GetChain(Indexes[f]);
    end
  else
    begin
    SetLength(Result,Length(IDs));
    for f:=0 to High(Result) do
      Result[f]:=FModelMan.LayerByIx(0).GetChain(IDs[f]);
    end;
end;

procedure TDescriptorGen.ClearAll;
begin
  FChains[TargetIndex]:=nil;
  FChains[ProbeIndex]:=nil;
  FSurfaceData[TargetIndex]:=EmptySurfaceReport;
  FSurfaceData[ProbeIndex]:=EmptySurfaceReport;
  FFullInnerContacts[TargetIndex]:=nil;
  FFullInnerContacts[ProbeIndex]:=nil;
  FSidechainInnerContacts[TargetIndex]:=nil;
  FSidechainInnerContacts[ProbeIndex]:=nil;
  FFullCrossContacts:=nil;
  FSideChainCrossContacts:=nil;
  FSelectedResidues[TargetIndex]:=nil;
  FSelectedResidues[ProbeIndex]:=nil;
  FNeighbourLists[TargetIndex]:=nil;
  FNeighbourLists[ProbeIndex]:=nil;
end;


constructor TDescriptorGen.Create;
begin
  inherited;
  FPDBSurf:=TPDBSurface.Create(NumSurfacePoints);
  { TODO : Set parameters }
  FModelMan:=TPDBModelMan.Create(Config.MonomersPath);
end;

destructor TDescriptorGen.Destroy;
begin
  FPDBSurf.Free;
  FModelMan.Free;
  inherited Destroy;
end;

procedure TDescriptorGen.GenerateAllDescriptors(const Parameters: TCDParameters);
begin
  ClearAll;
  LoadModels(Parameters);
  SetChainSurface(TargetIndex);
  SetChainSurface(ProbeIndex);
end;

end.

