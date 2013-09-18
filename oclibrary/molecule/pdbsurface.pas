{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain. Enjoy.
********************************************************************************
Author: Ludwig Krippahl
Date: 16.8.2013
Purpose:
  Interfaces pdb manager with surface calculations
  Allows selection of parts of pdb structures, neighbour calculations,
  inter-chain contacts and so forth
Requirements:
Revisions:
To do:
*******************************************************************************}

unit pdbsurface;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, surface, pdbmolecules,
  molecules, geomhash, quicksort, pdbparser, mmcifparser, oclconfiguration,
  stringutils;

type

  //Remove?

  TCuboid=array[0..1] of TCoord;
  TAtomRec=record
    Atom:TAtom;
    ASA:TFloat;
    ResidueRec:Integer;
  end;
  TAtomRecs=array of TAtomRec;

  TResidueRec=record
    Residue:TMolecule;
    ContactHull:TCuboid;
    ASA,IsolatedASA:TFloat;
    ContactSurfaces:TFloats;
    ContactResidues:TIntegers;
  end;
  TResidueRecs=array of TResidueRec;

  TGroup=record
    //group for surface characterization and contacts
    ASA:TFLoat;
    Chains:TMolecules;
    AtomRecs:TAtomRecs;
    ResidueRecs:TResidueRecs;
    Name:string;
  end;
  TGroups=array of TGroup;

  TSurfaceReport=record
    Residues:TMolecules;
    IsolatedSurface:TFloats;
    IsolatedSidechainSurface:TFLoats;
    SASurface:TFloats;
    SASidechainSurface:TFloats;
  end;



  { TPdbSurface }

  TPdbSurface = class
  protected
    FGoldenSphere:TCoords;
    FProbeRadius:TFloat;
    FMinHashCell:TFloat;     //cell width for geometric hashing
    FUseOnlyRes,FExcludeRes:TSimpleStrings;
                          //list of residue types to use or to exclude
                          //if FUseOnlyRes is nil, use all but those in FExcludeRes
                          //otherwise, use only those in FUseOnlyRes

  public
    property MinHashCell:TFloat read FMinHashCell write FMinHashCell;
    property ProbeRadius:TFloat read FProbeRadius write FProbeRadius;
    property UseOnlyRes:TSimpleStrings read FUseOnlyRes write FUseOnlyRes;
    property ExcludeRes:TSimpleStrings read FExcludeRes write FExcludeRes;

    constructor Create(NumPoints:Integer);
                      //number of points to use for generating each sphere
    destructor Destroy; override;
    function CalcASAs(Atoms:TAtoms):TFloats;
    function UseResidue(Residue:TMolecule):Boolean;
    function BuildIndexes(Chains:TMolecules):TSurfaceReport;
    procedure IsolatedResidueSurface(var Report:TSurfaceReport);
        //fills report with isolated surfaces for all residues in the report
    procedure ResiduesSAS(var Report:TSurfaceReport);
        ////fills report with ASA for full and sidechain residues in the report
    function SurfaceStatistics(Chains:TMolecules):TSurfaceReport;
    procedure InnerContactMatrixes(const Surf:TSurfaceReport; out Full,Sidechain:TMatrix);
    procedure OuterContactMatrixes(const Target,Probe:TSurfaceReport; out Full,Sidechain:TMatrix);

  end;

  function EmptySurfaceReport:TSurfaceReport;

implementation

function EmptySurfaceReport:TSurfaceReport;
begin
  with Result do
    begin
    Residues:=nil;
    IsolatedSurface:=nil;
    IsolatedSidechainSurface:=nil;
    SASurface:=nil;
    SASidechainSurface:=nil;
    end;
end;

{ TPdbSurface }

constructor TPdbSurface.Create(NumPoints:Integer);
begin
  inherited Create;
  FGoldenSphere:=GoldenSpiralPoints(NumPoints);
  FUseOnlyRes:=nil;
  FExcludeRes:=nil;
  FMinHashCell:=2;
  FProbeRadius:=1.4;
end;

destructor TPdbSurface.Destroy;
begin
  FGoldenSphere:=nil;
  FUseOnlyRes:=nil;
  FExcludeRes:=nil;
  inherited Destroy;
end;

function TPdbSurface.CalcASAs(Atoms: TAtoms): TFloats;

var
  f:Integer;
  points:TCoords;
  radii:TFloats;

begin
  SetLength(points,Length(Atoms));
  SetLength(radii,Length(Atoms));
  for f:=0 to High(Atoms) do
    begin
    points[f]:=Atoms[f].Coords;
    radii[f]:=Atoms[f].Radius+FProbeRadius;
    end;
    Result:=SRSurface(points,radii,FGoldenSphere,FMinHashCell);
end;


function TPdbSurface.UseResidue(Residue: TMolecule): Boolean;
begin
  if (FExcludeRes=nil) and (FUseOnlyRes=nil) then
    Result:=True
  else if (FUseOnlyRes=nil) and not IsInArray(Residue.Name,FExcludeRes) then
    Result:=True
  else if IsInArray(Residue.Name,FUseOnlyRes) then
    Result:=True
  else Result:=False;
end;

function TPdbSurface.BuildIndexes(Chains: TMolecules): TSurfaceReport;

var
  f,r,totres:Integer;

begin
  totres:=0;
  for f:=0 to High(Chains) do
    totres:=totres+Chains[f].GroupCount;
  SetLength(Result.Residues,totres);
  totres:=-1;
  for f:=0 to High(Chains) do
    for r:=0 to Chains[f].GroupCount-1 do
      if UseResidue(Chains[f].GetGroup(r)) then
        begin
        Inc(totres);
        Result.Residues[totres]:=Chains[f].GetGroup(r);
        end;
  SetLength(Result.Residues,totres);
end;

procedure TPdbSurface.IsolatedResidueSurface(var Report:TSurfaceReport);

var
  f,g:Integer;
  atoms:TAtoms;
  surfs:TFloats;
  tot:TFloat;

begin
  SetLength(Report.IsolatedSurface,Length(Report.Residues));
  SetLength(Report.IsolatedSidechainSurface,Length(Report.Residues));
  for f:=0 to High(Report.Residues) do
    begin
    atoms:=Report.Residues[f].AllAtoms;
    surfs:=CalcASAs(atoms);
    Report.IsolatedSurface[f]:=Sum(surfs);
    tot:=0;
    for g:=0 to High(atoms) do
        if not AtomIsAABackBone(atoms[g]) then
          tot:=tot+surfs[g];
    Report.IsolatedSidechainSurface[f]:=tot;
    end;
end;

procedure TPdbSurface.ResiduesSAS(var Report:TSurfaceReport);

var
  f,a,ar,totatoms:Integer;
  cchain,cres:TMolecule;
  allatoms,resatoms:TAtoms;
  starts,ends:TIntegers; //where each residue starts and ends
  surfvals:TFloats;
  tot,totsc:TFloat;

begin
  with Report do
    begin
    SetLength(SASurface,Length(Residues));
    SetLength(starts,Length(Residues));
    SetLength(ends,Length(Residues));
    SetLength(SASidechainSurface,Length(Residues));
    totatoms:=0;
    for f:=0 to High(Residues) do
    totatoms:=totatoms+Residues[f].AtomCount;
    SetLength(allatoms,totatoms);
    a:=-1;
    for f:=0 to High(Residues) do
      begin
      resatoms:=Residues[f].AllAtoms;
      starts[f]:=a+1;
      for ar:=0 to High(resatoms) do
        begin
        Inc(a);
        allatoms[a]:=resatoms[ar];
        end;
      ends[f]:=a;
      end;
    surfvals:=CalcASAS(allatoms);
    for f:=0 to High(Residues) do
      begin
      tot:=0;
      totsc:=0;
      for a:=starts[f] to ends[f] do
        begin
        tot:=tot+surfvals[a];
        if not AtomIsAABackbone(allatoms[a]) then
          totsc:=totsc+surfvals[a];
        end;
      SASurface[f]:=tot;
      SASidechainSurface[f]:=totsc;
      end;
    end;
end;

function TPdbSurface.SurfaceStatistics(Chains: TMolecules): TSurfaceReport;
begin
  Result:=BuildIndexes(Chains);
  IsolatedResidueSurface(Result);
  ResiduesSAS(Result);
end;

procedure TPdbSurface.InnerContactMatrixes(const Surf: TSurfaceReport; out
  Full, Sidechain: TMatrix);
begin

end;

procedure TPdbSurface.OuterContactMatrixes(const Target, Probe: TSurfaceReport;
  out Full, Sidechain: TMatrix);
begin

end;

end.

