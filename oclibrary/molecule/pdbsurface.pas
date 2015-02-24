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
  stringutils, molutils;

type

  TSurfaceReport=record
    PdbFile:string;
    Chains:TMolecules;
    Residues:TMolecules;
    IsolatedSurface:TFloats;
    IsolatedSidechainSurface:TFLoats;
    SASurface:TFloats;
    SASidechainSurface:TFloats;
    FullInnerContacts:TMatrix;
    SidechainInnerContacts:TMatrix;
  end;

  TSurfaceReports=array of TSurfaceReport;

  TCrossContactReport=record
    TargetIx,ProbeIx:Integer;
    FullCrossContacts:TMatrix;
    SidechainCrossContacts:TMatrix;
  end;

  TCrossContactReports=array of TCrossContactReport;

  { TPdbSurface }

  TPdbSurface = class
  protected
    FGoldenSphere:TCoords;
    FProbeRadius:TFloat;
    FMinHashCell:TFloat;  //cell width for geometric hashing
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
        //fills report with ASA for full and sidechain residues in the report
    function SurfaceStatistics(Chains:TMolecules):TSurfaceReport;

    procedure CalcInnerContacts(var Surf:TSurfaceReport);
        //the Sidechain matrix contains the surface contacts between each sidechain and any other atom
    procedure CalcInnerDistanceContacts(var Surf:TSurfaceReport;const Dist:TFloat);
        //same as CalcInnerContacts, but sets surface contact to 1 if any within Dist, 0 otherwise

    function CrossContactReport(Target,Probe:TSurfaceReport):TCrossContactReport;
    function CrossContactReports(const Surfaces:TSurfaceReports):TCrossContactReports;
        //the Sidechain matrix contains the surface contacts between each sidechain and any other atom

    function CrossDistanceContactReport(const Target,Probe:TSurfaceReport; const Dist:TFloat):TCrossContactReport;
    function CrossDistanceContactReports(const Surfaces:TSurfaceReports;const Dist:TFloat):TCrossContactReports;
        //same as cross contact report but using distance instead of surface

    //function SurfaceReport(PdbFile: string; Chains: TSimpleStrings): TSurfaceReport;
  end;

  function EmptySurfaceReport:TSurfaceReport;
  function EmptyCrossContactReport:TCrossContactReport;


implementation

function EmptySurfaceReport:TSurfaceReport;
begin
  with Result do
    begin
    PdbFile:='';
    Chains:=nil;
    Residues:=nil;
    IsolatedSurface:=nil;
    IsolatedSidechainSurface:=nil;
    SASurface:=nil;
    SASidechainSurface:=nil;
    FullInnerContacts:=nil;
    SidechainInnerContacts:=nil;
    end;
end;

function EmptyCrossContactReport: TCrossContactReport;
begin
 with Result do
   begin
   FullCrossContacts:=nil;
   SidechainCrossContacts:=nil;
   TargetIx:=-1;
   ProbeIx:=-1;
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
  Result:=ResidueIsAminoAcid(Residue);
  if Result then
    begin
    if (FExcludeRes=nil) and (FUseOnlyRes=nil) then
      Result:=True
    else if (FUseOnlyRes=nil) and not IsInArray(Residue.Name,FExcludeRes) then
      Result:=True
    else if IsInArray(Residue.Name,FUseOnlyRes) then
      Result:=True
    else Result:=False;
    end;
end;

function TPdbSurface.BuildIndexes(Chains: TMolecules): TSurfaceReport;

var
  f,r,totres:Integer;

begin
  totres:=0;
  Result.Chains:=Chains;
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
  SetLength(Result.Residues,totres+1);
end;

procedure TPdbSurface.IsolatedResidueSurface(var Report:TSurfaceReport);

var
  f,g:Integer;
  atoms,scatoms:TAtoms;
  surfs,scsurfs:TFloats;
  tot:TFloat;

begin
  SetLength(Report.IsolatedSurface,Length(Report.Residues));
  SetLength(Report.IsolatedSidechainSurface,Length(Report.Residues));
  for f:=0 to High(Report.Residues) do
    begin
    atoms:=Report.Residues[f].AllAtoms;
    scatoms:=NoBackbone(atoms);
    surfs:=CalcASAs(atoms);
    scsurfs:=CalcASAs(scatoms);
    Report.IsolatedSurface[f]:=Sum(surfs);
    Report.IsolatedSidechainSurface[f]:=Sum(scsurfs);
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



procedure TPdbSurface.CalcInnerContacts(var Surf:TSurfaceReport);

var
  f,r1,r2:Integer;
  hulls:array of TCuboid;
  atoms,scatoms:TAtoms;
  fullsurf,scsurf:TFloat;
  scasas,asas:TFloats;

begin
  with Surf do
    begin
    SetLength(hulls,Length(Residues));
    SetLength(FullInnerContacts,Length(Residues),Length(Residues));
    SetLength(SidechainInnerContacts,Length(Residues),Length(Residues));
    for r1:=0 to High(FullInnerContacts) do
      for r2:=0 to High(FullInnerContacts) do
        begin
        FullInnerContacts[r1,r2]:=0;
        SidechainInnerContacts[r1,r2]:=0;
        end;
    for f:=0 to High(Residues) do
      hulls[f]:=CalcHull(Residues[f].AllAtoms,2*FProbeRadius);
    for r1:=0 to High(Residues)-1 do
      for r2:=r1+1 to High(Residues) do
        if InContact(hulls[r1],hulls[r2]) then
          begin
          atoms:=AppendAtomsToArray(Residues[r1].AllAtoms,Residues[r2].AllAtoms);
          scatoms:=NoBackbone(atoms);
          asas:=CalcASAs(atoms);
          scasas:=CalcASAs(scatoms);
          fullsurf:=0.5*(Surf.IsolatedSurface[r1]+Surf.IsolatedSurface[r2]-Sum(asas));
          scsurf:=0.5*(Surf.IsolatedSidechainSurface[r1]+Surf.IsolatedSidechainSurface[r2]-Sum(scasas));;
          if fullsurf>Tiny then
            begin
            FullInnerContacts[r1,r2]:=fullsurf;
            FullInnerContacts[r2,r1]:=fullsurf;
            end;
          if scsurf>Tiny then
            begin
            SidechainInnerContacts[r1,r2]:=scsurf;
            SidechainInnerContacts[r2,r1]:=scsurf;
            end;
          end;
    end;
end;

procedure TPdbSurface.CalcInnerDistanceContacts(var Surf: TSurfaceReport;
  const Dist: TFloat);
var
  f,r1,r2:Integer;
  hulls:array of TCuboid;
  atomsr1,atomsr2,scatomsr1,scatomsr2:TAtoms;


begin
  with Surf do
    begin
    SetLength(hulls,Length(Residues));
    SetLength(FullInnerContacts,Length(Residues),Length(Residues));
    SetLength(SidechainInnerContacts,Length(Residues),Length(Residues));
    for r1:=0 to High(FullInnerContacts) do
      for r2:=0 to High(FullInnerContacts) do
        begin
        FullInnerContacts[r1,r2]:=0;
        SidechainInnerContacts[r1,r2]:=0;
        end;
    for f:=0 to High(Residues) do
      hulls[f]:=CalcCenterHull(Residues[f].AllAtoms,Dist);
    for r1:=0 to High(Residues)-1 do
      begin
      atomsr1:=Residues[r1].AllAtoms;
      scatomsr1:=NoBackbone(atomsr1);
      for r2:=r1+1 to High(Residues) do
        if InContact(hulls[r1],hulls[r2]) then
          begin
          atomsr2:=Residues[r2].AllAtoms;
          scatomsr2:=NoBackbone(atomsr2);
          if AtomsInContact(atomsr1,atomsr2,Dist) then
            begin
            FullInnerContacts[r1,r2]:=1;
            FullInnerContacts[r2,r1]:=1;
            end;
          if AtomsInContact(scatomsr1,scatomsr2,Dist) then
            begin
            SidechainInnerContacts[r1,r2]:=1;
            SidechainInnerContacts[r2,r1]:=1;
            end;
          end;
      end;
    end;
end;

function TPdbSurface.CrossContactReport(Target, Probe: TSurfaceReport
  ): TCrossContactReport;

var
  f,rt,rp:Integer;
  hullst,hullsp:array of TCuboid;
  atoms,scatoms:TAtoms;
  fullsurf,scsurf:TFloat;
  asas,scasas:TFloats;

begin
  SetLength(hullst,Length(Target.Residues));
  SetLength(hullsp,Length(Probe.Residues));
  with Result do
    begin
    SetLength(FullCrossContacts,Length(Target.Residues),Length(Probe.Residues));
    SetLength(SidechainCrossContacts,Length(Target.Residues),Length(Probe.Residues));
    for rt:=0 to High(FullCrossContacts) do
      for rp:=0 to High(FullCrossContacts[0]) do
        begin
        FullCrossContacts[rt,rp]:=0;
        SidechainCrossContacts[rt,rp]:=0;
        end;

    for f:=0 to High(hullst) do
      hullst[f]:=CalcHull(Target.Residues[f].AllAtoms,2*FProbeRadius);
    for f:=0 to High(hullsp) do
      hullsp[f]:=CalcHull(Probe.Residues[f].AllAtoms,2*FProbeRadius);

    for rt:=0 to High(hullst) do
      for rp:=0 to High(hullsp) do
        if InContact(hullst[rt],hullsp[rp]) then
          begin
          atoms:=AppendAtomsToArray(Target.Residues[rt].AllAtoms,Probe.Residues[rp].AllAtoms);
          scatoms:=NoBackbone(atoms);
          asas:=CalcASAs(atoms);
          scasas:=CalcASAs(scatoms);
          fullsurf:=0.5*(Target.IsolatedSurface[rt]+Probe.IsolatedSurface[rp]-Sum(asas));
          if fullsurf>Tiny then
            FullCrossContacts[rt,rp]:=fullsurf;
          scsurf:=0.5*(Target.IsolatedSidechainSurface[rt]+Probe.IsolatedSidechainSurface[rp]-Sum(scasas));;
          if scsurf>Tiny then
            SidechainCrossContacts[rt,rp]:=scsurf;
          end;
    end;
end;

function TPdbSurface.CrossContactReports(const Surfaces: TSurfaceReports
  ): TCrossContactReports;

var f,g,count:Integer;

begin
  count:=0;
  SetLength(Result,f*(f-1) div 2);
  for f:=0 to High(Surfaces)-1 do
    for g:=f+1 to High(Surfaces) do
      begin
      Result[count]:=CrossContactReport(Surfaces[f],Surfaces[g]);
      Result[count].TargetIx:=f;
      Result[count].ProbeIx:=g;
      Inc(count);
      end;
end;

function TPdbSurface.CrossDistanceContactReport(const Target,
  Probe: TSurfaceReport; const Dist: TFloat): TCrossContactReport;

var
  f,rt,rp:Integer;
  hullst,hullsp:array of TCuboid;
  atomsrt,atomsrp,scatomsrt,scatomsrp:TAtoms;

begin
  SetLength(hullst,Length(Target.Residues));
  SetLength(hullsp,Length(Probe.Residues));
  with Result do
    begin
    SetLength(FullCrossContacts,Length(Target.Residues),Length(Probe.Residues));
    SetLength(SidechainCrossContacts,Length(Target.Residues),Length(Probe.Residues));
    for rt:=0 to High(FullCrossContacts) do
      for rp:=0 to High(FullCrossContacts[0]) do
        begin
        FullCrossContacts[rt,rp]:=0;
        SidechainCrossContacts[rt,rp]:=0;
        end;

    for f:=0 to High(hullst) do
      hullst[f]:=CalcCenterHull(Target.Residues[f].AllAtoms,Dist);
    for f:=0 to High(hullsp) do
      hullsp[f]:=CalcCenterHull(Probe.Residues[f].AllAtoms,Dist);

    for rt:=0 to High(hullst) do
      begin
      atomsrt:=Target.Residues[rt].AllAtoms;
      scatomsrt:=NoBackbone(atomsrt);
      for rp:=0 to High(hullsp) do
        if InContact(hullst[rt],hullsp[rp]) then
          begin
          atomsrp:=Probe.Residues[rp].AllAtoms;
          scatomsrp:=NoBackbone(atomsrp);
          if AtomsInContact(atomsrt,atomsrp,Dist) then
            FullCrossContacts[rt,rp]:=1;
          if AtomsInContact(scatomsrt,scatomsrp,Dist) then
            SidechainCrossContacts[rt,rp]:=1;
          end;
      end
    end;
end;

function TPdbSurface.CrossDistanceContactReports(
  const Surfaces: TSurfaceReports; const Dist: TFloat): TCrossContactReports;

var f,g,count:Integer;

begin
  count:=0;
  SetLength(Result,f*(f-1) div 2);
  for f:=0 to High(Surfaces)-1 do
    for g:=f+1 to High(Surfaces) do
      begin
      Result[count]:=CrossDistanceContactReport(Surfaces[f],Surfaces[g],Dist);
      Result[count].TargetIx:=f;
      Result[count].ProbeIx:=g;
      Inc(count);
      end;
end;

{
function TPdbSurface.SurfaceReport(PdbFile: string; Chains: TSimpleStrings): TSurfaceReport;


begin

end;
 }
end.

