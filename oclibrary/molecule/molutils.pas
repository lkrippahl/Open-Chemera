{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 23.1.2014
Purpose:
  Utility functions for molecules and (selected) parts of molecules
Revisions:
To do:
*******************************************************************************}

unit molutils;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, geomutils, molecules, selections, geomhash;

procedure CenterMolecule(Molecule:TMolecule;Selection:TSelection=nil);
function FindCenter(Molecule:TMolecule;Selection:TSelection=nil):TCoord;
function ListCoords(Molecule:TMolecule;Selection:TSelection=nil):TCoords;overload;
function ListCoords(const Molecule:TMolecule;const AtomName:string):TCoords;overload;
function ListRadii(Molecule:TMolecule;Selection:TSelection=nil):TFloats;
function AtomsInContact(const Atoms1,Atoms2:TAtoms;const Dist:TFloat):Boolean;
  //contact is center to center distance, not surface to surface
function CalcHull(Atoms:TAtoms;Rad:TFloat):TCuboid;
  //Returns the cuboid region defined by all atom surfaces expanded by the Rad value
function CalcHulls(Groups:TMolecules;Rad:TFloat):TCuboids;
    //Returns the cuboid regions defined by all atom surfaces expanded by the Rad value

function CalcCenterHull(Atoms:TAtoms;Dist:TFloat):TCuboid;
  //Returns the cuboid region defined by all atom centers expanded by the Dist value
procedure GroupsInContact(const Groups1, Groups2: TMolecules; const Dist: TFloat;
            out Interface1,Interface2:TIntegers);
  //returns indexes of groups of each Mol1 and Mol2 within distance of the other
  //groups must be terminal (with only atoms, not groups)
  //TODO: make more efficient

procedure GroupContacts(const Groups1, Groups2: TMolecules; const Dist: TFloat;
            out Interface1,Interface2:TIntegers);
  //returns indexes of each pairwise contact from Mol1 to Mol2 within distance
  //groups must be terminal (with only atoms, not groups)
  //TODO: make more efficient



function NeighbourIndexes(FromCoords,ToCoords:TCoords;Dist:TFloat):TIntegers;
  //returns an array with the indexes of FromCoords that are within Dist of any ToCoords
  //uses geomhash for efficiency


implementation

procedure CenterMolecule(Molecule: TMolecule; Selection: TSelection);

var
  cent:TCoord;

begin
  cent:=FindCenter(Molecule,Selection);
  Molecule.Transform(Simmetric(cent));
end;

function FindCenter(Molecule: TMolecule; Selection: TSelection): TCoord;

var
  coords:TCoords;
  f:Integer;

begin
  coords:=Molecule.AllCoords;
  Result:=MidPoint(coords);
end;

function ListCoords(Molecule: TMolecule; Selection: TSelection): TCoords;

var
  atoms:TAtoms;
  f:Integer;

begin
  atoms:=Molecule.AllAtoms;
  SetLength(Result,Length(atoms));
  for f:=0 to High(atoms) do
    Result[f]:=atoms[f].Coords;
end;

function ListCoords(const Molecule: TMolecule; const AtomName: string): TCoords;

var
  atoms:TAtoms;
  f,count:Integer;

begin
  atoms:=Molecule.AllAtoms;
  SetLength(Result,Length(atoms));
  count:=0;
  for f:=0 to High(atoms) do
    if atoms[f].Name=AtomName then
      begin
      Result[count]:=atoms[f].Coords;
      Inc(count);
      end;
  SetLength(Result,count);
end;

function ListRadii(Molecule: TMolecule; Selection: TSelection): TFloats;
var
  atoms:TAtoms;
  f:Integer;

begin
  atoms:=Molecule.AllAtoms;
  SetLength(Result,Length(atoms));
  for f:=0 to High(atoms) do
    Result[f]:=atoms[f].Radius;
end;

function AtomsInContact(const Atoms1, Atoms2: TAtoms; const Dist: TFloat): Boolean;

{ TODO : Geometric hashing if many atoms }


var
  f,g:Integer;

begin
  Result:=False;
  for f:=0 to High(Atoms1) do
    for g:=0 to High(Atoms2) do
      if Distance(Atoms1[f].Coords,Atoms2[g].Coords)<=Dist then
        begin
        Result:=True;
        Break;
        end;
end;


function CalcHull(Atoms:TAtoms;Rad:TFloat):TCuboid;

var
  f:Integer;
  c1,c2:TCoord;

begin
  if Atoms=nil then
    begin
    c1:=NullVector;
    c2:=NullVector;
    end
  else
    begin
    c1:=Atoms[0].Coords;
    c2:=c1;
    for f:=1 to High(Atoms) do
      begin
      c1:=Min(c1,Add(Atoms[f].Coords,-Atoms[f].Radius));
      c2:=Max(c2,Add(Atoms[f].Coords,Atoms[f].Radius))
      end;
    for f:=0 to 2 do
      begin
      c1[f]:=c1[f]-Rad;
      c2[f]:=c2[f]+Rad;
      end;
    end;
  Result[0]:=c1;
  Result[1]:=c2;
end;

function CalcHulls(Groups: TMolecules; Rad: TFloat): TCuboids;

var f:Integer;

begin
  SetLength(Result,Length(Groups));
  for f:=0 to High(Groups) do
    Result[f]:=CalcHull(Groups[f].AllAtoms,Rad);
end;

function CalcCenterHull(Atoms: TAtoms; Dist: TFloat): TCuboid;

var
  f:Integer;
  c1,c2:TCoord;

begin
  if Atoms=nil then
    begin
    c1:=NullVector;
    c2:=NullVector;
    end
  else
    begin
    c1:=Atoms[0].Coords;
    c2:=c1;
    for f:=1 to High(Atoms) do
      begin
      c1:=Min(c1,Atoms[f].Coords);
      c2:=Max(c2,Atoms[f].Coords)
      end;
    c1:=Subtract(c1,Dist);
    c2:=Add(c2,Dist);
    end;
  Result[0]:=c1;
  Result[1]:=c2;
end;

procedure GroupsInContact(const Groups1, Groups2: TMolecules; const Dist: TFloat;
          out Interface1,Interface2:TIntegers);

var
  hulls1,hulls2:TCuboids;
  f,ix1,ix2:Integer;
  interfixs1,interfixs2:TIntegers;

begin
  hulls1:=CalcHulls(Groups1,Dist);
  hulls2:=CalcHulls(Groups2,Dist);
  interfixs1:=FilledInts(Length(hulls1),0);
  interfixs2:=FilledInts(Length(hulls2),0);
  for ix1:=0 to High(hulls1) do
    for ix2:=0 to High(hulls2) do
      if InContact(hulls1[ix1],hulls2[ix2]) and
        AtomsInContact(Groups1[ix1].GroupAtoms,Groups2[ix2].GroupAtoms,Dist) then
          begin
          interfixs1[ix1]:=1;
          interfixs2[ix2]:=1;
          end;
  Interface1:=nil;
  for f:=0 to High(interfixs1) do
    if interfixs1[f]>0 then
        AddToArray(f,Interface1);
  Interface2:=nil;
  for f:=0 to High(interfixs2) do
    if interfixs2[f]>0 then
        AddToArray(f,Interface2);
end;

procedure GroupContacts(const Groups1, Groups2: TMolecules; const Dist: TFloat;
  out Interface1, Interface2: TIntegers);

var
  hulls1,hulls2:TCuboids;
  f,ix1,ix2:Integer;
  interfixs1,interfixs2:TIntegers;

begin
  Interface1:=nil;
  Interface2:=nil;
  hulls1:=CalcHulls(Groups1,Dist);
  hulls2:=CalcHulls(Groups2,Dist);
  for ix1:=0 to High(hulls1) do
    for ix2:=0 to High(hulls2) do
      if InContact(hulls1[ix1],hulls2[ix2]) and
        AtomsInContact(Groups1[ix1].GroupAtoms,Groups2[ix2].GroupAtoms,Dist) then
          begin
          AddToArray(ix1,Interface1);
          AddToArray(ix2,Interface2);
          end;
end;

function NeighbourIndexes(FromCoords, ToCoords: TCoords; Dist: TFloat
  ): TIntegers;

var
  f,ix:Integer;
  hasher:TGeomHasher;
  tmprads:TFloats;

begin
  SetLength(Result,Length(FromCoords));
  tmprads:=FilledFloats(Length(ToCoords),Dist);
  hasher:=TGeomHasher.Create(ToCoords,Dist,tmprads);
  ix:=0;
  for f:=0 to High(FromCoords) do
    if hasher.IsInnerPoint(FromCoords[f]) then
        begin
        Result[ix]:=f;
        Inc(ix);
        end;
  SetLength(Result,ix);
  hasher.Free;
end;

end.

