{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain. Enjoy.
********************************************************************************
Author: Ludwig Krippahl
Date: 16.7.2011
Purpose:
  Surface calculations for generic sets of spheres
Requirements:
Revisions:
To do:
  Implement LCPO (http://en.wikipedia.org/wiki/Accessible_surface_area)
*******************************************************************************}

unit surface;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, geomhash, threedcalc;

function GoldenSpiralPoints(Num:Integer):TCoords;
  //Returns Num points distributed around a sphere centered at origin, radius 1
  //uses the Golden Section spiral algorithm, based on an implementation by
  //Patrick Boucher at http://www.softimageblog.com/archives/115

function SRSurface(Points:TCoords;Radii:TFloats;SpherePoints:TCoords;MinHashCell:Single=0):TFloats;
  //Based on the Shrake-Rupley ASA algorithm. SpherePoints is the base set of points
  //at the surface of a radius 1 sphere, centered at the origin, to estimate
  //surface area. Returns the surface area for each sphere defined by Points and Radii
  //MinHashCell is the minimum size of the geometric hash grid cells, which is also
  //at least as large as twice the largest atomic radius
  //NOTE: for ASA add probe radius to radius of each atom.


implementation

function GoldenSpiralPoints(Num:Integer):TCoords;

var
  inc,off,y,r,phi:TFloat;
  k:Integer;

begin
  SetLength(Result,Num);
  inc:=Pi * (3 - Sqrt(5));
  off:= 2 / Num;
  for k:=0 to Num-1 do
    begin
    y := k * off - 1 + (off / 2);
    r := Sqrt(1 - y*y);
    phi := k * inc;
    Result[k,0]:=Cos(phi)*r;
    Result[k,1]:=y;
    Result[k,2]:= Sin(phi)*r;
    end;
end;

function SRSurface(Points:TCoords;Radii:TFloats;SpherePoints:TCoords;MinHashCell:Single=0):TFloats;

function LargestRadius:TFloat;

var f:Integer;

begin
  Result:=0;
  for f:=0 to high(Radii) do
    if Radii[f]>Result then Result:=Radii[f];
end;

var
  ghash:TGeomHasher;
  sphere:TCoord;
  f,g,h,countouts:Integer;
  ixs:TIntegers;
  intersect:Boolean;

begin
  SetLength(Result,Length(Points));
  if Points=nil then Exit;
  ghash:=TGeomHasher.Create(Points,LargestRadius*2);
  for f:=0 to High(Points) do
    begin
    ixs:=ghash.ListNeighours(Points[f]);
    countouts:=0;
    for g:=0 to High(SpherePoints) do
      begin
      sphere:=AddVectors(ScaleVector(SpherePoints[g],Radii[f]),Points[f]);
      intersect:=False;
      for h:=0 to High(Ixs) do
        if (ixs[h]<>f) and (Distance(sphere,Points[ixs[h]])<Radii[ixs[h]]) then
          begin
          intersect:=True;
          Break;
          end;
      if not intersect then Inc(countouts);
      end;
    Result[f]:=4*Pi*Sqr(Radii[f])/Length(SpherePoints)*countouts;
    end;
  ghash.Free;
end;

end.

