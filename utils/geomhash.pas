{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 20.3.2011
Purpose:
  Geometric hashing class
Requirements:
Revisions:
To do:
*******************************************************************************}
unit geomhash;

{$mode objfpc}

interface

uses
  Classes, SysUtils, basetypes, threedcalc;

type

  { TGeomHasher }

  TGeomHasher=class
  //The size of the grid cell must be > 1e-6 and the points array must not
  //be empty
  private
    FHashGrid:array of array of array of TOCIntegers;  //indexes of hashed points
    FShiftToGrid:TOCCoord;
    FInvGridStep:TOCFloat;                             // 1/gridstep
    FHighX,FHighY,FHighZ:Integer;                      //max grid cells
  public
    constructor Create(Points:TOCCoords;GridStep:TOCFloat);
    function ListNeighours(C:TOCCoord):TOCIntegers;
      //Returns all indexes in the 27 grid cells belonging to the neighbourhood
      //of the cell corresponding to C
  end;

implementation

{ TGeomHasher }

constructor TGeomHasher.Create(Points:TOCCoords;GridStep:TOCFloat);

  procedure Setup;

  var
    f:Integer;
    minc,maxc:TOCCoord;

  begin
  maxc:=Points[0];
  minc:=maxc;
  for f:=1 to High(Points) do
    begin
    minc:=Min(minc,Points[f]);
    maxc:=Max(maxc,Points[f]);
    end;
  FShiftToGrid:=ScaleVector(minc,-1);
  FHighX:=Trunc((maxc[0]-minc[0])/GridStep);
  FHighY:=Trunc((maxc[1]-minc[1])/GridStep);
  FHighZ:=Trunc((maxc[2]-minc[2])/GridStep);
  SetLength(FHashGrid,FHighX+1,FHighY+1,FHighZ+1);
  end;

var
  f,x,y,z:Integer;
  c:TOCCoord;

begin
  inherited Create;
  Assert(Points<>nil,'Empty points array for hashing');
  Assert(GridStep>1e-6,'GridStep too small');
  GridStep:=1/GridStep;
  FInvGridStep:=GridStep;
  Setup;
  for x:=0 to FHighX do
    for y:=0 to FHighY do
      for z:=0 to FHighZ do
        FHashGrid[x,y,z]:=nil;
  for f:=0 to High(Points) do
    begin
    c:=AddVectors(Points[f],FShiftToGrid);
    x:=Trunc(c[0]*GridStep);
    y:=Trunc(c[1]*GridStep);
    z:=Trunc(c[2]*GridStep);
    AddToArray(f, FHashGrid[x,y,z]);
    end;
end;

function TGeomHasher.ListNeighours(C:TOCCoord):TOCIntegers;


var
  count,x,y,z,x1,x2,y1,y2,z1,z2,f:Integer;

procedure Bounds(out B1,B2:Integer; const Val:TOCFloat; const Hi:Integer);

begin
  B1:=Trunc(Val)-1;
  B2:=Trunc(Val)+1;
  if B1>Hi then B1:=Hi;
  if B1<0 then B1:=0;
  if B2>Hi then B2:=Hi;
  if B2<0 then B2:=0;
end;

begin
  C:=ScaleVector(AddVectors(C,FShiftToGrid),FInvGridStep);
  Result:=nil;
  Bounds(x1,x2,C[0],FHighX);
  Bounds(y1,y2,C[1],FHighY);
  Bounds(z1,z2,C[2],FHighZ);
  count:=0;
  for x:=x1 to x2 do
    for y:=y1 to y2 do
      for z:=z1 to z2 do
       count:=count+Length(FHashGrid[x,y,z]);
  SetLength(Result,count);
  count:=0;
  for x:=x1 to x2 do
    for y:=y1 to y2 do
      for z:=z1 to z2 do
        for f:=0 to High(FHashGrid[x,y,z]) do
          begin
          Result[count]:=FHashGrid[x,y,z,f];
          Inc(count);
          end;
end;

end.

