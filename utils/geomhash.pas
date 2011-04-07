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
  //TO USE: create, AssignCoords, Hash. Then use NeighbourWithin to query
  private
    FHashGrid:array of array of array of TOCIntegers; //indexes of hashed coords
    FCoords:TOCCoords;
    FMinC,FMaxC,FMinusMinc:TOCCoord;
    FGridStep,FOneOverGridStep:TOCFloat;
    FHighX,FHighY,FHighZ:Integer;                     //High values for hashgrid
    procedure FindMaxMin;
    procedure AddToGrid(C:TOCCoord;Ix:Integer);
  public
    property Coords:TOCCoords read FCoords;
    procedure AssignCoords(const ACoords:TOCCoords);   //copies coordinates array
    procedure Hash(AGridStep:TOCFloat);
    function NeighboursWithin(C:TOCCoord; Dist:TOCFloat):TOCIntegers;
  end;

implementation

{ TGeomHasher }

procedure TGeomHasher.FindMaxMin;

var f:Integer;

begin
  if FCoords<>nil then
    begin
    FMaxC:=FCoords[0];
    FMinC:=FMaxC;
    for f:=1 to High(FCoords) do
      begin
      FMinC:=Min(FMinC,FCoords[f]);
      FMaxC:=Max(FMaxC,FCoords[f]);
      end;
    FMinusMinC:=ScaleVector(FMinC,-1);
    end;
end;

procedure TGeomHasher.AddToGrid(C: TOCCoord; Ix: Integer);

var x,y,z:Integer;

begin
  C:=AddVectors(ScaleVector(C,FOneOverGridStep),FMinusMinC);
  x:=Trunc(C[0]);
  y:=Trunc(C[1]);
  z:=Trunc(C[2]);
  Assert((x>=0) and (x<=FHighX) and
         (y>=0) and (y<=FHighX) and
         (z>=0) and (z<=FHighX),'Out of hasgrid bounds');
  AddToArray(Ix, FHashGrid[x,y,z]);
end;

procedure TGeomHasher.AssignCoords(const ACoords: TOCCoords);
begin
  FCoords:=Copy(ACoords, 0, Length(ACoords));
  FHashGrid:=nil;
end;

procedure TGeomHasher.Hash(AGridStep: TOCFloat);

var f,g,h:Integer;

begin
  FindMaxMin;
  FGridStep:=AGridStep;
  FOneOverGridStep:=1/AGridStep;
  FHighX:=Trunc((FMaxC[0]-FMinC[0])/FGridStep);
  FHighY:=Trunc((FMaxC[1]-FMinC[1])/FGridStep);
  FHighZ:=Trunc((FMaxC[2]-FMinC[2])/FGridStep);
  SetLength(FHashGrid,FHighX+1,FHighY+1,FHighZ+1);
  for f:=0 to FHighX do
    for g:=0 to FHighY do
      for h:=0 to FHighZ do
        FHashGrid[f,g,h]:=nil;
  for f:=0 to High(FCoords) do
    AddToGrid(FCoords[f],f);
end;

function TGeomHasher.NeighboursWithin(C:TOCCoord; Dist: TOCFloat): TOCIntegers;

var
  x,y,z,x1,x2,y1,y2,z1,z2,f:Integer;

procedure Bounds(out B1,B2:Integer; const Val,Min:TOCFloat; const Hi:Integer);

begin
  B1:=Trunc(Val-Min-Dist);
  B2:=Trunc(Val-Min+Dist);
  if B1>Hi then B1:=Hi;
  if B1<0 then B1:=0;
  if B2>Hi then B2:=Hi;
  if B2<0 then B2:=0;
end;

begin
  Result:=nil;
  Bounds(x1,x2,C[0],FMinC[0],FHighX);
  Bounds(y1,y2,C[1],FMinC[1],FHighY);
  Bounds(z1,z2,C[2],FMinC[2],FHighZ);
  for x:=x1 to x2 do
    for y:=y1 to y2 do
      for z:=z1 to z2 do
       for f:=0 to High(FHashGrid[x,y,z]) do
        if Distance(FCoords[FHashGrid[x,y,z,f]],C)<=Dist then
          AddToArray(FHashGrid[x,y,z,f],Result);
end;

end.

