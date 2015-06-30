{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 1.6.2011
Purpose:
  Line grids for 3D shape representation, and computing surface and core regions
  NOTE: unlike original BiGGER implementation, linegrid segments are aligned
  with the Z coordinate, for a more intuitive [X,Y,Z] indexation of the matrix.
Requirements:
Revisions:
To do:

*******************************************************************************}

unit linegrids;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, geomutils, geomhash;

type

  TLineSegment=array [0..1] of Integer;  //start, end of segment
  TGridLine=array of TLineSegment;       //disjunt segments sorted from lower to high
  TLineArray=array of TGridLine;
  TGridPlane=array of array of TGridLine;

  TGridShape=record
      Grid:TGridPlane;
      NonEmpty:TLineArray;    //X oriented array of non-empty gridlines
      CellCounts:array of array of Integer;
      //CellsBelowMiddle,CellsAboveMiddle:array of array of Integer;
        //tested, and does not improve performance because when core splits domain there
        //is already too much overlap for more pruning with cellcounts
      TotalCount:Integer;
      ZMax:Integer;
  end;

  //For domain representation

  TDomainBlock=record
    //values for domain cuboid hull
    ProbeEndX,ProbeEndY,ProbeEndZ,TargetEndX,TargetEndY,TargetEndZ:Integer;
    XOffset,YOffset,ZOffset:Integer;
    DomainEndX,DomainEndY,DomainEndZ:Integer;
  end;


  TDomainGrid=record
    Block:TDomainBlock;
    Shape:TGridShape;                 //full domain
  end;




  { TDockingGrid }
  //grids associated with a structure. Includes base, core and surface
  TDockingGrid=class
    protected
      //atomic coordinates are shifted so that min-rad-FSurfThickness==0
      FCoords,FCoreCutCoords:TCoords;
      FRads,FCoreCutRads:TFloats;
      FResolution:TFloat;
      FTransVec:TCoord;
        //translation vector for adjusting coordinates
        //from coords to grid: add
        //from grid to coords: subtract
      FBase,FSurf,FCore:TGridShape;
      procedure SetSurfCoreLine(X,Y:Integer;var Neighs,Surfs,Cores:TIntegers);
      procedure SetLineVals(Line:TGridLine;var ValArray:TIntegers;Val:Integer);
      procedure UpdateNeighbours(Line:TGridLine;var Neighs:TIntegers;Mult:Integer);
      procedure ResetBaseLines(X,Y:Integer;var Neighs,Surfs,Cores:TIntegers);
        //recomputes neighbour counts and resets surf and core vals to zero
      procedure BaseLineIncY(X:Integer; var Y:Integer;var Neighs,Surfs,Cores:TIntegers);
        //incrementally changes Y to next value and updates
        //neighbour surf and core vals.

      procedure BuildBaseGrid;      //translates structure and creates FBase
      procedure BuildSurfCoreGrids; //builds from Base grid

    public
      property Base:TGridShape read FBase;
      property Core:TGridShape read FCore;
      property Surf:TGridShape read FSurf;
      property Resolution:TFloat read FResolution;
      property TransVec:TCoord read FTransVec;
      constructor Create(AResolution:TFloat);
      procedure BuildFromSpheres(Coords:TCoords;Rads:TFloats;
                CoreCuts:TCoords=nil;CoreCutRads:TFloats=nil);

    end;

  TDockingGrids=array of TDockingGrid;

  function Intersect(const Line1,Line2:TGridLine):TGridLine;overload;
  procedure Intersect(var Line:TGridLine; const NewMin,NewMax:Integer);
  procedure DisplaceLine(var Line:TGridLine;Displacement:Integer);
  function IntegersToLine(const Ints:TIntegers;const Threshold:Integer=0;
              Limit1:Integer=-1;Limit2:Integer=-1):TGridLine;overload;
  procedure SetIntegersLine(const Line:TGridLine;var Ints:TIntegers; Val:Integer);
  procedure ComputeShapeStats(var Shape:TGridShape);
  function CountSegments(const Grid:TGridPlane):Integer;
  function CountCells(const Grid:TGridPlane):Integer;overload;
  function CountCells(const Line:TGridLine):Integer;overload;
  function CountCellsBetween(const Line:TGridLine;const MinZ,MaxZ:Integer):Integer;
  procedure GetLineExtremes(const Line:TGridLine; out LineMin,LineMax:Integer);
  function NewLine(Ix1,Ix2:Integer):TGridLine;

implementation

function Intersect(const Line1, Line2: TGridLine): TGridLine;

var
  ix,ni1,ni2:Integer;
  i1,i2,ll1,ll2:Integer;
  top,bot:Integer;
begin
  SetLength(Result,Length(Line1)+Length(Line2));
  ix:=0;
  i1:=0;
  i2:=0;
  ll1:=Length(Line1);
  ll2:=Length(Line2);
  while (i1<ll1) and (i2<ll2) do
    begin
    top:=Min(Line1[i1,1],Line2[i2,1]);
    bot:=Max(Line1[i1,0],Line2[i2,0]);
    if top>=bot then
      begin
      Result[ix,0]:=bot;
      Result[ix,1]:=top;
      Inc(ix);
      end;
    if Line1[i1,1]>=Line2[i2,1] then ni2:=i2+1 else ni2:=i2;
    if Line1[i1,1]<=Line2[i2,1] then ni1:=i1+1 else ni1:=i1;
    i1:=ni1;
    i2:=ni2;
    end;
  SetLength(Result,ix);
end;

procedure Intersect(var Line: TGridLine; const NewMin, NewMax: Integer);

var
  f,lastvalid:Integer;

begin
  lastvalid:=-1;
  for f:=0 to High(Line) do
    begin
    if (Line[f,1]>=NewMin) and (Line[f,0]<=NewMax) then
      //line segment is inside new interval
      begin
      Line[f,0]:=Max(NewMin,Line[f,0]);
      Line[f,1]:=Min(NewMax,Line[f,1]);
      Inc(lastvalid);
      if lastvalid<f then
       Line[lastvalid]:=Line[f];
      end
    end;
  SetLength(Line,lastvalid+1);
end;

procedure DisplaceLine(var Line: TGridLine; Displacement: Integer);

var f:Integer;

begin
  for f:=0 to High(Line) do
    begin
    Line[f,0]:=Line[f,0]+Displacement;
    Line[f,1]:=Line[f,1]+Displacement;
    end;
end;

function IntegersToLine(const Ints: TIntegers; const Threshold:Integer; Limit1:Integer;
  Limit2: Integer): TGridLine;

var
  f,curr:Integer;
  isin:Boolean;

begin
  if Limit1<0 then Limit1:=0;
  if Limit2<0 then Limit2:=High(Ints);
  SetLength(Result,Limit2-Limit1+1);
  curr:=-1;
  isin:=False;
  for f:=Limit1 to Limit2 do
    begin
    if Ints[f]>Threshold then
      begin
        if not isin then
          begin
          isin:=True;
          Inc(curr);
          Result[curr,0]:=f;
          end;
      end
    else if isin then
      begin
      Result[curr,1]:=f-1;
      isin:=False;
      end;
    end;
  if isin then Result[curr,1]:=Limit2;
  SetLength(Result,curr+1);
end;

procedure SetIntegersLine(const Line: TGridLine; var Ints: TIntegers;
  Val: Integer);

var f,ff:Integer;

begin
  for f:=0 to High(Line) do
    for ff:=Line[f,0] to Line[f,1] do
      Ints[ff]:=Val;
end;

procedure ComputeShapeStats(var Shape: TGridShape);

var
  x,y:Integer;
  tmpline:TIntegers;

begin
  with Shape do
    begin
    SetLength(NonEmpty,Length(Grid));
    SetLength(CellCounts,Length(Grid),Length(Grid[0]));
    //SetLength(CellsBelowMiddle,Length(Grid),Length(Grid[0]));
    //SetLength(CellsAboveMiddle,Length(Grid),Length(Grid[0]));
    ZMax:=0;
    TotalCount:=0;
    for x:=0 to High(Grid) do
      begin
      tmpline:=FilledInts(Length(Grid[0]),0);
      for y:=0 to High(Grid[0]) do
        if Grid[x,y]<>nil then
          begin
          tmpline[y]:=1;
          if ZMax<Grid[x,y,High(Grid[x,y]),1] then
            ZMax:=Grid[x,y,High(Grid[x,y]),1];
          CellCounts[x,y]:=CountCells(Grid[x,y]);
          TotalCount:=TotalCount+CellCounts[x,y];
          //CellsBelowMiddle[x,y]:=CountCellsBetween(Grid[x,y],0,ZMax div 2);
          //CellsAboveMiddle[x,y]:=CountCellsBetween(Grid[x,y],ZMax div 2,ZMax);
          end;
      NonEmpty[x]:=IntegersToLine(tmpline,0);
      end;
    end;
end;

function CountSegments(const Grid:TGridPlane):Integer;

var x,y:Integer;

begin
  Result:=0;
  for x:=0 to High(Grid) do
    for y:=0 to High(Grid[x]) do
      Result:=Result+Length(Grid[x,y]);
end;

function CountCells(const Grid:TGridPlane):Integer;

var x,y,z:Integer;

begin
  Result:=0;
  for x:=0 to High(Grid) do
    for y:=0 to High(Grid[x]) do
      for z:=0 to High(Grid[x,y]) do
        Result:=Result+Grid[x,y,z,1]-Grid[x,y,z,0]+1;
end;

function CountCells(const Line: TGridLine): Integer;

var z:Integer;

begin
  Result:=0;
  for z:=0 to High(Line) do
    Result:=Result+Line[z,1]-Line[z,0]+1;
end;

function CountCellsBetween(const Line: TGridLine; const MinZ, MaxZ: Integer
  ): Integer;

var z:Integer;

begin
  Result:=0;
  for z:=0 to High(Line) do
    if (Line[z,0]<=MaxZ) and (Line[z,1]>=MinZ) then
      Result:=Result+Min(MaxZ,Line[z,1])-Max(Line[z,0],MinZ)+1;
end;

procedure GetLineExtremes(const Line: TGridLine; out LineMin, LineMax: Integer);

begin
  if Line=nil then
    begin
    LineMin:=-1;
    LineMax:=-2;
    end
  else
    begin
    LineMin:=Line[0,0];
    LineMax:=Line[High(Line),1];
    end;
end;

function NewLine(Ix1, Ix2: Integer): TGridLine;
begin
  SetLength(Result,1);
  Result[0,0]:=Ix1;
  Result[0,1]:=Ix2;
end;



{ TDockingGrid }

procedure TDockingGrid.SetSurfCoreLine(X, Y: Integer; var Neighs, Surfs,
  Cores: TIntegers);

var z,zz:Integer;

begin
  for z:=0 to High(FBase.Grid[X,Y]) do
    for zz:=FBase.Grid[X,Y,z,0] to FBase.Grid[X,Y,z,1] do
      if Neighs[zz]>=27 then
        Cores[zz]:=1
      else Surfs[zz]:=1;
  FCore.Grid[X,Y]:=IntegersToLine(Cores);
  FSurf.Grid[X,Y]:=IntegersToLine(Surfs);
end;

procedure TDockingGrid.SetLineVals(Line: TGridLine; var ValArray: TIntegers;
  Val: Integer);

var z,zz:Integer;

begin
  for z:=0 to High(Line) do
    for zz:=Line[z,0] to Line[z,1] do
      ValArray[zz]:=Val;
end;

procedure TDockingGrid.UpdateNeighbours(Line: TGridLine; var Neighs: TIntegers;
  Mult: Integer);

var z,zz,bot,top:Integer;

begin
  for z:=0 to High(Line) do
    begin
    bot:=Line[z,0];
    top:=Line[z,1];
    Neighs[bot-1]:=Neighs[bot-1]+Mult;
    Neighs[top+1]:=Neighs[top+1]+Mult;
    Neighs[bot]:=Neighs[bot]+Mult;
    if bot<top then
      begin
      Neighs[bot]:=Neighs[bot]+Mult;
      Neighs[top]:=Neighs[top]+2*Mult;
      end;
    for zz:=bot+1 to top-1 do
      Neighs[zz]:=Neighs[zz]+3*Mult;
    end;

end;

procedure TDockingGrid.ResetBaseLines(X, Y: Integer; var Neighs, Surfs,
  Cores: TIntegers);

var f,g:Integer;

begin
  for f:=0 to High(Neighs) do
    begin
    Neighs[f]:=0;
    Surfs[f]:=0;
    Cores[f]:=0;
    end;
  for f:=X-1 to X+1 do
    for g:=Y-1 to Y+1 do
      with FBase do
        UpdateNeighbours(Grid[f,g],Neighs,1);
end;

procedure TDockingGrid.BaseLineIncY(X: Integer; var Y: Integer; var Neighs,
  Surfs, Cores: TIntegers);

begin
  UpdateNeighbours(FBase.Grid[X-1,Y-1],Neighs,-1);
  UpdateNeighbours(FBase.Grid[X,Y-1],Neighs,-1);
  UpdateNeighbours(FBase.Grid[X+1,Y-1],Neighs,-1);

  SetLineVals(FBase.Grid[X,Y],Surfs,0);
  SetLineVals(FBase.Grid[X,Y],Cores,0);
  Inc(Y);


  UpdateNeighbours(FBase.Grid[X-1,Y+1],Neighs,1);
  UpdateNeighbours(FBase.Grid[X,Y+1],Neighs,1);
  UpdateNeighbours(FBase.Grid[X+1,Y+1],Neighs,1);

end;

procedure TDockingGrid.BuildBaseGrid;

var
  hash:TGeomHasher;
  halfres,maxrad:TFloat;
  top:TCoord;
  xyzpoint:TCoord;
  x,y,z:Integer;
  zline:TIntegers;


begin
  //Adjust coordinates
  maxrad:=Max(Max(FRads),FResolution);

  //Translate and size guaranteeing one layer of empty grid cells
  FTransVec:=Min(FCoords);
  FTransVec:=Add(Simmetric(FTransVec),maxrad+FResolution);
  FCoords:=Add(FTransvec,FCoords);
  top:=Add(Max(FCoords),maxrad+1.5*FResolution);
  SetLength(FBase.Grid,Round(top[0]/FResolution),Round(top[1]/FResolution));

  SetLength(zline,Round(top[2]/FResolution));
  FBase.ZMax:=High(zline);
  hash:=TGeomHasher.Create(FCoords,maxrad,FRads);
  halfres:=0.5*FResolution;
  for x:=0 to High(FBase.Grid) do
    begin
    xyzpoint[0]:=x*FResolution+halfres;
    for y:=0 to High(FBase.Grid[x]) do
      begin
      xyzpoint[1]:=y*FResolution+halfres;
      for z:=0 to High(zline) do
        begin
        xyzpoint[2]:=z*FResolution+halfres;
        if hash.IsInnerPoint(xyzpoint) then
          zline[z]:=1
        else zline[z]:=0;
        FBase.Grid[x,y]:=IntegersToLine(zline);
        end;
      end;
    end;
  hash.Free;
end;

procedure TDockingGrid.BuildSurfCoreGrids;

var
  x,y:Integer;
  neighs,surfs,cores:TIntegers;

begin
  SetLength(FCore.Grid,Length(FBase.Grid),Length(FBase.Grid[0]));
  SetLength(FSurf.Grid,Length(FBase.Grid),Length(FBase.Grid[0]));
  SetLength(neighs,FBase.ZMax+1);
  SetLength(surfs,FBase.ZMax+1);
  SetLength(cores,FBase.ZMax+1);
  for x:=1 to High(FBase.Grid)-1 do
    begin
    y:=1;
    ResetBaseLines(x,y,neighs,surfs,cores);
    repeat
      SetSurfCoreLine(x,y,neighs,surfs,cores);
      BaseLineIncY(x,y,neighs,surfs,cores);
    until y>=High(FBase.Grid[x])-1;
    SetSurfCoreLine(x,y,neighs,surfs,cores);
    end;
end;

constructor TDockingGrid.Create(AResolution:TFLoat);
begin
  inherited Create;
  FTransVec:=NullVector;
  FResolution:=AResolution;
end;

procedure TDockingGrid.BuildFromSpheres(Coords: TCoords; Rads: TFloats;
  CoreCuts: TCoords; CoreCutRads: TFloats);

begin
  Assert(Length(Coords)=Length(Rads),'Coordinates and radius values do not match');
  FCoords:=Copy(Coords,0,Length(Coords));
  FRads:=Copy(Rads,0,Length(Rads));
  FCoreCutCoords:=Copy(CoreCuts,0,Length(CoreCuts));
  FCoreCutRads:=Copy(CoreCutRads,0,Length(CoreCutRads));
  BuildBaseGrid;
  BuildSurfCoreGrids;
  ComputeShapeStats(FSurf);
  ComputeShapeStats(FCore);
  FCoords:=nil;
  FRads:=nil;
  FCoreCutCoords:=nil;
  FCoreCutRads:=nil;
end;


end.

