{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 23.12.2015
Purpose:
  Generate rotations for docking
Revisions:
TO DO:
  optimized rotations. Receives array of coords for protein and a cutoff value in
  angstrom
  Returns list of rotations guaranteeing distance
*******************************************************************************}

unit rotations;

{$mode objfpc}{$H+}

interface

uses
  SysUtils,basetypes,geomutils,surface,Classes; //classes for debug only

type

  { TRotationManager }

  TRotationManager = class
    private
    FPoints:TCoords;
       //points defining the shape of the protein
    FAxes:TCoords;
    FRotationsPerAxis:TIntegers;
    FSelectedAxes:TIntegers;
    FAxialDistance:TFloat;

    FParaMatrix,FOrtoMatrix:TMatrix;
     //these matrices, axes x points, have the distances
     //paralel and orthogonal to each axis for each point

    FAxisDistances:TMatrix;
     //squared distances between axes

    FForcedSingleRotation:TFloat;
     //Set >0 to force one single rotation per quaternion
    FRotations:TQuaternions;
    FRotationSelectedAxis:TIntegers;
    FRotationPoints:TCoordGroups;

    function AxisDistance(Axis1,Axis2:Integer):TFloat;
      //Maximum atom displacement for two axes

    function AxisDistances(Axis:Integer):TFloats;
      //Returns distance from one axis to all others

    function CoverByColumn(Distances:TMatrix;Displacement:TFloat):TIntegers;
      //Returns rows covered by each column
    function CountBelow(const Vec:TFloats; ColCoverCount,Cols:TIntegers; Thresh:TFloat):TFloat;
      //counts number of columns below threshold
    function ListBelow(const Vec:TFloats; Cols:TIntegers; Thresh:TFloat):TIntegers;
      //returns indexes of Cols for those below threshold
    function MaxCoverageRow(const Distances:TMatrix;const Displacement:TFloat;
                                const ColCoverCount,Rows,Cols:TIntegers):Integer;
      //returns index of Rows for row with mose coverage


    procedure PopMaxCoverageRow(const Distances:TMatrix;
                                const Displacement:TFloat;
                                const ColCoverCount:TIntegers;
                                var Rows,Cols,Selecters:TIntegers);
     //Finds row in Rows with best coverage of columns Cols
     //removes this from rows
     //adds to Cols
     //Due to optimization, rows and cols can be unsorted.

    function MinDistances(Distances: TMatrix; RowIxs: TIntegers): TFloats;overload;
     //Minimum per column of distances in rows RowIxs


    function CreateQuaternions(Axis:Integer;MinDist:TFloat):TQuaternions;
     //Creates a set of quaternions for the axis spaced according to the displacement
     //if FForcedSingleRotation>0 then only that rotation is created

    function IsRedundant(Points:TCoords;MinDist:TFloat):Boolean;
     //True if there is no set of points in FRotationPoints with no point outside mindist

    function AddQuaternions(Axis:Integer;MinDist,Filter:TFloat;SelectedIndex:Integer):Integer;
     //Adds a filtered set of quaternions to FRotations

    public
    property AxialDistance:TFloat read FAxialDistance;


    function Rotations:TQuaternions;
    function SelectedAxes:TCoords;
    function AxisIndexes:TIntegers;

    constructor Create(APoints:TCoords);
        //copies points
    procedure ClearMatrices;
    procedure InitializeAxes(BaseCount:Integer);
        // Create initial axes and distance matrix
        // clears
    procedure ComputeAxisDistances;
    procedure SelectAxesByCoverage(Displacement:TFloat);
    procedure SelectAxesByMaxmin(MinDist: TFloat;MaxAxes:Integer);
        //does not use FAxialDistances matrix
    procedure ForceSingleRotation(Monomers:Integer);
    procedure GenerateQuaternions(MinDist,Filter:TFloat);

    procedure MaxDistQuaternions(BaseCount:Integer;Displacement:TFloat;
                                Symmetry:Integer=-1;MaxAxes:Integer=-1);
        //Generate axes and quaterniosns using the maximum distance greedy optimization
        //does not create the distance matrix

    function GetSelectedAxes:TCoords;
    function MinDist(Axis:Integer):TFloat;
    function RotationCount(Axis:Integer):Integer;
    end;

function BaseSampleAxes(Count:Integer):TCoords;
   // Count axes spread over hemisphere, using the spiral algorithm

//these are deprecated (leave them?)
function UniformSampleAngles(ZSteps:Integer;out Axes:TCoords):TQuaternions;
function FixedZSampleAngles(ZSteps:Integer;ZRot:TFloat;out Axes:TCoords;Margin:Integer=0):TQuaternions;
function EulerSampleZYX(ZSteps:Integer):TQuaternions;

function SpacedPoints(Coords: TCoords; NumPoints: Integer): TIntegers;
  //returns a set of NumPoints integers with indexes for the most widely spaced points

function SpacedCoords(Coords: TCoords; NumPoints: Integer): TCoords;
  //returns an array of coords for the most widely spaced points


implementation


function TRotationManager.AxisDistance(Axis1, Axis2: Integer): TFloat;

var
  g:Integer;
  x,y,tmpa,tmpb:TFloat;
  p1x,p1y,p2x,p2y:TFloat;
  sp1x,sp2x:TFloat;
  a1x,a1y:TFloat;
  dist1,dist2,dist,tmp:TFloat;

begin
  if Axis1 = Axis2 then
    Result:=0
  else
    begin
    OrthogonalCoords(FAxes[Axis1],FAxes[Axis2],x,y);
    dist:=0;
    for g:=0 to High(FPoints) do
      begin
      a1x:=FParaMatrix[Axis1,g];
      a1y:=FOrtoMatrix[Axis1,g];

      //points are paralel*(x,y)+- ortho*(y,-x)
      tmpa:=FParaMatrix[Axis2,g]*x;
      tmpb:=FOrtoMatrix[Axis2,g]*y;
      p1x:=tmpa+tmpb;
      p2x:=tmpa-tmpb;
      tmpa:=FParaMatrix[Axis2,g]*y;
      tmpb:=-FOrtoMatrix[Axis2,g]*x;
      p1y:=tmpa+tmpb;
      p2y:=tmpa-tmpb;

      sp1x:=Sqr(p1x-a1x);
      sp2x:=Sqr(p2x-a1x);

      dist1:=Min(sp1x+Sqr(p1y-a1y),sp2x+Sqr(p2y-a1y));
      dist2:=Min(sp1x+Sqr(p1y+a1y),sp2x+Sqr(p2y+a1y));
      dist:=Max(dist,Max(dist1,dist2));

      dist1:=Min(sp1x+Sqr(p1y-a1y),sp1x+Sqr(p1y+a1y));
      dist2:=Min(sp2x+Sqr(p2y-a1y),sp2x+Sqr(p2y+a1y));
      dist:=Max(dist,Max(dist1,dist2));
      end;
    Result:=Sqrt(dist);
    end;
end;

function TRotationManager.AxisDistances(Axis: Integer): TFloats;

var f:Integer;

begin
  SetLength(Result,Length(FAxes));
  for f:=0 to High(Result) do
    Result[f]:=AxisDistance(Axis,f);
end;

function TRotationManager.CoverByColumn(Distances: TMatrix; Displacement: TFloat
  ): TIntegers;

var f,g:Integer;

begin
  SetLength(Result,Length(Distances));
  for f:=0 to High(Distances) do
    begin
    Result[f]:=0;
    for g:=0 to High(Distances) do
      if Distances[g,f]<Displacement then
        Result[f]:=Result[f]+1;
    end;
end;

function TRotationManager.CountBelow(const Vec:TFloats;
                            ColCoverCount,Cols:TIntegers; Thresh:TFloat):TFloat;

var f:Integer;

begin
  Result:=0;
  for f:=0 to High(Cols) do
    if Thresh>Vec[Cols[f]] then
      Result:=Result+1/(ColCoverCount[Cols[f]]+1);
end;

function TRotationManager.ListBelow(const Vec: TFloats; Cols: TIntegers;
  Thresh: TFloat): TIntegers;

var f,ix:Integer;

begin
  ix:=0;
  SetLength(Result,Length(Cols));
  for f:=0 to High(Cols) do
    if Thresh>Vec[Cols[f]] then
      begin
      Result[ix]:=f;
      Inc(ix);
      end;
  SetLength(Result,ix);
end;

function TRotationManager.MaxCoverageRow(const Distances: TMatrix;
  const Displacement: TFloat; const ColCoverCount, Rows, Cols: TIntegers
  ): Integer;

var
  f:Integer;
  count,maxcount:TFloat;

begin
  Result:=-1;
  maxcount:=-1;
  for f:=0 to High(Rows) do
    begin
    count := CountBelow(Distances[Rows[f]],ColCoverCount,Cols,Displacement);
    if count>maxcount then
      begin
      maxcount:=count;
      Result:=f;
      end;
    end;
end;

procedure TRotationManager.PopMaxCoverageRow(const Distances: TMatrix;
  const Displacement: TFloat; const ColCoverCount:TIntegers;
  var Rows, Cols, Selecters: TIntegers);

var
  f,maxix:Integer;
  coldrops:TIntegers;


begin
  maxix :=MaxCoverageRow(Distances,Displacement,ColCoverCount,Rows, Cols);
  coldrops:=ListBelow(Distances[Rows[maxix]],Cols,Displacement);
  //writeln('drops',length(coldrops));
  AddToArray(Rows[maxix],Selecters);
  RemoveFromArray(maxix,Rows);
  //writeln('rows',Length(Rows));
  RemoveFromArray(coldrops,Cols);
  //WriteLn('cols',Length(Cols));
end;

function TRotationManager.MinDistances(Distances: TMatrix; RowIxs: TIntegers
  ): TFloats;

var
  f,g:Integer;
  row:TFloats;

begin
  if RowIxs<>nil then
    begin
    Result:=Copy(Distances[RowIxs[0]],0,Length(Distances[RowIxs[0]]));
    for f:=1 to High(RowIxs) do
      begin
      row:=Distances[RowIxs[f]];
      for g:=0 to High(Result) do
        if Result[g]>row[g] then Result[g]:=row[g];
      end;
    end
  else Result:=nil
end;

function TRotationManager.CreateQuaternions(Axis: Integer; MinDist: TFloat): TQuaternions;

var
  f:Integer;
  step,maxdist:TFloat;

begin
  if FForcedSingleRotation>0 then
    begin
    SetLength(Result,1);
    Result[0]:=RotationQuaternion(FAxes[Axis],FForcedSingleRotation);
    end
  else
    begin
    maxdist:=Max(FOrtoMatrix[Axis]);
    SetLength(Result,Round(2*maxdist*Pi/MinDist));
    if Result<>nil then
      begin
      step:=2*maxdist*Pi/Length(Result);
      for f:=0 to High(Result) do
        Result[f]:=RotationQuaternion(FAxes[Axis],f*step);
      end;
    end;
end;

function TRotationManager.IsRedundant(Points: TCoords;
          MinDist: TFloat): Boolean;

var
  f,g:Integer;
  dist,maxdist:TFloat;

begin
  Result:=False;
  for f:=0 to High(FRotationPoints) do
    begin
    Result:=True;
    for g:=0 to High(Points) do
      if Distance(Points[g],FRotationPoints[f,g])>MinDist then
        begin
        Result:=False;
        Break;
        end;
    if Result then Break;
    end;
end;

function TRotationManager.AddQuaternions(Axis: Integer;
          MinDist,Filter: TFloat;SelectedIndex:Integer):Integer;

var
  quaternions:TQuaternions;
  f:Integer;
  toinsert:TBooleans;
  insertcount:Integer;
  points:TCoordGroups;
  oldlen:Integer;

begin
  quaternions:=CreateQuaternions(Axis,MinDist);
  SetLength(points,Length(quaternions));
  SetLength(toinsert,Length(quaternions));
  insertcount:=0;
  for f:=0 to High(quaternions) do
    begin
    points[f]:=Rotate(FPoints,quaternions[f]);
    toinsert[f]:= not IsRedundant(points[f],Filter);
    if toinsert[f] then Inc(insertcount);
    end;
  if insertcount>0 then
    begin
    oldlen:=Length(FRotations);
    SetLength(FRotations,oldlen+insertcount);
    SetLength(FRotationSelectedAxis,oldlen+insertcount);
    SetLength(FRotationPoints,oldlen+insertcount);
    for f:=0 to High(toinsert) do
      if toinsert[f] then
        begin
        FRotations[oldlen]:=quaternions[f];
        FRotationSelectedAxis[oldlen]:=SelectedIndex;
        FRotationPoints[oldlen]:=points[f];
        Inc(oldlen);
        end;
    end;
  Result:=insertcount;
end;

function TRotationManager.Rotations: TQuaternions;
begin
  Result:=Copy(FRotations,0,Length(FRotations));
end;

function TRotationManager.SelectedAxes: TCoords;

var f:Integer;

begin
  SetLength(Result,Length(FSelectedAxes));
  for f:=0 to High(Result) do
    Result[f]:=FAxes[FSelectedAxes[f]];
end;

function TRotationManager.AxisIndexes: TIntegers;
begin
  Result:=Copy(FRotationSelectedAxis,0,Length(FRotationSelectedAxis));
end;

constructor TRotationManager.Create(APoints: TCoords);
begin
  inherited Create;
  FPoints:=Copy(APoints,0,Length(APoints));
end;

procedure TRotationManager.ClearMatrices;
begin
  FAxes:=nil;
  FSelectedAxes:=nil;
  FRotations:=nil;
  FParaMatrix:=nil;
  FOrtoMatrix:=nil;
  FAxisDistances:=nil;
end;

procedure TRotationManager.InitializeAxes(BaseCount:Integer);

var
  f,g:Integer;
  x,y:TFloat;

begin
  FAxes:=BaseSampleAxes(BaseCount);
  SetLength(FParaMatrix,Length(FAxes),Length(FPoints));
  SetLength(FOrtoMatrix,Length(FAxes),Length(FPoints));
  for f:=0 to High(FAxes) do
    for g:=0 to High(FPoints) do
      begin
      OrthogonalCoords(FAxes[f],FPoints[g],x,y);
      FParaMatrix[f,g]:=x;
      FOrtoMatrix[f,g]:=y;
      end;
end;

procedure TRotationManager.ComputeAxisDistances;

var
  ax1,ax2:Integer;
  dist:TFloat;

  //debug
  sl:TStringList;
  s:string;

begin
  SetLength(FAxisDistances,Length(FAxes),Length(FAxes));
  FAxisDistances[High(FAxes),High(FAxes)]:=0;
  for ax1:=0 to High(FAxes)-1 do
    begin
    FAxisDistances[ax1,ax1]:=0;
    for ax2:=ax1+1 to High(FAxes) do
      begin
      dist:=AxisDistance(ax1,ax2);
      FAxisDistances[ax1,ax2]:=dist;
      FAxisDistances[ax2,ax1]:=dist;
      end;
    end;
    {
    sl:=TStringList.Create;
    for ax1:=0 to High(FAxes) do
      begin
      s:='';
      for ax2:=0 to High(FAxes) do
        s:=s+FLoatToStrF(Sqrt(FAxisDistances[ax1,ax2]),ffFixed,5,2)+#9;
      sl.add(s);
      end;
    sl.SaveToFile('H:/matrix.txt');
    sl.free;
    }
end;

procedure TRotationManager.SelectAxesByMaxmin(MinDist: TFloat;MaxAxes:Integer);

var
  maxd:TFloat;
  maxix:Integer;
  mins:TFLoats;
  rows:TIntegers;
  distances:TMatrix;

begin
  FSelectedAxes:=FilledInts(1,High(FAxes));
      //This selects the last axis as te first point. The last axis is the Z axis
  rows:=FilledInts(1,0);  //this is just a 0...N index to avoid using the FAxialDistances matrix
  SetLength(distances,MaxAxes);
  distances[0]:=AxisDistances(FSelectedAxes[0]);
  maxd := MinDist+1;
  while (maxd>MinDist) and (Length(FSelectedAxes)<MaxAxes) do
    begin
    mins:=MinDistances(distances,rows);
    MaxValIx(mins, maxd, maxix);
    AddToArray(maxix,FSelectedAxes);
    AddToArray(High(FSelectedAxes),rows);
    distances[High(FSelectedAxes)]:=AxisDistances(maxix);
    end;
  FAxialDistance:=maxd;
  //writeln(maxd);
  //writeln(Length(FSelectedAxes));
end;

procedure TRotationManager.ForceSingleRotation(Monomers: Integer);
begin
  if Monomers>0 then
    FForcedSingleRotation:=2*Pi/Monomers
  else FForcedSingleRotation:=-1;
end;

procedure TRotationManager.GenerateQuaternions(MinDist, Filter: TFloat);

var f:Integer;

begin
  FRotations:=nil;
  FRotationSelectedAxis:=nil;
  FRotationPoints:=nil;
  SetLength(FRotationsPerAxis,Length(FSelectedAxes));
  for f:=0 to High(FSelectedAxes) do
    FRotationsPerAxis[f]:=AddQuaternions(FSelectedAxes[f],MinDist,Filter,f);
end;

procedure TRotationManager.MaxDistQuaternions(BaseCount: Integer;
  Displacement: TFloat; Symmetry: Integer; MaxAxes: Integer);
begin
  if MaxAxes<1 then MaxAxes:=BaseCount;
  ForceSingleRotation(Symmetry);
  InitializeAxes(BaseCount);
  //writeln('Initialized');
  SelectAxesByMaxmin(Displacement,MaxAxes);
  //writeln('Selected');
  GenerateQuaternions(Displacement*2,Displacement);
  //writeln('Generated');
end;

procedure TRotationManager.SelectAxesByCoverage(Displacement: TFloat);

var
  rows,columns:TIntegers;
  f,g:Integer;
  sl:TStringList;
  s:string;
  colcovercount:TIntegers;

begin
  FSelectedAxes:=nil;
  colcovercount:=CoverByColumn(FAxisDistances,Displacement);
  SetLength(rows,Length(FAxes));
  SetLength(columns,Length(FAxes));
  for f:=0 to High(rows) do
    begin
    rows[f]:=f;
    columns[f]:=f;
    end;
  while columns<>nil do
    PopMaxCoverageRow(FAxisDistances, Displacement, colcovercount, rows, columns, FSelectedAxes);
  {
  sl:=TStringList.Create;
  for f:=0 to High(FSelectedAxes) do
      begin
      s:='';
      for g:=0 to High(FSelectedAxes) do
        s:=s+FLoatToStrF(Sqrt(FAxisDistances[FSelectedAxes[f],FSelectedAxes[g]]),ffFixed,5,2)+#9;
      sl.add(s);
      end;
    sl.SaveToFile('H:/matrix_reduced.txt');
  sl.Free;

  sl:=TStringList.Create;
  for f:=0 to High(FSelectedAxes) do
      sl.add(IntToStr(FSelectedAxes[f]));
  sl.SaveToFile('H:/list.txt');
  sl.Free;
  }
  //writeln(Length(FAxes),Length(FSelectedAxes));

end;

function TRotationManager.GetSelectedAxes: TCoords;

var f:Integer;

begin
  SetLength(Result,Length(FSelectedAxes));
  for f:=0 to High(FSelectedAxes) do
    Result[f]:=FAxes[FSelectedAxes[f]];
end;

function TRotationManager.MinDist(Axis: Integer): TFloat;

var f:Integer;

begin
  Result:=1e6;
  for f:=0 to High(FSelectedAxes) do
    if (f<>Axis) and (FAxisDistances[FSelectedAxes[Axis],FSelectedAxes[f]]<Result) then
      Result:=FAxisDistances[FSelectedAxes[Axis],FSelectedAxes[f]];
  Result:=Sqrt(Result);
end;

function TRotationManager.RotationCount(Axis: Integer): Integer;
begin
  Result:=FRotationsPerAxis[Axis];
end;

function BaseSampleAxes(Count: Integer): TCoords;
//returns approximately count axis with z>=0

var
  f,ix:Integer;

begin
  Result:=GoldenSpiralPoints(Count*2);
  ix:=0;
  for f:=0 to High(Result) do
    if Result[f,2]>=0 then
      begin
      Result[ix]:=Result[f];
      //writeln('*',f,result[ix,0],result[ix,1],result[ix,2]);
      Inc(ix);
      end;
  //WriteLn(Count,ix);
  SetLength(Result,ix)
end;


function UniformSampleAngles(ZSteps:Integer;Out Axes:TCoords):TQuaternions;

//Rotates around each axis by ZSteps. Axes are generated by points distributed
//around sphere (only for z>0, to eliminate symmetry redundancy)

var
sphere:TCoords;
f,g,ix:Integer;
rot:TFloat;
zaxis:TCoord;
begin

if ZSteps=1 then
  begin
  SetLength(Result,1);
  Result[0]:=IdentityQuaternion;
  end
else
  begin
  sphere:=GoldenSpiralPoints(ZSteps*ZSteps);
  SetLength(Axes,length(sphere)*ZSteps+1); //maximum lengths, to be cut
  SetLength(Result,Length(sphere)*ZSteps+1);
  Result[0]:=IdentityQuaternion;
  ix:=1;
  for f:=0 to High(sphere) do
    begin
    if sphere[f,2]>=0 then //only top hemisphere
      for g:=1 to ZSteps-1 do
        begin
        rot := 2*Pi/(ZSteps)*g;
        Result[ix]:=RotationQuaternion(sphere[f],rot);
        Inc(ix);
        end;
    end;
  SetLength(Axes,ix);
  SetLength(Result,ix);
  end;
end;

function FixedZSampleAngles(ZSteps: Integer; ZRot: TFloat; out Axes:TCoords;Margin:Integer=0): TQuaternions;
//Rotates around Z axis ZSteps, for each rotates pole to one point in sphere

var
ix,f,g:Integer;
sphere:TCoords;

procedure PushRotation(Axis:TCoord;Rot:TFloat);

begin
  Result[ix]:=RotationQuaternion(Axis,Rot);
  Axes[ix]:=Axis;
  Inc(ix);
end;

begin
if ZSteps=1 then
  begin
  SetLength(Result,1);
  Result[0]:=IdentityQuaternion;
  end
else
  begin
  //get ZSteps^2 points uniformly distributed around a sphere
  sphere:=GoldenSpiralPoints(ZSteps*ZSteps);
  SetLength(Result,Length(sphere)*(2*Margin+1)+1);
  SetLength(Axes,Length(sphere)*(2*Margin+1)+1);
  ix:=0;
  for f:=0 to High(sphere) do
    if sphere[f,2]>=0 then
      begin
      PushRotation(sphere[f],ZRot);
      for g:=1 to Margin do
        begin
        PushRotation(sphere[f],ZRot+g*2*Pi/ZSteps);
        PushRotation(sphere[f],ZRot-g*2*Pi/ZSteps);
        end;
      end;
  end;
  SetLength(Result,ix);
  SetLength(Axes,ix);
end;

function EulerSampleZYX(ZSteps: Integer): TQuaternions;


var
x,y,z,ix:Integer;
qz,qy,qx:TQuaternion;

begin
//get ZSteps^2 points uniformly distributed around a sphere
SetLength(Result,ZSteps*ZSteps*ZSteps div 2);
ix:=0;
for z:=0 to ZSteps-1 do
  begin
  qz:=RotationQuaternion(Coord(0,0,1),2*PI/ZSteps*z);
  for y:=0 to ZSteps-1 do
    begin
    qy:=RotationQuaternion(Coord(0,1,0),2*PI/ZSteps*y);
    for x:=0 to (ZSteps div 2)-1 do
      begin
      qx:=RotationQuaternion(Coord(1,0,0),2*PI/ZSteps*x);
      Result[ix]:=Multiply(Multiply(qz,qy),qx);
      Inc(ix);
      end;
    end;
  end;
end;

function FarthestPoint(const Points,Subset:TCoords;LastPoint:Integer=-1):Integer;
  //return point with largest minimum distance from subset of points

var
  f,g:Integer;
  best,tmp,d:TFloat;

begin
  if LastPoint<0 then
    LastPoint:= High(Subset);
  Result := -1;
  best := -1;
  for f:=0 to High(Points) do
    begin
    tmp := Distance(Points[f],Subset[0]);
    for g := 1 to LastPoint do
      begin
      d := Distance(Points[f],Subset[g]);
      if d<tmp then
        tmp := d;
      end;
    if tmp>best then
      begin
      best := tmp;
      Result := f;
      end;
    end;
end;

function SpacedPoints(Coords: TCoords; NumPoints: Integer): TIntegers;

var
  subset:TCoords;
  f,ix:Integer;

begin
  if NumPoints>Length(Coords) then NumPoints:=Length(Coords);
  SetLength(Result,NumPoints);
  SetLength(subset,NumPoints);
  subset[0]:=MidPoint(Coords);
  ix:=FarthestPoint(Coords,subset,0);
  subset[0]:=Coords[ix];
  Result[0]:=ix;
  for f:=1 to NumPoints-1 do
    begin
    Result[f]:=FarthestPoint(Coords,subset,f-1);
    subset[f]:=Coords[Result[f]];
    //writeln(result[f]);
    end;
end;

function SpacedCoords(Coords: TCoords; NumPoints: Integer): TCoords;

var
  selixs:TIntegers;
  f:Integer;

begin
  selixs:=SpacedPoints(Coords,NumPoints);
  SetLength(Result,Length(selixs));
  for f:=0 to High(selixs) do
    Result[f]:=Coords[selixs[f]];
end;


end.

