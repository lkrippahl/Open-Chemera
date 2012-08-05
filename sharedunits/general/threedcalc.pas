{*******************************************************************************
This source file is public domain. It is provided as is, with no explicit
or implied warranties regarding anything whatsoever.
********************************************************************************
Author: Ludwig Krippahl
Date: 31.3.2011
Purpose:
  Geometric docking module. Imported from original BiGGER.
Revisions:
TO DO:
  check function names and usefulness
*******************************************************************************}

{$mode objfpc}{$H+}

unit threedcalc;

interface

uses basetypes;

type
  TMatrix=array[0..2] of TCoord;
  TMatrixArray=array of TMatrix;
const
  MCIdentityMatrix:TMatrix=((1,0,0),(0,1,0),(0,0,1));
  MCXYZRotation=0;
  MCRotationXYAxis=1;

function AddToVector(const v1:TCoord;r:TFloat):TCoord;
function LargestVector(const v1,v2:TCoord):TCoord;
function SmallestVector(const v1,v2:TCoord):TCoord;
function DotProduct(const v1,v2:TCoord):TFloat;
function Normalise(const v1: TCoord): TCoord;
function CrossProduct(const v1, v2: TCoord): TCoord;
function BuildBase(const v1, v2: TCoord): TMatrix;
function RotMatrix(const m1, m2: TMatrix): TMatrix;
function InvertBase(const m: TMatrix): TMatrix;
function GetPlacementMatrix(f1, f2, m1, m2: TCoord): TMatrix; overload;
function GetPlacementMatrix(f1, m1: TCoord): TMatrix; overload;
function Distance(const Point1,Point2:TCoord):TFLoat;
function RotVector(const m: TMAtrix; const v: TCoord): TCoord;
function AddVectors(const v1,v2:TCoord):TCoord;
function PointsToVector(const Origin,Destination:TCoord):TCoord;
function SubtractVectors(const Initial,Subtract:TCoord):TCoord;
function ScaleVector(const v1:TCoord; scale:TFloat):TCoord;
function SimmetricVector(const v1:TCoord):TCoord;
function MidPoint(const v1,v2:TCoord):TCoord;
function BuildRotation(const R1:TCoord; RotType:Integer):TMatrix;
function OldBuildRotation(const R1:TCoord):TMatrix; //for old cdock compatibility only
function VectorSize(const v:TCoord):Double;
function EqualVectors(const v1,v2:TCoord):Boolean;
function RotateXToXY(const v1:TCoord):TMatrix;
function TransformVector(const Pivot:TCoord;Rotation:TMatrix;Place,Vector:TCoord):TCoord;
function RotAndPlace(const Rotation:TMatrix;const Place,Vector:TCoord):TCoord;overload;
procedure RotAndPlace(const Rotation:TMatrix;const Place:TCoord;var Coords:TCoords);overload;
function XRotation(Angle:TFloat):TMatrix;
function YRotation(Angle:TFloat):TMatrix;
function ZRotation(Angle:TFloat):TMatrix;
function DihedralAngle(const c1, c2, c3, c4: TCoord):TFloat;
function RotationAngle(Axis,Pivot,P1,P2:TCoord):TFloat;
function AxisRotationMatrix(Axis:TCoord;Angle:TFloat):TMatrix;
function AngularDifference(A1,A2:TFloat):TFloat;
function AngleIsBetween(Low,High,Angle:TFloat):Boolean;
function FixAngle(Angle:TFloat):TFloat;
function Foot(const LinePoint,LineVector,Point:TCoord):TCoord;
function TopCorner(const c1,c2:TCoord):TCoord;
function BottomCorner(const c1,c2:TCoord):TCoord;
function SmallestAxis(c:TCoord):Integer;
function LongestCoord(const c:TCoord):TFloat;

implementation

uses Math;

const TOL=0.00001;

function AddToVector(const v1:TCoord;r:TFLoat):TCoord;

begin
  Result[0]:=v1[0]+r;
  Result[1]:=v1[1]+r;
  Result[2]:=v1[2]+r;
end;

function LargestVector(const v1,v2:TCoord):TCoord;

begin
  if v1[0]>v2[0] then Result[0]:=v1[0] else Result[0]:=v2[0];
  if v1[1]>v2[1] then Result[1]:=v1[1] else Result[1]:=v2[1];
  if v1[2]>v2[2] then Result[2]:=v1[2] else Result[2]:=v2[2];
end;

function SmallestVector(const v1,v2:TCoord):TCoord;

begin
  if v1[0]<v2[0] then Result[0]:=v1[0] else Result[0]:=v2[0];
  if v1[1]<v2[1] then Result[1]:=v1[1] else Result[1]:=v2[1];
  if v1[2]<v2[2] then Result[2]:=v1[2] else Result[2]:=v2[2];
end;

function DotProduct(const v1,v2:TCoord):TFloat;

begin
  Result:=v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
end;

function XRotation(Angle:TFloat):TMatrix;

var s,c:TFloat;

begin
  s:=Sin(Angle);
  c:=Cos(Angle);
  Result[0,0]:=1;
  Result[0,1]:=0;
  Result[0,2]:=0;
  Result[1,0]:=0;
  Result[1,1]:=c;
  Result[1,2]:=-s;
  Result[2,0]:=0;
  Result[2,1]:=s;
  Result[2,2]:=c;
end;

function YRotation(Angle:TFloat):TMatrix;

var s,c:TFloat;

begin
  s:=Sin(Angle);
  c:=Cos(Angle);
  Result[0,0]:=c;
  Result[0,1]:=0;
  Result[0,2]:=s;
  Result[1,0]:=0;
  Result[1,1]:=1;
  Result[1,2]:=0;
  Result[2,0]:=-s;
  Result[2,1]:=0;
  Result[2,2]:=c;
end;

function ZRotation(Angle:TFloat):TMatrix;

var s,c:TFloat;

begin
  s:=Sin(Angle);
  c:=Cos(Angle);
  Result[0,0]:=c;
  Result[0,1]:=-s;
  Result[0,2]:=0;
  Result[1,0]:=s;
  Result[1,1]:=c;
  Result[1,2]:=0;
  Result[2,0]:=0;
  Result[2,1]:=0;
  Result[2,2]:=1;
end;


function EqualVectors(const v1,v2:TCoord):Boolean;

begin
  Result:=(v1[0]=v2[0]) and
          (v1[1]=v2[1]) and
          (v1[2]=v2[2]);
end;


function VectorSize(const v:TCoord):Double;

begin
  Result:=sqrt(sqr(v[0])+sqr(v[1])+sqr(v[2]));
end;

function OldBuildRotation(const R1:TCoord):TMatrix;


{Ú                                        ¿}
{³ -siné  -(cosé)(cosí)  -(cosé)(siní)  0 ³}
{³                                        ³}
{³  cosé  -(siné)(cosí)  -(siné)(siní)  0 ³}
{³                                        ³}
{³   0         siní          -cosí      0 ³}
{À                                        Ù}

{Ú                                                       ¿}
{³  c1c3 + s1s2s3        -c1s3 + c3s1s2          c2s1    ³}
{³                                                       ³}
{³      c2s3                  c2c3               -s2     ³}
{³                                                       ³}
{³ -c3s1 + c1s2s3         s1s3 + c1c3s2          c1c2    ³}
{À                                                       Ù}

{var c1,c2,c3,s1,s2,s3:TFloat48;
    m1,m2:matriz;
    f,g,t:integer;
begin
     c1:=cos(c[0]*pic);
     c2:=cos(c[1]*pic);
     c3:=cos(c[2]*pic);
     s1:=sin(c[0]*pic);
     s2:=sin(c[1]*pic);
     s3:=sin(c[2]*pic);
     m[0,0]:=c1*c3+s1*s2*s3;
     m[1,0]:=c3*s1*s2-c1*s3;
     m[2,0]:=c2*s1;
     m[0,1]:=c2*s3;
     m[1,1]:=c2*c3;
     m[2,1]:=-s2;
     m[0,2]:=c1*s2*s3-c3*s1;
     m[1,2]:=s1*s3+c1*c3*s2;
     m[2,2]:=c1*c2;
end;}


var c1,c2,c3,s1,s2,s3:TFloat;

begin
     c1:=Cos(R1[0]);
     c2:=Cos(R1[1]);
     c3:=Cos(R1[2]);
     s1:=Sin(R1[0]);
     s2:=Sin(R1[1]);
     s3:=Sin(R1[2]);
     Result[0,0]:=(c1*c3+s1*s2*s3);
     Result[0,1]:=(c3*s1*s2-c1*s3);
     Result[0,2]:=(c2*s1);
     Result[1,0]:=(c2*s3);
     Result[1,1]:=(c2*c3);
     Result[1,2]:=(-s2);
     Result[2,0]:=(c1*s2*s3-c3*s1);
     Result[2,1]:=(s1*s3+c1*c3*s2);
     Result[2,2]:=(c1*c2);
end;

function BuildRotation(const R1:TCoord;RotType:Integer):TMatrix;

begin
  case RotType of
    MCXYZRotation:Result:=RotMatrix(RotMatrix(ZRotation(R1[2]),YRotation(R1[1])),XRotation(R1[0]));
    MCRotationXYAxis:
      begin
      Result:=RotMatrix(YRotation(R1[1]),XRotation(R1[0]));
      Result:=RotMatrix(InvertBase(Result),RotMatrix(ZRotation(R1[2]),Result));
      end;
    else Result:=MCIdentityMatrix;
    end;
end;

function PointsToVector(const Origin,Destination:TCoord):TCoord;
begin
  Result[0]:=Destination[0]-Origin[0];
  Result[1]:=Destination[1]-Origin[1];
  Result[2]:=Destination[2]-Origin[2];
end;

function SubtractVectors(const Initial,Subtract:TCoord):TCoord;
begin
  Result[0]:=Initial[0]-Subtract[0];
  Result[1]:=Initial[1]-Subtract[1];
  Result[2]:=Initial[2]-Subtract[2];
end;

function RotVector(const  m: TMAtrix; const v: TCoord): TCoord;
begin
     result[0]:=m[0,0]*v[0]+m[0,1]*v[1]+m[0,2]*v[2];
     result[1]:=m[1,0]*v[0]+m[1,1]*v[1]+m[1,2]*v[2];
     result[2]:=m[2,0]*v[0]+m[2,1]*v[1]+m[2,2]*v[2];
end;

function AddVectors(const v1,v2:TCoord):TCoord;

begin
  Result[0]:=v1[0]+v2[0];
  Result[1]:=v1[1]+v2[1];
  Result[2]:=v1[2]+v2[2];
end;

function MidPoint(const v1,v2:TCoord):TCoord;

begin
  Result:=ScaleVector(AddVectors(v1,v2),0.5);
end;


function ScaleVector(const v1:TCoord; scale:TFloat):TCoord;

begin
  Result[0]:=v1[0]*scale;
  Result[1]:=v1[1]*scale;
  Result[2]:=v1[2]*scale;
end;

function SimmetricVector(const v1:TCoord):TCoord;

begin
  Result[0]:=-v1[0];
  Result[1]:=-v1[1];
  Result[2]:=-v1[2];
end;

function Distance(const Point1,Point2:TCoord):Double;

begin
  Result:=Sqrt(Sqr(Point1[0]-Point2[0])+
               Sqr(Point1[1]-Point2[1])+
               Sqr(Point1[2]-Point2[2]));
end;

function Normalise(const v1: TCoord): TCoord;

var r:double;
    f:integer;

begin
     r:=sqrt(sqr(v1[0])+sqr(v1[1])+sqr(v1[2]));
     if r>0 then
        for f:=1 to 3 do result[f]:=v1[f]/r
        else
            begin
            result[0]:=0;
            result[1]:=0;
            result[2]:=0;
            end;
end;

function CrossProduct(const v1, v2: TCoord): TCoord;
begin
     Result[0]:=v1[1]*v2[2]-v1[2]*v2[1];
     Result[1]:=v1[2]*v2[0]-v1[0]*v2[2];
     Result[2]:=v1[0]*v2[1]-v1[1]*v2[0];
end;

function BuildBase(const v1, v2: TCoord): TMatrix;
{v1 e v2 normalisados}
var t,t2:TCoord;
    f:integer;

begin
     t:=CrossProduct(v1,v2);
     t:=SimmetricVector(CrossProduct(t,v1));
     t2:=CrossProduct(v1,t);
     t:=Normalise(t);
     t2:=Normalise(t2);
     for f:=1 to 3 do
         begin
         Result[f,0]:=v1[f];
         Result[f,1]:=t[f];
         Result[f,2]:=t2[f];
         end;
end;

function RotMatrix(const m1, m2: TMatrix): TMatrix;

var f,g:integer;

begin
     for f:=1 to 3 do
     for g:=1 to 3 do
     Result[f,g]:=m1[f,0]*m2[0,g]+m1[f,1]*m2[1,g]+m1[f,2]*m2[2,g];
end;

function InvertBase(const m: TMatrix): TMatrix;

var f,g:integer;

begin
     for f:=1 to 3 do
         for g:=1 to 3 do
             Result[f,g]:=m[g,f];
end;

function GetPlacementMatrix(f1, f2, m1, m2: TCoord): TMatrix; overload;

var mm,mf:TMatrix;

begin
     f1:=Normalise(f1);
     m1:=Normalise(m1);
     f2:=Normalise(f2);
     m2:=Normalise(m2);
     mm:=BuildBase(m1,m2);
     mm:=InvertBase(mm);
     mf:=BuildBase(f1,f2);
     Result:=RotMatrix(mf,mm);
end;

function GetPlacementMatrix(f1, m1: TCoord): TMatrix; overload;

var f2:TCoord;

begin
     f2:=CrossProduct(f1,m1);
     if (Abs(f2[0])>TOL) or (Abs(f2[1])>TOL) or (Abs(f2[2])>TOL) then
       Result:=GetPlacementMatrix(f1,f2,m1,f2)
     else
      begin
      Result:=MCIdentityMatrix;
      if DotProduct(f1,m1)<0 then
        begin
        Result[0,0]:=-1;
        Result[1,1]:=-1;
        end;
      end;
end;

function RotateXToXY(const v1:TCoord):TMatrix;

var
  l,s,c:TFloat;

begin
  l:=Sqrt(Sqr(v1[1])+Sqr(v1[2]));
  s:=v1[2]/l;;
  c:=v1[1]/l;
  Result[0,0]:=1;
  Result[0,1]:=0;
  Result[0,2]:=0;
  Result[1,0]:=0;
  Result[1,1]:=c;
  Result[1,2]:=s;
  Result[2,0]:=0;
  Result[2,1]:=-s;
  Result[2,2]:=c;
end;

function TransformVector(const Pivot:TCoord;Rotation:TMatrix;Place,Vector:TCoord):TCoord;

begin
  Result:=SubtractVectors(Vector,Pivot);
  Result:=RotAndPlace(Rotation,Place,Result);
end;

function RotAndPlace(const Rotation:TMatrix;const Place,Vector:TCoord):TCoord;

begin
  Result:=RotVector(Rotation,Vector);
  Result:=AddVectors(Result,Place);
end;

procedure RotAndPlace(const Rotation:TMatrix;const Place:TCoord;var Coords:TCoords);

var f:Integer;

begin
  for f:=0 to High(Coords) do
    Coords[f]:=RotAndPlace(Rotation,Place,Coords[f]);
end;

function DihedralAngle(const c1, c2, c3, c4: TCoord):TFloat;

var m:TMatrix;
    v1,v2,v3:TCoord;
    f:integer;

begin
     for f:=1 to 3 do
         begin
         v2[f]:=c1[f]-c2[f];
         v1[f]:=c3[f]-c2[f];
         v3[f]:=c4[f]-c3[f];
         end;
     v1:=normalise(v1);
     v2:=normalise(v2);
     v3:=normalise(v3);
     m:=BuildBase(v1,v2);
     m:=InvertBase(m);
     v3:=RotVector(m,v3);
     result:=arccos(v3[1]);
     if v3[2]<0 then result:=-result;
end;

function RotationAngle(Axis,Pivot,P1,P2:TCoord):TFloat;

var
  t0:TFloat;
  vs:TFloat;
begin
  Axis:=Normalise(Axis);
  t0:=DotProduct(SubtractVectors(P1,Pivot),Axis);
  Pivot:=AddVectors(Pivot,ScaleVector(Axis,t0));
  P1:=SubtractVectors(P1,Pivot);
  P2:=SubtractVectors(P2,Pivot);
  vs:=VectorSize(P1)*VectorSize(P2);
  if vs>TinyFloat then
    begin
    Result:=DotProduct(P1,P2)/vs;
    Result:=ArcCos(Result);
    P2:=CrossProduct(P1,P2);
    if DotProduct(Axis,P2)<0 then
      Result:=2*PI-Result;
    end
  else Result:=0;
end;

function AxisRotationMatrix(Axis:TCoord;Angle:TFloat):TMatrix;

var Tmp:TMatrix;

begin
  Axis:=Normalise(Axis);
  Tmp:=GetPlacementMatrix(Coord(0,0,1),Axis);
  Result:=ZRotation(Angle);
  Result:=RotMatrix(InvertBase(Tmp),RotMatrix(Result,Tmp));
end;

function AngularDifference(A1,A2:TFloat):TFloat;

begin
  Result:=FixAngle(A2-A1);
end;

function AngleIsBetween(Low,High,Angle:TFloat):Boolean;

begin
  Result:=(High>=Angle) and (Low<=Angle);
end;

function FixAngle(Angle:TFloat):TFloat;

begin
  Result:=Angle-2*PI*Trunc(Angle/(2*PI));
  if Result<0 then
  Result:=Result+2*PI;
end;

function Foot(const LinePoint,LineVector,Point:TCoord):TCoord;

begin
  Result:=AddVectors(LinePoint,ScaleVector(LineVector,
    DotProduct(SubtractVectors(Point,LinePoint),LineVector)/
    Sqr(VectorSize(LineVector))));
end;

function TopCorner(const c1,c2:TCoord):TCoord;

begin
  if c1[0]>c2[0] then Result[0]:=c1[0] else Result[0]:=c2[0];
  if c1[1]>c2[1] then Result[1]:=c1[1] else Result[1]:=c2[1];
  if c1[2]>c2[2] then Result[2]:=c1[2] else Result[2]:=c2[2];
end;

function BottomCorner(const c1,c2:TCoord):TCoord;

begin
  if c1[0]<c2[0] then Result[0]:=c1[0] else Result[0]:=c2[0];
  if c1[1]<c2[1] then Result[1]:=c1[1] else Result[1]:=c2[1];
  if c1[2]<c2[2] then Result[2]:=c1[2] else Result[2]:=c2[2];
end;

function SmallestAxis(c:TCoord):Integer;
begin
  Result:=1;
  c[0]:=Abs(c[0]);
  c[1]:=Abs(c[1]);
  c[2]:=Abs(c[2]);
  if c[1]<c[0] then Result:=2;
  if c[2]<c[Result] then Result:=3;
end;

function LongestCoord(const c:TCoord):TFloat;

begin
  Result:=Abs(c[0]);
  if Abs(c[1])>Result then Result:=Abs(c[1]);
  if Abs(c[2])>Result then Result:=Abs(c[2]);
end;

end.

