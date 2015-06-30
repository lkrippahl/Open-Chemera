{*******************************************************************************
This source file is public domain. It is provided as is, with no explicit
or implied warranties regarding anything whatsoever.
********************************************************************************
Author: Ludwig Krippahl
Date: 18.2.2012
Purpose:
  Vector, matrix and other geometry utilities (rotation matrixes, etc)
  (replaces deprecated threedcalc)
Revisions:
*******************************************************************************}
unit geomutils;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes;

type
  TRotMatrix=array[0..2,0..2] of TFloat;
  TRotMatrixArray=array of TRotMatrix;
  TQuaternion=array[0..3] of TFloat;
  TQuaternions=array of TQuaternion;

const
  IdentityMatrix:TRotMatrix=((1,0,0),(0,1,0),(0,0,1));
  IdentityQuaternion:TQuaternion=(1,0,0,0);
  PI=3.1415926535897932;

  //rotation types
  MCXYZRotation=0;
  MCRotationXYAxis=1;


function Quaternion(r,i,j,k:TFloat):TQuaternion;

function Add(v1:TCoord;r:TFloat):TCoord;overload;
function Add(v1,v2:TCoord):TCoord;overload;
function Add(vec:TCoord;Coords:TCoords):TCoords;overload;
function Add(Floats:TFloats;Float:TFloat):TFloats;overload;
function Add(Coords:TCoords;Floats:TFloats):TCoords;overload;

function Subtract(v1:TCoord;r:TFloat):TCoord;overload;
function Subtract(v1,v2:TCoord):TCoord;overload;
function Subtract(Coords:TCoords;vec:TCoord):TCoords;overload;

function Multiply(Vec:TCoord;R:TFloat):TCoord;overload;
function Multiply(const Vals:TFloats;Mult:TFloat):TFloats;
function Multiply(const Cs:TCoords;Mult:TFloat):TCoords;
function Multiply(Q1,Q2:TQuaternion):TQuaternion;overload;
function Multiply(Mat1,Mat2:TRotMatrix):TRotMatrix;overload;

function DotProduct(v1,v2:TCoord):TFloat;
function CrossProduct(v1, v2: TCoord): TCoord;

function Norm(Vec:TCoord):TFloat;

function Conjugated(const Quat:TQuaternion):TQuaternion;

procedure Normalize(var Vec:TCoord);overload;
procedure Normalize(var Quat:TQuaternion);overload;

function Scaled(const Vec:TCoord; const Scale:TFloat):TCoord;

function Simmetric(Vec:TCoord):TCoord;overload;
function Simmetric(Coords:TCoords):TCoords;overload;

function Distance(Vec1,Vec2:TCoord):TFloat;overload;
function Distance(Vecs1,Vecs2:TCoords):TFloats;overload;
function Distance(Quart1,Quart2:TQuaternion):TFloat;overload; //euclidean

function Rotate(Vec:TCoord;RotMat:TRotMatrix):TCoord;overload;
function Rotate(Vecs:TCoords;RotMat:TRotMatrix):TCoords;overload;
function Rotate(Vec:TCoord;Quaternion:TQuaternion):TCoord;overload;
function Rotate(Vecs:TCoords;Quaternion:TQuaternion):TCoords;overload;

function RotationQuaternion(Axis:TCoord;Rotation:TFloat):TQuaternion;
//Axis must be normalized
//Returns Quaternion defining rotaion around axis

function RotationTo(const VFrom,VTo:TCoord):TQuaternion;
// Vectors are assumed to be normalized

function StaticRMSD(Coords1,Coords2:TCoords):TFloat;
// Assumes same number of coords in both arrays


//Old functions, with matrices
{ TODO : Replace these with quaternions? }
function RotAndPlace(Rotation:TRotMatrix;Place,Vector:TCoord):TCoord;overload;
procedure RotAndPlace(Rotation:TRotMatrix;Place:TCoord;var Coords:TCoords);overload;
function BuildRotation(R1:TCoord;RotType:Integer):TRotMatrix;
function XRotation(Angle:TFloat):TRotMatrix;
function YRotation(Angle:TFloat):TRotMatrix;
function ZRotation(Angle:TFloat):TRotMatrix;
function InvertBase(M: TRotMatrix): TRotMatrix;


{

function DotProduct(v1,v2:TCoord):TFloat;
function Normalise(v1: TCoord): TCoord;
function CrossProduct(v1, v2: TCoord): TCoord;
function BuildBase(v1, v2: TCoord): TRotMatrix;
function RotMatrix(m1, m2: TRotMatrix): TRotMatrix;
function InvertBase(m: TRotMatrix): TRotMatrix;
function GetPlacemenTRotMatrix(f1, f2, m1, m2: TCoord): TRotMatrix; overload;
function GetPlacemenTRotMatrix(f1, m1: TCoord): TRotMatrix; overload;
function Distance(Point1,Point2:TCoord):Double;
function RotVector(m: TRotMatrix; v: TCoord): TCoord;
function AddVectors(v1,v2:TCoord):TCoord;
function PointsToVector(Origin,Destination:TCoord):TCoord;
function SubtractVectors(Initial,Subtract:TCoord):TCoord;

function MidPoint(v1,v2:TCoord):TCoord;
function BuildRotation(R1:TCoord;RotType:Integer):TRotMatrix;
function OldBuildRotation(R1:TCoord):TRotMatrix; //for old cdock compatibility only... remove?
function VectorSize(v:TCoord):Double;
function EqualVectors(v1,v2:TCoord):Boolean;
function RotateXToXY(v1:TCoord):TRotMatrix;
function TransformVector(Pivot:TCoord;Rotation:TRotMatrix;Place,Vector:TCoord):TCoord;
function RotAndPlace(Rotation:TRotMatrix;Place,Vector:TCoord):TCoord;overload;
procedure RotAndPlace(Rotation:TRotMatrix;Place:TCoord;var Coords:TCoords);overload;
function XRotation(Angle:TFloat):TRotMatrix;
function YRotation(Angle:TFloat):TRotMatrix;
function ZRotation(Angle:TFloat):TRotMatrix;
function DihedralAngle(c1, c2, c3, c4: TCoord):TFloat;
function RotationAngle(Axis,Pivot,P1,P2:TCoord):TFloat;
function AxisRotationMatrix(Axis:TCoord;Angle:TFloat):TRotMatrix;
function AngularDifference(A1,A2:TFloat):TFloat;
function AngleIsBetween(Low,High,Angle:TFloat):Boolean;
function FixAngle(Angle:TFloat):TFloat;
function Foot(LinePoint,LineVector,Point:TCoord):TCoord;
function TopCorner(c1,c2:TCoord):TCoord;
function BottomCorner(c1,c2:TCoord):TCoord;
function SmallestAxis(c:TCoord):Integer;
function LongestCoord(c:TCoord):TFloat;
   }
implementation

uses Math;

function Quaternion(r, i, j, k: TFloat): TQuaternion;
begin
  Result[0]:=r;
  Result[1]:=i;
  Result[2]:=j;
  Result[3]:=k;
end;


function Add(v1:TCoord;r:TFloat):TCoord;

begin
  Result[0]:=v1[0]+r;
  Result[1]:=v1[1]+r;
  Result[2]:=v1[2]+r;
end;

function Add(v1,v2:TCoord):TCoord;

begin
  Result[0]:=v1[0]+v2[0];
  Result[1]:=v1[1]+v2[1];
  Result[2]:=v1[2]+v2[2];
end;


function Add(vec:TCoord;Coords:TCoords):TCoords;overload;

var f:Integer;

begin
  SetLength(Result,Length(Coords));
  for f:=0 to High(Result) do
    Result[f]:=Add(vec,Coords[f]);
end;

function Add(Floats: TFloats; Float: TFloat): TFloats;

var f:Integer;

begin
  SetLength(Result,Length(Floats));
  for f:=0 to High(Result) do
    Result[f]:=Floats[f]+Float;
end;

function Add(Coords: TCoords; Floats: TFloats): TCoords;

var f:Integer;

begin
  SetLength(Result,Length(Floats));
  for f:=0 to High(Result) do
    Result[f]:=Add(Coords[f],Floats[f]);
end;

function Subtract(v1:TCoord;r:TFloat):TCoord;

begin
  Result[0]:=v1[0]-r;
  Result[1]:=v1[1]-r;
  Result[2]:=v1[2]-r;
end;

function Subtract(v1,v2:TCoord):TCoord;

begin
  Result[0]:=v1[0]-v2[0];
  Result[1]:=v1[1]-v2[1];
  Result[2]:=v1[2]-v2[2];
end;


function Subtract(Coords:TCoords;vec:TCoord):TCoords;

var f:Integer;

begin
  SetLength(Result,Length(Coords));
  for f:=0 to High(Result) do
    Result[f]:=Subtract(Coords[f],vec);
end;

function Multiply(Vec: TCoord; R: TFloat): TCoord;
begin
  Result[0]:=Vec[0]*R;
  Result[1]:=Vec[1]*R;
  Result[2]:=Vec[2]*R;
end;

function Multiply(const Vals: TFloats; Mult: TFloat): TFloats;

var f:Integer;

begin
  SetLength(Result,Length(Vals));
  for f:=0 to High(Vals) do
    Result[f]:=Vals[f]*Mult;
end;

function Multiply(const Cs: TCoords; Mult: TFloat): TCoords;

var f:Integer;

begin
  SetLength(Result,Length(Cs));
  for f:=0 to High(Cs) do
    Result[f]:=Multiply(Cs[f],Mult);
end;


function Multiply(Q1, Q2: TQuaternion): TQuaternion;
begin
  Result[0]:= Q1[0]*Q2[0] - Q1[1]*Q2[1] - Q1[2]*Q2[2] - Q1[3]*Q2[3];
  Result[1]:= Q1[0]*Q2[1] + Q1[1]*Q2[0] + Q1[2]*Q2[3] - Q1[3]*Q2[2];
  Result[2]:= Q1[0]*Q2[2] - Q1[1]*Q2[3] + Q1[2]*Q2[0] + Q1[3]*Q2[1];
  Result[3]:= Q1[0]*Q2[3] + Q1[1]*Q2[2] - Q1[2]*Q2[1] + Q1[3]*Q2[0];
end;

function Multiply(Mat1, Mat2: TRotMatrix): TRotMatrix;

var f,g:integer;

begin
  for f:=0 to 2 do
    for g:=0 to 2 do
      Result[f,g]:=Mat1[f,0]*Mat2[0,g]+Mat1[f,1]*Mat2[1,g]+Mat1[f,2]*Mat2[2,g];
end;

function DotProduct(v1, v2: TCoord): TFloat;
begin
  Result:=v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
end;

function CrossProduct(v1, v2: TCoord): TCoord;
begin
  Result[0]:=v1[1]*v2[2]-v1[2]*v2[1];
  Result[1]:=v1[2]*v2[0]-v1[0]*v2[2];
  Result[2]:=v1[0]*v2[1]-v1[1]*v2[0];
end;

function Norm(Vec:TCoord):TFloat;

begin
  Result:=Sqrt(Sqr(Vec[0])+Sqr(Vec[1])+Sqr(Vec[2]));
end;

function Conjugated(const Quat: TQuaternion): TQuaternion;
begin
  Result[0]:=Quat[0];
  Result[1]:=-Quat[1];
  Result[2]:=-Quat[2];
  Result[3]:=-Quat[3];
end;

procedure Normalize(var Vec: TCoord);

var n:TFloat;

begin
  n:=Sqrt(Sqr(Vec[0])+Sqr(Vec[1])+Sqr(Vec[2]));
  if n>TINY then
    begin
    Vec[0]:=Vec[0]/n;
    Vec[1]:=Vec[1]/n;
    Vec[2]:=Vec[2]/n;
    end
  else //Return X vector
    begin
    Vec[0]:=1;
    Vec[1]:=0;
    Vec[2]:=0;
    end;
end;

procedure Normalize(var Quat: TQuaternion);

var n:TFloat;

begin
  n:=Sqrt(Sqr(Quat[0])+Sqr(Quat[1])+Sqr(Quat[2])+Sqr(Quat[3]));
  if n>TINY then
    begin
    Quat[0]:=Quat[0]/n;
    Quat[1]:=Quat[1]/n;
    Quat[2]:=Quat[2]/n;
    Quat[3]:=Quat[3]/n;
    end
  else
    Quat:=IdentityQuaternion;
end;

function Scaled(const Vec: TCoord; const Scale: TFloat): TCoord;
begin
  Result[0]:=Vec[0]*Scale;
  Result[1]:=Vec[1]*Scale;
  Result[2]:=Vec[2]*Scale;
end;

function Simmetric(Vec: TCoord): TCoord;
begin
  Result[0]:=-Vec[0];
  Result[1]:=-Vec[1];
  Result[2]:=-Vec[2];
end;

function Simmetric(Coords: TCoords): TCoords;

var f:Integer;

begin
  SetLength(Result,Length(Coords));
  for f:=0 to High(Coords) do Result[f]:=Simmetric(Coords[f]);
end;

function Distance(Vec1, Vec2: TCoord): TFloat;
begin
  Result:=Sqrt(Sqr(Vec1[0]-Vec2[0])+Sqr(Vec1[1]-Vec2[1])+Sqr(Vec1[2]-Vec2[2]));
end;

function Distance(Vecs1, Vecs2: TCoords): TFloats;

var f:Integer;

begin
  SetLength(Result,Min(Length(Vecs1),Length(Vecs2)));
  for f:=0 to High(Result) do
    Result[f]:=Distance(Vecs1[f],Vecs2[f]);
end;

function Distance(Quart1, Quart2: TQuaternion): TFloat;
begin
  Result:=Sqrt(Sqr(Quart1[0]-Quart2[0])+Sqr(Quart1[1]-Quart2[1])+
               Sqr(Quart1[2]-Quart2[2])+Sqr(Quart1[3]-Quart2[3]));
end;

function Rotate(Vec: TCoord; RotMat: TRotMatrix): TCoord;
begin
  Result[0]:=RotMat[0,0]*Vec[0]+RotMat[0,1]*Vec[1]+RotMat[0,2]*Vec[2];
  Result[1]:=RotMat[1,0]*Vec[0]+RotMat[1,1]*Vec[1]+RotMat[1,2]*Vec[2];
  Result[2]:=RotMat[2,0]*Vec[0]+RotMat[2,1]*Vec[1]+RotMat[2,2]*Vec[2];
end;

function Rotate(Vecs: TCoords; RotMat: TRotMatrix): TCoords;

var f:Integer;

begin
  SetLength(Result,Length(Vecs));
  for f:=0 to High(Result) do
    Result[f]:=Rotate(Vecs[f],RotMat);
end;

function Rotate(Vec: TCoord; Quaternion: TQuaternion): TCoord;

var
  conjugate:TQuaternion;
  tmp:TQuaternion;

begin
  conjugate:=Conjugated(Quaternion);
  tmp[0]:=0;
  tmp[1]:=Vec[0];
  tmp[2]:=Vec[1];
  tmp[3]:=Vec[2];
  tmp:=Multiply(Quaternion,tmp);
  tmp:=Multiply(tmp,conjugate);
  Result[0]:=tmp[1];
  Result[1]:=tmp[2];
  Result[2]:=tmp[3];
end;

function Rotate(Vecs: TCoords; Quaternion: TQuaternion): TCoords;

var
  conjugate:TQuaternion;
  tmp:TQuaternion;
  f:Integer;

begin
  conjugate:=Conjugated(Quaternion);
  SetLength(Result,Length(Vecs));
    for f:=0 to High(Result) do
      begin
      tmp[0]:=0;
      tmp[1]:=Vecs[f,0];
      tmp[2]:=Vecs[f,1];
      tmp[3]:=Vecs[f,2];
      tmp:=Multiply(Quaternion,tmp);
      tmp:=Multiply(tmp,conjugate);
      Result[f,0]:=tmp[1];
      Result[f,1]:=tmp[2];
      Result[f,2]:=tmp[3];
      end;
end;


function RotationQuaternion(Axis: TCoord; Rotation: TFloat): TQuaternion;

var
  sint,cost:TFloat;

begin
  cost:=Cos(Rotation/2);
  sint:=Sin(Rotation/2);
  Result[0]:=cost;
  Result[1]:=sint*Axis[0];
  Result[2]:=sint*Axis[1];
  Result[3]:=sint*Axis[2];
  Normalize(Result);
end;

function RotationTo(const VFrom,VTo: TCoord): TQuaternion;
// (from http://www.gamedev.net/topic/429507-finding-the-quaternion-betwee-two-vectors/?p=3856228#entry3856228


var
  tmp:TCoord;
  dot,n1,n2:TFloat;

begin
  tmp:=CrossProduct(VFrom,VTo);
  n1:=Sqr(VFrom[0])+Sqr(VFrom[1])+Sqr(VFrom[2]);
  n2:=Sqr(VTo[0])+Sqr(VTo[1])+Sqr(VTo[2]);
  dot:=DotProduct(VFrom,VTo);
  Result[0]:= Sqrt(n1*n2) + dot;
  Result[1]:=tmp[0];
  Result[2]:=tmp[1];
  Result[3]:=tmp[2];
  Normalize(Result);
end;

function StaticRMSD(Coords1, Coords2: TCoords): TFloat;

var f:Integer;

begin
  Result:=0;
  for f:=0 to High(Coords1) do
    Result:=Result+Sqr(Coords1[f,0]-Coords2[f,0])
                  +Sqr(Coords1[f,1]-Coords2[f,1])
                  +Sqr(Coords1[f,2]-Coords2[f,2]);
  Result:=Sqrt(Result/Length(Coords1));
end;

function RotAndPlace(Rotation: TRotMatrix; Place, Vector: TCoord): TCoord;

begin
  Result:=Rotate(Vector,Rotation);
  Result:=Add(Result,Place);
end;

procedure RotAndPlace(Rotation:TRotMatrix;Place:TCoord;var Coords:TCoords);

var f:Integer;

begin
  for f:=0 to High(Coords) do
    Coords[f]:=RotAndPlace(Rotation,Place,Coords[f]);
end;

function BuildRotation(R1: TCoord; RotType: Integer): TRotMatrix;

begin
  case RotType of
    MCXYZRotation:Result:=Multiply(Multiply(ZRotation(R1[2]),YRotation(R1[1])),XRotation(R1[0]));
    MCRotationXYAxis:
      begin
      Result:=Multiply(YRotation(R1[1]),XRotation(R1[0]));
      Result:=Multiply(InvertBase(Result),Multiply(ZRotation(R1[2]),Result));
      end;
    else Result:=IdentityMatrix;
    end;
end;

function XRotation(Angle: TFloat): TRotMatrix;

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


function YRotation(Angle: TFloat): TRotMatrix;

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

function ZRotation(Angle: TFloat): TRotMatrix;

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

function InvertBase(M: TRotMatrix): TRotMatrix;

var f,g:integer;

begin
  for f:=0 to 2 do
    for g:=0 to 2 do
      Result[f,g]:=M[g,f];
end;


end.



