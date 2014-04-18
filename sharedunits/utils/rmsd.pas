{*******************************************************************************
This source file is public domain. It is provided as is, with no explicit
or implied warranties regarding anything whatsoever.
********************************************************************************
Author: Ludwig Krippahl
Date: 2014 02 24
Purpose:
  Root mean square deviation computations
Revisions:
TODO:
  This was adapted from old Chemera code on Delphi. Probably needs a revision
  and some rewriting

  Also: see better algorithm in Coutsias, Seok, Dill,
    Using quaternions to calculate RMSD. J Comput Chem 2004
    http://dx.doi.org/10.1002/jcc.20110

*******************************************************************************}
unit rmsd;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, geomutils, cgradient;


type

  { TRMSDCalculator }

  TRMSDCalculator=class
  protected
    FFixed,FMobile,FPositions,FDerivs: TCoords;
    FRotation,FTranslation,FCenterVec:TCoord;
    FBaseRotation:TRotMatrix;
    FCGMinimiser:TCGMinimiser;
    FRmsd:TFloat;
    function EvaluateAndDerive(var Point:TFloats;var Deriv:TFloats):TFloat;dynamic;
    function Evaluate(const Point:TFloats):TFloat;dynamic;
    procedure GetDerivatives(var Deriv:TFloats);
    function CurrentRotation:TRotMatrix;
  public
    procedure Initialize;
    constructor Create;
    property CenterVec:TCoord read FCenterVec;
    property PlacedTranslation:TCoord read FTranslation;
    property PlacedRotation:TRotMatrix read FBaseRotation;
    procedure AddCoordinates(Fixed,Mobile:TCoord);overload;
    procedure AddCoordinates(Fixed,Mobile:TCoords);overload;
    property Rmsd:TFloat read FRmsd;
    function Minimise(MaxIterations:Integer;MinVariation:TFloat):TFloat;
    procedure SetCoords(Fixed,Mobile:TCoords);
    procedure BuildPositions;
    function GetCoord(Index:Integer):TCoord;
    function TransformedCoords(const Coords:TCoords):TCoords;
  end;

implementation

{ TRMSDCalculator }

procedure TRMSDCalculator.AddCoordinates(Fixed, Mobile: TCoord);
begin
  SetLength(FFixed,Length(FFixed)+1);
  SetLength(FMobile,Length(FFixed));
  FFixed[High(FFixed)]:=Fixed;
  FMobile[High(FMobile)]:=Mobile;
end;

procedure TRMSDCalculator.AddCoordinates(Fixed, Mobile: TCoords);

var
  f,len:Integer;

begin
  len:=Length(FFixed);
  SetLength(FFixed,Length(Fixed)+len);
  SetLength(FMobile,Length(Mobile)+len);
  for f:=0 to High(Fixed) do
    begin
    FFixed[f+len]:=Fixed[f];
    FMobile[f+len]:=Mobile[f];
    end;

end;

procedure TRMSDCalculator.BuildPositions;

var
  f:Integer;
  RotMatrix:TRotMatrix;

begin
  SetLength(FPositions,Length(FMobile));
  RotMatrix:=CurrentRotation;
  for f:=0 to High(FPositions) do
    FPositions[f]:=RotAndPlace(RotMatrix,FTranslation,FMobile[f]);
end;

constructor TRMSDCalculator.Create;
begin
  inherited Create;
  FBaseRotation:=BuildRotation(Coord(0,0,0),MCXYZRotation);
  FCGMinimiser:=TCGMinimiser.Create;
  FCGMinimiser.GetValueFunction:=@Evaluate;
  FCGMinimiser.ValueAndDerivativeFunction:=@EvaluateAndDerive;
end;

function TRMSDCalculator.CurrentRotation:TRotMatrix;
begin
  Result:=Multiply(BuildRotation(FRotation,MCXYZRotation),FBaseRotation);
end;

function TRMSDCalculator.Evaluate(const Point: TFloats): TFloat;

var
  f:Integer;
  RotMatrix:TRotMatrix;
  c:TCoord;

begin
  Result:=0;
  for f:=0 to 2 do
    begin
    FRotation[f]:=Point[f];
    FTranslation[f]:=Point[f+2];
    end;
  RotMatrix:=CurrentRotation;
  for f:=0 to High(FMobile) do
    begin
    c:=RotAndPlace(RotMatrix,FTranslation,FMobile[f]);
    Result:=Result+Sqr(c[0]-FFixed[f,0])+Sqr(c[1]-FFixed[f,1])+
      Sqr(c[2]-FFixed[f,2]);
    FPositions[f]:=c;
    end;
end;

function TRMSDCalculator.EvaluateAndDerive(var Point: TFloats;
  var Deriv: TFloats): TFloat;
var
  f:Integer;
  RotMatrix:TRotMatrix;
  Coords:TCoord;
  Dist:TFloat;

begin
  Result:=0;
  for f:=0 to 2 do
    begin
    FRotation[f]:=Point[f];
    Point[f]:=0;
    FTranslation[f]:=Point[f+3];
    end;
  RotMatrix:=CurrentRotation;
  for f:=0 to High(FMobile) do
    begin
    Coords:=RotAndPlace(RotMatrix,FTranslation,FMobile[f]);
    Dist:=Sqr(Coords[0]-FFixed[f,0])+
          Sqr(Coords[1]-FFixed[f,1])+
          Sqr(Coords[2]-FFixed[f,2]);
    Result:=Result+Dist;
    FDerivs[f]:=Subtract(FFixed[f],Coords);
    FPositions[f]:=Coords;
    end;
  GetDerivatives(Deriv);


end;

function TRMSDCalculator.GetCoord(Index: Integer): TCoord;
begin
  Result:=FPositions[Index];
end;

function TRMSDCalculator.TransformedCoords(const Coords: TCoords): TCoords;

//adapted from Chemera 2.0.

var f:Integer;

begin
  SetLength(Result,Length(Coords));
  for f:=0 to High(Result) do
    begin
    Result[f]:=Subtract(Coords[f],CenterVec);
    Result[f]:=RotAndPlace(PlacedRotation,PlacedTranslation,Result[f])
    end;
end;

procedure TRMSDCalculator.GetDerivatives(var Deriv: TFloats);

const
  XVector:TCoord=(1,0,0);
  YVector:TCoord=(0,1,0);
  ZVector:TCoord=(0,0,1);

var
  v1,AtDeriv,Pivot,Axis:TCoord;
  f:Integer;
  TDeriv,RDeriv:TCoord;

function Torsion:TFloat;

var f:Integer;

begin
  Result:=0;
  for f:=0 to High(FDerivs) do
    begin
    AtDeriv:=FDerivs[f];
    if (AtDeriv[0]<>0) or (AtDeriv[1]<>0) or (AtDeriv[2]<>0) then
      begin
      v1:=Subtract(FPositions[f],Pivot);
      v1:=CrossProduct(Axis,v1);
      Result:=Result+DotProduct(v1,AtDeriv);
      end;
    end;
end;

function Shear:TCoord;

var f:Integer;

begin
  Result:=NullVector;
  for f:=0 to High(FFixed) do
    Result:=Add(FDerivs[f],Result);
end;

begin
  FBaseRotation:=CurrentRotation;
  FRotation:=NullVector;
  TDeriv:=Shear;
  Pivot:=FTranslation;
  Axis:=XVector;
  RDeriv[0]:=Torsion;
  Axis:=YVector;
  RDeriv[1]:=Torsion;
  Axis:=ZVector;
  RDeriv[2]:=Torsion;
  for f:=0 to 2 do
    begin
    Deriv[f]:=RDeriv[f];
    Deriv[f+3]:=TDeriv[f];
    end;
end;

procedure TRMSDCalculator.Initialize;
begin
  FMobile:=nil;
  FFixed:=nil;
  FDerivs:=nil;
end;

function TRMSDCalculator.Minimise(MaxIterations: Integer;
  MinVariation: TFloat): TFloat;

var
  InitialPoint,InitialDeriv:TFloats;
  f:Integer;

procedure FirstTranslation;

var f:Integer;
    FF:TCoord;
begin
  FF:=NullVector;
  for f:=0 to High(FFixed) do
    FF:=Add(FF,FFixed[f]);
  if FMobile<>nil then
    begin
    FF:=Multiply(FF,1/Length(FFixed));
    for f:=0 to 2 do InitialPoint[f+3]:=FF[f];
    end;
end;

procedure Center;

var f:Integer;

begin
  FCenterVec:=NullVector;
  for f:=0 to High(FMobile) do
    FCenterVec:=Add(FCenterVec,FMobile[f]);
  if FMobile<>nil then FCenterVec:=Multiply(FCenterVec,1/Length(FMobile));
  for f:=0 to High(FMobile) do
    FMobile[f]:=Subtract(FMobile[f],FCenterVec);
end;

begin
  Center;
  SetLength(FDerivs,Length(FFixed));
  SetLength(FPositions,Length(FFixed));
  SetLength(InitialPoint,6);
  SetLength(InitialDeriv,6);
  for f:=0 to 5 do InitialPoint[f]:=0;
  FirstTranslation;
  EvaluateAndDerive(InitialPoint,InitialDeriv);
  Result:=FCGMinimiser.Minimize(MaxIterations,MinVariation,InitialPoint,InitialDeriv);
  if FMobile<>nil then
    Result:=Sqrt(Result/Length(FMobile));
  FRmsd:=Result;
  FDerivs:=nil;
  FPositions:=nil;
  InitialPoint:=nil;
  InitialDeriv:=nil;
  FBaseRotation:=CurrentRotation;
  FRotation:=NullVector;
end;

procedure TRMSDCalculator.SetCoords(Fixed, Mobile: TCoords);
begin
  FFixed:=Copy(Fixed,0,Length(Fixed));
  FMobile:=Copy(Mobile,0,Length(Mobile));
end;

end.

