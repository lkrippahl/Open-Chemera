{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain. Enjoy.
********************************************************************************
Author: Ludwig Krippahl
Date: 7-12-2011
Purpose:
  Calculates measures of coevolution, phylogenetics, ...
Requirements:
Revisions:
To do:
*******************************************************************************}

unit coevolution;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, statistics;

type

  //Groups organisms by the correlation of two protein distance matrixes

  { TTwoProteinCorrelation }

  TTwoProteinCorrelation=class
    //Note:all correlations are squared
  protected
    FMatrix1,FMatrix2:TMatrix;
    function WithGroupCorrelation(Indexes:TIntegers;OrganismIx:Integer):TFloat;
    function AllGroupCorrelations(Indexes:TIntegers):TFloats;
    function GetBestGroup(MandatoryIndex,Size:Integer):TIntegers;
      //not implemented. Should
  public
    constructor Create(AMatrix1,AMatrix2:TMatrix);
    procedure GroupByCorrelation(MandatoryIndex:Integer; out Indexes:TIntegers;
                                out Correlations:TFloats);
    procedure TheilSenToQuery(QueryIx:Integer;out Slope:TFloat;out Distances:TFloats);
      //calculates a Theil-Sen regression, outputs the slope and the distances of
      //all points in the QueryIx column to the Theil-Sen line
  end;

implementation

{ TTwoProteinCorrelation }

function TTwoProteinCorrelation.WithGroupCorrelation(Indexes: TIntegers;
  OrganismIx: Integer): TFloat;

var
  f,g:Integer;
  vec1,vec2:TFloats;

begin
  SetLength(vec1,Length(Indexes));
  SetLength(vec2,Length(Indexes));
  for f:=0 to High(Indexes) do
    begin
    vec1[f]:=FMatrix1[OrganismIx,Indexes[f]];
    vec2[f]:=FMatrix1[OrganismIx,Indexes[f]];
    end;
  Result:=Sqr(Pearson(vec1,vec2));
end;

function TTwoProteinCorrelation.AllGroupCorrelations(Indexes: TIntegers
  ): TFloats;

var f:Integer;

begin
  SetLength(Result,Length(Indexes));
  for f:=0 to High(Result) do Result[f]:=WithGroupCorrelation(Indexes,Indexes[f]);
end;


function TTwoProteinCorrelation.GetBestGroup(MandatoryIndex, Size: Integer
  ): TIntegers;

var
  currentcombination,currentgroup:TIntegers;
  bestcorrelation:TFloat;


procedure SetCurrentGroup;

var
  f:Integer;

begin
  for f:=0 to High(currentcombination) do
      if currentcombination[f]>=MandatoryIndex then
        currentgroup[f+1]:=currentcombination[f]+1
      else currentgroup[f+1]:=currentcombination[f];
end;

var
  f:Integer;

begin
  SetLength(currentgroup,Size);
  currentgroup[0]:=MandatoryIndex;
  SetLength(currentcombination,Size-1);
  for f:=0 to High(currentcombination) do currentcombination[f]:=f;
  bestcorrelation:=-1;
  Result:=nil;
  while currentcombination<>nil do
    begin
    SetCurrentGroup;
    //TODO: evaluate group
    //not finished yet!

    NextCombination(Length(FMatrix1)-1,currentcombination);
    end;
end;

constructor TTwoProteinCorrelation.Create(AMatrix1, AMatrix2: TMatrix);
begin
  inherited Create;
  FMatrix1:=AMatrix1;
  FMatrix2:=AMatrix2;
end;

procedure TTwoProteinCorrelation.GroupByCorrelation(MandatoryIndex: Integer;
  out Indexes: TIntegers; out Correlations: TFloats);
begin
  //TODO: implement
end;

procedure TTwoProteinCorrelation.TheilSenToQuery(QueryIx: Integer; out
  Slope: TFloat; out Distances: TFloats);

var
  xs,ys:TFloats;
  f:Integer;
  b,divisor:TFloat;

begin
  xs:=Copy(FMatrix1[QueryIx],0,Length(FMatrix1));
  ys:=Copy(FMatrix2[QueryIx],0,Length(FMatrix2));
  Slope:=TheilSenEstimator(Xs,Ys);
  //Formula: distance of P(m,n) to Ax+By+C=0 is
  //d=Abs(Am+Bn+C)/Sqrt(A^2+B^2)
  b:=-Slope; //TODO check if slope not 0
  divisor:=Sqrt(1+Sqr(b));
  SetLength(Distances,Length(xs));
  for f:=0 to High(Distances) do
    Distances[f]:=Abs(b*xs[f]+ys[f])/divisor;
end;

end.

