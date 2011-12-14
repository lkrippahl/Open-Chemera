{*******************************************************************************
This source file is public domain. It is provided as is, with no explicit
or implied warranties regarding anything whatsoever.
********************************************************************************
Author: Ludwig Krippahl
Date: 3.12.2011
Purpose:
  Statistical abd combinatorial functions.
Requirements:
Revisions:
To do:

*******************************************************************************}
unit statistics;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, quicksort;

function Mean(V:TFloats):TFloat;
function SVariance(V:TFLoats):TFloat;//Unbiased estimator based on sample
function Variance(V:TFLoats):TFloat;
function SCovariance(V1,V2:TFloats):TFloat;//Unbiased estimator based on sample
function Covariance(V1,V2:TFloats):TFloat;
function Pearson(V1,V2:TFLoats):TFloat;
procedure NextCombination(Total:Integer;var Current:TIntegers);
          //Current is a 0 based index of elements, out of Total, and is
          //incremented to find the next combination of Current elements out of
          //Total.
function TheilSenEstimator(Xs,Ys:TFloats):TFloat;
          //Median of the slopes of all pairs of points
          //See http://en.wikipedia.org/wiki/Theil–Sen_estimator

implementation

function Mean(V:TFloats):TFloat;

begin
  if V<>nil then
    Result:=Sum(V)/Length(V)
  else Result:=0;
end;

function Variance(V:TFLoats):TFloat;

var
  m:TFloat;
  f:Integer;

begin
  Assert(Length(V)>1,'Insufficient values');
  m:=Mean(V);
  Result:=0;
  for f:=0 to High(V) do
    Result:=Result+Sqr(V[f]-m);
  Result:=Result/Length(V);
end;

function SVariance(V:TFLoats):TFloat;

var
  m:TFloat;
  f:Integer;

begin
  Assert(Length(V)>1,'Insufficient values');
  m:=Mean(V);
  Result:=0;
  for f:=0 to High(V) do
    Result:=Result+Sqr(V[f]-m);
  Result:=Result/(Length(V)-1);
end;

function Covariance(V1,V2:TFloats):TFloat;

var
  m1,m2:TFloat;
  f:Integer;

begin
  Assert(Length(V1)=Length(V2),'Mismatching vectors');
  Assert(Length(V1)>1,'Insufficient values');
  m1:=Mean(V1);
  m2:=Mean(V2);
  Result:=0;
  for f:=0 to High(V1) do
    Result:=Result+(V1[f]-m1)*(V2[f]-m2);
  Result:=Result/Length(V1);
end;


function SCovariance(V1,V2:TFloats):TFloat;

var
  m1,m2:TFloat;
  f:Integer;

begin
  Assert(Length(V1)=Length(V2),'Mismatching vectors');
  Assert(Length(V1)>1,'Insufficient values');
  m1:=Mean(V1);
  m2:=Mean(V2);
  Result:=0;
  for f:=0 to High(V1) do
    Result:=Result+(V1[f]-m1)*(V2[f]-m2);
  Result:=Result/(Length(V1)-1);
end;

function Pearson(V1,V2:TFLoats):TFloat;

begin
  Result:=Sqrt(Variance(V1)*Variance(V2));
  if Abs(Result)<TinyFloat then
    Result:=0
  else
    Result:=Covariance(V1,V2)/Result;
end;

procedure NextCombination(Total:Integer;var Current:TIntegers);

var
  f,pivot,top:Integer;

begin
  //find where to change
  top:=Total-1; //0 is first element index
  pivot:=-1;
  for f:=High(Current) downto 0 do
    if Current[f]<top then
      begin
      pivot:=f;
      Break;
      end
    else
      Dec(top);

  //update from pivot onwards
  if pivot>=0 then
    begin
    Inc(Current[pivot]);
    top:=Current[pivot]+1;
    for f:=pivot+1 to High(Current) do
      begin
      Current[f]:=top;
      Inc(top);
      end;
    end
  else Current:=nil;
end;

function TheilSenEstimator(Xs,Ys:TFloats):TFloat;

//TO DO: «Once the slope m has been determined, one may determine a line through
//the sample points by setting the y-intercept b to be the median of the values yi − mxi»
//http://en.wikipedia.org/wiki/Theil–Sen_estimator

var
  f,g:Integer;
  x1,y1,x2:TFloat;
  sorted:TIntegers;
  slopes:TFLoats;
  slopeix:Integer;

begin
  Assert(Length(Xs)=Length(Ys),'Mismatching Coordinate Vectors');
  SetLength(slopes,Length(Xs)*Length(Ys));
  slopeix:=0;
  for f:=0 to High(Xs)-1 do
    begin
    x1:=Xs[f];
    y1:=Ys[f];
    for g:=f+1 to High(Ys) do
      begin
      x2:=Xs[g];
      if Abs(x1-x2)>TinyFloat then
        begin
        slopes[slopeix]:=(Ys[g]-y1)/(x2-x1);
        Inc(slopeix);
        end;
      end;
    end;
  SetLength(slopes,slopeix);
  sorted:=QSAscendingIndex(slopes);
  Result:=slopes[sorted[Length(sorted) div 2]];
end;

end.
