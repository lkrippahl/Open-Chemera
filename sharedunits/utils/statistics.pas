{*******************************************************************************
This source file is public domain. It is provided as is, with no explicit
or implied warranties regarding anything whatsoever.
********************************************************************************
Author: Ludwig Krippahl
Date: 3.12.2011
Purpose:
  Statistical functions.
Requirements:
Revisions:
To do:

*******************************************************************************}
unit statistics;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes;

function Mean(V:TFloats):TFloat;
function SVariance(V:TFLoats):TFloat;//Unbiased estimator based on sample
function Variance(V:TFLoats):TFloat;
function SCovariance(V1,V2:TFloats):TFloat;//Unbiased estimator based on sample
function Covariance(V1,V2:TFloats):TFloat;
function Pearson(V1,V2:TFLoats):TFloat;

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

end.

