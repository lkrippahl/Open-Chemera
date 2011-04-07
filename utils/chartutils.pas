{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 9.1.2011
Purpose:
 Types and functions for charts, like histograms.
 No user interface, but may create graphical plots in Bmps.
Requirements:

Revisions:
To do:
*******************************************************************************}
unit chartutils;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, graphics;

function Histogram(data:TOCFloats;lower,upper:TOCFloats):TOCIntegers;overload;
  // build histogram with bin i between lower[i] and upper[i], which must be of same length.
function Histogram(data:TOCMatrix;lower,upper:TOCFloats):TOCIntegers;overload;
  // same as above but for all rows of a matrix
procedure CreateBins(min,max:TOCFloat;bins:Integer;var lower,upper:TOCFloats);
  // create a set of regular bins. NOTE: last bin is expanded by 3E-7, roughly twice the
  // epsilon for singles. TODO: Improve this
function CountBetween(vals:TOCFloats;lo,hi:TOCFloat):Integer;
function CountBetween(vals:TOCMatrix;lo,hi:TOCFloat):Integer;
procedure PlotXY(xs,ys:TOCFloats;x1,x2,y1,y2:TOCFloat;bmp:TBitmap;r:TRect;col:TColor);



implementation

function Histogram(data: TOCFloats; lower, upper: TOCFloats): TOCIntegers;

var
  f,g:Integer;
  d:TOCFloat;

begin
  Assert(Length(lower)=Length(upper),'Lengths of bin bounds do not match');
  SetLength(Result,Length(lower));
  for f:=0 to High(Result) do Result[f]:=0;
  for f:=0 to High(data) do
    begin
    d:=data[f];
    for g:=0 to High(lower) do
      begin
      if (d>=lower[g]) and (d<upper[g]) then
        begin
        Inc(Result[g]);
        Break;
        end;
      end;
    end;
end;

function Histogram(data: TOCMatrix; lower, upper: TOCFloats): TOCIntegers;

var
  f,g:Integer;
  tmp:TOCIntegers;

begin
  Assert(Length(lower)=Length(upper),'Lengths of bin bounds do not match');
  SetLength(Result,Length(lower));
  for f:=0 to High(Result) do Result[f]:=0;
  for f:=0 to High(data) do
    begin
    tmp:=Histogram(data[f],lower,upper);
    for g:=0 to High(Result) do
      Result[g]:=Result[g]+tmp[g];
    end;
end;

procedure CreateBins(min, max: TOCFloat; bins: Integer; var lower,
  upper: TOCFloats);

var
  step,curr:TOCFLoat;
  f:Integer;

begin
  step:=(max-min)/(bins+1);
  SetLength(upper,bins);
  SetLength(lower,bins);
  curr:=min;
  lower[0]:=min;
  for f:=1 to High(lower) do
    begin
    curr:=curr+step;
    lower[f]:=curr;
    upper[f-1]:=curr;
    end;
  curr:=curr+step+step*1e-7;
  upper[High(upper)]:=curr;
end;

function CountBetween(vals: TOCFloats; lo, hi: TOCFloat): Integer;

var f:Integer;

begin
  Result:=0;
  for f:=0 to High(vals) do
    if (vals[f]<=hi) and (vals[f]>=lo) then Inc(Result);
end;

function CountBetween(vals: TOCMatrix; lo, hi: TOCFloat): Integer;

var f:Integer;

begin
  Result:=0;
  for f:=0 to High(vals) do
    Result:=Result+CountBetween(vals[f],lo,hi);
end;

procedure PlotXY(xs, ys: TOCFloats; x1, x2, y1, y2: TOCFloat; bmp: TBitmap;r:TRect;col:TColor);

var
  f,x,y:Integer;
  lx,ly:TOCFloat;

begin
  lx:=(r.Right-r.Left)/(x2-x1);
  ly:=(r.Bottom-r.Top)/(y2-y1);
  Assert(Length(xs)=Length(ys),'Vectors do not match');
  for f:=0 to High(xs) do
    begin
    x:=Round((xs[f]-x1)*lx)+r.Left;
    y:=r.Bottom-Round((ys[f]-x1)*ly);
    if (x>=r.Left) and (x<=r.Right)
    and (y>=r.Top) and (y<=r.Bottom) then
      bmp.Canvas.Pixels[x,y]:=col;
    end;
end;

end.

