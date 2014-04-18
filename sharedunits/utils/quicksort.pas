{*******************************************************************************
This source file is public domain. It is provided as is, with no explicit
or implied warranties regarding anything whatsoever.
********************************************************************************
Author: Ludwig Krippahl
Date: 31.3.2011
Purpose:
  Sorting
Revisions:
TODO: Only TFloats sorting is implemented in QSAscendingIndex; the others
  convert the data to sort as float. This should be improved by implementing
  QS on the other value types (? bad code reuse but good for large arrays)
*******************************************************************************}
unit quicksort;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes;

type
  TQSInt64Vals=array of Int64;

function QSAscendingIndex(const Values:TFloats):TIntegers;overload;
function QSAscendingIndex(const Values:TQSInt64Vals):TIntegers;overload;
function QSAscendingIndex(const Values:TIntegers):TIntegers;overload;
function QSIntsAsVals(Ints:array of Integer):TFloats;
function QSIndexAsRanking(Ix:TIntegers):TIntegers;
function QSZeroBasedIndex(Count:Integer):TIntegers;
function QSSorted(const Values:TIntegers):TIntegers;overload;
function QSSorted(const Values:TFloats):TFloats;overload;


implementation

function QSIndexAsRanking(Ix:TIntegers):TIntegers;

var f:Integer;

begin
  SetLength(Result,Length(Ix));
  for f:=0 to High(Ix) do
    Result[Ix[f]]:=f;
end;

function QSIntsAsVals(Ints:array of Integer):TFloats;

var f:Integer;

begin
  SetLength(Result,Length(Ints));
  for f:=0 to High(Ints) do
    Result[f]:=Ints[f]
end;

function QSAscendingIndex(const Values:TFloats):TIntegers;

var
  t:Integer;
  vpiv:Double;

procedure Quick(ini,fim:Integer);

var
  le,ri:Integer;

begin
  vpiv:=Values[Result[ini]];
  le:=ini;
  ri:=fim;
  repeat
    while Values[Result[le]]<vpiv do
      Inc(le);
    while Values[Result[ri]]>vpiv do
      Dec(ri);
    if (ri>=le) then
           begin
           t:=Result[ri];
           Result[ri]:=Result[le];
           Result[le]:=t;
           Inc(le);
           Dec(ri);
           end;
  until (ri<le);
  if ri>ini  then
     quick(ini,ri);
  if fim>le then
     quick(le,fim);
end;

begin
  SetLength(Result,Length(Values));
  for t:=0 to High(Result) do Result[t]:=t;
  if (High(Values)>0) then
    Quick(0,High(Values));
end;

function QSAscendingIndex(const Values:TQSInt64Vals):TIntegers;overload;

var
  t:Integer;
  vpiv:Double;

procedure Quick(ini,fim:Integer);

var
  le,ri:Integer;

begin
  vpiv:=Values[Result[ini]];
  le:=ini;
  ri:=fim;
  repeat
    while Values[Result[le]]<vpiv do
      Inc(le);
    while Values[Result[ri]]>vpiv do
      Dec(ri);
    if (ri>=le) then
           begin
           t:=Result[ri];
           Result[ri]:=Result[le];
           Result[le]:=t;
           Inc(le);
           Dec(ri);
           end;
  until (ri<le);
  if ri>ini  then
     quick(ini,ri);
  if fim>le then
     quick(le,fim);
end;

begin
  SetLength(Result,Length(Values));
  for t:=0 to High(Result) do Result[t]:=t;
  if (High(Values)>0) then
    Quick(0,High(Values));
end;

function QSZeroBasedIndex(Count:Integer):TIntegers;

var f:Integer;

begin
  SetLength(Result,Count);
  for f:=0 to Count-1 do
    Result[f]:=f;
end;

function QSSorted(const Values: TIntegers): TIntegers;

var
  ixs:TIntegers;
  f:Integer;

begin
  ixs:=QSAscendingIndex(Values);
  SetLength(Result,Length(ixs));
  for f:=0 to High(Result) do
    Result[f]:=Values[ixs[f]];
end;

function QSSorted(const Values: TFloats): TFloats;

var
  ixs:TIntegers;
  f:Integer;

begin
  ixs:=QSAscendingIndex(Values);
  SetLength(Result,Length(ixs));
  for f:=0 to High(Result) do
    Result[f]:=Values[ixs[f]];
end;

function QSAscendingIndex(const Values:TIntegers):TIntegers;overload;

var
  tmpvals:TFloats;
  f:Integer;

begin
  SetLength(tmpvals,Length(Values));
  for f:=0 to High(tmpvals) do
    tmpvals[f]:=Values[f];
  Result:=QSAscendingIndex(tmpvals);
end;

end.

