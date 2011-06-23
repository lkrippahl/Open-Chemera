{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 31.3.2011
Purpose:
  Sorting
Revisions:
*******************************************************************************}
unit quicksort;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes;

type
  TQSInt64Vals=array of Int64;

function QSAscendingIndex(const Values:TOCFloats):TOCIntegers;overload;
function QSAscendingIndex(const Values:TQSInt64Vals):TOCIntegers;overload;
function QSIntsAsVals(Ints:array of Integer):TOCFloats;
function QSIndexAsRanking(Ix:TOCIntegers):TOCIntegers;
function QSZeroBasedIndex(Count:Integer):TOCIntegers;

implementation

Uses forms;

function QSIndexAsRanking(Ix:TOCIntegers):TOCIntegers;

var f:Integer;

begin
  SetLength(Result,Length(Ix));
  for f:=0 to High(Ix) do
    Result[Ix[f]]:=f;
end;

function QSIntsAsVals(Ints:array of Integer):TOCFloats;

var f:Integer;

begin
  SetLength(Result,Length(Ints));
  for f:=0 to High(Ints) do
    Result[f]:=Ints[f]
end;

function QSAscendingIndex(const Values:TOCFloats):TOCIntegers;

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

function QSAscendingIndex(const Values:TQSInt64Vals):TOCIntegers;overload;

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

function QSZeroBasedIndex(Count:Integer):TOCIntegers;

var f:Integer;

begin
  SetLength(Result,Count);
  for f:=0 to Count-1 do
    Result[f]:=f;
end;

end.

