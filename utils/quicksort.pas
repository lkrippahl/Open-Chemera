unit quicksort;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils; 

type
  TQSIndex=array of Integer;
  TQSValues=array of Double;
  TQSInt64Vals=array of Int64;

function QSAscendingIndex(const Values:TQSValues):TQSIndex;overload;
function QSAscendingIndex(const Values:TQSInt64Vals):TQSIndex;overload;
procedure AddToIndex(var Index:TQSIndex;Value:Integer);
procedure AddToValues(var Values:TQSValues; Value:Double);
function QSIntsAsVals(Ints:array of Integer):TQSValues;
function QSIndexAsRanking(Ix:TQSIndex):TQSIndex;

implementation

Uses forms;

function QSIndexAsRanking(Ix:TQSIndex):TQSIndex;

var f:Integer;

begin
  SetLength(Result,Length(Ix));
  for f:=0 to High(Ix) do
    Result[Ix[f]]:=f;
end;

function QSIntsAsVals(Ints:array of Integer):TQSValues;

var f:Integer;

begin
  SetLength(Result,Length(Ints));
  for f:=0 to High(Ints) do
    Result[f]:=Ints[f]
end;

procedure AddToIndex(var Index:TQSIndex;Value:Integer);

begin
  SetLength(Index,Length(Index)+1);
  Index[High(Index)]:=Value;
end;

procedure AddToValues(var Values:TQSValues; Value:Double);

begin
  SetLength(Values,Length(Values)+1);
  Values[High(Values)]:=Value;
end;


function QSAscendingIndex(const Values:TQSValues):TQSIndex;

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

function QSAscendingIndex(const Values:TQSInt64Vals):TQSIndex;overload;

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


end.

