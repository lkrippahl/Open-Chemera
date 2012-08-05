{*******************************************************************************
This source file is public domain. It is provided as is, with no explicit
or implied warranties regarding anything whatsoever.
********************************************************************************
Author: Ludwig Krippahl
Date: 9.1.2011
Purpose:
  Base data types (string arrays, coords, etc) and utility functions.
Requirements:
Revisions:
To do:
  Fix caps on arguments
*******************************************************************************}

unit basetypes;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils;

type
  TSimpleStrings = array of string;
  {$IFDEF SINGLEPRECISION}
  TFloat = Single;
  {$ELSE}
  TFloat = Double; //This is the default
  {$ENDIF}
  TFloats = array of TFloat;
  TDoubles = array of Double;
  TSingles = array of Single;
  TMatrix= array of TFloats;
  TCoord = array [0..2] of TFloat; //3D coords ou 3D vector
  TCoords = array of TCoord;
  TIntegers = array of Integer;
  TCardinals = array of Cardinal;
  TBooleans = array of Boolean;

  procedure AddToArray(i:Integer; var a:TIntegers); overload;
  procedure AddToArray(s:string; var a:TSimpleStrings); overload;
  procedure AddToArray(c:TCoord; var a:TCoords); overload;
  procedure AddToArray(f:TFLoat; var a:TFloats); overload;
  procedure AddToArray(c:Cardinal; var a:TCardinals); overload;
  procedure AddToArray(b:Boolean; var a:TBooleans); overload;

  //Warning: removing does not preserve ordering (exchanges with last element)
  procedure RemoveFromArray(Ix:Integer;var A:TIntegers);overload;
  procedure RemoveFromArray(Ix:Integer;var A:TCoords);overload;
  procedure RemoveFromArray(Ix:Integer;var A:TFloats);overload;
  procedure RemoveFromArray(Ix:Integer;var A:TCardinals);overload;
  procedure RemoveFromArray(Ix:Integer;var A:TSimpleStrings);overload;

  //Return last index, or -1 if not found
  function IndexOf(const i:Integer; const a:TIntegers):Integer;overload;
  function IndexOf(const C:Cardinal; const A:TCardinals):Integer;overload;
  function IndexOf(const s:string; const a:TSimpleStrings):Integer;overload;

  function Concatenate(Coords1,Coords2:TCoords):TCoords;overload;
  function Concatenate(A1,A2:TSimpleStrings):TSimpleStrings;overload;

  procedure AddUniqueToArray(const Elm:string;var Arr:TSimpleStrings);overload;
  procedure AddUniqueToArray(const Elm:Cardinal;var Arr:TCardinals);overload;
  procedure AddUniqueToArray(const Elm:Integer;var Arr:TIntegers);overload;

  procedure AppendToArray(const Suffix:TCoords; var Arr:TCoords);overload;
  procedure AppendToArray(const Suffix:TSimpleStrings; var Arr:TSimpleStrings);overload;
  procedure AppendToArray(const Suffix:TCardinals; var Arr:TCardinals);overload;
  procedure AppendToArray(const Suffix:TIntegers; var Arr:TIntegers);overload;



  function Min(vals:TFLoats):TFLoat;overload;
  function Max(vals:TFLoats):TFLoat;overload;
  function Min(vals:TMatrix):TFLoat;overload;
  function Max(vals:TMatrix):TFLoat;overload;

  function MinIx(vals:TFLoats):Integer;overload;
  function MaxIx(vals:TFLoats):Integer;overload;

  function Sum(vals:TFloats):TFloat;overload;
  function Sum(vals:TIntegers):Integer;overload;


  function Min(const C1,C2:TCoord):TCoord;overload;
  function Max(const C1,C2:TCoord):TCoord;overload;
  function Coord(X,Y,Z:TFloat):TCoord;

  function StringToFloats(S:string):TFloats;
    //converts a string of numbers, separated by white spaces

  // Array generation utils
  function FilledInts(Len,Val: Integer): TIntegers;

  function StringToFloat(const S:String): TFloat;
    //tries comma and point as decimal separator



const
  NullVector:TCoord=(0,0,0);

  // approx 10x eps, for some safety checks on small numbers

  {$IFDEF SINGLEPRECISION}
  TinyFloat:TFloat = 1e-6;
  {$ELSE}
  TinyFloat:TFloat = 1e-14;
  {$ENDIF}

implementation

var
  PointSeparator, CommaSeparator: TFormatSettings;

procedure AddToArray(s:string; var a:TSimpleStrings); overload;

begin
  SetLength(a,Length(a)+1);
  a[High(a)]:=s;
end;

procedure AddToArray(c:TCoord; var a:TCoords); overload;

begin
  SetLength(a,Length(a)+1);
  a[High(a)]:=c;
end;

procedure AddToArray(i:Integer; var a:TIntegers); overload;
begin
  SetLength(a,Length(a)+1);
  a[High(a)]:=i;
end;

procedure AddToArray(f:TFloat; var a:TFloats); overload;
begin
  SetLength(a,Length(a)+1);
  a[High(a)]:=f;
end;

procedure AddToArray(c:Cardinal; var a:TCardinals); overload;

begin
  SetLength(a,Length(a)+1);
  a[High(a)]:=c;
end;

procedure AddToArray(b:Boolean; var a:TBooleans); overload;

begin
  SetLength(a,Length(a)+1);
  a[High(a)]:=b;
end;

procedure RemoveFromArray(Ix:Integer;var A:TIntegers);overload;

begin
  A[Ix]:=A[High(A)];
  SetLength(A,Length(A)-1);
end;

procedure RemoveFromArray(Ix:Integer;var A:TFloats);overload;

begin
  A[Ix]:=A[High(A)];
  SetLength(A,Length(A)-1);
end;

procedure RemoveFromArray(Ix:Integer;var A:TCardinals);overload;

begin
  A[Ix]:=A[High(A)];
  SetLength(A,Length(A)-1);
end;

procedure RemoveFromArray(Ix:Integer;var A:TCoords);overload;

begin
  A[Ix]:=A[High(A)];
  SetLength(A,Length(A)-1);
end;

procedure RemoveFromArray(Ix:Integer;var A:TSimpleStrings);overload;

begin
  A[Ix]:=A[High(A)];
  SetLength(A,Length(A)-1);
end;

function Concatenate(Coords1,Coords2:TCoords):TCoords;overload;

var f,len1:Integer;

begin
  len1:=Length(Coords1);
  SetLength(Result,len1+Length(Coords2));
  for f:=0 to High(Coords1) do
    Result[f]:=Coords1[f];
  for f:=len1 to High(Result) do
    Result[f]:=Coords2[f-len1];
end;

function Concatenate(A1,A2:TSimpleStrings):TSimpleStrings;overload;

var f,len1:Integer;

begin
  len1:=Length(A1);
  SetLength(Result,len1+Length(A2));
  for f:=0 to High(A1) do
    Result[f]:=A1[f];
  for f:=len1 to High(Result) do
    Result[f]:=A2[f-len1];
end;

function IndexOf(const i:Integer; const a:TIntegers):Integer;

begin
  Result:=High(a);
  while (Result>=0) and (a[Result]<>i) do
    Dec(Result);
end;

function IndexOf(const C:Cardinal; const A:TCardinals):Integer;overload;

begin
  Result:=High(A);
  while (Result>=0) and (A[Result]<>C) do
    Dec(Result);
end;


function IndexOf(const s:string; const a:TSimpleStrings):Integer;

begin
  Result:=High(a);
  while (Result>=0) and (a[Result]<>s) do
    Dec(Result);
end;


procedure AddUniqueToArray(const Elm:string;var Arr:TSimpleStrings);overload;

begin
  if IndexOf(Elm,Arr)<0 then
      AddToArray(Elm,Arr);
end;

procedure AddUniqueToArray(const Elm:Cardinal;var Arr:TCardinals);overload;

begin
  if IndexOf(Elm,Arr)<0 then
      AddToArray(Elm,Arr);
end;

procedure AddUniqueToArray(const Elm:Integer;var Arr:TIntegers);overload;

begin
  if IndexOf(Elm,Arr)<0 then
      AddToArray(Elm,Arr);
end;

procedure AppendToArray(const Suffix: TCoords; var Arr: TCoords);
var
  f,len:Integer;

begin
  if Suffix<>nil then
    begin
      len:=Length(Arr);
      SetLength(Arr,len+Length(Suffix));
      for f:=0 to High(Suffix) do
        Arr[f+len]:=Suffix[f];
    end;
end;

procedure AppendToArray(const Suffix:TSimpleStrings; var Arr:TSimpleStrings);overload;

var
  f,len:Integer;

begin
  if Suffix<>nil then
    begin
      len:=Length(Arr);
      SetLength(Arr,len+Length(Suffix));
      for f:=0 to High(Suffix) do
        Arr[f+len]:=Suffix[f];
    end;
end;

procedure AppendToArray(const Suffix:TCardinals; var Arr:TCardinals);overload;

var
  f,len:Integer;

begin
  if Suffix<>nil then
    begin
      len:=Length(Arr);
      SetLength(Arr,len+Length(Suffix));
      for f:=0 to High(Suffix) do
        Arr[f+len]:=Suffix[f];
    end;
end;

procedure AppendToArray(const Suffix:TIntegers; var Arr:TIntegers);overload;

var
  f,len:Integer;

begin
  if Suffix<>nil then
    begin
      len:=Length(Arr);
      SetLength(Arr,len+Length(Suffix));
      for f:=0 to High(Suffix) do
        Arr[f+len]:=Suffix[f];
    end;
end;


function Min(vals: TFLoats): TFLoat;

var f:Integer;

begin
  Result:=0;
  if vals<>nil then
    begin
    Result:=vals[0];
    for f:=1 to High(vals) do
      if vals[f]<Result then Result:=vals[f];
    end;
end;

function Max(vals: TFLoats): TFLoat;

var f:Integer;

begin
  Result:=0;
  if vals<>nil then
    begin
    Result:=vals[0];
    for f:=1 to High(vals) do
      if vals[f]>Result then Result:=vals[f];
    end;
end;

function Min(vals: TMatrix): TFLoat;

var f,g:Integer;

begin
  Result:=0;
  if (vals<>nil) and (vals[0]<>nil) then
    begin
    Result:=vals[0,0];
    for f:=0 to High(vals) do
      for g:=0 to High(vals[f]) do
        if vals[f,g]<Result then Result:=vals[f,g];
    end;
end;


function Max(vals: TMatrix): TFLoat;

var f,g:Integer;

begin
  Result:=0;
  if (vals<>nil) and (vals[0]<>nil) then
    begin
    Result:=vals[0,0];
    for f:=0 to High(vals) do
      for g:=0 to High(vals[f]) do
        if vals[f,g]>Result then Result:=vals[f,g];
    end;
end;

function Min(const C1,C2:TCoord):TCoord;overload;

begin
  if C1[0]<C2[0] then Result[0]:=C1[0] else Result[0]:=C2[0];
  if C1[1]<C2[1] then Result[1]:=C1[1] else Result[1]:=C2[1];
  if C1[2]<C2[2] then Result[2]:=C1[2] else Result[2]:=C2[2];
end;

function Max(const C1,C2:TCoord):TCoord;overload;

begin
  if C1[0]>C2[0] then Result[0]:=C1[0] else Result[0]:=C2[0];
  if C1[1]>C2[1] then Result[1]:=C1[1] else Result[1]:=C2[1];
  if C1[2]>C2[2] then Result[2]:=C1[2] else Result[2]:=C2[2];
end;

function MinIx(vals:TFLoats):Integer;overload;

var
  f:Integer;
  mi:TFloat;

begin
  Result:=-1;
  if vals<>nil then
    begin
    Result:=0;
    mi:=vals[0];
    for f:=1 to High(vals) do
      if vals[f]<mi then
        begin
        mi:=vals[f];
        Result:=f;
        end;
    end;
end;

function MaxIx(vals:TFLoats):Integer;overload;

var
  f:Integer;
  mi:TFloat;

begin
  Result:=-1;
  if vals<>nil then
    begin
    Result:=0;
    mi:=vals[0];
    for f:=1 to High(vals) do
      if vals[f]>mi then
        begin
        mi:=vals[f];
        Result:=f;
        end;
    end;
end;

function Sum(vals:TFloats):TFloat;overload;

var f:Integer;

begin
  Result:=0;
  for f:=0 to High(vals) do
    Result:=Result+vals[f];
end;

function Sum(vals:TIntegers):Integer;overload;

var f:Integer;

begin
  Result:=0;
  for f:=0 to High(vals) do
    Result:=Result+vals[f];
end;


function Coord(X,Y,Z:TFloat):TCoord;

begin
  Result[0]:=X;
  Result[1]:=Y;
  Result[2]:=Z;
end;

function FilledInts(Len,Val: Integer): TIntegers;

var f:Integer;

begin
  SetLength(Result,Len);
  for f:=0 to Len-1 do Result[f]:=Val
end;

function StringToFloat(const S: String): TFloat;
begin
  if Pos('.', S) > 0 then Result := StrToFloat(S,PointSeparator)
  else Result := StrToFloat(S,CommaSeparator);
end;

function StringToFloats(S:string):TFloats;

var
  f:Integer;
  t:string;

begin
  t:='';
  Result:=nil;
  for f:=1 to Length(S) do
    begin
    if S[f]>' ' then
      t:=t+s[f]
    else
      if t<>'' then
        begin
        AddToArray(StrToFloat(t),Result);
        t:='';
        end;
    end;
end;

procedure InitializeGlobals;

begin
  PointSeparator := DefaultFormatSettings;
  PointSeparator.DecimalSeparator := '.';
  PointSeparator.ThousandSeparator := '#';
  CommaSeparator := DefaultFormatSettings;
  CommaSeparator.DecimalSeparator := ',';
  CommaSeparator.ThousandSeparator := '#';
end;

initialization
  InitializeGlobals;
end.

