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
To do: Comments
*******************************************************************************}

unit basetypes;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils;

const
  OCTinyFloat=1e-10;  //for testing zero, differences, etc
                      //make it 1e6 larger than epsilon for TFloat

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

  procedure AddToArray(i:Integer; var a:TIntegers); overload;
  procedure AddToArray(s:string; var a:TSimpleStrings); overload;
  procedure AddToArray(c:TCoord; var a:TCoords); overload;
  procedure AddToArray(f:TFLoat; var a:TFloats); overload;

  function Concatenate(Coords1,Coords2:TCoords):TCoords;overload;
  function Concatenate(A1,A2:TSimpleStrings):TSimpleStrings;overload;

  procedure AddUniqueToArray(s:string; var a:TSimpleStrings); overload;

  function IndexOf(const i:Integer; const a:TIntegers):Integer;overload;
  function IndexOf(const s:string; const a:TSimpleStrings):Integer;overload;
  function Distance(c1,c2:TCoord):TFloat;
  function Min(vals:TFLoats):TFLoat;overload;
  function Max(vals:TFLoats):TFLoat;overload;
  function Min(vals:TMatrix):TFLoat;overload;
  function Max(vals:TMatrix):TFLoat;overload;
  function Min(const C1,C2:TCoord):TCoord;overload;
  function Max(const C1,C2:TCoord):TCoord;overload;
  function Coord(X,Y,Z:TFloat):TCoord;
  function StringToFloats(S:string):TFloats;
    //converts a string of numbers, separated by white spaces

  // Array generation utils
  function FilledInts(Len,Val: Integer): TIntegers;



const
  NullVector:TCoord=(0,0,0);

implementation

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

procedure AddUniqueToArray(s:string; var a:TSimpleStrings); overload;

begin
  if IndexOf(s,a)<0 then AddToArray(s,a);
end;

function IndexOf(const i:Integer; const a:TIntegers):Integer;

begin
  Result:=High(a);
  while (Result>=0) and (a[Result]<>i) do
    Dec(Result);
end;

function IndexOf(const s:string; const a:TSimpleStrings):Integer;

begin
  Result:=High(a);
  while (Result>=0) and (a[Result]<>s) do
    Dec(Result);
end;


function Distance(c1, c2: TCoord): TFloat;
begin
  Result:=Sqrt(Sqr(c1[1]-c2[1])+Sqr(c1[2]-c2[2])+Sqr(c1[0]-c2[0]));
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

end.

