{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 9.1.2011
Purpose:
  Base data types (string arrays, coords, etc) and utility functions

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
                      //make it 1e6 larger than epsilon for TOCFloat

type
  TOCStrings = array of string;
  TOCFloat = Double; //can be changed for different precision
  TOCFloats = array of TOCFloat;
  TOCDoubles = array of Double;
  TOCSingles = array of Single;
  TOCMatrix= array of TOCFloats;
  TOCCoord = array [0..2] of TOCFloat; //3D coords ou 3D vector
  TOCCoords = array of TOCCoord;
  TOCIntegers = array of Integer;

  procedure AddToArray(i:Integer; var a:TOCIntegers); overload;
  procedure AddToArray(s:string; var a:TOCStrings); overload;
  procedure AddToArray(c:TOCCoord; var a:TOCCoords); overload;
  procedure AddToArray(f:TOCFLoat; var a:TOCFloats); overload;

  function Concatenate(Coords1,Coords2:TOCCoords):TOCCoords;overload;
  function Concatenate(A1,A2:TOCStrings):TOCStrings;overload;

  procedure AddUniqueToArray(s:string; var a:TOCStrings); overload;

  function IndexOf(const i:Integer; const a:TOCIntegers):Integer;overload;
  function IndexOf(const s:string; const a:TOCStrings):Integer;overload;
  function Distance(c1,c2:TOCCoord):TOCFloat;
  function Min(vals:TOCFLoats):TOCFLoat;overload;
  function Max(vals:TOCFLoats):TOCFLoat;overload;
  function Min(vals:TOCMatrix):TOCFLoat;overload;
  function Max(vals:TOCMatrix):TOCFLoat;overload;
  function Min(const C1,C2:TOCCoord):TOCCoord;overload;
  function Max(const C1,C2:TOCCoord):TOCCoord;overload;
  function Coord(X,Y,Z:TOCFloat):TOCCoord;
  function StringToFloats(S:string):TOCFloats;
    //converts a string of numbers, separated by white spaces

  // Array generation utils
  function FilledInts(Len,Val: Integer): TOCIntegers;



const
  NullVector:TOCCoord=(0,0,0);

implementation

procedure AddToArray(s:string; var a:TOCStrings); overload;

begin
  SetLength(a,Length(a)+1);
  a[High(a)]:=s;
end;

procedure AddToArray(c:TOCCoord; var a:TOCCoords); overload;

begin
  SetLength(a,Length(a)+1);
  a[High(a)]:=c;
end;

procedure AddToArray(i:Integer; var a:TOCIntegers); overload;
begin
  SetLength(a,Length(a)+1);
  a[High(a)]:=i;
end;

procedure AddToArray(f:TOCFloat; var a:TOCFloats); overload;
begin
  SetLength(a,Length(a)+1);
  a[High(a)]:=f;
end;

function Concatenate(Coords1,Coords2:TOCCoords):TOCCoords;overload;

var f,len1:Integer;

begin
  len1:=Length(Coords1);
  SetLength(Result,len1+Length(Coords2));
  for f:=0 to High(Coords1) do
    Result[f]:=Coords1[f];
  for f:=len1 to High(Result) do
    Result[f]:=Coords2[f-len1];
end;

function Concatenate(A1,A2:TOCStrings):TOCStrings;overload;

var f,len1:Integer;

begin
  len1:=Length(A1);
  SetLength(Result,len1+Length(A2));
  for f:=0 to High(A1) do
    Result[f]:=A1[f];
  for f:=len1 to High(Result) do
    Result[f]:=A2[f-len1];
end;

procedure AddUniqueToArray(s:string; var a:TOCStrings); overload;

begin
  if IndexOf(s,a)<0 then AddToArray(s,a);
end;

function IndexOf(const i:Integer; const a:TOCIntegers):Integer;

begin
  Result:=High(a);
  while (Result>=0) and (a[Result]<>i) do
    Dec(Result);
end;

function IndexOf(const s:string; const a:TOCStrings):Integer;

begin
  Result:=High(a);
  while (Result>=0) and (a[Result]<>s) do
    Dec(Result);
end;


function Distance(c1, c2: TOCCoord): TOCFloat;
begin
  Result:=Sqrt(Sqr(c1[1]-c2[1])+Sqr(c1[2]-c2[2])+Sqr(c1[0]-c2[0]));
end;

function Min(vals: TOCFLoats): TOCFLoat;

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

function Max(vals: TOCFLoats): TOCFLoat;

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

function Min(vals: TOCMatrix): TOCFLoat;

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


function Max(vals: TOCMatrix): TOCFLoat;

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

function Min(const C1,C2:TOCCoord):TOCCoord;overload;

begin
  if C1[0]<C2[0] then Result[0]:=C1[0] else Result[0]:=C2[0];
  if C1[1]<C2[1] then Result[0]:=C1[1] else Result[1]:=C2[1];
  if C1[2]<C2[2] then Result[0]:=C1[2] else Result[2]:=C2[2];
end;

function Max(const C1,C2:TOCCoord):TOCCoord;overload;

begin
  if C1[0]>C2[0] then Result[0]:=C1[0] else Result[0]:=C2[0];
  if C1[1]>C2[1] then Result[0]:=C1[1] else Result[1]:=C2[1];
  if C1[2]>C2[2] then Result[0]:=C1[2] else Result[2]:=C2[2];
end;

function Coord(X,Y,Z:TOCFloat):TOCCoord;

begin
  Result[0]:=X;
  Result[1]:=Y;
  Result[2]:=Z;
end;

function FilledInts(Len,Val: Integer): TOCIntegers;

var f:Integer;

begin
  SetLength(Result,Len);
  for f:=0 to Len-1 do Result[f]:=Val
end;

function StringToFloats(S:string):TOCFloats;

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

