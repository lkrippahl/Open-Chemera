{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain. Enjoy.
********************************************************************************
Author: Ludwig Krippahl
Date: 9.1.2011
Purpose:
  Base data types (string arrays, coords, etc)

Requirements:
Revisions:
To do: Comments
*******************************************************************************}

unit basetypes;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils;

type
  TOCStrings = array of string;
  TOCFloat = Real; //can be changed for different precision
  TOCFloats = array of TOCFloat;
  TOCDoubles = array of Double;
  TOCSingles = array of Single;
  TOCMatrix= array of TOCFloats;
  TOCCoord = array [1..3] of TOCFloat; //3D coords ou 3D vector
  TOCCoords = array of TOCCoord;
  TOCIntegers = array of Integer;

  procedure AddString(s:string;var a:TOCStrings); // add a string to an array
  procedure AddCoord(c:TOCCoord; var a:TOCCoords); // add a coord to an array
  function AddVectors(c1,c2:TOCCoord):TOCCoord;  // sum of two vectors (as coord)
  function ScaleVector(c:TOCCoord;s:TOCFloat):TOCCoord; // multiply vector (as coord)
  procedure AddInteger(i:Integer; var a:TOCIntegers); // adds an integer to an array
  function Distance(c1,c2:TOCCoord):TOCFloat;
  function Min(vals:TOCFLoats):TOCFLoat;overload;
  function Max(vals:TOCFLoats):TOCFLoat;overload;
  function Min(vals:TOCMatrix):TOCFLoat;overload;
  function Max(vals:TOCMatrix):TOCFLoat;overload;

  // Array generation utils
  function FilledInts(Len,Val: Integer): TOCIntegers;


const
  NullVector:TOCCoord=(0,0,0);

implementation

procedure AddString(s:string;var a:TOCStrings);

begin
  SetLength(a,Length(a)+1);
  a[High(a)]:=s;
end;

procedure AddCoord(c:TOCCoord; var a:TOCCoords);

begin
  SetLength(a,Length(a)+1);
  a[High(a)]:=c;
end;

function AddVectors(c1, c2: TOCCoord): TOCCoord;
begin
  Result[1]:=c1[1]+c2[1];
  Result[2]:=c1[2]+c2[2];
  Result[3]:=c1[3]+c2[3];
end;

function ScaleVector(c: TOCCoord; s: TOCFloat): TOCCoord;
begin
  Result[1]:=c[1]*s;
  Result[2]:=c[2]*s;
  Result[3]:=c[3]*s;
end;

procedure AddInteger(i: Integer; var a: TOCIntegers);
begin
  SetLength(a,Length(a)+1);
  a[High(a)]:=i;
end;

function Distance(c1, c2: TOCCoord): TOCFloat;
begin
  Result:=Sqrt(Sqr(c1[1]-c2[1])+Sqr(c1[2]-c2[2])+Sqr(c1[3]-c2[3]));
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

function FilledInts(Len,Val: Integer): TOCIntegers;

var f:Integer;

begin
  SetLength(Result,Len);
  for f:=0 to Len-1 do Result[f]:=Val
end;

end.

