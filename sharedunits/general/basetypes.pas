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
  Classes, SysUtils,
  {$IFDEF MSWINDOWS}
    Windows;
  {$ELSE}
    LclIntf;
  {$ENDIF}

const
  //Tiny is a small number that generally should be rounded to zero,
  //as in most cases it is due to rounding errors.
  {$IFDEF SINGLEPRECISION}
  Tiny = 1e-6;
  {$ELSE}
  Tiny = 1e-12;
  {$ENDIF}

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
  TCoordGroups = array of TCoords;
  TIntegers = array of Integer;
  TIntegerTable = array of TIntegers;
  TCardinals = array of Cardinal;
  TBooleans = array of Boolean;

  TCuboid=array[0..1] of TCoord;
  TCuboids=array of TCuboid;

  procedure ForceNotZero(var Val:TFloat);
    //ensures Val<=-Tiny or Val>=Tiny

  //inserts at the beginning of the array
  procedure PushIntoArray(I:Integer; var A:TIntegers); overload;
  procedure PushIntoArray(I:string; var A:TSimpleStrings); overload;

  //Appends at the end of the array
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

  procedure RemoveFromArray(Ixs:TIntegers;var A:TIntegers);overload;

  //Return last index, or -1 if not found
  function IndexOf(const i:Integer; const a:TIntegers):Integer;overload;
  function IndexOf(const C:Cardinal; const A:TCardinals):Integer;overload;
  function IndexOf(const s:string; const a:TSimpleStrings):Integer;overload;

  function IsInArray(const i:Integer; const a:TIntegers):Boolean;overload;
  function IsInArray(const C:Cardinal; const A:TCardinals):Boolean;overload;
  function IsInArray(const s:string; const a:TSimpleStrings):Boolean;overload;

  function Concatenate(Coords1,Coords2:TCoords):TCoords;overload;
  function Concatenate(A1,A2:TSimpleStrings):TSimpleStrings;overload;

  function Slice(const Ints:TIntegers;Ix1,Ix2:Integer):TIntegers;overload;
  function Slice(const Floats:TFloats;Ix1,Ix2:Integer):TFloats;overload;
  function Slice(const Coords:TCoords;Ixs:TIntegers):TCoords;overload;

  function CountInArray(I:Integer; const Ints:TIntegers):Integer;overload;


  procedure AddUniqueToArray(const Elm:string;var Arr:TSimpleStrings);overload;
  procedure AddUniqueToArray(const Elm:Cardinal;var Arr:TCardinals);overload;
  procedure AddUniqueToArray(const Elm:Integer;var Arr:TIntegers);overload;

  procedure AppendToArray(const Suffix:TCoords; var Arr:TCoords);overload;
  procedure AppendToArray(const Suffix:TSimpleStrings; var Arr:TSimpleStrings);overload;
  procedure AppendToArray(const Suffix:TCardinals; var Arr:TCardinals);overload;
  procedure AppendToArray(const Suffix:TIntegers; var Arr:TIntegers);overload;
  procedure AppendToArray(const Suffix:TFloats; var Arr:TFloats);overload;

  function Min(vals:TFLoats):TFLoat;overload;
  function Max(vals:TFLoats):TFLoat;overload;
  function Min(vals:TMatrix):TFLoat;overload;
  function Max(vals:TMatrix):TFLoat;overload;
  function Min(vals:TCoords):TCoord;overload;
  function Max(vals:TCoords):TCoord;overload;
  function Min(vals:TIntegers):Integer;overload;
  function Max(vals:TIntegers):Integer;overload;

  function MinIx(vals:TFLoats):Integer;overload;
  function MaxIx(vals:TFLoats):Integer;overload;

  procedure MinValIx(Vals:TFloats;out MinVal:TFloat;out MinIx:Integer);overload;
  procedure MaxValIx(Vals: TFloats; out MinVal: TFloat; out MinIx: Integer);overload;

  function Sum(vals:TIntegers):Integer;overload;
  function Sum(vals:TFloats):TFloat;overload;
  function Sum(vals:TCoords):TCoord;overload;
  function Sum(const Vals1,Vals2:TFloats):TFloats;overload;

  function Min(const C1,C2:TCoord):TCoord;overload;
  function Max(const C1,C2:TCoord):TCoord;overload;
  function Min(const F1,F2:TFloat):TFloat;overload;
  function Max(const F1,F2:TFloat):TFloat;overload;
  function Min(const I1,I2:Integer):Integer;overload;
  function Max(const I1,I2:Integer):Integer;overload;
  function Coord(X,Y,Z:TFloat):TCoord;

  function IsEqual(A1,A2:TIntegers):Boolean;overload;

  function StringToFloat(const S:String): TFloat;
    //tries comma and point as decimal separator

  function StringToFloats(S:string):TFloats;
    //converts a string of numbers, separated by white spaces
    //tries comma and point as decimal separator (uses StringToFloat)

  function StringsToFloats(SS:TSimpleStrings):TFloats;
    //converts to float using StringToFloat

  // Array generation utils
  function FilledInts(Len,Val: Integer): TIntegers;
  function FilledFloats(Len:Integer;Val:TFloat):TFLoats;

  function ScaleMatrix(Mat:TMatrix;Scale:TFloat):TMatrix;
  function AddMatrices(Mat1,Mat2:TMatrix):TMatrix;
  function InContact(Cuboid1,Cuboid2:TCuboid):Boolean;
  function Average(Fs:TFloats):TFloat;overload;
  function Average(Ints:TIntegers):TFloat;overload;

  function Median(Fs:TFloats):TFloat;overload;
  function Median(Ints:TIntegers):Integer;overload;

  function Variance(Fs:TFloats;Avrg:TFloat):TFloat;

  function IsBetween(Val,Extreme1,Extreme2:TFloat):Boolean;

  //for time profiling
  function GetTickCount : DWORD;
  function GetTimeInteval(var StartTick:DWORD):Integer;

const
  NullVector:TCoord=(0,0,0);

  // approx 10x eps, for some safety checks on small numbers

  {$IFDEF SINGLEPRECISION}
  TinyFloat:TFloat = 1e-6;
  {$ELSE}
  TinyFloat:TFloat = 1e-14;
  {$ENDIF}

implementation

uses quicksort;

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

procedure ForceNotZero(var Val: TFloat);
begin
  if (Val<0) and (Val>-Tiny) then
    Val:=-Tiny
  else if (Val>=0) and (Val<Tiny) then
    Val:=Tiny;
end;

procedure PushIntoArray(I: Integer; var A: TIntegers);

var f:Integer;

begin
  SetLength(A,Length(A)+1);
  for f:=High(A) downto 1 do
    A[f]:=A[f-1];
  A[0]:=I;
end;

procedure PushIntoArray(I: string; var A: TSimpleStrings);

var f:Integer;

begin
  SetLength(A,Length(A)+1);
  for f:=High(A) downto 1 do
    A[f]:=A[f-1];
  A[0]:=I;
end;

procedure AddToArray(i:Integer; var a:TIntegers); overload;
begin
  SetLength(a,Length(a)+1);
  a[High(a)]:=i;
end;

procedure AddToArray(f: TFLoat; var a: TFloats);
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

function Slice(const Ints: TIntegers; Ix1, Ix2: Integer): TIntegers;

var f:Integer;

begin
  SetLength(Result,Ix2-Ix1+1);
  for f:=Ix1 to Ix2 do
    Result[f-Ix1]:=Ints[f];
end;

function Slice(const Floats: TFloats; Ix1, Ix2: Integer): TFloats;

var f:Integer;

begin
  SetLength(Result,Ix2-Ix1+1);
  for f:=Ix1 to Ix2 do
    Result[f-Ix1]:=Floats[f];
end;

function Slice(const Coords: TCoords; Ixs: TIntegers): TCoords;

var f:Integer;

begin
  SetLength(Result,Length(Ixs));
  for f:=0 to High(Ixs) do
    Result[f]:=Coords[Ixs[f]];
end;


function CountInArray(I: Integer; const Ints: TIntegers): Integer;

var f:Integer;

begin
  Result:=0;
  for f:=0 to High(Ints) do
    if Ints[f]=I then Inc(Result);
end;

procedure RemoveFromArray(Ixs: TIntegers; var A: TIntegers);

var f,top:Integer;

begin
  top:=Length(A);
  f:=0;
  while f<top do
    begin
    if IsInArray(f,Ixs) then
      begin
      Dec(top);
      A[f]:=A[top];
      end;
    Inc(f);
    end;
  SetLength(A,top);
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

function IsInArray(const i:Integer; const a:TIntegers):Boolean;overload;

begin
  Result:=IndexOf(i,a)>=0;
end;

function IsInArray(const C:Cardinal; const A:TCardinals):Boolean;overload;

begin
  Result:=IndexOf(C,A)>=0;
end;

function IsInArray(const s:string; const a:TSimpleStrings):Boolean;overload;

begin
  Result:=IndexOf(s,a)>=0;
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

procedure AppendToArray(const Suffix: TFloats; var Arr: TFloats);

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

function Min(vals: TCoords): TCoord;

var f:Integer;

begin
  if Length(vals)>0 then
    begin
    Result:=vals[0];
    for f:=1 to High(vals) do
      Result:=Min(Result,vals[f]);
    end
  else Result:=NullVector;
end;

function Max(vals: TCoords): TCoord;

var f:Integer;

begin
  if Length(vals)>0 then
    begin
    Result:=vals[0];
    for f:=1 to High(vals) do
      Result:=Max(Result,vals[f]);
    end
  else Result:=NullVector;
end;

function Min(vals: TIntegers): Integer;

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

function Max(vals: TIntegers): Integer;

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

function Min(const F1,F2:TFloat):TFloat;overload;
begin
  if F1<F2 then Result:=F1 else Result:=F2;
end;

function Max(const F1,F2:TFloat):TFloat;overload;

begin
  if F1>F2 then Result:=F1 else Result:=F2;
end;

function Min(const I1,I2:Integer):Integer;overload;

begin
  if I1<I2 then Result:=I1 else Result:=I2;
end;


function Max(const I1,I2:Integer):Integer;overload;

begin
  if I1>I2 then Result:=I1 else Result:=I2;
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

function Sum(vals: TCoords): TCoord;

var f:Integer;

begin
  Result:=NullVector;
  if Length(vals)=0 then
    for f:=0 to High(vals) do
      begin
      Result[0]:=Result[0]+vals[f,0];
      Result[1]:=Result[1]+vals[f,1];
      Result[2]:=Result[2]+vals[f,2];
      end;
end;

function Sum(const Vals1, Vals2: TFloats): TFloats;

var f,len:Integer;

begin
  len:=Min(Length(Vals1),Length(Vals2));
  SetLength(Result,len);
  for f:=0 to High(Result) do
    Result[f]:=Vals1[f]+Vals2[f];
end;

procedure MinValIx(Vals: TFloats; out MinVal: TFloat; out MinIx: Integer);

var
  f:Integer;

begin
  MinIx:=-1;
  MinVal:=0;
  if vals<>nil then
    begin
    MinIx:=0;
    MinVal:=Vals[0];
    for f:=1 to High(Vals) do
      if vals[f]<MinVal then
        begin
        MinVal:=Vals[f];
        MinIx:=f;
        end;
    end;
end;


procedure MaxValIx(Vals: TFloats; out MinVal: TFloat; out MinIx: Integer);

var
  f:Integer;

begin
  MinIx:=-1;
  MinVal:=0;
  if vals<>nil then
    begin
    MinIx:=0;
    MinVal:=Vals[0];
    for f:=1 to High(Vals) do
      if Vals[f]>MinVal then
        begin
        MinVal:=Vals[f];
        MinIx:=f;
        end;
    end;
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

function StringsToFloats(SS: TSimpleStrings): TFloats;

var f:Integer;

begin
  SetLength(Result,Length(SS));
  for f:=0 to High(SS) do
    Result[f]:=StringToFloat(SS[f]);
end;

function FilledInts(Len,Val: Integer): TIntegers;

var f:Integer;

begin
  SetLength(Result,Len);
  for f:=0 to Len-1 do Result[f]:=Val
end;

function FilledFloats(Len: Integer; Val: TFloat): TFLoats;

var f:Integer;

begin
  SetLength(Result,Len);
  for f:=0 to Len-1 do Result[f]:=Val
end;

function IsEqual(A1, A2: TIntegers): Boolean;

var f:Integer;

begin
  Result:=Length(A1)=Length(A2);
  if Result then
    for f:=0 to High(A1) do
      if A1[f]<>A2[f] then
        begin
        Result:=False;
        Break;
        end;
end;

function StringToFloat(const S: String): TFloat;
begin
  if Pos('.', S) > 0 then Result := StrToFloat(S,PointSeparator)
  else Result := StrToFloat(S,CommaSeparator);
end;

function ScaleMatrix(Mat: TMatrix; Scale: TFloat): TMatrix;

var x,y:Integer;

begin
  SetLength(Result,Length(Mat),Length(Mat[0]));
  for x:=0 to High(Mat) do
    for y:=0 to High(Mat[x]) do
      Result[x,y]:=Mat[x,y]*Scale;
end;

function AddMatrices(Mat1, Mat2: TMatrix): TMatrix;

var x,y:Integer;

begin
  if (Mat1=nil) or (Mat2=nil) or (Length(Mat1)<>Length(Mat2))
  or (Length(Mat1[0])<>Length(Mat2[0])) then
    Result:=nil
  else
    begin
    SetLength(Result,Length(Mat1),Length(Mat1[0]));
    for x:=0 to High(Mat1) do
      for y:=0 to High(Mat1[x]) do
        Result[x,y]:=Mat1[x,y]+Mat2[x,y];
    end;
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
        AddToArray(StringToFloat(t),Result);
        t:='';
        end;
    end;
  if t<>'' then
     AddToArray(StringToFloat(t),Result);
end;

function InContact(Cuboid1,Cuboid2:TCuboid):Boolean;

var c1,c2:TCoord;

begin
  c1:=Max(Cuboid1[0],Cuboid2[0]);
  c2:=Min(Cuboid1[1],Cuboid2[1]);
  Result:=(c1[0]<=c2[0]) and (c1[1]<=c2[1]) and (c1[2]<=c2[2]);
end;

procedure InitializeGlobals;

begin
  PointSeparator := DefaultFormatSettings;
  PointSeparator.DecimalSeparator := '.';
  PointSeparator.ThousandSeparator := ',';
  CommaSeparator := DefaultFormatSettings;
  CommaSeparator.DecimalSeparator := ',';
  CommaSeparator.ThousandSeparator := '.';
end;

function Average(Fs:TFloats):TFloat;

var f:Integer;

begin
  Result:=0;
  if Length(Fs)>0 then
    begin
    for f:=0 to High(Fs) do
      Result:=Result+Fs[f];
    Result:=Result/Length(Fs);
    end;
end;

function Average(Ints:TIntegers):TFloat;

var f:Integer;

begin
  Result:=0;
  if Length(Ints)>0 then
    begin
    for f:=0 to High(Ints) do
      Result:=Result+Ints[f];
    Result:=Result/Length(Ints);
    end;
end;

function Median(Fs: TFloats): TFloat;

var
  ixs:TIntegers;
  ix:Integer;

begin
  Result:=0;
  if Length(Fs)>0 then
    begin
    ixs:=QSAscendingIndex(Fs);
    ix:=Trunc(Length(Fs)/2);
    Result:=Fs[ixs[ix]];
    end;
end;

function Median(Ints: TIntegers):Integer;

var
  ixs:TIntegers;
  ix:Integer;

begin
  Result:=0;
  if Length(Ints)>0 then
    begin
    ixs:=QSAscendingIndex(Ints);
    ix:=Trunc(Length(Ints)/2);
    Result:=Ints[ixs[ix]];
    end;
end;

function Variance(Fs:TFloats;Avrg:TFloat):TFloat;

var f:Integer;

begin
  Result:=0;
  for f:=0 to High(Fs) do
    Result:=Result+Sqr(Fs[f]-Avrg);
  Result:=Result/Length(Fs);
end;

function IsBetween(Val, Extreme1, Extreme2: TFloat): Boolean;
begin
  Result := ((Val>=Extreme1) and (Val<=Extreme2)) or
            ((Val<=Extreme1) and (Val>=Extreme2));
end;

function GetTickCount: DWORD;
  //(from OrphPort package)
   {On Windows, this is number of milliseconds since Windows was
   started. On non-Windows platforms, LCL returns number of
   milliseconds since Dec. 30, 1899, wrapped by size of DWORD.
   This value can overflow LongInt variable when checks turned on,
   so "wrap" value here so it fits within LongInt.
  Also, since same thing could happen with Windows that has been
   running for at least approx. 25 days, override it too.}
begin
{$IFDEF MSWINDOWS}
  Result := Windows.GetTickCount mod High(LongInt);
{$ELSE}
  Result := LclIntf.GetTickCount mod High(LongInt);
{$ENDIF}
end;

function GetTimeInteval(var StartTick: DWORD): Integer;
//returns time interval and resets timer

var endtick:DWORD;

begin
  endtick:=GetTickCount;
  Result:=endtick-StartTick;
  StartTick:=endtick;
end;


initialization
  InitializeGlobals;
end.

