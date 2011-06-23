{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 7.4.2011
Purpose:
  Utilities for handling strings and string arrays
Requirements:
Revisions:
To do:
*******************************************************************************}

unit ocstringutils;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes;

function SplitString(AString:string;Separator:string):TOCStrings;
function GetInteger(AString:string;Start,Finish:Integer; out Val:Integer):Boolean;overload;
function GetInteger(AString:string;Start,Finish:Integer):Integer;overload;
function GetFloat(AString:string;Start,Finish:Integer; out Val:Double):Boolean;overload;
function GetFloat(AString:string;Start,Finish:Integer):Double;overload;
function GetString(AString:string;Start,Finish:Integer):string;
function LastIndexOf(s:string;a:TOCStrings):Integer;
  // returns last index of string, -1 if not found
function Deblank(s:string):string;
  // returns string minus blank spaces (only >32)
function TrimmedBlanks(const S:string):string;
  //removes all <=' ' characters from the extremities of the string;
function StringAsArray(s:string):TOCStrings;
procedure SaveToFile(ss:TOCStrings;filename:string);
function FlattenStrings(SS:TOCStrings;Sep:string):string;
function SnipString(var S:string;const Spacer:string):string;
  //cuts S to first spacer, removing spacer, returns part before spacer
  //returns empty string and leaves S if spacer is not found
function CountInString(const S:string; const C:Char):Integer;
function CleanString(const S:string; const C:Char):string;
  //returns S minus all instances of C
function AppendToLength(const S:string; const C:Char; const Len:Integer):string;
  //returns S filled with C until length Len


implementation

function SplitString(AString:string;Separator:string):TOCStrings;

var
  p,len:Integer;

begin
  Result:=nil;
  len:=Length(Separator)-1;
  repeat
    p:=Pos(Separator,AString);
    if p>0 then
      begin
      if p>1 then AddToArray(Copy(AString,1,p-1),Result);
      Delete(AString,1,p+len);
      end;
  until p<1;
  if AString<>'' then AddToArray(AString,Result);
end;

function GetInteger(AString:string;Start,Finish:Integer; out Val:Integer):Boolean;

var
  f:Integer;
  s:string;

begin
  s:='';
  f:=Start;
  while (f<=Length(AString)) and (f<=Finish) do
    begin
    if AString[f]<>' ' then s:=s+AString[f];
    Inc(f);
    end;
  try
    Val:=StrToInt(s);
    Result:=True;
  except
    Val:=0;
    Result:=False;
  end;
end;

function GetFloat(AString:string;Start,Finish:Integer; out Val:Double):Boolean;

var
  f:Integer;
  s:string;

begin
  s:='';
  f:=Start;
  while (f<=Length(AString)) and (f<=Finish) do
    begin
    if AString[f] in [',','.'] then s:=s+DecimalSeparator
      else if AString[f]<>' ' then s:=s+AString[f];
    Inc(f);
    end;
  try
    Val:=StrToFloat(s);
    Result:=True;
  except
    Result:=False;
    Val:=0;
  end;
end;

function GetInteger(AString:string;Start,Finish:Integer):Integer;

begin
  if not GetInteger(AString,Start,Finish,Result) then
    Result:=-1;
end;

function GetFloat(AString:string;Start,Finish:Integer):Double;
begin
  if not GetFloat(AString,Start,Finish,Result) then
    Result:=-1;
end;

function GetString(AString:string;Start,Finish:Integer):string;

var
  f:Integer;

begin
  Result:='';
  f:=Start;
  while (f<=Length(AString)) and (f<=Finish) do
    begin
    if AString[f]<>' ' then Result:=Result+AString[f];
    Inc(f);
    end;
end;

function LastIndexOf(s: string; a: TOCStrings): Integer;
begin
  Result:=High(a);
  while (Result>=0) and (a[Result]<>s) do
    Dec(Result);
end;

function Deblank(s: string): string;

var f:Integer;

begin
  Result:='';
  for f:=1 to Length(s) do
    if s[f]>' ' then Result:=Result+s[f];
end;

function StringAsArray(s: string): TOCStrings;
begin
  SetLength(Result,1);
  Result[0]:=s;
end;

procedure SaveToFile(ss: TOCStrings; filename: string);

var
  sl:TStringList;
  f:Integer;

begin
  sl:=TStringList.Create;
  for f:=0 to High(ss) do
    sl.Add(ss[f]);
  sl.SaveToFile(filename);
  sl.Free;
end;

function FlattenStrings(SS: TOCStrings; Sep: string): string;

var f:Integer;

begin
  if Length(SS)>0 then
    begin
    Result:=SS[0];
    for f:=1 to High(SS) do Result:=Result+Sep+SS[f];
    end
  else Result:='';
end;

function SnipString(var S:string;const Spacer:string):string;

var
  ix:Integer;
begin
  Result:='';
  ix:=Pos(Spacer,S);
  if ix>0 then
    begin
    Result:=Copy(S,1,Ix-1);
    Delete(S,1,Ix+Length(Spacer)-1);
    end;
end;

function CountInString(const S: string; const C: Char): Integer;

var f:Integer;

begin
  Result:=0;
  for f:=1 to Length(S) do
    if S[f]=C then Inc(Result);
end;

function TrimmedBlanks(const S:string):string;

var
  ix1,ix2:Integer;
begin
  ix2:=Length(S);
  ix1:=1;

  //find first non space;
  while (ix1<ix2) and (S[ix1]<=' ') do
    Inc(ix1);

  //find last non space;
  while (ix2>0) and (S[ix2]<=' ') do
    Dec(ix2);

  if ix2>=ix1 then Result:=Copy(s,ix1,ix2-ix1+1)
  else Result:='';
end;

function CleanString(const S:string; const C:Char):string;

var f,count:Integer;

begin
  SetLength(Result,Length(S));
  count:=0;
  for f:=1 to Length(S) do
    if S[f]<>C then
      begin
      Inc(count);
      Result[count]:=S[f];
      end;
  SetLength(Result,count);
end;

function AppendToLength(const S:string; const C:Char; const Len:Integer):string;

var f,lens:Integer;

begin
  lens:=Length(S);
  if lens>=Len then Exit(S); //nothing to append
  SetLength(Result,Len);
  for f:=1 to lens do
    Result[f]:=S[f];
  for f:=lens+1 to Len do
    Result[f]:=C;
end;

end.

