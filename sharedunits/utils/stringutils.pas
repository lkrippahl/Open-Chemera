{*******************************************************************************
This source file is public domain. It is provided as is, with no explicit
or implied warranties regarding anything whatsoever.
********************************************************************************
Author: Ludwig Krippahl
Date: 7.4.2011
Purpose:
  Utilities for handling strings and string arrays  (TSimpleStrings)
Requirements:
Revisions:
To do:
*******************************************************************************}

unit stringutils;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes,strutils;




function GrabWord(var Text:string;const Sep:string=' '):string;
  //deletes separator if found at beginning
  //Returns substring to separator deleting that from Text

function GrabBetween(var Text:string;const Sep1,Sep2:string):string;
  //returns substring between separators deleting that from Text

function SplitString(Text:string;const Sep:string):TSimpleStrings;overload;
  //splits using grabword (i.e. skipping repeated separators)
procedure SplitString(Text:string;Words:TStringList;const Sep:string=' ');overload;
  //same as above, but adds result to the Words TStringList, which must have been
  //created by caller

function GetInteger(AString:string;Start,Finish:Integer; out Val:Integer):Boolean;overload;
function GetInteger(AString:string;Start,Finish:Integer):Integer;overload;
function GetFloat(AString:string;Start,Finish:Integer; out Val:Double):Boolean;overload;
function GetFloat(AString:string;Start,Finish:Integer):Double;overload;
function GetString(AString:string;Start,Finish:Integer):string;
function LastIndexOf(s:string;a:TSimpleStrings):Integer;
  // returns last index of string, -1 if not found
function Deblank(s:string):string;
  // returns string minus blank spaces (only >32)
function TrimmedBlanks(const S:string):string;
  //removes all <=' ' characters from the extremities of the string;
function StringAsArray(s:string):TSimpleStrings;
procedure SaveToFile(ss:TSimpleStrings;filename:string);
function FlattenStrings(SS:TSimpleStrings;Sep:string):string;
function SnipString(var S:string;const Spacer:string):string;
  //cuts S to first spacer, removing spacer, returns part before spacer
  //returns empty string and leaves S if spacer is not found
function CountInString(const S:string; const C:Char):Integer;
function CleanString(const S:string; const C:Char):string;overload;
  //returns S minus all instances of C
function CleanString(const s:string):string;overload;
  //returns all characters >= space
function AppendToLength(const S:string; const C:Char; const Len:Integer):string;
  //returns S filled with C until length Len
function FixLineBreaks(s:string):string;
  //converts to current OS
function UncasedCompare(S1,S2:string):Boolean;
  // string comparison indifferent to case
function AsSimpleStrings(const Sl:TStrings):TSimpleStrings;
procedure AppendToStringList(const Ss:TSimpleStrings;const Sl:TStrings);
function CleanFileName(FileName:string):string;
  //removes illegal characters from filename.
  //TO DO: check if illegal character list is correct
function ReplaceHexCodes(S:string):string;
  //replace "=HH" hex codes with character of same hex code

implementation

function GrabWord(var Text:string;const Sep:string=' '):string;
//consumes the original string

var p:Integer;

begin
  while Pos(Sep,Text)=1 do
    Text:=Copy(Text,Length(Sep)+1,Length(Text));
  p:=Pos(Sep,Text);
  if p<=0 then
    begin
    Result:=Text;
    Text:='';
    end
  else
    begin
    Result:=Copy(Text,1,p-1);
    Text:=Copy(Text,p+Length(Sep),Length(Text));
    end;
end;


function GrabBetween(var Text:string;const Sep1,Sep2:string):string;
//consumes the original string, Text

var p1,p2:Integer;

begin
  p1:=Pos(Sep1,Text);
  p2:=Pos(Sep2,Text);
  if p1<p2 then
    begin
    Result:=Copy(Text,p1+Length(Sep1),p2-p1-Length(Sep1));
    Text:=Copy(Text,1,p1-1)+Copy(Text,p2+Length(Sep2),Length(Text));
    end;
end;

function SplitString(Text: string; const Sep: string): TSimpleStrings;overload;
begin
  Result:=nil;
  while Text<>'' do
    AddToArray(GrabWord(Text,Sep),Result);
end;

procedure SplitString(Text:string;Words:TStringList;const Sep:string=' ');overload;
begin
  while Text<>'' do Words.Add(GrabWord(Text,Sep));
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

function LastIndexOf(s: string; a: TSimpleStrings): Integer;
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

function StringAsArray(s: string): TSimpleStrings;
begin
  SetLength(Result,1);
  Result[0]:=s;
end;

procedure SaveToFile(ss: TSimpleStrings; filename: string);

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

function FlattenStrings(SS: TSimpleStrings; Sep: string): string;

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

function CleanString(const s:string):string;overload;

var
  f:Integer;

begin
  Result:='';
  for f:=1 to Length(s) do
    if s[f]>=' ' then Result:=Result+s[f];
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

function FixLineBreaks(s:string):string;

var f:Integer;

begin
  Result:='';
  for f:=1 to Length(s) do
    if s[f]=#10 then Result:=Result+LineEnding
end;

function UncasedCompare(S1,S2:string):Boolean;
  // string comparison indifferent to case
begin
  Result:=UpperCase(s1)=UpperCase(s2);
end;

function AsSimpleStrings(const Sl:TStrings):TSimpleStrings;

var f:Integer;

begin
  SetLength(Result,Sl.Count);
  for f:=0 to Sl.Count-1 do
    Result[f]:=Sl.Strings[f];
end;

procedure AppendToStringList(const Ss:TSimpleStrings;const Sl:TStrings);

var f:Integer;

begin
  for f:=0 to High(Ss) do Sl.Add(Ss[f]);
end;

function CleanFileName(FileName:string):string;

var
  f,c:Integer;

begin
  SetLength(Result,Length(FileName));
  c:=0;
  for f:=1 to Length(FileName) do
    if (FileName[f]>=' ') and (Pos(FileName[f], '/?<>\:*|"^')<1) then
      begin
      Inc(c);
      Result[c]:=FileName[f];
      end;
  SetLength(Result,c);
end;

function ReplaceHexCodes(S:string):string;
  //replace "=HH" hex codes with character of same hex code

var f:Integer;

begin
  Result:='';
  f:=1;
  while f<=Length(S) do
    begin
    if (f<=Length(S)-2) and
        (S[f]='=') and
        (Pos(S[f+1],'0123456789ABCDEF')>0) and
        (Pos(S[f+2],'0123456789ABCDEF')>0) then
        begin
          Result:=Result+Char(Hex2Dec(S[f+1]+S[f+2]));
          Inc(f,3);
        end
    else
      begin
      Result:=Result+S[f];
      Inc(f,1);
      end;
    end;
end;

end.

