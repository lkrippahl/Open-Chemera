unit ocstringutils;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes,LCLProc;

function SplitString(AString:string;Separator:string):TOCStrings;
function GetInteger(AString:string;Start,Finish:Integer; out Val:Integer):Boolean;overload;
function GetInteger(AString:string;Start,Finish:Integer):Integer;overload;
function GetFloat(AString:string;Start,Finish:Integer; out Val:Double):Boolean;overload;
function GetFloat(AString:string;Start,Finish:Integer):Double;overload;
function GetString(AString:string;Start,Finish:Integer):string;
function LastIndexOf(s:string;a:TOCStrings):Integer; // returns last index of string, -1 if not found
function Deblank(s:string):string; // returns string minus blank spaces (only >32)
function StringAsArray(s:string):TOCStrings;
procedure SaveToFile(ss:TOCStrings;filename:string);
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
      if p>1 then AddString(Copy(AString,1,p-1),Result);
      Delete(AString,1,p+len);
      end;
  until p<1;
  if AString<>'' then AddString(AString,Result);
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


end.

