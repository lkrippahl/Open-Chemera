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
function LastIndexOf(s:string;a:TSimpleStrings):Integer;overload;
function LastIndexOf(I:Integer;a:TIntegers):Integer;overload;
  // returns last index of string, -1 if not found
function FirstIndexOf(s:string;a:TSimpleStrings):Integer;overload;
function FirstIndexOf(I:Integer;a:TIntegers):Integer;overload;
  // -1 if not found

function FirstByPrefix(Prefix:string;SStrings:TSimpleStrings):Integer;
  // first index starting with Prefix. Case sensitive
function FirstContaining(Substr:string;SStrings:TSimpleStrings):Integer;
  // first index of string containing Substr


function LookupByPrefix(Prefix:string;SStrings:TSimpleStrings):string;
  // returns string after first prefix, or '' if not found

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
function AsSimpleStrings(const Sl:TStrings):TSimpleStrings; overload;
function AsSimpleStrings(const Sl:TStrings; FirstIx,LastIx:Integer):TSimpleStrings; overload;

function CopyStrings(SS:TSimpleStrings; FirstIx,LastIx:Integer):TSimpleStrings;

procedure AppendToStringList(const Ss:TSimpleStrings;const Sl:TStrings);
function CleanFileName(FileName:string):string;
  //removes illegal characters from filename.
  //TO DO: check if illegal character list is correct
function ReplaceHexCodes(S:string):string;
  //replace "=HH" hex codes with character of same hex code
function PosUncased(Substr,Str:string; Start:Integer=1):Integer;
  //returns position of Substr in Str, ignoring case

//simple XML-like functions
function TagString(Tag,Text:string):string;
  //returns the string with <tag>... </tag>
function GetTaggedField(Tag:string; var Text:string; ForceEnd:Boolean=False):string;
  //returns the string between tags, deleting from text. If ForceEnd then
  //will force end of field at next Tag
function GetTaggedFields(Tag:string; var Text:string; ForceEnd:Boolean=False):TSimpleStrings;
  //returns all fields with this tag
function DecodeQP(const Text:string):string;
  //Decodes  quoted-printable text
function HtmltoAscii(const Text:string):string;
  //converts special & html chars to ascii
function BasicASCII(const Text:string):string;
function ReadAsSimpleStrings(const FileName:string):TSimpleStrings;

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

function LastIndexOf(I:Integer;a:TIntegers):Integer;overload;
begin
  Result:=High(a);
  while (Result>=0) and (a[Result]<>I) do
    Dec(Result);
end;

function FirstIndexOf(s:string;a:TSimpleStrings):Integer;

var f:Integer;

begin
  Result:=-1;
  for f:=0 to High(a) do
    if s=a[f] then
      begin
      Result:=f;
      Break;
      end;
end;

function FirstIndexOf(I:Integer;a:TIntegers):Integer;

var f:Integer;

begin
  Result:=-1;
  for f:=0 to High(a) do
    if I=a[f] then
      begin
      Result:=f;
      Break;
      end;
end;


function FirstByPrefix(Prefix:string;SStrings:TSimpleStrings):Integer;

var f:Integer;

begin
  Result:=-1;
  for f:=0 to High(SStrings) do
    if Pos(Prefix,SStrings[f])=1 then
      begin
      Result:=f;
      Break;
      end;
end;

function FirstContaining(Substr:string;SStrings:TSimpleStrings):Integer;

var f:Integer;

begin
  Result:=-1;
  for f:=0 to High(SStrings) do
    if Pos(Substr,SStrings[f])>0 then
      begin
      Result:=f;
      Break;
      end;
end;


function LookupByPrefix(Prefix:string;SStrings:TSimpleStrings):string;

var ix:Integer;

begin
  ix:=FirstByPrefix(Prefix,SStrings);
  if ix<0 then Result:=''
  else Result:=Copy(SStrings[ix],Length(Prefix)+1,Length(SStrings));
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

function AsSimpleStrings(const Sl:TStrings; FirstIx,LastIx:Integer):TSimpleStrings;

var f:Integer;

begin
  SetLength(Result,LastIx-FirstIx+1);
  for f:=FirstIx to LastIx do
    Result[f-FirstIx]:=Sl.Strings[f];
end;

function CopyStrings(SS:TSimpleStrings; FirstIx,LastIx:Integer):TSimpleStrings;

var f:Integer;

begin
  SetLength(Result,LastIx-FirstIx+1);
  for f:=FirstIx to LastIx do
    Result[f-FirstIx]:=SS[f];
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

function PosUncased(Substr,Str:string; Start:Integer=1):Integer;

var
  f,g:Integer;
  ucsubstr:string;
  last,max:Integer;

begin
  ucsubstr:=UpperCase(Substr);
  Result:=-1;
  last:=Length(Substr);
  max:=Length(Str)-last;
  f:=Start-1;
  while (Result<0) and (f<max) do
    begin
    Result:=f+1;
    for g:=last downto 1 do
      if (Substr[g]<>Str[f+g]) and (ucsubstr[g]<>Str[f+g]) then
        begin
        Result:=-1;
        Break;
        end;
    Inc(f);
    end;
end;

function TagString(Tag,Text:string):string;

begin
  Result:='<'+Tag+'>'+Text+'</'+Tag+'>';
end;

function GetTaggedField(Tag:string; var Text:string; ForceEnd:Boolean=False):string;
  //returns the string between tags, deleting from text. If ForceEnd then
  //will force end of field at next Tag

  //TODO: allow properties

var
  pstart, pend, pnext:Integer;
  staglen,etaglen:Integer;

begin
  pstart:=PosUncased('<'+Tag+'>',Text);
  staglen:=Length(Tag)+2; // change this for properties
  pend:=PosUncased('</'+Tag+'>',Text,pstart+staglen+1);
  etaglen:=Length(Tag)+3; // change this for properties
  if ForceEnd then
    begin
    pnext:=PosUncased('<'+Tag+'>',Text,pstart+3+Length(Tag));
    if (pend<1) or (pnext<pend) then
      begin
      pend:=pnext;
      etaglen:=Length(Tag)+2;
      end;
    if pend<1 then
      begin
      pend:=Length(Text);
      etaglen:=0;
      end;
    end;
  Result:='';
  if (pstart>0) and (pend>0) and (pend>pstart) then
    begin
    Result:=Copy(Text,pstart+staglen,pend-pstart-staglen);
    Delete(Text,pstart,pend-pstart+etaglen-1);
    end;
end;

function GetTaggedFields(Tag:string; var Text:string; ForceEnd:Boolean=False):TSimpleStrings;

var
  fld:string;

begin
  Result:=nil;
  repeat
    fld:=GetTaggedField(Tag,Text,ForceEnd);
    if fld<>'' then AddToArray(fld,Result);
  until fld='';
end;

function DecodeQP(const Text:string):string;

var
  totlen:Integer;
  f,current:Integer;

begin
  totlen:=Length(Text);
  SetLength(Result,totlen);
  current:=0;
  f:=1;
  while f<=totlen do
    begin
    if (Text[f]='=') and (f<totlen-1) then
      begin
      if Text[f+1]<' ' then
        begin
        Inc(f,2);
        if Text[f]<' ' then
          Inc(f); // skip lines with cr+nl
        end
      else
        begin
        Inc(current);
        try
          Result[current]:=Chr(Hex2Dec(Text[f+1]+Text[f+2]));
          Inc(f,3);
          //check if line break was replaced by space, but still broken in file
          if Result[current]=' ' then
            begin
            if Text[f]<' ' then Inc(f);
            if Text[f]<' ' then Inc(f);
            end;
        except
          Inc(current);
          Result[current]:=Text[f];
          Inc(f);
        end;
        end;
      end
    else
      begin
      Inc(current);
      Result[current]:=Text[f];
      Inc(f);
      end;
    end;
  SetLength(Result,current);
end;

function HtmltoAscii(const Text:string):string;

var
  totlen:Integer;
  f,current:Integer;
  tmp:string;

procedure ConvertHtml(var s:string);

begin
    s:=StringReplace(s,'&lt;','<',[rfReplaceAll, rfIgnoreCase]);
    s:=StringReplace(s,'&gt;','>',[rfReplaceAll, rfIgnoreCase]);
    //s:=StringReplace(s,#13+#10,'',[rfReplaceAll, rfIgnoreCase]);
    //s:=StringReplace(s,#13,'',[rfReplaceAll, rfIgnoreCase]);
    s:=StringReplace(s,'<br/>',#13+#10,[rfReplaceAll, rfIgnoreCase]);
    s:=StringReplace(s,'<br>',#13+#10,[rfReplaceAll, rfIgnoreCase]);
    s:=StringReplace(s,'&quot;','"',[rfReplaceAll, rfIgnoreCase]);
    s:=StringReplace(s,'&amp;','&',[rfReplaceAll, rfIgnoreCase]);
    s:=StringReplace(s,'&apos;','''',[rfReplaceAll, rfIgnoreCase]);
    s:=StringReplace(s,'&nbsp;',' ',[rfReplaceAll, rfIgnoreCase]);
end;


begin
  totlen:=Length(Text);
  SetLength(Result,totlen);
  current:=0;
  f:=1;
  while f<totlen-1 do
    begin
    if (Text[f]='&') and (Text[f+1]='#') then
      begin
      f:=f+2;
      tmp:='';
      while (f<totlen) and (Text[f]<>';') do
        begin
        tmp:=tmp+Text[f];
        Inc(f);
        end;
      Inc(f);
      try
        Result[current]:=Chr(StrToInt(tmp));
        Inc(current);
      except
      end;
      end
    else
      begin
      Result[current]:=Text[f];
      Inc(f);
      Inc(current);
      end;
    end;
  SetLength(Result,current);
  ConvertHtml(Result);
end;

function BasicASCII(const Text: string): string;

var
  f,c:Integer;

begin
  c:=0;
  SetLength(Result,Length(Text));
  for f:=1 to Length(Text) do
    if Ord(Text[f])<=127 then
      begin
      Inc(c);
      Result[c]:=Text[f];
      end;
  SetLength(Result,c);
end;

function ReadAsSimpleStrings(const FileName: string): TSimpleStrings;

var tmp:TStringList;

begin
  tmp:=TStringList.Create;
  tmp.LoadFromFile(FileName);
  Result:=AsSimpleStrings(tmp);
  tmp.Free;
end;

end.

