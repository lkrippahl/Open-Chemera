unit miscutils;

{$mode objfpc}{$H+}

interface

uses basetypes;


function Max(s1,s2:Single):Single;overload;
function Min(s1,s2:Single):Single;overload;
function Max(s1,s2:Double):Double;overload;
function Min(s1,s2:Double):Double;overload;
function Max(s1,s2:Integer):Integer;overload;
function Min(s1,s2:Integer):Integer;overload;
function Min(ss:TSingles):Single;overload;
function Max(ss:TSingles):Single;overload;
function Max(ss:TIntegers):Integer;overload;
function IsBetween(db,d1,d2:Single):Boolean;overload;
function IsBetween(db,d1,d2:Integer):Boolean;overload;
procedure Switch(var i1,i2:Integer);overload;
procedure Switch(var i1,i2:Single);overload;
procedure ForceBetween(var i:Integer; low, upper:Integer);overload;
procedure ForceBetween(var i:Single; low, upper:Single);overload;


implementation


function Min(ss:TSingles):Single;overload;

var f:Integer;

begin
  if ss=nil then
    Result:=0
  else
    begin
    Result:=ss[0];
    for f:=1 to High(ss) do
      if Result>ss[f] then Result:=ss[f];
    end;
end;

function Max(ss:TIntegers):Integer;overload;

var f:Integer;

begin
  if ss=nil then
    Result:=0
  else
    begin
    Result:=ss[0];
    for f:=1 to High(ss) do
      if Result<ss[f] then Result:=ss[f];
    end;
end;


function Max(ss:TSingles):Single;overload;

var f:Integer;

begin
  if ss=nil then
    Result:=0
  else
    begin
    Result:=ss[0];
    for f:=1 to High(ss) do
      if Result<ss[f] then Result:=ss[f];
    end;
end;


procedure ForceBetween(var i:Integer; low, upper:Integer);

begin
  if i<low then i:=low;
  if i>upper then i:=upper;
end;

procedure ForceBetween(var i:Single; low, upper:Single);

begin
  if i<low then i:=low;
  if i>upper then i:=upper;
end;


procedure Switch(var i1,i2:Integer);

var i:Integer;

begin
  i:=i1;i1:=i2;i2:=i;
end;

procedure Switch(var i1,i2:Single);

var i:Single;

begin
  i:=i1;i1:=i2;i2:=i;
end;

function Max(s1,s2:Integer):Integer;

begin
  if s1>s2 then Result:=s1 else Result:=s2;
end;

function Min(s1,s2:Integer):Integer;

begin
  if s1<s2 then Result:=s1 else Result:=s2;
end;

function Max(s1,s2:Single):Single;

begin
  if s1>s2 then Result:=s1 else Result:=s2;
end;

function Min(s1,s2:Single):Single;

begin
  if s1<s2 then Result:=s1 else Result:=s2;
end;


function Max(s1,s2:Double):Double;overload;
begin
  if s1>s2 then Result:=s1 else Result:=s2;
end;

function Min(s1,s2:Double):Double;overload;
begin
  if s1<s2 then Result:=s1 else Result:=s2;
end;

function IsBetween(db,d1,d2:Single):Boolean;

begin
  Result:=((db<=d2) and (db>=d1)) or
          ((db>=d2) and (db<=d1));
end;

function IsBetween(db,d1,d2:Integer):Boolean;

begin
  Result:=((db<=d2) and (db>=d1)) or
          ((db>=d2) and (db<=d1));
end;


end.


end.

