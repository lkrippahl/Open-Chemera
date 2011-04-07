{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 7.4.2011
Purpose:
  Utilities for handling simple strings
Requirements:
Revisions:
To do:
*******************************************************************************}

unit stringutils;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes;

function SnipString(var S:string;const Spacer:string):string;
  //cuts S to first spacer, removing spacer, returns part before spacer
  //returns empty string and leaves S if spacer is not found

function TrimmedBlanks(const S:string):string;
  //removes all <=' ' characters from the extremities of the string;

implementation

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

end.

