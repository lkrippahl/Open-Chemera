{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 30.5.2011
Purpose:
  Parser for CAPS results files
  (CAPS: Co-Evolution Analysis using Protein Sequences, Fares et al)
Requirements:
Revisions:
To do: Mostly everything. Still at quick fix stage...
*******************************************************************************}

unit capsparser;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, ocstringutils, debugutils;

type

  //TODO: this to be the real interface to the caps file
  TCAPSFile=class
  protected

  public

  end;

  //TODO: this is a quick-fix function, remove
  procedure QuickCAPSParser(FileName:string;out Q1,Q2:string;out IXs1,IXs2:TOCIntegers);


implementation

//TODO: this is a quick-fix function, remove
procedure QuickCAPSParser(FileName:string;out Q1,Q2:string;out IXs1,IXs2:TOCIntegers);

procedure ProcessIxs(S:string);

var
  f,ix:Integer;
  issecond:Boolean;

begin
  issecond:=False;
  while S<>'' do
    begin
    while (S<>'') and (S[1]<>'(') do
      begin
      if S[1]='[' then issecond:=True;
      Delete(S,1,1);
      end;
    Delete(S,1,1);
    if S<>'' then
      begin
      ix:=StrToInt(Copy(S,1,Pos(')',s)-1))-1;
      if issecond then AddToArray(ix,IXs2)
      else AddToArray(ix,IXs1);
      end;
    end;
end;

var
  sl:TStringList;
  f:Integer;

begin
  sl:=TStringList.Create;
  sl.LoadFromFile(FileName);
  f:=0;
  repeat Inc(f) until Pos('Amino acid sequence alignment',sl.Strings[f])=1;
  Inc(f,3);
  Q1:=CleanString(sl.Strings[f],'-');
  repeat Inc(f) until Pos('Amino acid sequence alignment',sl.Strings[f])=1;
  Inc(f,3);
  Q2:=CleanString(sl.Strings[f],'-');
  repeat Inc(f) until Pos('Groups of inter-molecular Co-evolving amino acids:',sl.Strings[f])=1;
  Inc(f,2);
  while f<sl.Count do
    begin
    ProcessIxs(sl.Strings[f]);
    Inc(f);
    end;
  sl.Free;
end;

end.

