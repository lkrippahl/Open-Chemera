{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 24.4.2011
Purpose:
  residue interactions (amino acids, for now). Includes Mutual Information
  and other measures for coevolution
Requirements:
  Uses alignment because interaction matrices are the same as substitution
  matrices (TSubMatrix)
Revisions:
To do:
*******************************************************************************}

unit resinteraction;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, alignment, basetypes;

function AverageInteraction(const Res1,Res2:string; const IntMat:TSubMatrix;GapMarker:Char):TOCFloat;
function MutualInformation(const s1, s2: string; ScaleMi:Boolean=False): TOCFloat;



implementation

function AverageInteraction(const Res1,Res2:string; const IntMat:TSubMatrix; GapMarker:Char):TOCFloat;
//assumes Res strings have the same length

var f,c:Integer;

begin
  Result:=0;
  c:=0;
  for f:=1 to Length(Res1) do
    if (Res1[f]<>GapMarker) and (Res2[f]<>GapMarker) then
      begin
      Result:=Result+GetSubScore(IntMat,Res1[f],Res2[f]);
      Inc(c);
      end;
  if c>0 then Result:=Result/c;
end;

function MutualInformation(const s1, s2: string; ScaleMi:Boolean=False): TOCFloat;

const
  Ln2=0.69314718056; //ln(2) to convert to log2

var
  f,g:Integer;
  symb1,symb2:string;
  step,min,mult:TOCFloat;
  marg1,marg2:array of TOCFloat;
  cross:array of array of TOCFLoat;

procedure BuildSymbols;

var f,g:Integer;

begin
  symb1:='';
  symb2:='';
  for f:=1 to Length(s1) do
    if Pos(s1[f],symb1)<1 then symb1:=symb1+s1[f];
  for f:=1 to Length(s2) do
    if Pos(s2[f],symb2)<1 then symb2:=symb2+s2[f];
  SetLength(marg1,Length(symb1));
  SetLength(marg2,Length(symb2));
  SetLength(cross,Length(symb1),Length(symb2));
  for f:=0 to Length(symb1)-1 do
    marg1[f]:=0;
  for f:=0 to Length(symb2)-1 do
    marg2[f]:=0;
  for f:=0 to High(cross) do
    for g:=0 to High(cross[0]) do
      cross[f,g]:=0;
end;

procedure Count;

var f,i1,i2:Integer;

begin
  step:=1/Length(s1);
  for f:=1 to Length(s1) do
    begin
    i1:=Pos(s1[f],symb1)-1;
    i2:=Pos(s2[f],symb2)-1;
    marg1[i1]:=marg1[i1]+step;
    marg2[i2]:=marg2[i2]+step;
    cross[i1,i2]:=cross[i1,i2]+step;
    end;
end;

begin
  Result:=0;
  BuildSymbols;
  Count;
  min:=step/2;
  for f:=0 to High(cross) do
    for g:=0 to High(cross[0]) do
      begin
      mult:=marg1[f]*marg2[g];
      if cross[f,g]>min then
        Result:=Result+cross[f,g]*Ln(cross[f,g]/mult);
      end;
  Result:=Result/Ln2;
  if ScaleMI then Result:=Result*Ln(Sqrt(Length(symb1)*Length(symb2)));
end;

end.

