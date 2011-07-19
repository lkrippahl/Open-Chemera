{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 23.6.2011
Purpose:
  Parser for MSA blocks in refinement.
Requirements:
Revisions:
To do:
  This implementation is still provisional, just a quick fix for current block
  formats.
  We sould settle on block file format first...
*******************************************************************************}

unit msablocksparser;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, alignment, basetypes;

type
  TBlock=record
    Seqs:TOCStrings;
  end;
  TBlockSet=record
    Blocks:array of TBlock;
    Positions:TOCIntegers;
  end;

function LoadBlocks(const BlockFile:string):TBlockSet;
procedure MergeBlocks(var MSA:TMSA; BlockSet:TBlockSet);
function ShiftToSequence(const MSA:TMSA; const BlockSet:TBlockSet; const SeqIx:Integer):TBlockSet;
  //Returns a blockset shifted to start on positions relative to the sequence SeqIx,
  //with 1 corresponding to the first elemento of the sequence
  //note that gaps are ignored



implementation

function EmptyBlockSet:TBlockSet;

begin
  Result.Blocks:=nil;
  Result.Positions:=nil;
end;

function CopyBlockSet(FromBlockSet:TBlockSet):TBlockSet;

var f:Integer;

begin
  SetLength(Result.Blocks,Length(FromBlockSet.Blocks));
  for f:=0 to High(Result.Blocks) do
    Result.Blocks[f].Seqs:=Copy(FromBlockSet.Blocks[f].Seqs,0,Length(FromBlockSet.Blocks[f].Seqs));
  Result.Positions:=Copy(FromBlockSet.Positions,0,Length(FromBlockSet.Positions));
end;

function LoadBlocks(const BlockFile:string):TBlockSet;

var
  curblock,f:Integer;
  sl:TStringList;
  s:string;

procedure FixBlockFormat;

var
  g:Integer;

begin
  g:=0;
  while g<sl.Count do
    begin
    if pos('#',sl.Strings[g])=1 then sl.Delete(g)
    else if sl.Strings[g]='***Block***' then
      begin
      sl.Strings[g]:='BLOCK:'+Copy(sl.Strings[g+1],Pos(':',sl.Strings[g+1])+1,Length(sl.Strings[g+1]));
      sl.Delete(g+1);
      sl.Delete(g+2);
      Inc(g);
      end
    else if Pos(';',sl.Strings[g])>0 then
      sl.Strings[g]:=Copy(sl.Strings[g],1,Pos(';',sl.Strings[g])-1)
    else Inc(g);
    end;
end;

begin
  Result:=EmptyBlockSet;
  curblock:=-1;
  sl:=TStringList.Create;
  sl.LoadFromFile(BlockFile);
  FixBlockFormat;
  for f:=0 to sl.Count-1 do
    begin
    s:=sl.Strings[f];
    if Pos('BLOCK:',s)=1 then
      begin
      Inc(curblock);
      SetLength(Result.Blocks,curblock+1);
      Result.Blocks[curblock].Seqs:=nil;
      SetLength(Result.Positions,curblock+1);
      Result.Positions[curblock]:=StrToInt(Copy(s,7,Length(s)));
      end
    else if s='END' then
      Break              //there should be no blocks only 3 wide
    else
      if (s<>'') and (Pos('No solution!',s)<1) then
        AddToArray(s,Result.Blocks[curblock].Seqs);
    end;
  sl.Free;
end;

procedure MergeBlocks(var MSA:TMSA; BlockSet:TBlockSet);
//TO DO: error checking!
// also, perhaps best to merge the blockset and not from the file name (or overload)

var
  f,g,h,p0:Integer;
  s,t:string;
begin
  for f:=0 to High(MSA.Alignment) do
    begin
    s:=MSA.Alignment[f];
    for g:=0 to High(BlockSet.Blocks) do
      begin
      p0:=BlockSet.Positions[g]-1;
      t:=BlockSet.Blocks[g].Seqs[f];
      for h:=1 to Length(t) do
        s[h+p0]:=t[h];
      end;
    MSA.Alignment[f]:=s;
    end;
end;

function ShiftToSequence(const MSA:TMSA; const BlockSet:TBlockSet; const SeqIx:Integer):TBlockSet;

var
  curix,f:Integer;
  s:string;
  index:TOCIntegers;

begin
  Result:=CopyBlockSet(BlockSet);
  s:=MSA.Alignment[SeqIx];
  //index the elements of the string in the alignment
  curix:=0;
  SetLength(index,Length(s));
  for f:=0 to High(Index) do
    begin
    if s[f+1]<>MSA.GapMarker then Inc(curix);
    index[f]:=curix;
    end;
  for f:=0 to High(Result.Positions) do
    Result.Positions[f]:=index[Result.Positions[f]];
end;

end.

