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
  This is preliminary work. Block files should include the sequence IDs
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
procedure MergeBlocks(var MSA:TMSA; BlockFile:string);

implementation

function EmptyBlockSet:TBlockSet;

begin
  Result.Blocks:=nil;
  Result.Positions:=nil;
end;

function LoadBlocks(const BlockFile:string):TBlockSet;

var
  curblock,f:Integer;
  sl:TStringList;
  s:string;

begin
  Result:=EmptyBlockSet;
  curblock:=-1;
  sl:=TStringList.Create;
  sl.LoadFromFile(BlockFile);
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
    else if Pos('END',s)=1 then
      Break
    else
      begin
      AddToArray(s,Result.Blocks[curblock].Seqs);
      end;
    end;
  sl.Free;
end;

procedure MergeBlocks(var MSA:TMSA; BlockFile:string);
//TO DO: error checking!
// also, perhaps best to merge the blockset and not from the file name (or overload)

var
  blocks:TBlockSet;
  f,g,h,p0:Integer;
  s,t:string;
begin
  blocks:=LoadBlocks(BlockFile);
  for f:=0 to High(MSA.Alignment) do
    begin
    s:=MSA.Alignment[f];
    for g:=0 to High(blocks.Blocks) do
      begin
      p0:=blocks.Positions[g]-1;
      t:=blocks.Blocks[g].Seqs[f];
      for h:=1 to Length(t) do
        s[h+p0]:=t[h];
      end;
    MSA.Alignment[f]:=s;
    end;
end;

end.

