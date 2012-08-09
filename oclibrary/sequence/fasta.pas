{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 9.1.2011
Purpose:
  Parser for FASTA
Requirements:
Revisions:
To do:
  Merge with alignment unit?
*******************************************************************************}
unit fasta;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils;

type
  TFastaData=record
    ID:string;
    Sequence:string;
  end;
  TFastaSeqs = array of TFastaData;

  { TFastaReader }

  TFastaReader=class
  private
    FSeqs:TFastaSeqs;
    procedure AddSeq(seq:TFastaData);
  public
    property Seqs:TFastaSeqs read FSeqs;
    constructor Create(fromfile:string);
    procedure Clear;
    procedure Free;
    procedure Load(fromfile:string);
  end;

implementation

{ TFastaReader }

procedure TFastaReader.AddSeq(seq: TFastaData);
begin
  SetLength(FSeqs,Length(FSeqs)+1);
  FSeqs[High(FSeqs)]:=seq;
end;

constructor TFastaReader.Create(fromfile: string);
begin
  inherited Create;
  Load(fromfile);
end;

procedure TFastaReader.Clear;
begin
  FSeqs:=nil;
end;

procedure TFastaReader.Free;
begin
  if Self<>nil then
    begin
    Clear;
    inherited Free;
    end;
end;

procedure TFastaReader.Load(fromfile: string);

var
  buf:TStringList;
  cur:TFastaData;
  f:Integer;
  s:string;

begin
  cur.ID:='';
  cur.Sequence:='';
  buf:=TStringList.Create;
  buf.LoadFromFile(fromfile);
  for f:=0 to buf.Count-1 do
    begin
    s:=buf.Strings[f];
    if Pos('>',s)=1 then
      begin
      if cur.ID<>'' then
        AddSeq(cur);
      cur.ID:=Copy(s,2,Length(s));
      cur.Sequence:='';
      end
    else cur.Sequence:=cur.Sequence+s;
    end;
  if cur.Sequence<>'' then AddSeq(cur);
  buf.Free;
end;

end.

