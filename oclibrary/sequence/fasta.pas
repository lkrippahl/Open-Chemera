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
  Classes, SysUtils, sequence;

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
    function ExtractID(IdLine:string):string;
    function ExtractOrganism(IdLine:string):string;

  public
    property Seqs:TFastaSeqs read FSeqs;
    constructor Create(fromfile:string);
    procedure Clear;
    procedure Free;
    procedure Load(fromfile:string);

    function AsOCSequences: TOCSequences;
    function SequenceByID(ID:string):string;
  end;

  { TFastaWriter }

  TFastaWriter=class
  private
    FSeqs:TFastaSeqs;
  public
    procedure AddSeq(Seq:TFastaData);overload;
    procedure AddSeq(ID,Sequence:string);overload;
    procedure AddSeq(Seq:TOCSequence;IncludeOrganism:Boolean=False);overload;
    procedure AppendSeqs(Seqs:TOCSequences;IncludeOrganism:Boolean=False);
    procedure SaveToFile(FileName:string;MaxCharsPerLine:Integer=70);
    procedure Clear;
    procedure Free;
  end;

implementation

{ TFastaWriter }

procedure TFastaWriter.AddSeq(Seq: TFastaData);
begin
  SetLength(FSeqs,Length(FSeqs)+1);
  FSeqs[High(FSeqs)]:=Seq;
end;

procedure TFastaWriter.AddSeq(ID, Sequence: string);
begin
  SetLength(FSeqs,Length(FSeqs)+1);
  FSeqs[High(FSeqs)].ID:=ID;
  FSeqs[High(FSeqs)].Sequence:=Sequence;
end;

procedure TFastaWriter.AddSeq(Seq: TOCSequence;IncludeOrganism:Boolean=False);
begin
  SetLength(FSeqs,Length(FSeqs)+1);
  FSeqs[High(FSeqs)].ID:=Seq.ID;
  if IncludeOrganism then
    FSeqs[High(FSeqs)].ID:=FSeqs[High(FSeqs)].ID+'|os='+Seq.Organism;
  FSeqs[High(FSeqs)].Sequence:=Seq.Sequence;
end;

procedure TFastaWriter.AppendSeqs(Seqs: TOCSequences;IncludeOrganism:Boolean=False);

var ol,f:Integer;

begin
  ol:=Length(FSeqs);
  SetLength(FSeqs,ol+Length(Seqs));
  for f:=0 to High(Seqs) do
    begin
    FSeqs[ol+f].ID:=Seqs[f].ID;
    if IncludeOrganism then
      FSeqs[ol+f].ID:=FSeqs[ol+f].ID+'|os='+Seqs[f].Organism;
    FSeqs[ol+f].Sequence:=Seqs[f].Sequence;
    end;
end;

procedure TFastaWriter.SaveToFile(FileName: string; MaxCharsPerLine:Integer=70);

var
  sl:TStringList;
  f:Integer;
  tmp:string;

begin
  sl:=TStringList.Create;
  for f:=0 to High(FSeqs) do
    begin
    tmp:=Copy(FSeqs[f].ID,1,MaxCharsPerLine-1);
    sl.Add('>'+tmp);
    tmp:=FSeqs[f].Sequence;
    while tmp<>'' do
      begin
      sl.Add(Copy(tmp,1,MaxCharsPerLine));
      Delete(tmp,1,MaxCharsPerLine);
      end;
    end;
  sl.SaveToFile(FileName);
  sl.Free;
end;

procedure TFastaWriter.Clear;
begin
  FSeqs:=nil;
end;

procedure TFastaWriter.Free;
begin
  Clear;
  inherited;
end;

{ TFastaReader }

procedure TFastaReader.AddSeq(seq: TFastaData);
begin
  SetLength(FSeqs,Length(FSeqs)+1);
  FSeqs[High(FSeqs)]:=seq;
end;

function TFastaReader.ExtractID(IdLine: string): string;

var ix:Integer;

begin
  ix:=Pos(' ',IdLine);
  if ix>0 then Result:=Copy(IdLine,1,ix-1)
  else Result:=IdLine;
end;

function TFastaReader.ExtractOrganism(IdLine: string): string;

var ix:Integer;

begin
  ix:=Pos('OS="',UpperCase(IdLine));
  if ix>0 then
      begin
      Result:=Copy(IdLine,ix+4,Length(IdLine));
      Delete(Result,Pos('"',Result),Length(Result));
      end
  else
    begin
    ix:=Pos('[ORGANISM=',UpperCase(IdLine));
    if ix>0 then
      begin
      Result:=Copy(IdLine,ix+10,Length(IdLine));
      Delete(Result,Pos(']',Result),Length(Result));
      end
    else Result:='';
    end
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

function TFastaReader.AsOCSequences: TOCSequences;

var f:Integer;

begin
  SetLength(Result,Length(FSeqs));
  for f:=0 to High(Result) do
    begin
    Result[f]:=EmptySequence;
    Result[f].ID:=ExtractId(FSeqs[f].ID);
    Result[f].Organism:=ExtractOrganism(FSeqs[f].ID);
    Result[f].Sequence:=FSeqs[f].Sequence;
    end;
end;

function TFastaReader.SequenceByID(ID: string): string;

var
  f:Integer;

begin
  Result:='';
  for f:=0 to High(FSeqs) do
    if FSeqs[f].ID=ID then
      begin
      Result:=FSeqs[f].Sequence;
      Break;
      end;
end;

end.

