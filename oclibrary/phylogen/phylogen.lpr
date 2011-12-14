{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain. Enjoy.
********************************************************************************
Author: Ludwig Krippahl
Date: 3.12.2011
Purpose:
  Sequence comparison and phylogenetic analysis. Currently only compares trees
  from matching organisms
Requirements:
  LCL for some units (e.g. FileUtil)
Revisions:
To do:

*******************************************************************************}
program phylogen;

{$mode objfpc}{$H+}

uses
  {$IFDEF UNIX}{$IFDEF UseCThreads}
  cthreads,
  {$ENDIF}{$ENDIF}
  Classes, SysUtils, CustApp, statistics, basetypes, stringutils, quicksort,
  coevolution;

type

  { TPhylogen }

  TPhylogen = class(TCustomApplication)
  protected
    procedure DoRun; override;
    procedure LoadMatrices(FileName:string; out M1, M2: TMatrix);
    function CorrelationSqr(M1,M2:TMatrix;Columns:TIntegers):TFLoats;
    function FindMinimum(M1,M2:TMatrix;Columns:TIntegers; out CVec:TFloats):Integer;
    procedure CompareTrees(FileName:string);
    procedure TheilSen(FileName:string);
  public
    constructor Create(TheOwner: TComponent); override;
    destructor Destroy; override;
    procedure WriteHelp; virtual;
  end;

{ TPhylogen }

procedure TPhylogen.DoRun;
var
  fn, ErrorMsg: String;

begin
  // quick check parameters
  ErrorMsg:=CheckOptions('h f t','help file theilesen');
  if ErrorMsg<>'' then begin
    ShowException(Exception.Create(ErrorMsg));
    Terminate;
    Exit;
  end;

  // parse parameters
  if HasOption('h','help') then begin
    WriteHelp;
    Terminate;
    Exit;
  end;

  { add your program here }
  if HasOption('f','file') then
    begin
    fn:=GetOptionValue('f','file');
    CompareTrees(fn);
    end;

  if HasOption('t','theilsen') then
    begin
    fn:=GetOptionValue('t','theilsen');
    TheilSen(fn);
    end;


  // stop program loop
  Terminate;
end;

procedure TPhylogen.LoadMatrices(FileName: string; out M1, M2: TMatrix);
//Assumes text files with IX1, tab, ix2, tab, val1, tab, val2

procedure DecodeLine(Line:string;out Ix1,Ix2:Integer; out F1,F2:TFloat);

var
  ss:TSimpleStrings;
  f:Integer;

begin
  DecimalSeparator:='.';
  ss:=SplitString(Line,#9);
  Ix1:=StrToInt(ss[0]);
  Ix2:=StrToInt(ss[1]);
  F1:=StrToFloat(ss[2]);
  F2:=StrToFloat(ss[3]);
end;


var
  f:Integer;
  sl:TStringList;
  s:string;
  i1,i2:Integer;
  val1,val2:TFloat;

begin
  sl:=TStringList.Create;
  sl.LoadFromFile(FileName);
  while (sl.Count>0) and (sl.Strings[sl.Count-1]='') do
    sl.Delete(sl.Count-1);
  s:=sl.Strings[sl.Count-1];
  DecodeLine(s,i1,i2,val1,val2);
  SetLength(M1,i1,i2);
  SetLength(M2,i1,i2);
  for f:=0 to sl.Count-1 do
    begin
    s:=sl.Strings[f];
    DecodeLine(s,i1,i2,val1,val2);
    Dec(i1);
    Dec(i2);
    M1[i1,i2]:=val1;
    M2[i1,i2]:=val2;
    end;
end;

function TPhylogen.CorrelationSqr(M1, M2: TMatrix; Columns: TIntegers): TFLoats;

var
  f,g:Integer;
  vec1,vec2:TFloats;

begin
  SetLength(vec1,Length(Columns));
  SetLength(vec2,Length(Columns));
  SetLength(Result,Length(Columns));
  for f:=0 to High(Columns) do
    begin
    for g:=0 to High(Columns) do
      begin
      vec1[g]:=M1[Columns[f],Columns[g]];
      vec2[g]:=M2[Columns[f],Columns[g]];
      end;
    Result[f]:=Sqr(Pearson(vec1,vec2));
    end;
end;

function TPhylogen.FindMinimum(M1, M2: TMatrix; Columns: TIntegers;
  out CVec: TFloats): Integer;
begin
  CVec:=CorrelationSqr(M1,M2,Columns);
  Result:=MinIx(CVec);
end;

procedure TPhylogen.CompareTrees(FileName: string);

var
  mat1,mat2:TMatrix;
  current:TIntegers;
  f:Integer;
  correlvec:TFloats;

begin
  LoadMatrices(FileName,mat1,mat2);
  //set initial list, with all
  SetLength(current,Length(mat1));
  for f:=0 to High(Current) do current[f]:=f;

  while current<>nil do
    begin
    f:=FindMinimum(mat1,mat2,current,correlvec);
    WriteLn(current[f]+1,#9,correlvec[f]);
    RemoveFromArray(f,current);
    end;
end;

procedure TPhylogen.TheilSen(FileName: string);

var
  mat1,mat2:TMatrix;
  correlator:TTwoProteinCorrelation;
  distances:TFloats;
  slope:TFloat;
  f:Integer;

begin
  LoadMatrices(FileName,mat1,mat2);
  correlator:=TTwoProteinCorrelation.Create(mat1,mat2);
  correlator.TheilSenToQuery(0,slope,distances);
  WriteLn('Slope:'#9,slope);
  WriteLn('Distances:');
  for f:=0 to High(distances) do
      WriteLn(f,#9,distances[f]);
end;

constructor TPhylogen.Create(TheOwner: TComponent);
begin
  inherited Create(TheOwner);
  StopOnException:=True;
end;

destructor TPhylogen.Destroy;
begin
  inherited Destroy;
end;

procedure TPhylogen.WriteHelp;
begin
  { add your help code here }
  WriteLn('Usage: ',ExeName,' -f matrixfile');
end;

var
  Application: TPhylogen;

{$R *.res}

begin
  Application:=TPhylogen.Create(nil);
  Application.Title:='Phylogen';
  Application.Run;
  Application.Free;
end.
