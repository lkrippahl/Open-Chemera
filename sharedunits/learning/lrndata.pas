{*******************************************************************************
This source file is public domain. It is provided as is, with no explicit
or implied warranties regarding anything whatsoever.
********************************************************************************
Author: Ludwig Krippahl
Date: 23.12.2013
Purpose:
  Data management for machine learning tasks
Requirements:
Revisions:
To do: Still in experimental stage.
*******************************************************************************}

unit lrndata;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, stringutils, lrnconfig, filebuffer;

const
  //column types
  ColTypeFloat=0;
  ColTypeText=1;
  ColTypeInt=2;

  //column file versions
  ColNoneVersion=0;
  ColFLoatVersion=0;
  ColTextVersion=0;
  ColIntVersion=0;

  //table file versions
  BaseTableVersion=0;
  LearnTableVersion=0;
  DataManVersion=0;



type

{ TDataColumn }

TDataColumn=class
  protected
    FName:string;
  public
    property Name:string read FName write FName;
    function FloatVal(Ix:Integer):TFloat;virtual;abstract;
    function FloatVals:TFloats;virtual;abstract;
    function TextVal(Ix:Integer):string;virtual;abstract;
    function LastIx:Integer;virtual;abstract;
    function ColType:Integer;virtual;abstract;
    procedure SetColLength(Len: Integer);virtual;abstract;
    procedure WriteToBuffer(AFileBuffer:TFileBuffer);virtual;
    procedure ReadFromBuffer(AFileBuffer:TFileBuffer);virtual;
end;

TDataColumns=array of TDataColumn;

{ TTextDataColumn }

TTextDataColumn=class(TDataColumn)
private
  protected
    FVals:TSimpleStrings;
  public
    property Name:string read FName write FName;
    property Vals:TSimpleStrings read FVals write FVals;
    function FloatVals:TFloats;override;
    function FloatVal(Ix:Integer):TFloat;override;
    function TextVal(Ix:Integer):string;override;
    function LastIx:Integer;override;
    function ColType:Integer;override;
    procedure SetColLength(Len: Integer);override;
    procedure WriteToBuffer(AFileBuffer:TFileBuffer);override;
    procedure ReadFromBuffer(AFileBuffer:TFileBuffer);override;
end;

{ TFloatDataColumn }

TFloatDataColumn=class(TDataColumn)
private
  protected
    FVals:TFloats;
  public
    property Name:string read FName write FName;
    property Vals:TFloats read FVals write FVals;
    function FloatVal(Ix:Integer):TFloat;override;
    function FloatVals:TFloats;override;
    function TextVal(Ix:Integer):string;override;
    function LastIx:Integer;override;
    function ColType:Integer;override;
    procedure SetColLength(Len: Integer);override;
    procedure WriteToBuffer(AFileBuffer:TFileBuffer);override;
    procedure ReadFromBuffer(AFileBuffer:TFileBuffer);override;
end;

{ TIntegerDataColumn }

TIntegerDataColumn=class(TDataColumn)
private
  protected
    FVals:TIntegers;
  public
    property Name:string read FName write FName;
    property Vals:TIntegers read FVals write FVals;
    function FloatVal(Ix:Integer):TFloat;override;
    function FloatVals:TFloats;override;
    function TextVal(Ix:Integer):string;override;
    function LastIx:Integer;override;
    function ColType:Integer;override;
    procedure SetColLength(Len: Integer);override;
    procedure WriteToBuffer(AFileBuffer:TFileBuffer);override;
    procedure ReadFromBuffer(AFileBuffer:TFileBuffer);override;
end;



{ TDataTable }

TDataTable=class
  protected
    FColumns:TDataColumns;
    FLabels:TTextDataColumn;
    procedure ClearData;
  public
    Name:string;
    property Columns:TDataColumns read FColumns;
    constructor CreateFromFile(AName,FileName,Separator:string);
    procedure Free;virtual;
    procedure LoadSepFile(FileName,Separator:string;SkipLines:Integer=0;HasHeader:Boolean=True);virtual;
      //Loads file using Separator as column separator, but first skipping SkipLines
      //Uses TStrings as buffer; not appropriate for large files
      //Loads all data as text
      //if HasHeader uses first line for column names (after Skiplines)
    procedure ConvertToFloat(ColIx:Integer);overload;
    procedure ConvertToFloat(FromColIx,ToColIx:Integer);overload;
    procedure ConvertToFloat(ColIxs:TIntegers);overload;
    procedure DropColumns(ColNames:TSimpleStrings);overload;
    procedure DropColumns(ColIxs:TIntegers);overload;
    function ColumnIndex(ColName:string):Integer;
    function GetColumnIndexes(ColNames:TSimpleStrings):TIntegers;
    function Rows:Integer;
    procedure SetLabels(ColName:string);overload;
    procedure SetLabels(ColIx:Integer);overload;
    procedure WriteToBuffer(AFileBuffer:TFileBuffer);virtual;
    procedure ReadFromBuffer(AFileBuffer:TFileBuffer);virtual;
    function GetColumnNames:TSimpleStrings;overload;
    function GetColumnNames(Ixs:TIntegers):TSimpleStrings;overload;
    function GetLabels:TSimpleStrings;

end;

{ TLearnData }
//Assumes that all data columns are TFloatDataColumn

TLearnData=class(TDataTable)
  protected
    //Normalization values
    FSubtract,FDivide:TFloats;
  public
    //Classes for supervised learning
    //must start with 0 at least (non-negative)
    DataClasses:TIntegers;
    property Subtract:TFloats read FSubtract write FSubtract;
    property Divide:TFloats read FDivide write FDivide;
    procedure CalcAverageStDev;
      //Set FSubtract and FDivide constants to averages and standard devs
    procedure CalcMinusPlusOne;
      //Set FSubtract and FDivide constants to Min-1 and (Max-Min)/2, for scaling into [-1,1]
    procedure ScaleData;
      //Applies precomputed FSubtract and FDivide constants
    function GetColumn(ColName:string):TFloats;
      //Returns the pointer to the column values (not a copy)
    function CountClassPoints(DataClass:Integer):Integer;
      //DataClass <0 is all data points, same as Rows
    procedure WriteToBuffer(AFileBuffer:TFileBuffer);override;
    procedure ReadFromBuffer(AFileBuffer:TFileBuffer);override;
    procedure AddColumn(Col:TFloatDataColumn;SubVal,DivVal:TFloat);
end;

TLearnDataTables=array of TLearnData;

{ TDataManager }
//For most functions, tables need to have the same columns
//although not necessarily in the same order (the names must match)

TDataManager=class
  protected
    FTables:TLearnDataTables;
  public
    property Tables:TLearnDataTables read FTables;
    procedure ClearTables;
    procedure LoadTables(FileNames:TSimpleStrings;LabelCol:string='');
      //Loads tables, sets LabelCol column as labels, converts remainder to float
      { TODO : Only loads float data tables from tab separated files
        Generalize for other files? Or specific functions for each? }

    //procedure AddTable...

    procedure DropColumns(ColNames:TSimpleStrings);overload;
    procedure DropColumns(ColIxs:TIntegers);overload;
    procedure ScaleToMinusPlusOne;
    procedure ScaleToVectors(Subtract,Divide:TFloats);

    function SelectData(ColNames:TSimpleStrings;TableIxs:TIntegers;DataClass:Integer=-1):TMatrix;
      //Returns a copy of the selected data as a matrix
      //If DataClass<0 then all data is returned
      //first dimension of the matrix determines the data column
      //second dimension is the set of values for each colum
      //(this means that each column corresponds to a TFloats)
    procedure LoadFromFile(FileName:string);
    procedure SaveToFile(FileName:string);
      //uses one file for each table, FileName.ix with ix from zero to high(tables)
      //creates one catalog FileName with all file names

  end;


function CreateColumn(ColType:Integer):TDataColumn;


implementation

function CreateColumn(ColType: Integer): TDataColumn;
begin
  case ColType of
     ColTypeFloat:Result:=TFloatDataColumn.Create;
     ColTypeText:Result:=TTextDataColumn.Create;
     ColTypeInt:Result:=TIntegerDataColumn.Create;
     else Result:=nil;
  end;
end;

{ TDataColumn }

procedure TDataColumn.WriteToBuffer(AFileBuffer: TFileBuffer);
begin
  AFileBuffer.WriteInteger(ColNoneVersion);
  AFileBuffer.WriteString(Name);
end;

procedure TDataColumn.ReadFromBuffer(AFileBuffer: TFileBuffer);
begin
  AFileBuffer.GetInteger;//skiop ColNoneVersion;
  AFileBuffer.ReadString(FName);
end;

{ TDataManager }

procedure TDataManager.ClearTables;

var f:Integer;

begin
  for f:=0 to High(FTables) do FTables[f].Free;
  FTables:=nil;
end;

procedure TDataManager.LoadTables(FileNames: TSimpleStrings;LabelCol:string);

var f:Integer;

begin
  ClearTables;
  SetLength(FTables,Length(FileNames));
  for f:=0 to High(FileNames) do
    begin
    FTables[f]:=TLearnData.CreateFromFile(FileNames[f],FileNames[f],#9);
    if LabelCol<>'' then
      FTables[f].SetLabels(LabelCol);
    FTables[f].ConvertToFloat(0,High(FTables[f].Columns));
    end;
end;

procedure TDataManager.DropColumns(ColNames: TSimpleStrings);

var
  f:Integer;

begin
  for f:=0 to High(FTables) do
    FTables[f].DropColumns(ColNames);
end;

procedure TDataManager.DropColumns(ColIxs: TIntegers);

var
  f:Integer;

begin
  for f:=0 to High(FTables) do
    FTables[f].DropColumns(ColIxs);
end;

procedure TDataManager.ScaleToMinusPlusOne;

var
  mi,ma,s,d,tmp:TFloat;
  tmps:TFloats;
  colname:string;
  f,g,ix:Integer;


begin
  //Reset all table scale factors
  for f:=0 to High(FTables) do
    begin
    FTables[f].Subtract:=FilledFloats(Length(FTables[0].Columns),0);
    FTables[f].Divide:=FilledFloats(Length(FTables[0].Columns),1);
    end;
  //Compute scale factors
  for f:=0 to High(FTables[0].Columns) do
    begin
    colname:=FTables[0].Columns[f].Name;
    mi:=Min(TFloatDataColumn(FTables[0].Columns[f]).Vals);
    ma:=Max(TFloatDataColumn(FTables[0].Columns[f]).Vals);
    //Compute maximum and minimum
    for g:=1 to High(FTables) do
      begin
      tmps:=FTables[g].GetColumn(colname);
      tmp:=Min(tmps);
      if tmp<mi then mi:=tmp;
      tmp:=Max(tmps);
      if tmp>ma then ma:=tmp;
      end;
    d:=0.5*(ma-mi);
    s:=mi+d;

    //place in scaling factors
    FTables[0].Subtract[f]:=s;
    FTables[0].Divide[f]:=d;
    for g:=1 to High(FTables) do
      begin
      ix:=FTables[g].ColumnIndex(colname);
      FTables[g].Subtract[ix]:=s;
      FTables[g].Divide[ix]:=d;
      end;
    end;
  //Rescale
  for f:=0 to High(FTables) do
    FTables[f].ScaleData;

end;

procedure TDataManager.ScaleToVectors(Subtract, Divide: TFloats);
var
  mi,ma,s,d,tmp:TFloat;
  tmps:TFloats;
  colname:string;
  f,g,ix:Integer;


begin
  for f:=0 to High(FTables) do
    begin
    FTables[f].Subtract:=Copy(Subtract,0,Length(Subtract));
    FTables[f].Divide:=Copy(Divide,0,Length(Divide));
    FTables[f].ScaleData;
    end;
end;

function TDataManager.SelectData(ColNames: TSimpleStrings; TableIxs: TIntegers;
  DataClass: Integer): TMatrix;

var
  numrows,f,g,h,rix,cix:Integer;

begin
  numrows:=0;
  for f:=0 to High(TableIxs) do
    numrows:=numrows+FTables[TableIxs[f]].CountClassPoints(DataClass);
  SetLength(Result,Length(ColNames),numrows);

  for f:=0 to High(ColNames) do
    begin
    rix:=0;
    for g:=0 to High(TableIxs) do
      begin
      cix:=FTables[TableIxs[g]].ColumnIndex(ColNames[f]);
      for h:=0 to FTables[TableIxs[g]].Rows-1 do
         if (DataClass<0) or (DataClass=FTables[TableIxs[g]].DataClasses[h]) then
           begin
           Result[f,rix]:=TFloatDataColumn(FTables[TableIxs[g]].Columns[cix]).Vals[h];
           Inc(rix);
           end;
      end;
    end;
end;

procedure TDataManager.LoadFromFile(FileName: string);

var
  buffer:TFileBuffer;
  f:Integer;
  sl:TStringList;

begin
  ClearTables;
  sl:=TStringList.Create;
  sl.LoadFromFile(FileName);
  //skip DataManVersion on first line
  SetLength(FTables,sl.Count-1);
  for f:=0 to sl.Count-2 do
    begin
    buffer:=TFileBuffer.Create;
    if FileExists(sl.Strings[f+1]) then
       buffer.LoadFromFile(sl.Strings[f+1])
    else buffer.LoadFromFile(ExtractFilePath(FileName)+ExtractFileName(sl.Strings[f+1]));
    FTables[f]:=TLearnData.Create;
    FTables[f].ReadFromBuffer(buffer);
    buffer.Free;
    end;
  sl.Free;
end;


procedure TDataManager.SaveToFile(FileName: string);

var
  buffer:TFileBuffer;
  sl:TStringList;
  f:Integer;

begin
  sl:=TStringList.Create;
  sl.Add(IntToStr(DataManVersion));
  for f:=0 to High(FTables) do
    begin
    buffer:=TFileBuffer.Create;
    FTables[f].WriteToBuffer(buffer);
    sl.Add(FileName+'.'+IntToStr(f));
    buffer.FlushToFile(FileName+'.'+IntToStr(f));
    buffer.Free
    end;
  sl.SaveToFile(FileName);
  sl.Free;
end;

{ TLearnData }

procedure TLearnData.CalcAverageStDev;

var f:Integer;

begin
  SetLength(FSubtract,Length(FColumns));
  SetLength(FDivide,Length(FColumns));
  for f:=0 to High(FSubtract) do
    begin
    FSubtract[f]:=Average(FColumns[f].FloatVals);
    FDivide[f]:=Sqrt(Variance(FColumns[f].FloatVals,FSubtract[f]));
    end;

end;

procedure TLearnData.CalcMinusPlusOne;

var
  f:Integer;
  ma,mi:TFloat;

begin
  SetLength(FSubtract,Length(FColumns));
  SetLength(FDivide,Length(FColumns));
  for f:=0 to High(FSubtract) do
    begin
    mi:=Min(FColumns[f].FloatVals);
    ma:=Max(FColumns[f].FloatVals);
    FDivide[f]:=0.5*(ma-mi);
    FSubtract[f]:=mi+FDivide[f];
    end;
end;
procedure TLearnData.ScaleData;

var
  f,g:Integer;

begin
  for f:=0 to High(FColumns) do
    if FColumns[f].ColType=ColTypeFloat then
      for g:=0 to High(TFloatDataColumn(FColumns[f]).Vals) do
        TFloatDataColumn(FColumns[f]).Vals[g]:=
          (TFloatDataColumn(FColumns[f]).Vals[g]-FSubtract[f])/FDivide[f];
end;

function TLearnData.GetColumn(ColName:string): TFloats;

var ix:Integer;

begin
  ix:=ColumnIndex(ColName);
  Result:=TFloatDataColumn(FColumns[ix]).Vals;
end;

function TLearnData.CountClassPoints(DataClass: Integer): Integer;

var f:Integer;

begin
  if DataClass<0 then Result:=Rows
  else
    begin
    Result:=0;
    for f:=0 to High(DataClasses) do
      if DataClasses[f]=DataClass then
        Inc(Result);
    end;
end;

procedure TLearnData.WriteToBuffer(AFileBuffer: TFileBuffer);
begin
  inherited WriteToBuffer(AFileBuffer);
  AFileBuffer.WriteInteger(LearnTableVersion);
  AFileBuffer.WriteFloats(FSubtract);
  AFileBuffer.WriteFloats(FDivide);
  AFileBuffer.WriteIntegers(DataClasses);
end;

procedure TLearnData.ReadFromBuffer(AFileBuffer: TFileBuffer);
begin
  inherited ReadFromBuffer(AFileBuffer);
  AFileBuffer.GetInteger;//skip LearnTableVersion
  AFileBuffer.ReadFloats(FSubtract);
  AFileBuffer.ReadFloats(FDivide);
  AFileBuffer.ReadIntegers(DataClasses);
end;

procedure TLearnData.AddColumn(Col: TFloatDataColumn; SubVal, DivVal: TFloat);
begin
  SetLength(FColumns,Length(FColumns)+1);
  FColumns[High(FColumns)]:=Col;
  AddToArray(SubVal,FSubtract);
  AddToArray(DivVal,FDivide);
end;

{ TIntegerDataColumn }

function TIntegerDataColumn.FloatVal(Ix: Integer): TFloat;
begin
  Result:=FVals[Ix];
end;

function TIntegerDataColumn.FloatVals: TFloats;

var f:Integer;

begin
  SetLength(Result,Length(FVals));
  for f:=0 to High(Result) do
    Result[f]:=FVals[f];
end;

function TIntegerDataColumn.TextVal(Ix: Integer): string;
begin
  Result:=IntToStr(FVals[Ix]);
end;

function TIntegerDataColumn.LastIx: Integer;
begin
  Result:=High(Vals);
end;

function TIntegerDataColumn.ColType: Integer;
begin
  Result:=ColTypeInt;
end;

procedure TIntegerDataColumn.SetColLength(Len: Integer);
begin
  SetLength(FVals,Len);
end;

procedure TIntegerDataColumn.WriteToBuffer(AFileBuffer: TFileBuffer);
begin
  inherited WriteToBuffer(AFileBuffer);
  AFileBuffer.WriteInteger(ColIntVersion);
  AFileBuffer.WriteIntegers(FVals);

end;

procedure TIntegerDataColumn.ReadFromBuffer(AFileBuffer: TFileBuffer);
begin
  inherited ReadFromBuffer(AFileBuffer);
  AFileBuffer.GetInteger;//skip ColIntVersion;
  AFileBuffer.ReadIntegers(FVals);
end;

{ TDataTable }

procedure TDataTable.ClearData;

var f:Integer;

begin
  for f:=0 to High(FColumns) do
    FColumns[f].Free;
  FColumns:=nil;
  if FLabels<>nil then
    begin
    FLabels.Free;
    FLabels:=nil;
    end;
end;

constructor TDataTable.CreateFromFile(AName, FileName, Separator: string);
begin
  inherited Create();
  ClearData;
  Name:=AName;
  LoadSepFile(FileName,Separator);
end;

procedure TDataTable.Free;
begin
  ClearData;
  inherited Free;
end;

procedure TDataTable.LoadSepFile(FileName,Separator:string;SkipLines:Integer=0;HasHeader:Boolean=True);

var
  sl:TStringList;
  f,g,rix:Integer;
  ss:TSimpleStrings;

begin
  sl:=TStringList.Create;
  LrnReport('Loading '+FileName);
  sl.LoadFromFile(FileName);

  //Create columns
  ss:=SplitString(sl.Strings[SkipLines],Separator);
  if HasHeader then
    Inc(SkipLines); //skip header
  SetLength(FColumns,Length(ss));
  for f:=0 to High(FColumns) do
    begin
    FColumns[f]:=TTextDataColumn.Create;
    if HasHeader then
      FColumns[f].Name:=ss[f]
    else
      FColumns[f].Name:='Column '+IntToStr(f);
    FColumns[f].SetColLength(sl.Count-SkipLines);
    end;
  rix:=0;
  for f:=SkipLines to sl.Count-1 do
    begin
    ss:=SplitString(sl.Strings[f],Separator);
    if Length(ss)<>Length(FColumns) then
      raise Exception.Create('Mismatching rows on '+FileName);
    for g:=0 to High(FColumns) do
      TTextDataColumn(FColumns[g]).Vals[rix]:=ss[g];
    Inc(rix);
    end;
  sl.Free;
end;

procedure TDataTable.ConvertToFloat(ColIx: Integer);

var
  floatcol:TFloatDataColumn;
  floats:TFloats;
  f:Integer;
begin
  if FColumns[ColIx].ColType=ColTypeText then
    begin
    floats:=StringsToFloats(TTextDataColumn(FColumns[ColIx]).Vals);
    floatcol:=TFloatDataColumn.Create;
    floatcol.Name:=FColumns[ColIx].Name;
    floatcol.Vals:=floats;
    FColumns[ColIx].Free;
    FColumns[ColIx]:=floatcol;
    end
  else if FColumns[ColIx].ColType=ColTypeInt then
    begin
    floatcol:=TFloatDataColumn.Create;
    floatcol.Name:=FColumns[ColIx].Name;
    floatcol.SetColLength(FColumns[ColIx].LastIx+1);
    for f:=0 to High(floatcol.Vals) do
      floatcol.Vals[f]:=TIntegerDataColumn(FColumns[ColIx]).Vals[f];
    FColumns[ColIx].Free;
    FColumns[ColIx]:=floatcol;
    end
  else if FColumns[ColIx].ColType<>ColTypeFloat then
    raise Exception.Create('Cannot convert column '+FColumns[Colix].Name);
end;

procedure TDataTable.ConvertToFloat(FromColIx, ToColIx: Integer);

var f:Integer;

begin
  for f:=FromColIx to ToColIx do
    ConvertToFloat(f);
end;

procedure TDataTable.ConvertToFloat(ColIxs: TIntegers);

var f:Integer;

begin
  for f:=0 to High(ColIxs) do
    ConvertToFloat(ColIxs[f]);
end;

procedure TDataTable.DropColumns(ColNames: TSimpleStrings);

var
  ixs:TIntegers;
  f:Integer;

begin
  SetLength(ixs,Length(ColNames));
  for f:=0 to High(ixs) do
    ixs[f]:=ColumnIndex(ColNames[f]);
  DropColumns(ixs);
end;

procedure TDataTable.DropColumns(ColIxs: TIntegers);

var
  f,ix:Integer;

begin
  ix:=0;
  for f:=0 to High(FColumns) do
    begin
    if IsInArray(f,ColIxs) then
      FColumns[f].Free
    else
      begin
      if ix<f then
        FColumns[ix]:=FColumns[f];
      Inc(ix);
      end;
    end;
  SetLength(FColumns,ix);
end;

function TDataTable.ColumnIndex(ColName: string): Integer;

begin
  Result:=High(FColumns);
  while (Result>=0) and (FColumns[Result].Name<>ColName) do
    Dec(Result);
end;

function TDataTable.GetColumnIndexes(ColNames: TSimpleStrings): TIntegers;

var f:Integer;

begin
  SetLength(Result,Length(ColNames));
  for f:=0 to High(Result) do
    Result[f]:=ColumnIndex(ColNames[f]);
end;

function TDataTable.Rows: Integer;
begin
  if FColumns=nil then
    Result:=0
  else Result:=FColumns[0].LastIx+1;
end;

procedure TDataTable.SetLabels(ColName: string);
begin
  SetLabels(ColumnIndex(ColName));
end;

procedure TDataTable.SetLabels(ColIx: Integer);

var
  f:Integer;
  tmp:TIntegers;

begin
  if FLabels=nil then
    FLabels:=TTextDataColumn.Create;
  FLabels.SetColLength(FColumns[ColIx].LastIx+1);
  for f:=0 to High(FLabels.Vals) do
    FLabels.Vals[f]:=FColumns[ColIx].TextVal(f);
  tmp:=FilledInts(1,ColIx);
  DropColumns(tmp);
end;


procedure TDataTable.WriteToBuffer(AFileBuffer: TFileBuffer);

var f:Integer;

begin
  AFileBuffer.WriteInteger(BaseTableVersion);
  AFileBuffer.WriteString(Name);
  AFileBuffer.WriteInteger(Length(FColumns));
  for f:=0 to High(FColumns) do
    begin
    AFileBuffer.WriteInteger(FColumns[f].ColType);
    FColumns[f].WriteToBuffer(AFileBuffer);
    end;
  if FLabels=nil then
    AFileBuffer.WriteBoolean(False)
  else
    begin
    AFileBuffer.WriteBoolean(True);
    FLabels.WriteToBuffer(AFileBuffer);
    end;
end;

procedure TDataTable.ReadFromBuffer(AFileBuffer: TFileBuffer);

var f:Integer;

begin
  ClearData;
  AFileBuffer.GetInteger;//skip BaseTableVersion);
  AFileBuffer.ReadString(Name);
  SetLength(FColumns,AFileBuffer.GetInteger);
  for f:=0 to High(FColumns) do
    begin
    FColumns[f]:=CreateColumn(AFileBuffer.GetInteger);
    FColumns[f].ReadFromBuffer(AFileBuffer);
    end;
  if AFileBuffer.GetBoolean then
    begin
    FLabels:=TTextDataColumn.Create;
    FLabels.ReadFromBuffer(AFileBuffer);
    end;
end;

function TDataTable.GetColumnNames: TSimpleStrings;

var f:Integer;

begin
  SetLength(Result,Length(FColumns));
  for f:=0 to High(Result) do Result[f]:=FColumns[f].Name;
end;

function TDataTable.GetColumnNames(Ixs: TIntegers): TSimpleStrings;

var f:Integer;

begin
  SetLength(Result,Length(Ixs));
  for f:=0 to High(Result) do Result[f]:=FColumns[Ixs[f]].Name;
end;

function TDataTable.GetLabels: TSimpleStrings;
begin
  Result:=Copy(FLabels.Vals,0,Length(FLabels.Vals));
end;


{ TFloatDataColumn }

function TFloatDataColumn.FloatVals: TFloats;
begin
  Result:=FVals;
end;

function TFloatDataColumn.FloatVal(Ix: Integer): TFloat;
begin
  Result:=FVals[Ix];
end;

function TFloatDataColumn.TextVal(Ix: Integer): string;
begin
  Result:=FloatToStr(FVals[Ix]);
end;

function TFloatDataColumn.LastIx: Integer;
begin
  Result:=High(FVals);
end;

function TFloatDataColumn.ColType: Integer;
begin
  Result:=ColTypeFloat;
end;

procedure TFloatDataColumn.SetColLength(Len: Integer);
begin
  SetLength(FVals,Len);
end;

procedure TFloatDataColumn.WriteToBuffer(AFileBuffer: TFileBuffer);
begin
  inherited WriteToBuffer(AFileBuffer);
  AFileBuffer.WriteInteger(ColFloatVersion);
  AFileBuffer.WriteFloats(FVals);
end;

procedure TFloatDataColumn.ReadFromBuffer(AFileBuffer: TFileBuffer);
begin
  inherited ReadFromBuffer(AFileBuffer);
  AFileBuffer.GetInteger;//skip ColFloatVersion;
  AFileBuffer.ReadFloats(FVals);

end;

{ TTextDataColumn }

function TTextDataColumn.FloatVals: TFloats;
begin
  Result:=StringsToFloats(FVals);
end;

function TTextDataColumn.FloatVal(Ix: Integer): TFloat;
begin
  Result:=StrToFloat(FVals[Ix]);
end;

function TTextDataColumn.TextVal(Ix: Integer): string;
begin
  Result:=FVals[Ix];
end;

function TTextDataColumn.LastIx: Integer;
begin
  Result:=High(FVals);
end;

function TTextDataColumn.ColType: Integer;
begin
  Result:=ColTypeText;
end;

procedure TTextDataColumn.SetColLength(Len: Integer);
begin
  SetLength(FVals,Len);
end;

procedure TTextDataColumn.WriteToBuffer(AFileBuffer: TFileBuffer);
begin
  inherited WriteToBuffer(AFileBuffer);
  AFileBuffer.WriteInteger(ColTextVersion);
  AFileBuffer.WriteStrings(FVals);

end;

procedure TTextDataColumn.ReadFromBuffer(AFileBuffer: TFileBuffer);
begin
  inherited ReadFromBuffer(AFileBuffer);
  AFileBuffer.GetInteger;//skip ColTextVersion;
  AFileBuffer.ReadStrings(FVals);
end;

end.

