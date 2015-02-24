{*******************************************************************************
This source file is public domain. It is provided as is, with no explicit
or implied warranties regarding anything whatsoever.
********************************************************************************
Author: Ludwig Krippahl
Date: 26.12.2013
Purpose:
  Histogram management
Requirements:
Revisions:
To do:
  KDE on TFloatHistogram
  see http://www.mathworks.com/matlabcentral/fileexchange/14034-kernel-density-estimator
*******************************************************************************}
unit lrnhistogram;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, filebuffer;

const
  HistogramVersion=0;
  FloatHistogramVersion=0;

  HistTypeNone=0;
  HistTypeFloat=1;

type
  THistogram=class;
  THistograms=array of THistogram;
  THistMatrix=array of THistograms;

  { THistogram }

  THistogram=class
    protected
      FHistogram:TFloats;
      FLabels:TSimpleStrings;
      FDataCount:Integer; //Data is not loaded nor saved, but the number of points may be useful
      FSmoothing:TFloat;  //Smoothing factor used in CalcHistogram
      FIsLog:Boolean;
      FIsNormalized:Boolean;
    public
      Name:string;
      function Count:Integer;
      procedure CalcHistogram(NumBins:Integer=-1;Smoothing:TFloat=0);virtual;
        //-1 bins means auto binning
        //currently uses Rice rule Ceil(2*n^(1/3))
        { TODO : Implement   "A method for selecting the bin size of a time histogram"
          Shimazaki, Hideaki and Shinomoto, Shigeru, Neural Computatiom 2007 }
      procedure Normalize;
        //makes the sum of the bins equal to 1
      procedure ConvertToLog;
      procedure ConvertToP;
      procedure ReportHistogram(Sl:TStringList;WriteLabels:Boolean=False);virtual;abstract;
      procedure WriteToBuffer(ABuffer:TFileBuffer);virtual;
      procedure ReadFromBuffer(ABuffer:TFileBuffer);virtual;
      function HistogramType:Integer;virtual;
      function BinCount:Integer;
    end;



  { TFloatHistogram }
  TFloatHistogram=class(THistogram)
    protected
      FData:TFloats;
      FMin,FMax,FRange:TFloat;
      procedure SetFMin(Val:TFloat);
      procedure SetFMax(Val:TFloat);
    public
      property MinVal:TFloat read FMin write SetFMin;
      property MaxVal:TFloat read FMax write SetFMax;
      property Data:TFloats read FData write FData;
        //ATTENTION: does not manage data and data is not saved nor loaded
        //Only necessary to set up before computing histogram
      constructor Create(SourceData:TFloats);
      procedure CalcHistogram(NumBins:Integer=100;Smoothing:TFloat=0);override;
      procedure WriteToBuffer(ABuffer:TFileBuffer);override;
      procedure ReadFromBuffer(ABuffer:TFileBuffer);override;
      function HistogramType:Integer;override;
      function BinVal(X:TFLoat):TFloat;
      function BinVals(X:TFLoats):TFloats;
      procedure ReportHistogram(Sl:TStringList;WriteLabels:Boolean=False);override;
  end;


function CreateHistogram(HistType:Integer):THistogram;
procedure SaveHistMatrix(Hists:THistMatrix;FileName:string);
procedure LoadHistMatrix(out Hists:THistMatrix;FileName:string);
procedure FreeHistMatrix(var Hists:THistMatrix);


implementation

function CreateHistogram(HistType: Integer): THistogram;
begin
  case HistTYpe of
    HistTypeNone: Result:=THistogram.Create();
    HistTypeFloat: Result:=TFloatHistogram.Create(nil);
    else Result:=nil;
  end;
end;

procedure SaveHistMatrix(Hists: THistMatrix; FileName: string);

var
  f,g:Integer;
  buffer:TFileBuffer;

begin
  buffer:=TFileBuffer.Create;
  buffer.WriteInteger(Length(Hists));
  if Hists<>nil then
    for f:=0 to High(Hists) do
      begin
      buffer.WriteInteger(Length(Hists[f]));
      for g:=0 to High(Hists[f]) do
        begin
        buffer.WriteInteger(Hists[f,g].HistogramType);
        Hists[f,g].WriteToBuffer(buffer);
        end;
      end;
  buffer.FlushToFile(FileName);
  buffer.Free;
end;

procedure LoadHistMatrix(out Hists: THistMatrix; FileName: string);

var
  f,g:Integer;
  buffer:TFileBuffer;

begin
  buffer:=TFileBuffer.Create;
  WriteLn(FileName);
  buffer.LoadFromFile(FileName);
  SetLength(Hists,buffer.GetInteger);
  if Hists<>nil then
    for f:=0 to High(Hists) do
      begin
      SetLength(Hists[f],buffer.GetInteger);
      for g:=0 to High(Hists[f]) do
        begin
        Hists[f,g]:=CreateHistogram(buffer.GetInteger);
        Hists[f,g].ReadFromBuffer(buffer);
        end;
      end;
  buffer.FlushToFile(FileName);
  buffer.Free;
end;

procedure FreeHistMatrix(var Hists: THistMatrix);

var f,g,h:Integer;

begin
  for f:=0 to High(Hists) do
    for g:=0 to High(Hists[f]) do
      Hists[f,g].Free;
  Hists:=nil;
end;

{ THistogram }

function THistogram.Count: Integer;
begin
  Result:=FDataCount;
end;

procedure THistogram.CalcHistogram(NumBins: Integer; Smoothing: TFloat);
begin
  if NumBins<0 then
    NumBins:=Trunc(2*Exp(Ln(FDataCount)/3)+1);
  FHistogram:=FilledFloats(NumBins,Smoothing);
  FSmoothing:=Smoothing;
  FIsNormalized:=False;
  FIsLog:=False;
end;

procedure THistogram.Normalize;

var
  f:Integer;
  tot:TFloat;

begin
  tot:=Sum(FHistogram);
  for f:=0 to High(FHistogram) do
    FHistogram[f]:=FHistogram[f]/tot;
  FIsNormalized:=True;
end;

procedure THistogram.ConvertToLog;

var f:Integer;

begin
  for f:=0 to High(FHistogram) do
    FHistogram[f]:=Ln(FHistogram[f]);
  FIsLog:=True;
end;

procedure THistogram.ConvertToP;

var f:Integer;

begin
  for f:=0 to High(FHistogram) do
    FHistogram[f]:=Exp(FHistogram[f]);
  FIsLog:=False;
end;

procedure THistogram.WriteToBuffer(ABuffer: TFileBuffer);
begin
  ABuffer.WriteInteger(HistogramVersion);
  ABuffer.WriteInteger(FDataCount);
  ABuffer.WriteFloat(FSmoothing);
  ABuffer.WriteBoolean(FIsLog);
  ABuffer.WriteBoolean(FIsNormalized);
  ABuffer.WriteFloats(FHistogram);
  ABuffer.WriteStrings(FLabels);
  ABuffer.WriteString(Name);
end;

procedure THistogram.ReadFromBuffer(ABuffer: TFileBuffer);

begin
  ABuffer.GetInteger; //skip version;
  ABuffer.ReadInteger(FDataCount);
  ABuffer.ReadFloat(FSmoothing);
  ABuffer.ReadBoolean(FIsLog);
  ABuffer.ReadBoolean(FIsNormalized);
  ABuffer.ReadFloats(FHistogram);
  ABuffer.ReadStrings(FLabels);
  ABuffer.ReadString(Name);
end;

function THistogram.HistogramType: Integer;
begin
  Result:=HistTypeNone;
end;

function THistogram.BinCount: Integer;
begin
  Result:=Length(FHistogram);
end;

{ TFloatHistogram }

procedure TFloatHistogram.SetFMin(Val: TFloat);
begin
  FMin:=Val;
  FRange:=FMax-FMin;
end;

procedure TFloatHistogram.SetFMax(Val: TFloat);
begin
  FMax:=Val;
  FRange:=FMax-FMin;
end;

constructor TFloatHistogram.Create(SourceData: TFloats);
begin
  inherited Create;
  FData:=SourceData;
  FDataCount:=Length(FData);
end;

procedure TFloatHistogram.CalcHistogram(NumBins: Integer;Smoothing: TFloat);

var
  f,ix:Integer;
  total:TFloat;

begin
  inherited;
  NumBins:=Length(FHistogram);
    //if NumBins<0 the length is computed in the parent class
  FDataCount:=Length(FData);
  for f:=0 to High(FData) do
    begin
    ix:=Trunc((FData[f]-FMin)/FRange*NumBins);
    if ix<0 then ix:=0;
    if ix>=NumBins then ix:=High(FHistogram);
    FHistogram[ix]:=FHistogram[ix]+1;
    end;
end;

procedure TFloatHistogram.WriteToBuffer(ABuffer: TFileBuffer);
begin
  inherited WriteToBuffer(ABuffer);
  ABuffer.WriteInteger(FloatHistogramVersion);
  ABuffer.WriteFloat(FMin);
  ABuffer.WriteFloat(FMax);
  ABuffer.WriteFloat(FRange);
end;

procedure TFloatHistogram.ReadFromBuffer(ABuffer: TFileBuffer);
begin
  inherited ReadFromBuffer(ABuffer);
  ABuffer.GetInteger;//skip FloatHistogramVersion;
  ABuffer.ReadFloat(FMin);
  ABuffer.ReadFloat(FMax);
  ABuffer.ReadFloat(FRange);
end;

function TFloatHistogram.HistogramType: Integer;
begin
  Result:=HistTypeFloat;
end;

function TFloatHistogram.BinVal(X: TFLoat): TFloat;

begin
  if X<=FMin then Result:=FHistogram[0]
  else if X>=FMax then Result:=FHistogram[High(FHistogram)]
  else Result:=FHistogram[Trunc((X-FMin)/FRange*Length(FHistogram))]
end;

function TFloatHistogram.BinVals(X: TFLoats): TFloats;

var
  mult:TFloat;
  f:Integer;

begin
  mult:=1/FRange*Length(FHistogram);
  SetLength(Result,Length(X));
  for f:=0 to High(X) do
    if X[f]<=FMin then Result[f]:=FHistogram[0]
    else if X[f]>=FMax then Result[f]:=FHistogram[High(FHistogram)]
    else Result[f]:=FHistogram[Trunc((X[f]-FMin)*mult)];
end;

procedure TFloatHistogram.ReportHistogram(Sl: TStringList; WriteLabels: Boolean);

var
  f:Integer;
  s:string;

begin
  for f:=0 to High(FHistogram) do
    begin
    if WriteLabels and (FLabels<>nil) then
      s:=FLabels[f]+#9
    else
      s:='';
    s:=s+FloatToStr(FHistogram[f]);
    Sl.Add(s);
    end;
end;

end.

