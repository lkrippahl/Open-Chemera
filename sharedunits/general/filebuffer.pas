{*******************************************************************************
This source file is public domain. It is provided as is, with no explicit
or implied warranties regarding anything whatsoever.
********************************************************************************
Author: Ludwig Krippahl
Date: 9.1.2011 (from older version)
Purpose:
  buffer for storing binary data files in memory before writing or after reading
Requirements:
Revisions:
To do:
  Would a memory stream be better?
  Add a maximum memory footprint, read in chunks, to process very large files.
  Implement compression/decompression.
*******************************************************************************}

unit filebuffer;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes;

type

  { TFileBuffer }

  TFileBuffer=class
    private
    FBuffer:array of Byte;
    FPosition:Integer;
    FEndOfBuffer:Boolean;
    FCompressed,FEncrypted:Boolean;
    FFileName:string;
    procedure WritePointer(p:Pointer;Size:Integer);
    procedure ReadPointer(p:Pointer;Size:Integer);
    function GetEndOfBuffer: Boolean;
    function GetSize: Integer;
    public
    property BufferFileName:string read FFileName;
    property Position:Integer read FPosition;
    property Size:Integer read GetSize;
    procedure Compress;
    procedure Decompress;
    procedure WriteBoolean(B:Boolean);
    procedure ReadBoolean(out B:Boolean);
    procedure WriteInteger(I:Integer);
    procedure WriteSingle(S:Single);
    procedure WriteString(s:string);
    procedure WriteDouble(D:double);
    procedure ReadInteger(out I:Integer);
    procedure ReadSingle(out S:Single);
    procedure ReadString(out S:string);
    procedure ReadDouble(out D:Double);


    procedure ReadByte(out B:Byte);
    procedure WriteByte(B:Byte);
    function GetSingle:Single;
    function GetString:string;
    function GetBoolean:Boolean;
    function GetInteger:Integer;
    function GetDouble:Double;
    function GetChar:Char;
    function GetWord:Word;
    function GetLongWord:LongWord;
    function GetByte:Byte;
    function GetMemoryStream:TMemoryStream;
    procedure AddBuffer(Buffer:TFileBuffer);
    property EndOfBuffer:Boolean read GetEndOfBuffer;
    procedure ClearBuffer;
    procedure ResetBuffer;
    procedure Append;
    procedure FlushToFile(FileName:string);
    procedure LoadFromFile(FileName:string;ReadState:Boolean=True);
    procedure LoadCompressedFile(FileName:string);
    constructor Create;
    procedure WriteBlock(P:Pointer;Count:Integer);
    procedure ReadBlock(P:Pointer;Count:Integer);
    function DataLength:Integer;
    procedure WriteIntegers(const Ints:TIntegers);
    procedure ReadIntegers(out Ints:TIntegers);
    function GetIntegers:TIntegers;
    procedure WriteDoubles(const Ds:TDoubles);overload;
    procedure ReadDoubles(out Ds:TDoubles);overload;
    procedure WriteSingles(const Ds:TSingles);
    procedure ReadSingles(out Ds:TSingles);
    procedure WriteStrings(const ss:TSimpleStrings);
    procedure ReadStrings(out ss:TSimpleStrings);

    procedure ReadFloat(out F:TFloat);
    procedure WriteFloat(F:TFloat);
    procedure WriteFloats(const Ds:TFloats);overload;
    procedure ReadFloats(out Ds:TFloats);overload;
      //Floats are written and loaded as doubles and then converted
      //This is to prevent incompatibility with file formats if precision is changed
      //on compilation

    procedure WriteMemoryStream(Ms:TMemoryStream);
    procedure ReadMemoryStream(out Ms:TMemoryStream);
    end;

implementation

const
  BufferNormal=0;
  BufferCompressed=1;
  BufferEncrypted=2;

type
  TDummyByteArray=packed array[0..10000000] of Byte;
  PDummyByteArray=^TDummyByteArray;

{ TFileBuffer }

procedure TFileBuffer.AddBuffer(Buffer: TFileBuffer);

var
  b:Byte;
  f:Integer;

begin
  if Buffer.FBuffer<>nil then
    begin
    FPosition:=High(FBuffer);
    SetLength(FBuffer,Length(FBuffer)+Length(Buffer.FBuffer));
    for f:=0 to High(Buffer.FBuffer) do
      begin
      b:=Buffer.FBuffer[f];
      WritePointer(@b,1);
      end;
    end;
end;

procedure TFileBuffer.Append;
begin
  FPosition:=High(FBuffer);
  FEndOfBuffer:=True;
end;

procedure TFileBuffer.ClearBuffer;
begin
  FPosition:=-1;
  FBuffer:=nil;
  FEndOfBuffer:=True;
end;

procedure TFileBuffer.Compress;
begin
  FCompressed:=True;
end;

constructor TFileBuffer.Create;
begin
  inherited Create;
  ClearBuffer;
end;

procedure TFileBuffer.Decompress;
begin
  FCompressed:=False;
end;


procedure TFileBuffer.FlushToFile(FileName:string);

var
  Fil:file;
  b:Byte;

begin
  FFileName:=FileName;
  if FBuffer<>nil then
    begin
    ResetBuffer;
    try
    AssignFile(Fil,FileName);
    Rewrite(Fil,1);
    BlockWrite(Fil,FBuffer[0],Length(FBuffer));
    finally
    CloseFile(Fil);
    end;
    end;
end;

function TFileBuffer.GetBoolean: Boolean;
begin
  ReadBoolean(Result);
end;

function TFileBuffer.GetDouble: Double;
begin
  ReadDouble(Result);
end;

function TFileBuffer.GetEndOfBuffer: Boolean;
begin
  Result:=FEndOfBuffer;
end;

function TFileBuffer.GetInteger: Integer;
begin
  ReadInteger(Result);
end;

function TFileBuffer.GetSingle: Single;
begin
  ReadSingle(Result);
end;

function TFileBuffer.GetString: string;
begin
  ReadString(Result);
end;

procedure TFileBuffer.LoadFromFile(FileName: string;ReadState:Boolean=True);

var
  Fil:file;
  len:Integer;
  b:byte;

begin
  FFileName:=FileName;
  ClearBuffer;
  AssignFile(Fil,FileName);
  FileMode := 0;
  Reset(Fil,1);
  len:=FileSize(Fil);
  SetLength(FBuffer,Len);
  BlockRead(Fil,FBuffer[0],Len);
  ResetBuffer;
  CloseFile(Fil);
end;

procedure TFileBuffer.ReadBlock(P: Pointer; Count: Integer);
begin
  ReadPointer(P,Count);
end;

procedure TFileBuffer.ReadBoolean(out B: Boolean);

var i:Integer;

begin
  ReadInteger(I);
  B:=I<>0;
end;

procedure TFileBuffer.ReadByte(out B: Byte);
begin
  ReadPointer(@B,SizeOf(Byte));
end;

procedure TFileBuffer.ReadDouble(out D: Double);
begin
  ReadPointer(@D,SizeOf(Double));
end;

procedure TFileBuffer.ReadInteger(out I: Integer);
begin
  ReadPointer(@I,SizeOf(Integer));
end;

procedure TFileBuffer.ReadPointer(p: Pointer; Size: Integer);

var f:Integer;

begin
  if (Size>0) and (p<>nil) and (Length(FBuffer)>=FPosition+Size) then
    for f:=0 to Size-1 do
      begin
      PDummyByteArray(p)^[f]:=FBuffer[FPosition];
      Inc(FPosition);
      end
    else FEndOfBuffer:=True;
  if FPosition>High(FBuffer) then
    FEndOfBuffer:=True;
end;

procedure TFileBuffer.ReadSingle(out S: Single);
begin
  ReadPointer(@S,SizeOf(Single));
end;

procedure TFileBuffer.ReadString(out S: string);

var b:Byte;

begin
  s:='';
  repeat
    ReadPointer(@b,SizeOf(b));
    if b>0 then s:=s+Chr(b);
  until FEndOfBuffer or (b=0);
end;

procedure TFileBuffer.ResetBuffer;
begin
  FPosition:=0;
  FEndOfBuffer:=FPosition>=High(FBuffer);
end;

procedure TFileBuffer.WriteBlock(P: Pointer; Count: Integer);
begin
  WritePointer(p,Count);
end;

procedure TFileBuffer.WriteBoolean(B: Boolean);
begin
  if B then WriteInteger(1)
  else WriteInteger(0);
end;

procedure TFileBuffer.WriteByte(B: Byte);
begin
  WritePointer(@B,SizeOf(Byte));
end;

procedure TFileBuffer.WriteDouble(D: double);
begin
  WritePointer(@D,SizeOf(Double));
end;

procedure TFileBuffer.WriteInteger(I: Integer);
begin
  WritePointer(@I,SizeOf(Integer));
end;

procedure TFileBuffer.WritePointer(p: Pointer; Size: Integer);

var f:Integer;

begin
  if (Size>0) and (p<>nil) then
    begin
    if Length(FBuffer)<FPosition+Size+1 then
      SetLength(FBuffer,FPosition+Size+1);
    for f:=0 to Size-1 do
      begin
      Inc(FPosition);
      FBuffer[FPosition]:=PDummyByteArray(p)^[f];
      end;
    end;
end;

procedure TFileBuffer.WriteSingle(S: Single);
begin
  WritePointer(@S,SizeOf(S));
end;

procedure TFileBuffer.WriteString(s: string);

var f:Integer;
    b:Byte;

begin
  if s<>'' then
    for f:=1 to Length(s) do
      begin
      b:=Ord(s[f]);
      WritePointer(@b,SizeOf(Byte));
      end;
    b:=0;
    WritePointer(@b,SizeOf(Byte));
end;

function TFileBuffer.DataLength: Integer;
begin
  Result:=Length(FBuffer);
end;

function TFileBuffer.GetChar: Char;

var b:Byte;

begin
  ReadByte(b);
  Result:=Char(b);
end;

function TFileBuffer.GetWord: Word;
begin
  ReadBlock(@Result,2);
end;

function TFileBuffer.GetLongWord: LongWord;
begin
  ReadBlock(@Result,4);
end;

function TFileBuffer.GetByte: Byte;
begin
  ReadBlock(@Result,1);
end;

procedure TFileBuffer.ReadIntegers(out Ints:TIntegers);

var f:Integer;

begin
  SetLength(Ints,GetInteger);
  for f:=0 to High(Ints) do readInteger(Ints[f]);
end;

procedure TFileBuffer.WriteIntegers(const Ints:TIntegers);

var f:Integer;

begin
  WriteInteger(Length(Ints));
  for f:=0 to High(Ints) do WriteInteger(Ints[f]);
end;

procedure TFileBuffer.ReadDoubles(out Ds: TDoubles);
var f:Integer;

begin
  SetLength(Ds,GetInteger);
  for f:=0 to High(Ds) do ReadDouble(Ds[f]);
end;

procedure TFileBuffer.ReadSingles(out Ds: TSingles);
var f:Integer;

begin
  SetLength(Ds,GetInteger);
  for f:=0 to High(Ds) do ReadSingle(Ds[f]);
end;

procedure TFileBuffer.WriteDoubles(const Ds: TDoubles);
var f:Integer;

begin
  WriteInteger(Length(Ds));
  for f:=0 to High(Ds) do WriteDouble(Ds[f]);
end;

procedure TFileBuffer.WriteSingles(const Ds: TSingles);
var f:Integer;

begin
  WriteInteger(Length(Ds));
  for f:=0 to High(Ds) do WriteSingle(Ds[f]);
end;

procedure TFileBuffer.ReadStrings(out ss: TSimpleStrings);
var f,l:Integer;

begin
  l:=GetInteger;
  SetLength(ss,l);
  for f:=0 to High(ss) do ReadString(ss[f]);
end;

procedure TFileBuffer.ReadFloat(out F: TFloat);
begin
  F:=GetDouble;
end;

procedure TFileBuffer.WriteFloat(F: TFloat);
begin
  WriteDouble(F);
end;

procedure TFileBuffer.WriteFloats(const Ds: TFloats);

var f:Integer;

begin
  WriteInteger(Length(Ds));
  for f:=0 to High(Ds) do WriteDouble(Ds[f]);
end;

procedure TFileBuffer.ReadFloats(out Ds: TFloats);

var f:Integer;

begin
  SetLength(Ds,GetInteger);
  for f:=0 to High(Ds) do
    Ds[f]:=GetDouble;
end;

procedure TFileBuffer.WriteStrings(const ss: TSimpleStrings);
var f:Integer;

begin
  WriteInteger(Length(ss));
  for f:=0 to High(ss) do WriteString(ss[f]);
end;

procedure TFileBuffer.LoadCompressedFile(FileName: string);

{var InputStream:TFileStream;
    DecompStr:TDecompressionStream;
    Temp:array[0..999] of Byte;
    Ol,f,Count:Integer;  }

begin
  {InputStream:=TFileStream.Create(FileName,fmOpenRead);
  DecompStr:=TDecompressionStream.Create(InputStream);
  FBuffer:=nil;
  repeat
    Ol:=Length(FBuffer);
    Count:=DecompStr.Read(Temp,1000);
    SetLength(FBuffer,Length(FBuffer)+Count);
    for f:=0 to Count-1 do
      FBuffer[Ol+f]:=Temp[f];
  until Count<1000;
  ResetBuffer;
  DecompStr.Free;
  InputStream.Free;  }
end;

function TFileBuffer.GetIntegers: TIntegers;
  var f:Integer;
begin
  SetLength(Result,GetInteger);
  for f:=0 to High(Result) do ReadInteger(Result[f]);
end;

function TFileBuffer.GetMemoryStream: TMemoryStream;

var MSize:Int64;

begin
  Result:=TMemoryStream.Create;
  ReadBlock(@MSize,SizeOf(Size));
  Result.SetSize(MSize);
  Result.Seek(0,soFromBeginning);
  ReadBlock(Result.Memory,MSize);
end;

procedure TFileBuffer.ReadMemoryStream(out Ms: TMemoryStream);
begin
  Ms:=GetMemoryStream;
end;

procedure TFileBuffer.WriteMemoryStream(Ms: TMemoryStream);

var
  MSize:Int64;
  ByteBuffer:^Byte;

begin
  MSize:=Ms.Size;
  WritePointer(@MSize,SizeOf(Int64));
  Ms.Seek(0,soFromBeginning);
  WritePointer(Ms.Memory,MSize);
end;

function TFileBuffer.GetSize: Integer;
begin
  Result:=Length(FBuffer)
end;

end.
