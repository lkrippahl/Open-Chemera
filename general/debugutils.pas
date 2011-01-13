unit debugutils;

{$mode objfpc}{$H+}

{
Debug utilities

For debugging large variables, mostly.
use SetDebugFile first

TODO: Conditional compiling?

}

interface

uses
  Classes, SysUtils, basetypes;

procedure SetDebugFile(FileName:string);
procedure DebugReport(s:string);overload;
procedure DebugReport(ss:TOCStrings);overload;
procedure DebugReport(ints:TOCIntegers);overload;
procedure DebugReport(Coords:TOCCoords);overload;
procedure DebugReport(Mat:TOCMatrix);overload;

implementation

var DebugFileName:string='G:\log.txt';

procedure OpenFile(out Fil:TextFile);

begin
  AssignFile(Fil,DebugFileName);
  if FileExists(DebugFileName) then
    Append(Fil)
  else Rewrite(Fil);
end;

procedure LogHeader(var Fil:TextFile);

begin
  WriteLn(Fil,'Log entry: '+DateTimeToStr(Now));
end;

procedure DebugReport(s:string);overload;

var fil:TextFile;

begin
  OpenFile(fil);
  LogHeader(fil);
  WriteLn(fil,s);
  CloseFile(fil);
end;

procedure DebugReport(ss:TOCStrings);overload;

var
  fil:TextFile;
  f:Integer;

begin
  OpenFile(fil);
  LogHeader(fil);
  for f:=0 to High(ss) do
    WriteLn(fil,ss[f]);
  CloseFile(fil);
end;

procedure DebugReport(Ints:TOCIntegers);overload;

var
  fil:TextFile;
  f:Integer;
  s:string;

begin
  OpenFile(fil);
  LogHeader(fil);
  s:='';
  for f:=0 to High(Ints) do
    s:=s+IntToStr(Ints[f])+#9;
  WriteLn(fil,s);
  CloseFile(fil);
end;

procedure DebugReport(Coords:TOCCoords);overload;

var
  fil:TextFile;
  f:Integer;
  s:string;

begin
  OpenFile(fil);
  LogHeader(fil);
  for f:=0 to High(Coords) do
    WriteLn(fil,FloatToStr(Coords[f,1])+#9+FloatToStr(Coords[f,2])+#9+FloatToStr(Coords[f,3]));
  CloseFile(fil);
end;

procedure DebugReport(Mat:TOCMatrix);overload;

var
  fil:TextFile;
  f,g:Integer;
  s:string;

begin
  OpenFile(fil);
  LogHeader(fil);
  for f:=0 to High(Mat) do
    begin
    s:='';
    for g:=0 to High(Mat[f]) do
      s:=s+FloatToStr(Mat[f,g])+#9;
    WriteLn(fil,s);
    end;
  CloseFile(fil);
end;


procedure SetDebugFile(FileName:string);

begin
  DebugFileName:=FileName;
end;

end.

