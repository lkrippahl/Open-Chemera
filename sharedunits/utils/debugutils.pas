{*******************************************************************************
This source file is public domain. It is provided as is, with no explicit
or implied warranties regarding anything whatsoever.
********************************************************************************
Author: Ludwig Krippahl
Date: 9.1.2011
Purpose:
  Debug utilities
  For debugging large variables, mostly.
  use SetDebugFile to set a file, otherwise will debug to the console
Requirements:
Revisions:
To do:
  Conditional compiling, debug flag
*******************************************************************************}
unit debugutils;

{$mode objfpc}{$H+}


interface

uses
  Classes, SysUtils, basetypes;

procedure LogHeader;
procedure SetDebugFile(FileName:string);
procedure DebugReport(s:string);overload;
procedure DebugReport(ss:TSimpleStrings);overload;
procedure DebugReport(ints:TIntegers);overload;
procedure DebugReport(Coords:TCoords);overload;
procedure DebugReport(Coord:TCoord);overload;
procedure DebugReport(Mat:TMatrix);overload;
procedure DebugReport(I:Integer);overload;
procedure DebugReport(F:TFloat);overload;

implementation

var DebugFileName:string='';//g:\debug.txt';    //if no file, debugs to console

procedure OpenFile(out Fil:TextFile);

begin
  AssignFile(Fil,DebugFileName);
  if FileExists(DebugFileName) then
    Append(Fil)
  else Rewrite(Fil);
end;

procedure LogHeader;

var fil:TextFile;

begin
  if DebugFileName<>'' then
    begin
    OpenFile(fil);
    WriteLn(fil,'Log entry: '+DateTimeToStr(Now));
    CloseFile(fil);
    end
  else WriteLn('Log entry: '+DateTimeToStr(Now));
end;

procedure LogString(S:string);

var fil:TextFile;

begin
  if DebugFileName<>'' then
    begin
    OpenFile(fil);
    WriteLn(fil,s);
    CloseFile(fil);
    end
  else WriteLn(s);
end;

procedure DebugReport(s:string);overload;

begin
  LogString(s);
end;

procedure DebugReport(ss:TSimpleStrings);overload;

var
  f:Integer;

begin
  for f:=0 to High(ss) do
    LogString(ss[f]);

end;

procedure DebugReport(Ints:TIntegers);overload;

var
  f:Integer;
  s:string;

begin
  s:='';
  for f:=0 to High(Ints) do
    s:=s+IntToStr(Ints[f])+#9;
  LogString(s);
end;

procedure DebugReport(Coords:TCoords);overload;

var
  f:Integer;
  s:string;

begin
  for f:=0 to High(Coords) do
    LogString(FloatToStr(Coords[f,0])+#9+FloatToStr(Coords[f,1])+#9+FloatToStr(Coords[f,2]));

end;

procedure DebugReport(Coord:TCoord);overload;
begin
  LogString(FloatToStr(Coord[0])+#9+FloatToStr(Coord[1])+#9+FloatToStr(Coord[2]));
end;

procedure DebugReport(Mat:TMatrix);overload;

var
  f,g:Integer;
  s:string;

begin
  for f:=0 to High(Mat) do
    begin
    s:='';
    for g:=0 to High(Mat[f]) do
      s:=s+FloatToStr(Mat[f,g])+#9;
    LogString(s);
    end;
end;

procedure DebugReport(I:Integer);overload;

begin
  LogString(IntToStr(I));
end;

procedure DebugReport(F:TFloat);overload;

begin
  LogString(FloatToStr(F));
end;


procedure SetDebugFile(FileName:string);

begin
  DebugFileName:=FileName;
end;

end.

