unit clustalparser;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, alignment, ocstringutils;


function ReadClustal(fromfile:string):TMSA;

implementation

function ReadClustal(fromfile: string): TMSA;

var
  sl:TStringList;
  f,Ix:Integer;
  nam,seq,s:string;

begin
  sl:=TStringList.Create;
  sl.LoadFromFile(fromfile);
  Result.GapMarker:='-';
  for f:=1 to sl.Count-1 do // first line is header
    begin
    s:=sl.Strings[f];
    nam:=Deblank(Copy(s,1,18));
    if nam<>'' then
      begin
      seq:=Copy(s,19,Length(s));
      Delete(seq,Pos(' ',seq),Length(seq));
      end;
    Ix:=LastIndexOf(nam,Result.SequenceIDs);
    if Ix<0 then
      begin
      AddString(nam,Result.SequenceIDs);
      AddString(seq,Result.Alignment);
      end
    else
      Result.Alignment[Ix]:=Result.Alignment[Ix]+seq;
    end;
  sl.Free;
end;

end.

