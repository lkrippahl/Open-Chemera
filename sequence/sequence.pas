unit sequence;

{

Base unit for sequence data

}

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils;

type
  TOCSequence=record
    Sequence:string;
    Organism:string;
    Database:string;
    ID:string;
    Code:string;
    Description:string;
    end;
  TOCSequences=array of TOCSequence;

procedure AddSequence(se:TOCSequence; var ss:TOCSequences);

implementation

procedure AddSequence(se: TOCSequence; var ss: TOCSequences);

var l:Integer;

begin
  l:=Length(ss);
  SetLength(ss,l+1);
  ss[l]:=se;
end;

end.

