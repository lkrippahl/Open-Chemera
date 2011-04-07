{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 9.1.2011
Purpose:
  Basic sequence data structures and utils
  Based on UniProtKB and EBI results
Requirements:
Revisions:
To do:
*******************************************************************************}

unit sequence;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils;


type
  //these fields are based on the EBI results
  TOCSequence=record
    Sequence:string;
    Organism:string;
    Database:string;    // e.g. uniprot
    ID:string;          //protein identifier in the database. E.g. Q6AW43_9ANNE
    Code:string;         //protein accession code, e.g. Q6AW43

    Evidence:string;    //UniProtKB, evidence level, can be a number 1-5
                        //see http://www.uniprot.org/manual/protein_existence
                        //TODO: should it always be a number?
    Gene:string;
    Description:string;
    end;
  TOCSequences=array of TOCSequence;

procedure AddSequence(se:TOCSequence; var ss:TOCSequences);
function EmptySequence:TOCSequence;

implementation

procedure AddSequence(se: TOCSequence; var ss: TOCSequences);

var l:Integer;

begin
  l:=Length(ss);
  SetLength(ss,l+1);
  ss[l]:=se;
end;

function EmptySequence:TOCSequence;

begin
  with Result do
    begin
    Sequence:='';
    Organism:='';
    Database:='';
    ID:='';
    Evidence:='';
    Code:='';
    Gene:='';
    Description:='';
    end;
end;

end.

