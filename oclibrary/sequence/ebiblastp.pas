{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 21.09.2013
Purpose:
  Parser for EBI BLAST results files (in XML)
Requirements:
Revisions:
To do:
  Still very basic
  Assumes only one query sequence
*******************************************************************************}
unit ebiblastp;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, sequence, XMLRead, DOM;

  function ReadEBIBlast(FileName:string):TOCSequences;

implementation

function ReadEBIBlast(FileName:string):TOCSequences;

var
  doc:TXMLDocument;
  node:TDOMNode;
  hits,seqs:TDOMNOdeList;
  tmp,seq,id,organism:string;
  tot,f:Integer;

procedure AddToResult(IX:Integer;Sequence,Identifier,Organism:string);

begin
  Result[IX]:=EmptySequence;
  Result[IX].ID:=Identifier;
  Result[IX].Sequence:=Sequence;
  Result[IX].Organism:=Organism;
end;

begin
  Result:=nil;
  doc:=nil;
  hits:=nil;
  try
    ReadXMLFile(doc, FileName);

    hits:=doc.DocumentElement.GetElementsByTagName('hit');

    SetLength(Result,hits.Count+1);

    seqs:=doc.DocumentElement.GetElementsByTagName('sequence');
    id:=seqs[0].Attributes.GetNamedItem('name').NodeValue;
    seqs.free;

    node:=hits[0].FindNode('alignments');
    node:=node.FindNode('alignment');
    node:=node.FindNode('querySeq');
    seq:=node.TextContent;
    AddToResult(0,seq,id,'query');


    for f:=0 to hits.count-1 do
      begin
      id:=hits[f].Attributes.GetNamedItem('id').NodeValue;
      tmp:=hits[f].Attributes.GetNamedItem('description').NodeValue;

      //find organism description. May be OS= or Tax=
      if pos('OS=',tmp)>0 then
        Delete(tmp,1,Pos('OS=',tmp)+2)
      else
        Delete(tmp,1,Pos('Tax=',tmp)+3);

      //Delete next tag, if any
      if Pos('=',tmp)>0 then
        begin
        Delete(tmp,Pos('=',tmp),Length(tmp));
        while tmp[Length(tmp)]<>' ' do
          Delete(tmp,Length(tmp),1);
        Delete(tmp,Length(tmp),1);
        end;

      organism:=tmp;

      node:=hits[f].FindNode('alignments');
      node:=node.FindNode('alignment');
      node:=node.FindNode('matchSeq');
      seq:=node.TextContent;
      AddToResult(f+1,seq,id,organism);
      end;
    hits.Free;
    doc.Free;

  except
  end;

end;

end.

