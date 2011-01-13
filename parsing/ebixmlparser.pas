unit ebixmlparser;

{
  Parser for XML data from the European Bioinformatics Institute
    http://www.ebi.ac.uk/
}

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, sequence, alignment, XMLRead, DOM, xmltools;

type
  TBlastPData=record
    SequenceInfo:TOCSequences;
    //All sequences are matches found. The query sequence information is not provided
    Alignments:TSingleAlignments;
    //To do: parameters for the calculation.
    end;

function ReadBlastPFile(filename:string):TBlastPData;
  //Parser for single alignents,  WU-blastp  version 2.0
implementation

function ReadBlastPFile(filename:string): TBlastPData;

// To do: only storing alignment data. Should grab all the rest too.

var
  alnode,node:TDOMNode;
  doc:TXMLDocument;
  qid,currmid:string;

function GetSequenceInfo(node:TDOMNode):TOCSequence;

var
  desc,os:string;

begin
  desc:=GetNodeAttribute(node,'description');
  if Pos('OS=',desc)>0 then
    begin
    os:=Copy(desc,Pos('OS=',desc)+3,Length(desc));
    Delete(desc,Pos('OS=',desc),Length(desc));
    if Pos('GN=',os)>0 then
      Delete(os,Pos('GN=',os)-1,length(os));
    end
  else os:='';
  with Result do
    begin
    Description:=desc;
    Sequence:='';// alignment only returns the matching part of the sequence
    Organism:=os;
    Database:=GetNodeAttribute(node,'database');
    ID:=GetNodeAttribute(node,'id');
    currmid:=ID;
    Code:=GetNodeAttribute(node,'ac');
    end;
end;

function GetAlignment(node:TDOMNode):TSingleAlignment;

var nd:TDOMNode;

begin
  with Result do
    begin
    QueryId:=qid;
    MatchID:=currmid;
    Score:=StrToFloat(GetChildText(node,'score'));
    Bits:=StrToFloat(GetChildText(node,'bits'));
    Expectation:=StrToFloat(GetChildText(node,'expectation'));
    Probability:=StrToFloat(GetChildText(node,'probability'));
    Identity:=StrToInt(GetChildText(node,'identity'));
    Positives:=StrToInt(GetChildText(node,'positives'));
    nd:=node.FindNode('querySeq');
    QueryStart:=StrToInt(GetNodeAttribute(nd,'start'));
    QueryEnd:=StrToInt(GetNodeAttribute(nd,'end'));
    nd:=node.FindNode('matchSeq');
    MatchStart:=StrToInt(GetNodeAttribute(nd,'start'));
    MatchEnd:=StrToInt(GetNodeAttribute(nd,'end'));
    AlignedQuery:=GetChildText(node,'querySeq');
    AlignedMatch:=GetChildText(node,'matchSeq');
    end;
end;

begin
  Result.Alignments:=nil;
  Result.SequenceInfo:=nil;
  ReadXMLFile(doc, filename);
  node:=doc.DocumentElement.FindNode('Header');
  node:=node.FindNode('parameters');
  node:=node.FindNode('sequences');
  node:=node.FindNode('sequence');
  qid:=GetNodeAttribute(node,'name');
  node:=doc.DocumentElement.FindNode('SequenceSimilaritySearchResult');
  node:=node.FindNode('hits');
  node:=node.FindNode('hit');
  while node<>nil do
    begin
    AddSequence(GetSequenceInfo(node),Result.SequenceInfo);
    alnode:=node.FindNode('alignments');
    alnode:=alnode.FindNode('alignment');
    AddAlignment(GetAlignment(alnode),Result.Alignments);
    node:=node.NextSibling
    end;
  doc.Free;
end;

end.

