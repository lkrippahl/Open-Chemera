{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 9.1.2011
Purpose:
  Parser for PDB files
Requirements:
  Record supported:
  ATOM / HETATM
  CONNECT
  MODEL: increments model counter only
  TER: increments chain counter.

Revisions:
To do: Read pdb info and comments
*******************************************************************************}

unit pdbparser;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, ocstringutils, zstream;

type
  TPDBAtom = record
    IsHet:Boolean;      //heteroatom
    Serial:Integer;
    AtomName:string[4];
    AltLoc:string[1];
    ResName:string[3];
    ChainID:string[1];
    ResSeq:Integer;
    ICode:string[1];
    Coords:TOCCoord;
    Occupancy:Single;
    Temp:Single;
    Element:string[2];
    Charge:string[2];
    ModelNum:Integer;   //continuous count incremented on MODEL keyword, -1 if no MODEL
    ChainNum:Integer;
      //count incremented on TER, starting at 0. Can be used to index the chainIDs
    end;

  TPDBAtoms=array of TPDBAtom;

  TPDBConnection=record
    AtomID:Integer;
    Connects:TOCIntegers;
  end;

  TPDBConnections=array of TPDBConnection;

  TPDBLayerInfo=record
    Header,Title,Compound,Source,
    KeyWords,ExpTechnique,Author,Journal,Remarks:TOCStrings;
    UserComments:TOCStrings;
    end;

  { TPDBReader }

  TPDBReader=class
  private
    FChainIDs:TOCStrings;
    FAtoms:TPDBAtoms;
    FConnections:TPDBConnections;
    FAtomCount,FModelCount,FChainCount:Integer;
    FInfo:TPDBLayerInfo;
  public
    property Info:TPDBLayerInfo read FInfo;
    property Atoms:TPDBAtoms read FAtoms;
    property Connections:TPDBConnections read FConnections;
    property AtomCount:Integer read FAtomCount;
    property ChainCount:Integer read FChainCount;
    property ChainIDs:TOCStrings read FChainIDs;
    property ModelCount:Integer read FModelCount;
    constructor Create(FromFile:string = '');
    procedure IndexChains;                                    //creates FChainIDs
    procedure Load(FromFile:string);
    procedure ReadFromList(Buf: TStringList);
    procedure Clear;
  end;


implementation

{ TPDBReader }

constructor TPDBReader.Create(fromfile: string = '');
begin
  inherited Create;
  if fromfile<>'' then Load(fromfile);
end;

procedure TPDBReader.IndexChains;

var f:Integer;

begin
  SetLength(FChainIDs,FChainCount);
  for f:=0 to High(FAtoms) do
    FChainIDs[FAtoms[f].ChainNum]:=FAtoms[f].ChainID;
end;

procedure TPDBReader.Load(FromFile: string);
// if file extension is .gz assumes a gzip file

var
  buf:TStringList;
  zfs:TGZFileStream;

begin
  buf:=TStringList.Create;
  if Uppercase(ExtractFileExt(FromFile))='.GZ' then
    begin
    zfs:=TGZFileStream.Create(FromFile,gzopenread);
    buf.LoadFromStream(zfs);
    zfs.Free;
    end
    else buf.LoadFromFile(FromFile);
  ReadFromList(buf);
  buf.Free;
end;

procedure TPDBReader.ReadFromList(Buf: TStringList);

{
1 - 6          Record name     "ATOM"
7 - 11         Integer         serial        Atom serial number
13 - 16        Atom            name          Atom name
17             Character       altLoc        Alternate location indicator
18 - 20        Residue name    resName       Residue name
22             Character       chainID       Chain identifier
23 - 26        Integer         resSeq        Residue sequence number
27             Character       iCode         Code for insertion of residues
31 - 38        Real(8.3)       x             Orthogonal coordinates for X
39 - 46        Real(8.3)       y             Orthogonal coordinates for Y
47 - 54        Real(8.3)       z             Orthogonal coordinates for Z
55 - 60        Real(6.2)       occupancy     Occupancy
61 - 66        Real(6.2)       tempFactor    Temperature factor
68 - 70        Integer         ftNote        Footnote number, being deprecated
73 - 76        LString(4)      segID         Segment identifier, left-justified
77 - 78        LString(2)      element       Element symbol, right-justified
79 - 80        LString(2)      charge        Charge on the atom, IUPAC form}

function GetPdbAtom(s:string):TPDBAtom;

begin
  with Result do
      begin
      IsHet:=False;
      Serial:=GetInteger(s,7,11);
      AtomName:=GetString(s,13,16);
      AltLoc:=GetString(s,17,17);
      ResName:=GetString(s,18,20);
      ChainID:=GetString(s,22,22);
      ResSeq:=GetInteger(s,23,26);
      ICode:=GetString(s,27,27);
      Coords[0]:=GetFloat(s,31,38);
      Coords[1]:=GetFloat(s,39,46);
      Coords[2]:=GetFloat(s,47,54);
      Occupancy:=GetFloat(s,55,60);
      Temp:=GetFloat(s,61,66);
      Element:=GetString(s,77,78);
      ModelNum:=FModelCount-1; //-1 indicates ocurred before MODEL keyword
      ChainNum:=FChainCount;
      end
end;


function GetPdbConnect(s:string):TPDBConnection;

var i:Integer;

{
1 - 6 Record name "CONECT"
7 - 11 Integer serial Atom serial number
12 - 16 Integer serial Serial number of bonded atom
17 - 21 Integer serial Serial number of bonded atom
22 - 26 Integer serial Serial number of bonded atom
27 - 31 Integer serial Serial number of bonded atom
}

begin
  Result.AtomId:=GetInteger(s,7,11);
  Result.Connects:=nil;
  if GetInteger(s,12,16,i) then AddToArray(i,Result.Connects);
  if GetInteger(s,17,21,i) then AddToArray(i,Result.Connects);
  if GetInteger(s,22,26,i) then AddToArray(i,Result.Connects);
  if GetInteger(s,27,31,i) then AddToArray(i,Result.Connects);
end;

procedure ReadLines;

var
  f:Integer;
  s:string;

begin
  for f:=0 to Buf.Count-1 do
    begin
    s:=buf.Strings[f];
    if (Pos('ATOM ',s)=1) or (Pos('HETATM',s)=1) then
      begin
      FAtoms[FAtomCount]:=GetPDBAtom(s);
      FAtoms[FAtomCount].IsHet:=(Pos('HETATM',s)=1);
      Inc(FAtomCount);
      end
    else if Pos('CONNECT',s)=1 then
      begin
      SetLength(FConnections,Length(FConnections)+1);
      FConnections[High(FConnections)]:=GetPdbConnect(s);
      end
    else if Pos('TER',s)=1 then Inc(FChainCount)
    else if Pos('MODEL',s)=1 then Inc(FModelCount);
    end;
end;

begin
  Clear;
  SetLength(FAtoms,Buf.Count);//maximum number of atoms is number of lines in buffer;
  ReadLines;
  SetLength(FAtoms,FAtomCount);
  FChainCount:=FAtoms[High(FAtoms)].ChainNum+1; // in case one TER was missing
  FModelCount:=FAtoms[High(FAtoms)].ModelNum+1; // 0 if no models. This may be redundant, but with PDB files one never knows...
  IndexChains;
end;

procedure TPDBReader.Clear;
begin
  FAtoms:=nil;
  FConnections:=nil;
  FAtomCount:=0;
  FModelCount:=0;
  FChainCount:=0;
end;

end.

