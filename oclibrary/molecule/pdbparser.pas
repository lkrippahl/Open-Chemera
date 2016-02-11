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
  MODEL: assigns number to model
  ENDMDL: Increments model number
  TER: increments chain counter.

Revisions:
To do:
  Read pdb info and comments
  create multiple models when reading a file with several models
*******************************************************************************}

unit pdbparser;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, stringutils, zstream;

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
    Coords:TCoord;
    Occupancy:Single;
    OccTemp:Single;     //merge of occupation and temp, as used in some cases
    Temp:Single;
    Element:string[2];
    Charge:string[2];
    ModelNum:Integer;   //given by MODEL or 1 if no MODEL; incremented by ENDMDL
    ChainNum:Integer;
      //count incremented on TER, starting at 0. Can be used to index the chainIDs
    end;

  TPDBAtoms=array of TPDBAtom;

  TPDBConnection=record
    AtomID:Integer;
    Connects:TIntegers;
  end;

  TPDBConnections=array of TPDBConnection;

  TPDBInfo=record
    Header,Title,Compound,Source,
    KeyWords,ExpTechnique,Author,Journal,Remarks:TStrings;
    UserComments:TStrings;
    end;

  { TPDBReader }

  TPDBReader=class
  private
    FChainIDs:TSimpleStrings;
    FAtoms:TPDBAtoms;
    FConnections:TPDBConnections;
    FAtomCount,FModelCount,FChainCount:Integer;
    FInfo:TPDBInfo;
  public
    property Info:TPDBInfo read FInfo;
    property Atoms:TPDBAtoms read FAtoms;
    property Connections:TPDBConnections read FConnections;
    property AtomCount:Integer read FAtomCount;
    property ChainCount:Integer read FChainCount;
    property ChainIDs:TSimpleStrings read FChainIDs;
    property ModelCount:Integer read FModelCount;
    constructor Create(FromFile:string = '');
    procedure IndexChains;                                    //creates FChainIDs
    procedure Load(FromFile:string);
    procedure ReadFromList(Buf: TStringList);
    procedure Clear;
  end;

  function AtomRecord(AtomName,ResName,ChainName:string;AtomId,
                      ResId:Integer;Position:TCoord;Element:string;
                      Occupancy:TFloat=1.0;
                      Temperature:TFloat=1.0):string;

implementation

function AtomRecord(AtomName,ResName,ChainName:string;AtomId,
                    ResId:Integer;Position:TCoord;Element:string;
                    Occupancy:TFloat=1.0;
                    Temperature:TFloat=1.0):string;

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

begin
  if Length(AtomName)<4 then
    AtomName:=' '+AtomName;
    //this seems to be the rule on pdb files, although that is not what the documentation states.
  Result:='ATOM  '+
          RightJustify(AtomId,5)+' '+
          LeftJustify(AtomName,4)+' '+
          LeftJustify(ResName,3)+' '+
          LeftJustify(ChainName,1)+
          RightJustify(ResId,4)+'    '+
          RightJustify(Position[0],8,3)+
          RightJustify(Position[1],8,3)+
          RightJustify(Position[2],8,3)+
          RightJustify(Occupancy,6,2)+
          RightJustify(Temperature,6,2)+
          '           '+
          LeftJustify(Element,2);
end;

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

function GetPdbAtom(s:string;var OldChain:string):TPDBAtom;

begin
  with Result do
      begin
      IsHet:=False;
      Serial:=GetInteger(s,7,11);
      AtomName:=Deblank(GetString(s,13,16));

      AltLoc:=GetString(s,17,17);
      ResName:=GetString(s,18,20);
      ChainID:=GetString(s,22,22);
      if ChainId<>OldChain then
        begin
        Inc(FChainCount);
        OldChain:=ChainId;
        end;
      ResSeq:=GetInteger(s,23,26);
      ICode:=GetString(s,27,27);
      Coords[0]:=GetFloat(s,31,38);
      Coords[1]:=GetFloat(s,39,46);
      Coords[2]:=GetFloat(s,47,54);
      Occupancy:=GetFloat(s,55,60);
      Temp:=GetFloat(s,61,66);
      OccTemp:=GetFloat(s,55,66);
      Element:=GetString(s,77,78);
       // sometimes charge comes right after the element name
      if (Length(Element)>1) and
         (not (Element[2] in ['a'..'z'])) then
        Element:=Element[1];
      ModelNum:=FModelCount; //-1 indicates ocurred before MODEL keyword
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
  s,oldchain:string;

begin
  FChainCount:=-1;
  oldchain:='***';
  //FChainCount is increased by GetPDBAtom whenever the chain changes
  for f:=0 to Buf.Count-1 do
    begin
    s:=buf.Strings[f];
    if (Pos('ATOM ',s)=1) or (Pos('HETATM',s)=1) then
      begin
      FAtoms[FAtomCount]:=GetPDBAtom(s,oldchain);
      FAtoms[FAtomCount].IsHet:=(Pos('HETATM',s)=1);
      Inc(FAtomCount);
      end
    else if Pos('CONNECT',s)=1 then
      begin
      SetLength(FConnections,Length(FConnections)+1);
      FConnections[High(FConnections)]:=GetPdbConnect(s);
      end
    else if Pos('TER',s)=1 then
        oldchain:='***'
          //This forces an FChainCount increase
          //even if the chain ID remains the same
    else if Pos('MODEL',s)=1 then
        try
          oldchain:='***';
          FModelCount:=StrToInt(Deblank(Copy(s,8,Length(s))));
        except
          Inc(FModelCount);
        end
    else if Pos('ENDMDL',s)=1 then
        Inc(FModelCount);
    end;
end;

begin
  Clear;
  SetLength(FAtoms,Buf.Count);//maximum number of atoms is number of lines in buffer;
  ReadLines;
  SetLength(FAtoms,FAtomCount);
  if FAtomCount>0 then
    begin
    FChainCount:=FAtoms[High(FAtoms)].ChainNum+1; // in case one TER was missing
    FModelCount:=FAtoms[High(FAtoms)].ModelNum; // 0 if no models. This may be redundant, but with PDB files one never knows...
    IndexChains;
    end;
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
