{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 9.1.2011
Purpose:
  Class for storing PDB-like file data.
  TPDBModel corresponds to a model. Contains chains, residues, etc.
  TPDBModelSet contains several models.
  This is the base organization for molecules. The assumption is that a pdb-like
  organization is sufficient for all (proteins, DNA, ligands, etc).
Requirements:
  Requires PDBeChem ligand dictionary files. See
  http://www.ebi.ac.uk/msd-srv/chempdb
  for more information.

Revisions:
To do:
  Comments

  Adding H atoms

  add data from CCD
  http://www.wwpdb.org/ccd.html
  (need mmCIF first)

  Implement gromos force field (a separate unit for the force fields?)
  http://www.gromacs.org/Downloads/User_contributions/Force_fields

  TPDBInfo should belong to a TPDBFile (or TPDBLayer) intermediate class above
  the TPDBModel, as one file can contain multiple models


*******************************************************************************}

unit pdbmolecules;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, molecules, pdbparser, FileUtil, oclconfiguration,
  stringutils,LCLProc;

type

  TTemplate = record                    //residue templates
    Name: string;
    Atoms: TSimpleStrings;
    AtomicNumbers:TIntegers;
    AtomIDs: TIntegers;
    Coords: TCoords;
    Connects: array of TIntegers;
  end;

  TTemplates = array of TTemplate;

  { TPDBModel }

  TPDBModel = class
  protected
    FTemplates: TTemplates;
    FProtein: TMolecule;
    FFileName: string;
    FInfo: TPDBInfo;
  public
    property Info: TPDBInfo read FInfo write FInfo;
    property Molecule:TMolecule read FProtein;
    constructor Create(Templates: TTemplates;AID:Integer);
    procedure Free;
    procedure ClearChains;
    procedure CreateChains(const IDs: TSimpleStrings);
    function NewEmptyChain(ChainName: string; ChainID: integer): TMolecule;
    function NewChain(ChainName: string; ChainID, Size: integer): TMolecule;
    function GetChain(ChainIx: integer): TMolecule;overload;
    function GetChain(ChainID:string):TMolecule;overload;
    function GetResidue(ChainIx, ResIx: integer): TMolecule;
    function GetAtom(ChainIx, ResIx, AtomIx: integer): TAtom;
    function LoadPDB(FileName: string):TMolecule;
    procedure ResetTemplates(ATemplates:TTemplates);
    function ChainCount:Integer;
    function ResidueCount(ChainIx:Integer):Integer;
    procedure AssignAtomicData; //atomic number, vdW radius
    function TemplateIx(Name:string):Integer;
  end;

  TPDBModels = array of TPDBModel;

  { TPDBModelSet }

  { TPDBModelMan }

  TPDBModelMan = class
  protected
    FLayers: TPDBModels;
    FTemplates: TTemplates;
    procedure LoadTemplates(Path: string);
  public
    constructor Create(molCIFPath: string);
    function Count: integer;
    function LayerByIx(Ix: integer): TPDBModel;
    function AddNewLayer: TPDBModel;
    function LoadLayer(PdbFileName: string):TMolecule;
    procedure ClearLayers;
  end;

  function AtomIsAABackbone(Atom:TAtom):Boolean;

implementation

function AtomIsAABackbone(Atom:TAtom):Boolean;

begin
  with Atom do
    Result:= (Name='N') or (Name='O') or (Name='CA') or (Name='C');
end;


{ TPDBModel }

constructor TPDBModel.Create(Templates:TTemplates;AID:Integer);
begin
  inherited Create;
  FTemplates:=Templates;
  FProtein:=TMolecule.Create('',AID,nil);
end;

procedure TPDBModel.Free;

begin
  if Self <> nil then
  begin
    ClearChains;
    inherited;
  end;
end;

procedure TPDBModel.ClearChains;

begin
  //TODO Add on delete callbacks
  FProtein.Cleanup;
end;

procedure TPDBModel.CreateChains(const IDs: TSimpleStrings);

var
  f:integer;

begin
  for f := 0 to High(IDs) do
    FProtein.NewGroup(IDs[f], f + 1);
end;

function TPDBModel.NewEmptyChain(ChainName: string; ChainID: integer): TMolecule;
begin
  Result:=FProtein.NewGroup(ChainName, ChainId);
end;

function TPDBModel.NewChain(ChainName: string; ChainID, Size: integer): TMolecule;

begin
  Result := NewEmptyChain(ChainName, ChainId);
  Result.CreateEmptyGroups(Size);
end;


function TPDBModel.GetChain(ChainIx: integer): TMolecule;
begin
  Result := FProtein.GetGroup(ChainIx);
end;

function TPDBModel.GetChain(ChainID: string): TMolecule;

var f:Integer;

begin
  Result:=nil;
  for f:=0 to High(FProtein.Groups) do
    if FProtein.Groups[f].Name=ChainId then
      begin
      Result:=FProtein.Groups[f];
      Break;
      end;
end;

function TPDBModel.GetResidue(ChainIx, ResIx: integer): TMolecule;
begin
  Result := FProtein.GetGroup(ChainIx);
  if Result<>nil then Result:=Result.GetGroup(ResIx);
end;

function TPDBModel.GetAtom(ChainIx, ResIx, AtomIx: integer): TAtom;

var mol:TMolecule;

begin
  Result:=nil;
  mol := FProtein.GetGroup(ChainIx);
  if mol<>nil then mol:=mol.GetGroup(ResIx);
  if mol <>nil then Result := mol.GetAtom(AtomIx);
end;

function TPDBModel.LoadPDB(FileName: string):TMolecule;

var
  parser: TPDBReader;
  cr,f,cc: integer;
  cres:TMolecule;
  atom:TAtom;

begin
  ClearChains;
  parser := TPDBReader.Create(FileName);
  FProtein.Name:=ChangeFileExt(ExtractFileName(FileName),'');
  CreateChains(parser.ChainIDs);

  //read atoms

  cc:=-1;                                 //current chain
  for f := 0 to High(parser.Atoms) do
    with Parser.Atoms[f] do
      begin
      if cc<>ChainNum then                //ChainNumber is always >=0
        begin
        cc:=ChainNum;
        cr:=-1;                           //current residue number, also >=0
        end;
      if ResSeq<>cr then                  //get new residue
        begin
        cres:=FProtein.GetGroup(cc).NewGroup(ResName,ResSeq);
        cr:=ResSeq;
        end;
      atom:=cres.NewAtom(AtomName,Serial);
      atom.Coords:=Coords;
      end;

  //make extra connections
  //TO DO:
  for f:=0 to High(parser.Connections) do
    begin

    end;
  AssignAtomicData;
  Result:=FProtein;
end;

procedure TPDBModel.ResetTemplates(ATemplates: TTemplates);
begin
  FTemplates:=ATemplates;
end;

function TPDBModel.ChainCount: Integer;
begin
  Result:=FProtein.GroupCount;
end;

function TPDBModel.ResidueCount(ChainIx: Integer): Integer;
begin
  Result:=FProtein.GetGroup(ChainIx).GroupCount;
end;

procedure TPDBModel.AssignAtomicData;

var
  atoms:TAtoms;
  f,tempix,atomix:Integer;

begin
  atoms:=FProtein.AllAtoms;
  for f:=0 to High(atoms) do
    begin
    atoms[f].AtomicNumber:=-1;
    atoms[f].Radius:=Config.DefaultAtomicRadius;
    tempix:=TemplateIx(atoms[f].Parent.Name);
    if tempix>=0 then
      begin
      atomix:=LastIndexOf(atoms[f].Name,FTemplates[tempix].Atoms);
      if atomix>=0 then
        begin
        atoms[f].AtomicNumber:=FTemplates[tempix].AtomicNumbers[atomix];
        if Atoms[f].AtomicNumber>0 then
          atoms[f].Radius:=AtomData[atoms[f].AtomicNumber-1].VdWradius;
        end;
      end;
    end;
end;

function TPDBModel.TemplateIx(Name: string): Integer;
begin
  Result:=High(FTemplates);
  while (Result>=0) and (Name<>FTemplates[Result].Name) do
    Dec(Result);
end;

{ TPDBModelMan }

procedure TPDBModelMan.LoadTemplates(Path: string);
//TODO: using pdb templates, change to molCIF


var
  srec:TSearchRec;
  pdbparser:TPdbReader;

procedure SetTemplate;

var f,ix:Integer;

begin
  with FTemplates[High(FTemplates)] do
    begin
    Name:=ChangeFileExt(srec.Name,'');
    SetLength(Atoms,pdbparser.AtomCount);
    SetLength(AtomicNumbers,pdbparser.AtomCount);
    SetLength(AtomIds,pdbparser.AtomCount);
    SetLength(Coords,pdbparser.AtomCount);
    SetLength(Connects,pdbparser.AtomCount);
    //fill atom data
    for f:=0 to High(Atoms) do
      begin
      Atoms[f]:=pdbparser.Atoms[f].AtomName;
      AtomIDs[f]:=pdbparser.Atoms[f].Serial;
      Coords[f]:=pdbparser.Atoms[f].Coords;
      AtomicNumbers[f]:=AtomicNumber(pdbparser.Atoms[f].Element);
      //empty connections
      Connects[f]:=nil;
      end;
    for f:=0 to High(pdbparser.Connections) do
      begin
      ix:=IndexOf(pdbparser.Connections[f].AtomId,AtomIDs);
      Connects[ix]:=Copy(pdbparser.Connections[f].Connects,0,
        Length(pdbparser.Connections[f].Connects));
      end;
  end;
end;

var f:Integer;

begin
  //remove references to templates in layers
  for f:=0 to High(FLayers) do
    FLayers[f].ResetTemplates(FTemplates);

  //clear and load
  FTemplates:=nil;
  pdbparser:=TPdbReader.Create;
  if FindFirst(Path+'*.pdb',faAnyFile,srec)=0 then
    repeat
    SetLength(FTemplates,Length(FTemplates)+1);
    pdbparser.Load(Path+srec.Name);
    SetTemplate;
    until FindNext(srec)<>0;
  FindClose(srec);
  pdbparser.Free;

  //reassign templates to layers
  for f:=0 to High(FLayers) do
    FLayers[f].ResetTemplates(FTemplates);
end;

constructor TPDBModelMan.Create(molCIFPath:string);
begin
  inherited Create;
  LoadTemplates(molCIFPath);
end;

function TPDBModelMan.Count: integer;
begin
  Result := Length(FLayers);
end;

function TPDBModelMan.LayerByIx(Ix: integer): TPDBModel;
begin
  Result := FLayers[Ix];
end;

function TPDBModelMan.AddNewLayer: TPDBModel;
begin
  Result := TPDBModel.Create(FTemplates,Length(FLayers));
  SetLength(FLayers, Length(FLayers) + 1);
  FLayers[High(FLayers)] := Result;
end;

function TPDBModelMan.LoadLayer(PdbFileName: string):TMolecule;

begin
  Result:=AddNewLayer.LoadPDB(PdbFileName);
end;

procedure TPDBModelMan.ClearLayers;

var f:Integer;

begin
  for f:=0 to High(FLayers) do
    FLayers[f].Free;
  FLayers:=nil;
end;

end.


