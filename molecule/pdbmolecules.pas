{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 9.1.2011
Purpose:
  Class for storing PDB-like file data.
  TPDBLayer corresponds to a file. Contains chains, residues, etc.
  TPDBLayerSet contains several layers.
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

*******************************************************************************}

unit pdbmolecules;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, molecules, pdbparser, FileUtil;

type

  TTemplate = record                    //residue templates
    Name: string;
    Atoms: array of string[4];
    Elements:array of string[2];
    AtomIDs: TOCIntegers;
    Coords: TOCCoords;
    Connects: array of TOCIntegers;
  end;

  TTemplates = array of TTemplate;

  { TPDBLayer }

  TPDBLayer = class
  protected
    FTemplates: TTemplates;
    FProtein: TMolecule;
    FFileName: string;
    FInfo: TPDBLayerInfo;
  public
    property Info: TPDBLayerInfo read FInfo write FInfo;
    property Molecule:TMolecule read FProtein;
    constructor Create(Templates: TTemplates;AID:Integer);
    procedure Free;
    procedure ClearChains;
    procedure CreateChains(const IDs: TOCStrings);
    function NewEmptyChain(ChainName: string; ChainID: integer): TMolecule;
    function NewChain(ChainName: string; ChainID, Size: integer): TMolecule;
    function GetChain(ChainIx: integer): TMolecule;
    function GetResidue(ChainIx, ResIx: integer): TMolecule;
    function GetAtom(ChainIx, ResIx, AtomIx: integer): TAtom;
    function LoadPDB(FileName: string):TMolecule;
    procedure ResetTemplates(ATemplates:TTemplates);
  end;

  TPDBLayers = array of TPDBLayer;

  { TPDBLayerSet }

  { TPDBLayerMan }

  TPDBLayerMan = class
  protected
    FLayers: TPDBLayers;
    FTemplates: TTemplates;
    procedure LoadTemplates(Folder: string);
  public
    constructor Create(TemplateFolder: string);
    function Count: integer;
    function LayerByIx(Ix: integer): TPDBLayer;
    function AddNewLayer: TPDBLayer;
    function LoadLayer(PdbFileName: string):TMolecule;
    procedure ClearLayers;
  end;

implementation

{ TPDBLayer }

constructor TPDBLayer.Create(Templates:TTemplates;AID:Integer);
begin
  inherited Create;
  FTemplates:=Templates;
  FProtein:=TMolecule.Create('',AID,nil);
end;

procedure TPDBLayer.Free;

begin
  if Self <> nil then
  begin
    ClearChains;
    inherited;
  end;
end;

procedure TPDBLayer.ClearChains;

begin
  //TODO Add on delete callbacks
  FProtein.Cleanup;
end;

procedure TPDBLayer.CreateChains(const IDs: TOCStrings);

var
  f:integer;

begin
  for f := 0 to High(IDs) do
    FProtein.NewGroup(IDs[f], f + 1);
end;

function TPDBLayer.NewEmptyChain(ChainName: string; ChainID: integer): TMolecule;
begin
  Result:=FProtein.NewGroup(ChainName, ChainId);
end;

function TPDBLayer.NewChain(ChainName: string; ChainID, Size: integer): TMolecule;

begin
  Result := NewEmptyChain(ChainName, ChainId);
  Result.CreateEmptyGroups(Size);
end;


function TPDBLayer.GetChain(ChainIx: integer): TMolecule;
begin
  Result := FProtein.GetGroup(ChainIx);
end;

function TPDBLayer.GetResidue(ChainIx, ResIx: integer): TMolecule;
begin
  Result := FProtein.GetGroup(ChainIx);
  if Result<>nil then Result:=Result.GetGroup(ResIx);
end;

function TPDBLayer.GetAtom(ChainIx, ResIx, AtomIx: integer): TAtom;

var mol:TMolecule;

begin
  Result:=nil;
  mol := FProtein.GetGroup(ChainIx);
  if mol<>nil then mol:=mol.GetGroup(ResIx);
  if mol <>nil then Result := mol.GetAtom(AtomIx);
end;

function TPDBLayer.LoadPDB(FileName: string):TMolecule;

var
  parser: TPDBReader;
  cr,f,cc: integer;
  cres:TMolecule;
  atom,atom2:TAtom;

begin
  ClearChains;
  parser := TPDBReader.Create(FileName);
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
  Result:=FProtein;
end;

procedure TPDBLayer.ResetTemplates(ATemplates: TTemplates);
begin
  FTemplates:=ATemplates;
end;

{ TPDBLayerMan }

procedure TPDBLayerMan.LoadTemplates(Folder: string);

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
    SetLength(Elements,pdbparser.AtomCount);
    SetLength(AtomIds,pdbparser.AtomCount);
    SetLength(Coords,pdbparser.AtomCount);
    SetLength(Connects,pdbparser.AtomCount);
    //fill atom data
    for f:=0 to High(Atoms) do
      begin
      Atoms[f]:=pdbparser.Atoms[f].AtomName;
      AtomIDs[f]:=pdbparser.Atoms[f].Serial;
      Coords[f]:=pdbparser.Atoms[f].Coords;
      Elements[f]:=pdbparser.Atoms[f].Element;
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
  Folder:=AppendPathDelim(Folder);
  if FindFirst(Folder+'*.pdb',faAnyFile,srec)=0 then
    repeat
    SetLength(FTemplates,Length(FTemplates)+1);
    pdbparser.Load(Folder+srec.Name);
    SetTemplate;
    until FindNext(srec)<>0;
  FindClose(srec);
  pdbparser.Free;

  //reassign templates to layers
  for f:=0 to High(FLayers) do
    FLayers[f].ResetTemplates(FTemplates);
end;

constructor TPDBLayerMan.Create(TemplateFolder:string);
begin
  inherited Create;
  LoadTemplates(TemplateFolder);
end;

function TPDBLayerMan.Count: integer;
begin
  Result := Length(FLayers);
end;

function TPDBLayerMan.LayerByIx(Ix: integer): TPDBLayer;
begin
  Result := FLayers[Ix];
end;

function TPDBLayerMan.AddNewLayer: TPDBLayer;
begin
  Result := TPDBLayer.Create(FTemplates,Length(FLayers));
  SetLength(FLayers, Length(FLayers) + 1);
  FLayers[High(FLayers)] := Result;
end;

function TPDBLayerMan.LoadLayer(PdbFileName: string):TMolecule;

begin
  Result:=AddNewLayer.LoadPDB(PdbFileName);
end;

procedure TPDBLayerMan.ClearLayers;

var f:Integer;

begin
  for f:=0 to High(FLayers) do
    FLayers[f].Free;
  FLayers:=nil;
end;

end.


