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

  PDB data is being replaced with template data (charges, elements, ...). This
  may not be desireable (in LoadPDB)

*******************************************************************************}

unit pdbmolecules;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, molecules, pdbparser, FileUtil, oclconfiguration,
  stringutils,LCLProc;

type

  PDBChargeOrigin=(pdbNone,pdbOccTemp,pdbOccupancy,pdbCharge);
  PDBResidueTypes=set of (resAA,resNonAA);

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
    property FileName:string read FFileName;
    constructor Create(Templates: TTemplates;AID:Integer);
    procedure Free;
    procedure ClearChains;
    procedure CreateChains(const IDs: TSimpleStrings);
    procedure DeleteResidues(Options:PDBResidueTypes;
                             ResTypes:TSimpleStrings=nil;
                             OnDelete:TOnDeleteCallback=nil);overload;
    procedure DeleteResidues(Name: string; OnDelete: TOnDeleteCallback=nil);overload;
    procedure DeleteNonAAResidues(OnDelete: TOnDeleteCallback=nil);
    procedure DeleteWater(OnDelete: TOnDeleteCallback=nil);
    procedure DeleteChainsByName(ChainName:string;OnDelete:TOnDeleteCallback=nil);
    function NewEmptyChain(ChainName: string; ChainID: integer): TMolecule;
    function NewChain(ChainName: string; ChainID, Size: integer): TMolecule;
    function GetChain(ChainIx: integer): TMolecule;overload;
    function GetChain(ChainID:string):TMolecule;overload;
    procedure AppendChain(Chain:TMolecule);
    function GetResidue(ChainIx, ResIx: integer): TMolecule;
    function GetResidue(ChainName:string; ResId:Integer):TMolecule;
    function GetAtom(ChainIx, ResIx, AtomIx: integer): TAtom;
    function LoadPDB(PdbFileName: string; ChargeFrom:PDBChargeOrigin=pdbNone):TMolecule;
    procedure ResetTemplates(ATemplates:TTemplates);
    function ChainCount:Integer;
    function ListChains:TSimpleStrings;
    procedure RenameChains(FirstChar:Char='A';LastChar:Char='Z');
    function ResidueCount(ChainIx:Integer):Integer;
    procedure AssignAtomicData; //atomic number, vdW radius
    function TemplateIx(Name:string):Integer;
    function CopyChains(Chains:TSimpleStrings):TMolecule;
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
    function LayerByFileName(FileName: string): TPDBModel;
    function AddNewLayer: TPDBModel;
    function LoadLayer(PdbFileName: string; ChargeFrom:PDBChargeOrigin=pdbNone):TMolecule;
    function GetChains(Layer: Integer; Indexes: TIntegers): TMolecules;overload;
    function GetChains(Layer: Integer; IDs: TSimpleStrings): TMolecules;overload;
    procedure ClearLayers;
    procedure Free;

  end;

  function AtomIsAABackbone(Atom:TAtom):Boolean;
  function BackboneOnly(Atoms:TAtoms):TAtoms;
  function NoBackBone(Atoms:TAtoms):TAtoms;
  function ResidueIsAminoAcid(Residue:TMolecule):Boolean;
  function ChainSequence(Chain:TMolecule;MissingMarker:string='X'):string;
    //returns string with one letter aa code for each residue. Non AA residues
    //are and gaps in the sequence of residue IDs are filled with the MissingMarker
  procedure SaveToPDB(Molecule:TMolecule;FileName:string);
  function GetResidue(const Protein:TMolecule;const ChainName:string;const ResId:Integer):TMolecule;
    //assumes Protein is a protein, with chains and residues

implementation

function AtomIsAABackbone(Atom:TAtom):Boolean;

begin
  with Atom do
    Result:= (Name='N') or (Name='O') or (Name='CA') or (Name='C');
end;

function BackboneOnly(Atoms: TAtoms): TAtoms;

var f,c:Integer;

begin
  SetLength(Result,Length(Atoms));
  c:=0;
  for f:=0 to High(Atoms) do
    if AtomIsAABackbone(Atoms[f]) then
      begin
      Result[c]:=Atoms[f];
      Inc(c);
      end;
  SetLength(Result,c);
end;

function NoBackBone(Atoms: TAtoms): TAtoms;

var f,c:Integer;

begin
  SetLength(Result,Length(Atoms));
  c:=0;
  for f:=0 to High(Atoms) do
    if not AtomIsAABackbone(Atoms[f]) then
      begin
      Result[c]:=Atoms[f];
      Inc(c);
      end;
  SetLength(Result,c);
end;

function ResidueIsAminoAcid(Residue: TMolecule): Boolean;
begin
  Result:=(Residue.GroupCount=0) and (AAOneLetterCode(Residue.Name)<>'');
end;

function ChainSequence(Chain:TMolecule;MissingMarker:string='X'):string;

var
  f,previd:Integer;
  res:TMolecule;
  tmp:string;
begin
  Result:='';
  if Chain.GroupCount>0 then
    begin
    previd:=Chain.GetGroup(0).ID-1;
    for f:=0 to Chain.GroupCount-1 do
      begin
      Inc(previd);
      res:=Chain.GetGroup(f);
      while previd<res.ID do
        begin
        Inc(previd);
        Result:=Result+MissingMarker;{ TODO : this is very inneficient. However, not sure if worth improving, since this should never be called many times for the same pdb file }
        end;
      tmp:=AAOneLetterCode(res.Name);
      if tmp='' then
        Result:=Result+MissingMarker
      else Result:=Result+tmp;
      end;
    end;
end;

procedure SaveToPDB(Molecule: TMolecule; FileName: string);

var
  sl:TStringList;
  atoms:TAtoms;
  f:Integer;
  res,chain:TMolecule;
  oldchain:TMolecule;
  rname,chname:string;
  rid:Integer;
begin
  oldchain:=nil;
  sl:=TStringList.Create;
  atoms:=Molecule.AllAtoms;
  for f:=0 to High(atoms) do
    begin
    res:=atoms[f].Parent;
    rname:='';
    rid:=0;
    chname:='';
    if res<>nil then
      begin
      rname:=res.Name;
      rid:=res.ID;
      chain:=res.Parent;
      if chain<>nil then
        begin
        chname:=chain.Name;
        if (oldchain<>nil) and (oldchain<>chain) then
          sl.Add('TER       ');
        oldchain:=chain;
        end;
      end;
    sl.Add(AtomRecord(atoms[f].Name,rname,chname,atoms[f].ID,rid,atoms[f].Coords,
    Element(atoms[f].AtomicNumber)));
    end;
  sl.SaveToFile(FileName);
  sl.Free;

end;

function GetResidue(const Protein: TMolecule; const ChainName: string;
  const ResId: Integer): TMolecule;
begin
  Result:=Protein.GetGroup(ChainName);
  if Result<>nil then Result:=Result.GetGroupById(ResId);
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

procedure TPDBModel.DeleteResidues(Options: PDBResidueTypes;
  ResTypes: TSimpleStrings;OnDelete:TOnDeleteCallback=nil);

var
  c,r:Integer;
  chain,res:TMolecule;

  function ResidueToDelete(Residue:TMolecule):Boolean;

  begin
    Result:= ((resAA in Options) and (AAOneLetterCode(Residue.Name)<>'')) or
             ((resNonAA in Options) and (AAOneLetterCode(Residue.Name)='')) or
             (LastIndexOf(Residue.Name,ResTypes)>=0);
  end;

begin
  for c:=0 to FProtein.GroupCount-1 do
    begin
    chain:=FProtein.GetGroup(c);
    chain.TagAllAtoms(0);
    for r:=0 to chain.GroupCount-1 do
      begin
      res:=chain.GetGroup(r);
      if ResidueToDelete(res) then
        res.TagAllAtoms(1);
      end;
    chain.DeleteTaggedAtoms(1,OnDelete);
    end;
  FProtein.DeleteEmptyGroups;
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

begin
  Result:=FProtein.GetGroup(ChainId);
end;

procedure TPDBModel.AppendChain(Chain: TMolecule);
begin
  FProtein.AddGroup(Chain);
end;

function TPDBModel.GetResidue(ChainIx, ResIx: integer): TMolecule;
begin
  Result := FProtein.GetGroup(ChainIx);
  if Result<>nil then Result:=Result.GetGroup(ResIx);
end;

function TPDBModel.GetResidue(ChainName: string; ResId: Integer): TMolecule;
begin
  Result := GetChain(ChainName);
  if Result<>nil then Result:=Result.GetGroupById(ResId);
end;

function TPDBModel.GetAtom(ChainIx, ResIx, AtomIx: integer): TAtom;

var mol:TMolecule;

begin
  Result:=nil;
  mol := FProtein.GetGroup(ChainIx);
  if mol<>nil then mol:=mol.GetGroup(ResIx);
  if mol <>nil then Result := mol.GetAtom(AtomIx);
end;

function TPDBModel.LoadPDB(PdbFileName: string; ChargeFrom: PDBChargeOrigin
  ): TMolecule;

var
  parser: TPDBReader;
  cr,f,cc: integer;
  cres:TMolecule;
  atom:TAtom;

begin
  ClearChains;
  FFileName:=PdbFileName;
  parser := TPDBReader.Create(PdbFileName);
  FProtein.Name:=ChangeFileExt(ExtractFileName(PdbFileName),'');
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
      //element may be present in the PDB file
      //However, this is superseded if there is monomer template data
      atom.AtomicNumber:=AtomicNumber(Element);
      atom.Radius:=VdWRadius(atom.AtomicNumber);
      if ChargeFrom=pdbOccupancy then
        atom.Charge:=Occupancy
      else if ChargeFrom=pdbOccTemp then
        atom.Charge:=OccTemp
      else if ChargeFrom=pdbCharge then
        begin
        if Length(Charge)=2 then
          begin
          atom.Charge:=StrToFloat(Charge[1]);
          if Charge[2]='-' then atom.Charge:=-atom.Charge;
          end;
        end;
      end;

  { TODO : make extra connections }
  for f:=0 to High(parser.Connections) do
    begin

    end;
  //assign element, radius, charge, et al, from the monomer templates
  //will replace information on PDB
  { TODO : It may not be desireable to discard charges and element data on PDB in favour of templates. Check this... }
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

function TPDBModel.ListChains: TSimpleStrings;

begin
  Result:=FProtein.ListGroupNames;
end;

procedure TPDBModel.RenameChains(FirstChar: Char; LastChar: Char);

var
  f:Integer;
  name:Char;

begin
  name:=FirstChar;
  for f:=0 to FProtein.GroupCount-1 do
    begin
    FProtein.GetGroup(f).Name:=name;
    name:=Succ(name);
    if name>LastChar then
      name:=FirstChar;
    end;
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

function TPDBModel.CopyChains(Chains: TSimpleStrings): TMolecule;

var f:Integer;

begin
  if Chains=nil then
    Result:=nil
  else
    begin
    Result:=TMolecule.Create(FProtein.Name,FProtein.ID,nil);
    for f:=0 to High(Chains) do
      Result.AddGroup(TMolecule.CopyFrom(FProtein.GetGroup(Chains[f]),Result));
    end;
end;

procedure TPDBModel.DeleteResidues(Name: string;OnDelete:TOnDeleteCallback=nil);

var
  names:TSimpleStrings;

begin
  SetLength(names,1);
  names[0]:=Name;
  DeleteResidues([],names,OnDelete);
end;

procedure TPDBModel.DeleteNonAAResidues(OnDelete:TOnDeleteCallback=nil);

begin
  DeleteResidues([resNonAA],nil,OnDelete);
end;

procedure TPDBModel.DeleteWater(OnDelete:TOnDeleteCallback=nil);
begin
  DeleteResidues('HOH',OnDelete);
end;

procedure TPDBModel.DeleteChainsByName(ChainName: string;OnDelete:TOnDeleteCallback=nil);
var
  c:Integer;
  chain:TMolecule;

begin
  FProtein.TagAllAtoms(0);
  for c:=0 to FProtein.GroupCount-1 do
    begin
    chain:=FProtein.GetGroup(c);
    if chain.Name=ChainName then
      chain.TagAllAtoms(1);
    end;
  FProtein.DeleteTaggedAtoms(1,OnDelete);
  FProtein.DeleteEmptyGroups;
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

function TPDBModelMan.LayerByFileName(FileName: string): TPDBModel;

var f:Integer;

begin
  Result:=nil;
  for f:=0 to High(FLayers) do
    if FLayers[f].FileName=FileName then
      begin
      Result:=FLayers[f];
      Break;
      end;
end;

function TPDBModelMan.AddNewLayer: TPDBModel;
begin
  Result := TPDBModel.Create(FTemplates,Length(FLayers));
  SetLength(FLayers, Length(FLayers) + 1);
  FLayers[High(FLayers)] := Result;
end;

function TPDBModelMan.LoadLayer(PdbFileName: string; ChargeFrom: PDBChargeOrigin
  ): TMolecule;

begin
  Result:=AddNewLayer.LoadPDB(PdbFileName, ChargeFrom);
end;

function TPDBModelMan.GetChains(Layer: Integer; Indexes: TIntegers): TMolecules;

var f:Integer;

begin
  SetLength(Result,Length(Indexes));
  for f:=0 to High(Result) do
    Result[f]:=LayerByIx(Layer).GetChain(Indexes[f]);
end;

function TPDBModelMan.GetChains(Layer: Integer; IDs: TSimpleStrings
  ): TMolecules;

var f:Integer;

begin
  SetLength(Result,Length(IDs));
  for f:=0 to High(Result) do
    Result[f]:=LayerByIx(Layer).GetChain(IDs[f]);
end;

procedure TPDBModelMan.ClearLayers;

var f:Integer;

begin
  for f:=0 to High(FLayers) do
    FLayers[f].Free;
  FLayers:=nil;
end;

procedure TPDBModelMan.Free;
begin
  if Self<>nil then
    begin
    ClearLayers;
    inherited;
    end;
end;

end.


