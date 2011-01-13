{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain. Enjoy.
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
Revisions:
To do: Comments
*******************************************************************************}

unit pdbmolecules;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, molecules;

type

  TPDBLayerInfo=record
    Header,Title,Compound,Source,
    KeyWords,ExpTechnique,Author,Journal,Remarks:TOCStrings;
    UserComments:TOCStrings;
    end;

  { TPDBLayer }

  TPDBLayer=class
  protected
    FChains:TMolecules;
    FFileName:string;
    FInfo:TPDBLayerInfo;
  public
    property Info:TPDBLayerInfo read FInfo write FInfo;
    constructor Create;
    procedure Free;
    procedure ClearChains;
    function NewEmptyChain(ChainName:string;ChainID:Integer):TMolecule;
    function NewChain(ChainName:string;ChainID,Size:Integer):TMolecule;
    function ChainIndex(Chain:TMolecule):Integer;
    function GetChain(ChainIx:Integer):TMolecule;
    function GetResidue(ChainIx,ResIx:Integer):TMolecule;
    function GetAtom(ChainIx,ResIx,AtomIx:Integer):TAtom;
  end;

  TPDBLayers=array of TPDBLayer;

  { TPDBLayerSet }

  TPDBLayerSet=class
  protected
    FLayers:TPDBLayers;

  public
    constructor Create;
    function Count:Integer;
    function LayerByIx(Ix:Integer):TPDBLayer;
    function AddNewLayer:TPDBLayer;

  end;

implementation

{ TPDBLayer }

constructor TPDBLayer.Create;
begin
  inherited;
end;

procedure TPDBLayer.Free;

begin
  if Self<>nil then
    begin
    ClearChains;
    inherited;
    end;
end;

procedure TPDBLayer.ClearChains;

var f:Integer;

begin
  for f:=0 to High(FChains) do FChains[f].Free;
end;

function TPDBLayer.NewEmptyChain(ChainName:string;ChainID:Integer): TMolecule;
begin
  SetLength(FChains, Length(FChains)+1);
  Result:=TMolecule.Create(ChainName,ChainId,nil);
  FChains[High(FChains)]:=Result;
end;

function TPDBLayer.NewChain(ChainName:string;ChainID,Size: Integer): TMolecule;

var chain:TMolecule;

begin
  chain:=NewEmptyChain(ChainName,ChainId);
  chain.CreateEmptyGroups(Size);
end;

function TPDBLayer.ChainIndex(Chain: TMolecule): Integer;
begin
  Result:=High(FChains);
  while (Result>=0) and (FChains[Result]<>Chain) do
    Dec(Result);
end;

function TPDBLayer.GetChain(ChainIx: Integer): TMolecule;
begin
  Result:=FChains[ChainIx];
end;

function TPDBLayer.GetResidue(ChainIx, ResIx: Integer): TMolecule;
begin
  Result:=FChains[ChainIx].GetGroup(ResIx);
end;

function TPDBLayer.GetAtom(ChainIx, ResIx, AtomIx: Integer): TAtom;
begin
  Result:=FChains[ChainIx].GetGroup(ResIx).GetAtom(AtomIx);
end;

{ TPDBLayerSet }

constructor TPDBLayerSet.Create;
begin
  inherited;
end;

function TPDBLayerSet.Count: Integer;
begin
  Result:=Length(FLayers);
end;

function TPDBLayerSet.LayerByIx(Ix: Integer): TPDBLayer;
begin
  Result:=FLayers[Ix];
end;

function TPDBLayerSet.AddNewLayer: TPDBLayer;
begin
  Result:=TPDBLayer.Create;
  SetLength(FLayers,Length(FLayers)+1);
  FLayers[High(FLayers)]:=Result;
end;

end.

