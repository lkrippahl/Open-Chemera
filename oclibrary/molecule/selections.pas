{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 9.1.2011
Purpose:
  Stores selection lists (current selected atoms and molecules,
  and stored selections) and atom display options
Requirements:
Revisions:
To do:   everything
*******************************************************************************}

unit selections;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, molecules, base3ddisplay, basetypes,
  oclconfiguration;

const
  CPKMapName='CPK';

type

  TSelection = TAtoms;

  TBondSettings=record
    MaterialIx:Integer;
    Rad:TFloat;
    IsVisible:Boolean;
    Atom1Ix,Atom2Ix:Integer; //Indices of atoms in display settings
    Atom1,Atom2:TAtom;
  end;

  TBondSettingsArray=array of TBondSettings;

  TAtomSettings=record
    Atom:TAtom;
    AtomRad:TFloat;
    Highlight:Integer;
    IsVisible:Boolean;
    MaterialIX:Integer; //Index to material in the display manager
  end;

  TAtomSettingsArray=array of TAtomSettings;

  TMoleculeSettings=record
    AtomSettings:TAtomSettingsArray;
    BondSettings:TBondSettingsArray;
    Molecule:TMolecule;
    AtomTable,BondTable:TIntegerTable;
  end;

  TMoleculeSettingsArray=array of TMoleculeSettings;

  //to map atom properties to material indexes
  TMaterialMap=record
    MapName:string;
    MaterialIxs:TIntegers;
    MinVal,Mult:TFloat;
  end;

  TMaterialMaps=array of TMaterialMap;

  { TSettingsManager }

  TSettingsManager=class
  protected
    FMolSettings:TMoleculeSettingsArray;
    FMaterials:T3DMaterials;
    FMaps:TMaterialMaps;
    function MaterialIndex(Mat:T3DMaterial):Integer;overload;
    function MaterialIndex(const Map:TMaterialMap;Value:TFloat):Integer;overload;
    function MapIndex(Name:string):Integer;
    procedure AddMaterial(const Mat:T3DMaterial);
  public
    property Maps:TMaterialMaps read FMaps;
    property MolSettings:TMoleculeSettingsArray read FMolSettings;
    property Materials:T3DMaterials read FMaterials;
    constructor Create;
    procedure Free;
    function IndexFromMolecule(Mol:TMolecule):Integer;
    procedure AttachMolecule(AMolecule:TMolecule);
    procedure DetachMolecule(AMolecule:TMolecule);
    procedure CreateCPKMaterials;
    procedure AssignCPK(Molecule:TMolecule);
    procedure CreateMaterialTable(Molecule:TMolecule);
  end;

function DefaultAtomSettings(AAtom:TAtom):TAtomSettings;
function DefaultBondSettings(ABond:TAtomBond):TBondSettings;
  //Assumes atoms in bond are tagged with their indeces on the settings
function DefaultMolSettings:TMoleculeSettings;

implementation

function DefaultAtomSettings(AAtom:TAtom): TAtomSettings;
begin
  with Result do
    begin
    Highlight:=-1;
    MaterialIx:=-1;
    AtomRad:=0;
    Atom:=AAtom;
    IsVisible:=True;
    if AAtom<>nil then
      begin
      AtomRad:=AAtom.Radius;
      end;
    end;
end;

function DefaultBondSettings(ABond:TAtomBond): TBondSettings;

//Assumes atoms in bond are tagged with their indeces on the settings

begin
  with Result do
    begin
    MaterialIx:=-1;
    Rad:=0;
    IsVisible:=True;
    Atom1:=ABond.Atom1;
    Atom2:=ABond.Atom2;
    Atom1Ix:=Atom1.Tag;
    Atom2Ix:=Atom2.Tag;
    end;
end;

function DefaultMolSettings: TMoleculeSettings;
begin
  with Result do
    begin
    AtomSettings:=nil;
    BondSettings:=nil;
    Molecule:=nil;
    AtomTable:=nil;
    BondTable:=nil;
    end;
end;

{ TSettingsManager }

function TSettingsManager.IndexFromMolecule(Mol: TMolecule): Integer;
begin
  Result:=High(FMolSettings);
  while (Result>=0) and (FMolSettings[Result].Molecule<>Mol) do
    Dec(Result);
end;

function TSettingsManager.MaterialIndex(Mat: T3DMaterial): Integer;
begin
  Result:=High(FMaterials);
  while (Result>=0) and not IsSameMaterial(FMaterials[Result],Mat) do
    Dec(Result);
end;

function TSettingsManager.MaterialIndex(const Map: TMaterialMap; Value: TFloat): Integer;

var ix:Integer;

begin
  with Map do
    begin
    ix:=Round((Value-MinVal)*Mult);
    if ix<0 then Result:=0
    else if ix>High(MaterialIxs) then Result:=High(MaterialIxs)
    else Result:=Ix;
    end;
end;

function TSettingsManager.MapIndex(Name: string): Integer;
begin
  Result:=High(FMaps);
  while (Result>=0) and (FMaps[Result].MapName<>Name) do
    Dec(Result);
end;

procedure TSettingsManager.AddMaterial(const Mat: T3DMaterial);
begin
  SetLength(FMaterials,Length(FMaterials)+1);
  FMaterials[High(FMaterials)]:=Mat;
end;

constructor TSettingsManager.Create;
begin
  inherited;
  FMolSettings:=nil;
end;

procedure TSettingsManager.Free;
begin
  if Self<>nil then
    begin
    //cleanup stuff here
    inherited Free;
    end;
end;

procedure TSettingsManager.AttachMolecule(AMolecule: TMolecule);

var
  atoms:TAtoms;
  bonds:TAtomBonds;
  f:Integer;

begin
  SetLength(FMolSettings,Length(FMolSettings)+1);
  FMolSettings[High(FMolSettings)]:=DefaultMolSettings;
  with FMolSettings[High(FMolSettings)] do
    begin
    Molecule:=AMolecule;
    atoms:=AMolecule.AllAtoms;
    SetLength(AtomSettings,Length(atoms));

    for f:=0 to High(atoms) do
      begin
      AtomSettings[f]:=DefaultAtomSettings(atoms[f]);
      atoms[f].Tag:=f;      //use tags to index atoms
      end;

    bonds:=AMolecule.AllBonds;
    SetLength(BondSettings,Length(bonds));
    for f:=0 to High(bonds) do
      BondSettings[f]:=DefaultBondSettings(bonds[f]);

    end;
end;

procedure TSettingsManager.DetachMolecule(AMolecule: TMolecule);

var
  f, i, ix:Integer;
  tmp:TMoleculeSettingsArray;

begin
  ix:=IndexFromMolecule(AMolecule);
  if ix>=0 then
    begin
    i:=0;
    SetLength(tmp,Length(FMolSettings)-1);
    for f:=0 to High(FMolSettings) do
      if f<>ix then
        begin
        tmp[i]:=FMolSettings[f];
        Inc(i);
        end;
    FMolSettings[ix].AtomSettings:=nil;
    FMolSettings[ix].BondSettings:=nil;
    FMolSettings:=tmp;
    end;
end;

procedure TSettingsManager.CreateCPKMaterials;

var
  f,ix,matix:Integer;
  mat:T3DMaterial;

begin
  ix:=MapIndex(CPKMapName);
  if ix<0 then
    begin
    SetLength(FMaps,Length(FMaps)+1);
    ix:=High(FMaps);
    FMaps[ix].MapName:=CPKMapName;
    FMaps[ix].MinVal:=1;    //To map directly to atomic numbers
    FMaps[ix].Mult:=1;
    SetLength(FMaps[ix].MaterialIxs,Length(AtomData));
    for f:=0 to High(AtomData) do
      begin
      with AtomData[f] do
        begin
        mat.Ambient:=RGBAColor(0.1,0.1,0.1,1);
        mat.Specular:=RGBAColor(1,1,1,1); //RGBAColor(CPKColor[0],CPKColor[1],CPKColor[2],1);   //RGBAColor(1,1,1,1);
        mat.Diffuse:=RGBAColor(CPKColor[0],CPKColor[1],CPKColor[2],1);
        mat.Shine:=128;
        mat.Texture:='';
        end;
      matix:=MaterialIndex(mat);
      if matix<0 then
        begin
        AddMaterial(mat);
        matix:=High(FMaterials);
        end;
      FMaps[ix].MaterialIxs[f]:=matix;
      end;
    end;
end;

procedure TSettingsManager.AssignCPK(Molecule: TMolecule);

var
  f,ix,mapix:Integer;
  atomicnumber:Integer;

begin
  ix:=IndexFromMolecule(Molecule);
  if ix>=0 then
    begin
    mapix:=MapIndex(CPKMapName);
    if mapix<0 then
      begin
      CreateCPKMaterials;
      mapix:=MapIndex(CPKMapName);
      end;
    with FMolSettings[ix] do
      with FMaps[mapix] do
        for f:=0 to High(AtomSettings) do
          with AtomSettings[f] do
            begin
            atomicnumber:=Atom.AtomicNumber;
            if atomicnumber<0 then atomicnumber:=1;
            MaterialIX:=MaterialIxs[Round(Mult*(atomicnumber-MinVal))];
            end;
    end;
end;

procedure TSettingsManager.CreateMaterialTable(Molecule: TMolecule);

var
  f,ix,matix:Integer;
  count:TIntegers;

begin
  ix:=IndexFromMolecule(Molecule);
  if ix>=0 then
    with FMolSettings[ix] do
      begin
      SetLength(AtomTable,Length(FMaterials),Length(AtomSettings));
      SetLength(count,Length(FMaterials));
      for f:=0 to High(count) do
        count[f]:=0;
      for f:=0 to High(AtomSettings) do
        begin
        matix:=AtomSettings[f].MaterialIX;
        AtomTable[matix,count[matix]]:=f;
        Inc(count[matix]);
        end;
      for f:=0 to High(count) do
        SetLength(AtomTable[f],count[f]);
      end;
end;

end.

