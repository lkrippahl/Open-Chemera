{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 9.1.2011
Purpose:
  Manages the display of molecules and other objects. Stores all display
  settings.

Requirements:
Revisions:
To do:
*******************************************************************************}

unit displayobjects;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, base3ddisplay,LCLProc, molecules, oclconfiguration;

type
  TAtomSettings=record
    Atom:TAtom;
  end;
  TAtomSettingsArray=array of TAtomSettings;
  TBondSettings=record
    BondIx:Integer;
    Molecule:TMolecule;
  end;
  TBondSettingsArray=array of TBondSettings;

  { TDisplayManager }

  TDisplayManager=class
  private
     //These objects are provided externally.
    //This form neither creates nor destroys them
    FDisplay:IBase3DDisplay;        //Interface to display window
    FLayers:TMolecules;

    //SettingsLists
    FAtoms:TAtomSettingsArray;
    FBonds:TBondSettingsArray;

  public
    constructor Create(Display:IBase3DDisplay);
    procedure SetDisplay(Display:IBase3DDisplay);

    procedure Attach(Molecule:TMolecule);
      //add a molecule to rendering
    procedure Detach(Molecule:TMolecule);
      //remove one molecule

    function LayerIndex(Molecule:TMolecule):Integer;
    procedure Render;

    procedure OnDeleteAtoms;
      //called from molecules in order to remove tagged atoms and affected bonds
      //the atoms to be deleted are tagged with molecules.TagToDelete

    procedure Test;    //this is for testing and debugging displays
  end;


implementation

{ TDisplayManager }

constructor TDisplayManager.Create(Display: IBase3DDisplay);
begin
  inherited Create;
  SetDisplay(Display);
end;

procedure TDisplayManager.SetDisplay(Display: IBase3DDisplay);
begin
  FDisplay:=Display;
end;


procedure TDisplayManager.Attach(Molecule: TMolecule);
//TO DO: rendering, display settings, etc

begin
  SetLength(FLayers,Length(FLayers)+1);
  FLayers[High(FLayers)]:=Molecule;
  DebugLn('Attached');
end;

procedure TDisplayManager.Detach(Molecule: TMolecule);
// recreates FLayers without Molecule
// TO DO: Must change the settings arrays too
// TO DO: must delete object lists on display

var
  ix,f,c:Integer;
  tmp:TMolecules;

begin
  ix:=LayerIndex(Molecule);           //Detaching is costly, so check first
  if ix>=0 then
    begin
    SetLength(tmp,Length(FLayers)-1);
    c:=0;
    for f:=0 to High(FLayers) do
      if f<>ix then
        begin
        tmp[c]:=FLayers[f];
        Inc(c);
        end;
    //TO DO: Remove display settings, display lists, etc
    FLayers:=tmp;
    end;
end;

function TDisplayManager.LayerIndex(Molecule: TMolecule): Integer;
begin
  Result:=High(FLayers);
  while (Result>=0) and (FLayers[Result]<>Molecule) do
    Dec(Result);
end;

procedure TDisplayManager.Render;
//TO DO: mostly everything. Still just for testing
//TO DO: should check which have changed
//TO DO: group objects by material

var
  f,g:Integer;
  al:TAtoms;
  ol:T3DObjectList;

begin
  FDisplay.ClearObjectLists;
  for f:=0 to High(FLayers) do
    begin
    al:=FLayers[f].AllAtoms;
    SetLength(ol.Objects,1);
    ol.Material:=DefaultMaterial;
    for g:=0 to High(al) do
        begin
        if al[g].AtomicNumber>=0 then
          with AtomData[al[g].AtomicNumber-1] do
            begin
            ol.Material.Ambient:=RGBAColor(CPKColor[0],CPKColor[1],CPKColor[2],1);
            end
          else
            ol.Material:=DefaultMaterial;
        ol.Objects[0].ObjectType:=otSphere;
        ol.Objects[0].sphC:=al[g].Coords;
        ol.Objects[0].Rad:=al[g].Radius;
        FDisplay.AddObjectList(ol);
        end;
    end;
  FDisplay.Compile;
end;

procedure TDisplayManager.OnDeleteAtoms;
begin

end;


procedure TDisplayManager.Test;

const
  count=3;
  dist=10;

var
  ol:T3DObjectList;
  x,y,z:Integer;

begin
  ol.Material:=ColorMaterial(1,0,0,1);
  ol.Objects:=nil;
  for x:=-count to count do
    for y:=-count to count do
      for z:=-count to count do
        begin
        SetLength(ol.Objects,Length(ol.Objects)+1);
        with ol.Objects[High(ol.Objects)] do
          begin
          ObjectType:=otSphere;
          Rad:=2;
          sphC:=Coord(x*dist,y*dist,z*dist);
          end;
        end;
  FDisplay.AddObjectList(ol);
  FDisplay.Compile;
  DebugLn('Test');
end;

end.

