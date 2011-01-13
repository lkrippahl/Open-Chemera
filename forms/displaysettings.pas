{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain. Enjoy.
********************************************************************************
Author: Ludwig Krippahl
Date: 10.1.2011
Purpose:
  This form stores all display options for showing atoms and bonds.
Requirements:
Revisions:
To do:
  Add materials, display options
  Add support for other objects and molecules.
*******************************************************************************}

unit displaysettings;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, FileUtil, Forms, Controls, Graphics, Dialogs,
  basetypes, molecules, pdbmolecules,base3ddisplay;

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

  { TDisplaySettingsForm }

  TDisplaySettingsForm = class(TForm)
  private
    { private declarations }

    //These objects are provided externally.
    //This form neither creates nor destroys them
    FDisplay:IBase3DDisplay;        //Interface to display window
    FLayerSet:TPDBLayerSet;         //Where molecules are stored

    //SettingsLists
    FAtoms:TAtomSettingsArray;
    FBonds:TBondSettingsArray;

  public
    { public declarations }
    constructor Create(Display:IBase3DDisplay;Layers:TPDBLayerSet);
    procedure Render;
    //called from molecules in order to remove tagged atoms and affected bonds
    //the atoms to be deleted are tagged with molecules.TagToDelete
    procedure OnDeleteAtoms;
  end; 

var
  DisplaySettingsForm: TDisplaySettingsForm;

implementation

{$R *.lfm}

{ TDisplaySettingsForm }

constructor TDisplaySettingsForm.Create(Display: IBase3DDisplay;
  Layers: TPDBLayerSet);
begin

end;

procedure TDisplaySettingsForm.Render;
begin

end;

procedure TDisplaySettingsForm.OnDeleteAtoms;
begin

end;

end.

