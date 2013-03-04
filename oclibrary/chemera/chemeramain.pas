{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 9.1.2011
Purpose:
  This is the main form for Chemera.
Requirements:
Revisions:
To do:
*******************************************************************************}

unit chemeramain;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, FileUtil, Forms, Controls, Graphics, Dialogs, Menus,
  pdbmolecules, oglform, displayobjects, LCLProc, oclconfiguration;

type

  { TCmMainForm }

  TCmMainForm = class(TForm)
    MainMenu: TMainMenu;
    FileMn: TMenuItem;
    OpenDlg: TOpenDialog;
    OpenMn: TMenuItem;
    procedure FormActivate(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure OpenMnClick(Sender: TObject);
  private
    { private declarations }
    FInitializing:Boolean;
    FDisplay:TOpenGlForm;
    FDispMan:TDisplayManager;
    FMolecules:TPdbLayerMan;
    procedure InitChemera;
  public
    { public declarations }
  end; 

var
  CmMainForm: TCmMainForm;

implementation

{$R *.lfm}

{ TCmMainForm }

procedure TCmMainForm.FormCreate(Sender: TObject);
begin
  FInitializing:=True;
end;

procedure TCmMainForm.OpenMnClick(Sender: TObject);


begin
  OpenDlg.Filter:='';//Pdb file|*.pdb;GZipped PDB file|*.gz';
  if OpenDlg.Execute then
    begin
    FDispMan.Attach(
      FMolecules.LoadLayer(OpenDlg.FileName));
    FDispMan.Render;
    end;
end;

procedure TCmMainForm.InitChemera;
begin
  FInitializing:=False;

  //Load configuration data
  LoadAtomData;
  LoadAAData;
  FDisplay:=TOpenGLForm.Create(Self);
  FDisplay.Show;
  FDispMan:=TDisplayManager.Create(FDisplay);

  //DEBUG
  //FDispMan.Test;

  FMolecules:=TPdbLayerMan.Create(Config.MonomersPath);
end;

procedure TCmMainForm.FormActivate(Sender: TObject);
begin
  if FInitializing then InitChemera;
end;

end.

