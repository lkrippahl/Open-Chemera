{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain. Enjoy.
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
  Classes, SysUtils, FileUtil, Forms, Controls, Graphics, Dialogs, Menus, pdbmolecules,glutwindow;

type

  { TCmMainForm }

  TCmMainForm = class(TForm)
    MainMenu1: TMainMenu;
    MenuItem1: TMenuItem;
    procedure FormActivate(Sender: TObject);
    procedure FormCreate(Sender: TObject);
  private
    { private declarations }
    FInitializing:Boolean;
    FDisplay:TGlutWindow;
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

procedure TCmMainForm.InitChemera;
begin
  FDisplay:=TGlutWindow.Create(Top+Height,Left,600,Width);
end;

procedure TCmMainForm.FormActivate(Sender: TObject);
begin
  if FInitializing then InitChemera;
end;

end.

