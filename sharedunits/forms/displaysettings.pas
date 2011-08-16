{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 10.1.2011
Purpose:
  This form interfaces with the DisplayManager settings
Requirements:
Revisions:
To do:
  Add materials, display options
  Add support for other objects and molecules.
  Organize links to molecules to allow adding or removing one withou messing
  the display of the rest
*******************************************************************************}

unit displaysettings;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, FileUtil, Forms, Controls, Graphics, Dialogs,
  basetypes, molecules,base3ddisplay, displayobjects;

type

  { TDisplaySettingsForm }

  TDisplaySettingsForm = class(TForm)
  private
    { private declarations }
    FDispMan:TDisplayManager;


  public
    { public declarations }
    constructor Create(TheOwner:TComponent;DisplayManager: TDisplayManager);reintroduce;

  end;

implementation

{$R *.lfm}

{ TDisplaySettingsForm }

constructor TDisplaySettingsForm.Create(TheOwner:TComponent;
  DisplayManager: TDisplayManager);
begin
  inherited Create(TheOwner);
  FDispMan:=DisplayManager;
end;


end.

