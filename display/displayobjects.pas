{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain. Enjoy.
********************************************************************************
Author: Ludwig Krippahl
Date: 9.1.2011
Purpose:
  Manages objects to be displayed in the 3D renderer. A TDisplayObjects object
  is given to the 3D window, which then renders the objects.
  This class is used by the display managers; it does not read directly from
  molecules or other data.
Requirements:
Revisions:
To do:
*******************************************************************************}

unit displayobjects;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, base3ddisplay,LCLProc;

type

  { TDisplayManager }

  TDisplayManager=class
  private
    FDisplay:IBase3DDisplay;
  public
    constructor Create(Display:IBase3DDisplay);
    procedure SetDisplay(Display:IBase3DDisplay);
    //this is for testing and debugging displays
    procedure Test;
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

procedure TDisplayManager.Test;

var ol:T3DObjectList;

begin
  ol.Material:=ColorMaterial(1,0.5,0.1,1);
  SetLength(ol.Objects,1);
  with ol.Objects[0] do
    begin
    ObjecType:=otSphere;
    Rad:=2;
    sphC:=Coord(0,0,0);
    end;
  FDisplay.AddObjectList(ol);
  FDisplay.Compile;
  DebugLn('Test');
end;

end.

