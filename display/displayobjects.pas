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
          ObjecType:=otSphere;
          Rad:=2;
          sphC:=Coord(x*dist,y*dist,z*dist);
          end;
        end;
  FDisplay.AddObjectList(ol);
  FDisplay.Compile;
  DebugLn('Test');
end;

end.

