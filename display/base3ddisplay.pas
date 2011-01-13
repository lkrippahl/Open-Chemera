{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain. Enjoy.
********************************************************************************
Author: Ludwig Krippahl
Date: 9.1.2011
Purpose:
  This is the abstract base interface for the Chemera 3D display.
  All such forms, windows, etc must implement this interface.
  Also includes utility functions (eg for colors, coords, etc)
Requirements:
  Must provide the basic functions for display
    Cylinder
    Sphere
    Colors and stuff
Revisions:
To do:
  Material specs (only color at present)
  Add textures
  Add lights
  Add zoom, clip, etc
  Add motion controls
*******************************************************************************}

unit base3ddisplay;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, gl;


type
  TRGBAColor=array[0..3] of GLFLoat;

  //TODO: decide on callback parameters (Clicked, etc)
  TBase3DOnMouse=procedure(ObjectId:Integer; Action:Integer);

  //Material specifications. Should have the most features implemented
  T3DMaterial=record
    R,G,B,A:Single;
  end;

  //Display interface
  //TODO: fill in functions in paralel with GLUT window development
  IBase3DDisplay=interface

    //These functions return the object IDs
    //Quality is for the rendering of each object (e.g. number of polygons)
    //-1 in quality means to use the default, as set by each renderer
    function AddSphere(Coords:TOCCoords;Rad:TOCFloat;MatId:Integer;Quality:Integer=-1):Integer;
    function AddCylinder(C1,C2:TOCCoords;R1,R2:TOCFloat;MatId:Integer;Quality:Integer=-1):Integer;
    function AddLine(C1,C2:TOCCoords;R1,R2:TOCFloat;MatId:Integer):Integer;

    //Material IDs are independent of the object ids
    function AddMaterial(Material:T3DMaterial):Integer;

    procedure ClearObjects;
    procedure ClearMaterials;
    procedure Refresh;
    procedure Resize(Width,Height:Integer);
  end;

function RGBAColor(R,G,B,A:GLFloat):TRGBAColor;

implementation

function RGBAColor(R,G,B,A:GLFloat):TRGBAColor;

begin
  Result[0]:=R;
  Result[1]:=G;
  Result[2]:=B;
  Result[3]:=A;
end;


end.

