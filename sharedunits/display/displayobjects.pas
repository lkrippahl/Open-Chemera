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
  Display settings for atoms
  Display options for grids (colors, layers, which grids)
*******************************************************************************}

unit displayobjects;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, base3ddisplay,LCLProc, molecules, oclconfiguration,
  linegrids, geomutils;

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

  TGridSettings=record
    Grid:TGridPlane;
    TransVec:TCoord;
    Resolution:TFloat;
    Ambient,Diffuse,Specular:TRGBAColor;
  end;

  TGridSettingsArray=array of TGridSettings;

  { TDisplayManager }

  TDisplayManager=class
  private
     //These objects are provided externally.
    //This form neither creates nor destroys them
    FDisplay:IBase3DDisplay;        //Interface to display window
    FLayers:TMolecules;
    FGrids:TGridSettingsArray;

    //SettingsLists
    FAtoms:TAtomSettingsArray;
    FBonds:TBondSettingsArray;
    function GridAsObjects(Grid: TGridPlane;TransVec:TCoord;Res:TFloat;
                  Cubes:Boolean=True):T3DObjects;

  public
    constructor Create(Display:IBase3DDisplay);
    procedure SetDisplay(Display:IBase3DDisplay);

    procedure Attach(Molecule:TMolecule);overload;
    procedure Attach(AGrid:TGridPlane;ATransVec:TCoord;AResolution:TFloat;AColor:TRGBAColor);overload;

    procedure Detach(Molecule:TMolecule);
      //remove one molecule

    function LayerIndex(Molecule:TMolecule):Integer;
    procedure RenderLayers;
    procedure RenderGrids;
    procedure Render(ShowLayers:Boolean=True;ShowGrids:Boolean=True);

    procedure OnDeleteAtoms;
      //called from molecules in order to remove tagged atoms and affected bonds
      //the atoms to be deleted are tagged with molecules.TagToDelete

    procedure Test;    //this is for testing and debugging displays
  end;


implementation

{ TDisplayManager }

function TDisplayManager.GridAsObjects(Grid: TGridPlane; TransVec: TCoord;
  Res: TFloat; Cubes: Boolean): T3DObjects;
var
  x,y,z,f,zz:Integer;
  count,curr:Integer;

  procedure AddCuboid(GridX,GridY,GridZ1,GridZ2:Integer);

  var tl,br:TCoord;

  begin
    br:=Coord((GridX+0.1)*res,(GridY+0.1)*res,(GridZ1+0.1)*res);
    tl:=Coord((GridX+0.9)*res,(GridY+0.9)*res,(GridZ2+0.9)*res);
    br:=Add(br,TransVec);
    tl:=Add(tl,TransVec);
    Inc(curr);
    Result[curr].ObjectType:=otCuboid;
    Result[curr].cubTopLeft:=tl;
    Result[curr].cubBotRight:=br;
  end;

begin
  if Cubes then count:=CountCells(Grid)
  else count:=CountSegments(Grid);
  SetLength(Result,count);
  curr:=-1;
  for x:=0 to High(Grid) do
   for y:=0 to High(Grid[x]) do
     for z:=0 to High(Grid[x,y]) do
       if Cubes then
         for zz:=Grid[x,y,z,0] to Grid[x,y,z,1] do
           AddCuboid(x,y,zz,zz)
       else
          AddCuboid(x,y,Grid[x,y,z,0],Grid[x,y,z,1]);
  if curr<High(Result) then
    SetLength(Result,curr+1);
end;

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
  DebugLn('Attached:'+Molecule.Name);
end;

procedure TDisplayManager.Attach(AGrid:TGridPlane;ATransVec:TCoord;AResolution:TFloat;AColor:TRGBAColor);
begin
  SetLength(FGrids,Length(FGrids)+1);
  with FGrids[High(FGrids)] do
    begin
    Grid:=AGrid;
    TransVec:=ATransVec;
    Diffuse:=AColor;
    Specular:=AColor;
    Ambient:=AColor;
    Resolution:=AResolution;
    end;
  DebugLn('Attached Docking grid');
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

procedure TDisplayManager.RenderLayers;
//TO DO: mostly everything. Still just for testing
//TO DO: should check which have changed
//TO DO: group objects by material

var
  f,g:Integer;
  al:TAtoms;
  ol:T3DObjectList;

begin
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
//        if f=0 then
//          ol.Material.Ambient:=RGBAColor(1,0,0,1);
        ol.Objects[0].ObjectType:=otSphere;
        ol.Objects[0].sphC:=al[g].Coords;
        ol.Objects[0].Rad:=al[g].Radius;
        if f>0 then ol.Objects[0].Rad:=al[g].Radius/4;
        FDisplay.AddObjectList(ol);
        end;
    end;
end;

procedure TDisplayManager.RenderGrids;
//TO DO: mostly everything. Still just for testing
//TO DO: should check which have changed
//TO DO: group objects by material

var
  ol:T3DObjectList;
  gr:Integer;


begin
  ol.Material:=DefaultMaterial;
  for gr:=0 to High(FGrids) do
    begin
    ol.Material.Ambient:=FGrids[gr].Ambient;
    ol.Material.Diffuse:=FGrids[gr].Diffuse;
    ol.Material.Specular:=FGrids[gr].Specular;
    ol.Objects:=GridAsObjects(FGrids[gr].Grid,FGrids[gr].TransVec,FGrids[gr].Resolution,False);
    FDisplay.AddObjectList(ol);
    end;
end;
procedure TDisplayManager.Render(ShowLayers: Boolean; ShowGrids: Boolean);

begin
  FDisplay.ClearObjectLists;
  if ShowLayers then RenderLayers;
  if ShowGrids then RenderGrids;
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

