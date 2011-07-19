{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain. Enjoy.
********************************************************************************
Author: Ludwig Krippahl
Date: 14.7.2011
Purpose:
  Calculates exposed surface for pdb files
Requirements:
  LCL for some units (e.g. FileUtil)
Revisions:
To do:
  Add hidrogens (in pdb molecule)
  Use correct atomic radii
*******************************************************************************}
program pdbsurface;

{$mode objfpc}{$H+}

uses
  {$IFDEF UNIX}{$IFDEF UseCThreads}
  cthreads,
  {$ENDIF}{$ENDIF}
  Classes, SysUtils, basetypes, CustApp, surface, pdbmolecules, molecules, geomhash,
  pdbparser
  { you can add units after this };

type

  { TPdbSurface }
                      TPdbSurface = class(TCustomApplication)
  protected
    procedure DoRun; override;
  public
    constructor Create(TheOwner: TComponent); override;
    destructor Destroy; override;
    procedure WriteHelp; virtual;
  end;

{ TPdbSurface }

procedure TPdbSurface.DoRun;
var
  ErrorMsg: String;
  pdblayerman:TPDBLayerMan;
  pdblayer:TPDBLayer;
  srec:TSearchRec;
  atoms:TAtoms;
  points,sphere:TOCCoords;
  radii,asas:TOCFloats;

procedure CalcASAs;

var f:Integer;

begin
  SetLength(points,Length(atoms));
  SetLength(radii,Length(atoms));
  for f:=0 to High(atoms) do
    begin
    points[f]:=atoms[f].Coords;
    radii[f]:=2.4; //TODO: use correct atomic radii + 1.4
    end;
  asas:=SRSurface(points,radii,sphere);
end;

procedure SaveResidueASAs(FileName:string);
//NOTE: assumes atoms are sorted by residue

var
  res,lastres:TMolecule;
  sl:TStringList;
  s:string;
  asa:TOCFloat;
  f:Integer;

begin
  lastres:=nil;
  sl:=TStringList.Create;
  s:='';
  for f:=0 to High(atoms) do
    begin
    res:=TMolecule(atoms[f].Parent);
    if res<>lastres then
      begin
      if s<>'' then sl.add(s+#9+FloatToStr(asa));
      asa:=0;
      lastres:=res;
      s:=res.Name+' '+IntToStr(Res.Id);
      end;
    asa:=asa+asas[f];
    end;
  sl.SaveToFile(FileName);
  sl.Free;
end;

begin
  // quick check parameters
  if HasOption('h','help') or (ParamCount<1) then begin
    WriteHelp;
    Terminate;
    Exit;
  end;

  { add your program here }

  pdblayerman:=TPdbLayerMan.Create(ParamStr(1));
  sphere:=GoldenSpiralPoints(100);
  if FindFirst('*.pdb',faAnyFile,srec)=0 then
    repeat
    pdblayerman.LoadLayer(srec.Name);
    atoms:=pdblayerman.LayerByIx(0).Molecule.AllAtoms;
    CalcASAs;
    SaveResidueASAs(srec.Name+'.txt');
    pdblayerman.ClearLayers;
    until FindNext(srec)<>0;
  FindClose(srec);
  pdblayerman.Free;
  // stop program loop
  Terminate;

end;

constructor TPdbSurface.Create(TheOwner: TComponent);
begin
  inherited Create(TheOwner);
  StopOnException:=True;
end;

destructor TPdbSurface.Destroy;
begin
  inherited Destroy;
end;

procedure TPdbSurface.WriteHelp;
begin
  { add your help code here }
  writeln('Usage: pdbsurface ligands');
  writeln('Where ligands is the path to the ligands folder containing the pdb');
  writeln('files with the PDBeChem data for residues');
end;

var
  Application: TPdbSurface;

{$R *.res}

begin
  Application:=TPdbSurface.Create(nil);
  Application.Title:='PdbSurface';
  Application.Run;
  Application.Free;
end.

