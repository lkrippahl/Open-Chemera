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
  Use correct atomic radii. Currently using 1.6 averaga united atom radius
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
  sphere:TCoords;
  rescontact:Boolean;

function CalcASAs(Atoms:TAtoms):TFloats;

var
  f:Integer;
  points:TCoords;
  radii:TFloats;

begin
  SetLength(points,Length(Atoms));
  SetLength(radii,Length(Atoms));
  for f:=0 to High(Atoms) do
    begin
    points[f]:=Atoms[f].Coords;
    radii[f]:=3; //TODO: use correct atomic radii + 1.4
    end;
  Result:=SRSurface(points,radii,sphere);
end;

procedure SaveResidueASAs(Atoms:TAtoms; FileName:string);
//NOTE: assumes atoms are sorted by residue

var
  res,lastres:TMolecule;
  sl:TStringList;
  s:string;
  asa:TFloat;
  asas:TFLoats;
  f:Integer;

begin
  lastres:=nil;
  sl:=TStringList.Create;
  s:='';
  asas:=CalcASAs(Atoms);
  for f:=0 to High(Atoms) do
    begin
    res:=TMolecule(Atoms[f].Parent);
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

procedure ContactTable(FileName:string);

//TODO: inefficient to repeat surface calculation for each individual residue
// ?? but may not amount to much compared to the quadratic combinations...

var
  atoms1,atoms2:TAtoms;
  mol:TMolecule;
  f,len2,chainc,resc:Integer;
  c1,c2,r1,r2,contactpercentage:Integer;
  chain1,chain2:TMolecule;
  sr1,sr2,srboth:Single;
  sl:TStringList;

begin
  sl:=TStringList.Create;
  mol:=pdblayerman.LayerByIx(0).Molecule;
  chainc:=mol.GroupCount-1;
  for c1:=0 to chainc-1 do
    begin
    chain1:=mol.GetGroup(c1);
    for c2:=c1+1 to chainc do
      begin
      chain2:=mol.GetGroup(c2);
      for r1:=0 to chain1.GroupCount-1 do
        if chain2.GetGroup(r1).Name<>'HOH' then
        begin
        atoms1:=chain1.GetGroup(r1).AllAtoms;
        sr1:=Sum(CalcAsas(atoms1));
        for r2:=0 to chain2.GroupCount-1 do
          if chain2.GetGroup(r2).Name<>'HOH' then
          begin
          atoms2:=chain2.GetGroup(r2).AllAtoms;
          sr2:=Sum(CalcAsas(atoms2));

          //calculate both residues together
          len2:=Length(atoms2);
          SetLength(atoms2,len2+Length(atoms1));
          for f:=0 to High(atoms1) do
            atoms2[f+len2]:=atoms1[f];
          srboth:=Sum(CalcAsas(atoms2));
          contactpercentage:=100-Round(100*srboth/(sr2+sr1));
          if contactpercentage>1 then
            begin
            sl.Add(chain1.Name+#9+IntToStr(chain1.GetGroup(r1).ID)+#9+chain1.GetGroup(r1).Name+
                #9+chain2.Name+#9+IntToStr(chain2.GetGroup(r2).ID)+#9+chain2.GetGroup(r2).Name+
                #9+IntToStr(contactpercentage)+#9+FloatToStr(sr2+sr1-srboth));
            WriteLn(sl.Strings[sl.Count-1]);
            end;

          end;
        end;
      end;
    end;
  sl.SaveToFile(FileName);
  sl.Free;
end;

procedure SurfaceAreas;

var
  atoms:TAtoms;

begin
  atoms:=pdblayerman.LayerByIx(0).Molecule.AllAtoms;
  CalcASAs(atoms);
  SaveResidueASAs(atoms,srec.Name+'.txt');
end;


begin
  // quick check parameters
  if HasOption('h','help') or (ParamCount<0) then begin //TODO: check parameters
    WriteHelp;
    Terminate;
    Exit;
  end;

  if HasOption('c','contacts') then rescontact:=True
    else rescontact:=False;

  pdblayerman:=TPdbLayerMan.Create(ParamStr(1));
  sphere:=GoldenSpiralPoints(100);
  if FindFirst('*.pdb',faAnyFile,srec)=0 then
    repeat
    pdblayerman.LoadLayer(srec.Name);
    if rescontact then ContactTable(ChangeFileExt(srec.Name,'.contacts'))
    else SurfaceAreas;
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
  Application.Run;
  Application.Free;
end.

