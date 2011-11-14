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
  minhashcell:Single;

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
  Result:=SRSurface(points,radii,sphere,minhashcell);
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

procedure CalcHull(Atoms:TAtoms;out C1,C2:TCoord;Rad:TFloat);

var f:Integer;

begin
  if Atoms=nil then
    begin
    C1:=NullVector;
    C2:=NullVector;
    end
  else
    begin
    C1:=Atoms[0].Coords;
    C2:=C1;
    for f:=1 to High(Atoms) do
      begin
      C1:=Min(C1,Atoms[f].Coords);
      C2:=Max(C2,Atoms[f].Coords);
      end;
    for f:=0 to 2 do
      begin
      C1[f]:=C1[f]-Rad;
      C2[f]:=C2[f]+Rad;
      end;
    end;
end;

procedure ContactTable(FileName:string);

type
  TCuboid=array[0..1] of TCoord;
  TResidueRec=record
    Residue:TMolecule;
    ContactHull:TCuboid;
    Surface:TFloat;
  end;
  TChain=array of TResidueRec;

var
  chains:array of TChain;

procedure BuildChains;

var
  c,r:Integer;
  mol,chain:TMolecule;

begin
  mol:=pdblayerman.LayerByIx(0).Molecule;
  SetLength(chains,mol.GroupCount);
  for c:=0 to High(chains) do
    begin
    chain:=mol.GetGroup(c);
    SetLength(chains[c],chain.GroupCount);
    for r:=0 to High(chains[c]) do
      begin
      chains[c,r].Residue:=chain.GetGroup(r);
      //TODO: radius added to contact hull should not be fixed
      CalcHull(chains[c,r].Residue.AllAtoms,chains[c,r].ContactHull[0],
                chains[c,r].ContactHull[1],1.6);
      chains[c,r].Surface:=Sum(CalcASAS(chains[c,r].Residue.AllAtoms));
      end;
    end;
end;

function Intersect(Cuboid1,Cuboid2:TCuboid):Boolean;

var c1,c2:TCoord;

begin
  c1:=Max(Cuboid1[0],Cuboid2[0]);
  c2:=Min(Cuboid1[1],Cuboid2[1]);
  Result:=(c1[0]<=c2[0]) and (c1[1]<=c2[1]) and (c1[2]<=c2[2]);
end;

var
  f,c1,c2,r1,r2:Integer;
  mol:TMolecule;
  sr,sr1,sr2:TFloat;
  len1,contactpercentage:Integer;
  atoms1,atoms2,allatoms:TAtoms;
  sl:TStringList;

begin
  sl:=TStringList.Create;
  BuildChains;
  mol:=pdblayerman.LayerByIx(0).Molecule;
  for c1:=0 to High(chains)-1 do
    for c2:=c1+1 to High(chains) do
      for r1:=0 to High(chains[c1]) do
        for r2:=0 to High(chains[c2]) do
          if Intersect(chains[c1,r1].ContactHull,chains[c2,r2].ContactHull) then
          begin
          atoms1:=chains[c1,r1].Residue.AllAtoms;
          atoms2:=chains[c2,r2].Residue.AllAtoms;
          SetLength(allatoms,Length(atoms1)+Length(atoms2));
          for f:=0 to High(atoms1) do allatoms[f]:=atoms1[f];
          len1:=Length(atoms1);
          for f:=0 to High(atoms2) do allatoms[f+len1]:=atoms2[f];
          sr:=Sum(CalcASAs(allatoms));
          sr1:=chains[c1,r1].Surface;
          sr2:=chains[c2,r2].Surface;
          contactpercentage:=100-Round(100*sr/(sr1+sr2));
          if contactpercentage>1 then
            begin
            sl.Add(mol.GetGroup(c1).Name+#9+
                   IntToStr(chains[c1,r1].Residue.ID)+#9+
                   chains[c1,r1].Residue.Name+#9+
                   mol.GetGroup(c2).Name+#9+
                   IntToStr(chains[c2,r2].Residue.ID)+#9+
                   chains[c2,r2].Residue.Name+#9+
                   #9+IntToStr(contactpercentage)+#9+FloatToStr(sr2+sr1-sr));
            WriteLn(sl.Strings[sl.Count-1]);
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
  if HasOption('w','width') then
    begin
    WriteLn('Hashing cell width of ',GetOptionValue('w','width'),'A cells');
    minhashcell:=StrToFloat(GetOptionValue('w','width'));
    end
  else minhashcell:=2;
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
  writeln('Arguments (optional)');
  writeln('-c --contacts: report inter-chain residue contact surfaces instead of SAS');
  writeln('-w --width: set the minimum width of hashing grid cells (e.g. -w 5 for a 5A width).');
end;

var
  Application: TPdbSurface;

{$R *.res}

begin
  Application:=TPdbSurface.Create(nil);
  Application.Run;
  Application.Free;
end.

