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
  interfaces, Classes, SysUtils, basetypes, CustApp, surface, pdbmolecules, molecules,
  geomhash, quicksort, pdbparser, mmcifparser, oclconfiguration, stringutils
  { you can add units after this };

type
  TCuboid=array[0..1] of TCoord;
  TAtomRec=record
    Atom:TAtom;
    ASA:TFloat;
    ResidueRec:Integer;
  end;
  TAtomRecs=array of TAtomRec;

  TResidueRec=record
    Residue:TMolecule;
    ContactHull:TCuboid;
    ASA,IsolatedASA:TFloat;
    ContactSurfaces:TFloats;
    ContactResidues:TIntegers;
  end;
  TResidueRecs=array of TResidueRec;

  TGroup=record
    //group for surface characterization and contacts
    ASA:TFLoat;
    Chains:TMolecules;
    AtomRecs:TAtomRecs;
    ResidueRecs:TResidueRecs;
    Name:string;
  end;
  TGroups=array of TGroup;

  { TPdbSurface }
  TPdbSurface = class(TCustomApplication)
  protected
    procedure DoRun; override;
  public
    constructor Create(TheOwner: TComponent); override;
    destructor Destroy; override;
    procedure WriteHelp; virtual;
  end;

var
  GoldenSphere:TCoords;
  MinHashCell:Single;     //cell width for geometric hashing
  SortResidues:Boolean;   //sort residue contacts and surface, largest first
  MinSurf:Single;         //ASA cutoff for reporting residues at surface
  UseOnlyRes,ExcludeRes:TSimpleStrings;
                          //list of residue types to use or to exclude
                          //(exclude is ignored unles useonly is nil

function CalcASAs(Atoms:TAtoms;Rad:TFloat):TFloats;overload;

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
    radii[f]:=Atoms[f].Radius+Rad;
    end;
  Result:=SRSurface(points,radii,GoldenSphere,MinHashCell);

end;

procedure CalcASAs(var AtomRecs:TAtomRecs;Rad:TFloat);overload;

var
  asas:TFLoats;
  atoms:TAtoms;
  f:Integer;

begin
  SetLength(atoms,Length(AtomRecs));
  for f:=0 to High(atoms) do atoms[f]:=AtomRecs[f].Atom;
  asas:=CalcAsas(atoms,Rad);
  for f:=0 to High(AtomRecs) do
    AtomRecs[f].ASA:=asas[f];
end;

function UseResidue(ResName:string):Boolean;

begin
  if UseOnlyRes<>nil then
    Result:=IndexOf(ResName,UseOnlyRes)>=0
  else Result:=IndexOf(ResName,ExcludeRes)<0;
end;

procedure PopulateGroup(var Group:TGroup);//from chains creates atomrecs and residuerecs

var
  tmpatoms:TAtoms;
  tmpresidues:TMolecules;
  f,r,a,alen:Integer;


begin
  with Group do
    begin
    AtomRecs:=nil;
    ResidueRecs:=nil;
    ASA:=-1;
    for f:=0 to High(Chains) do
      begin
      tmpresidues:=Chains[f].Groups;
      for r:=0 to High(tmpresidues) do
        if UseResidue(tmpresidues[r].Name) then
          begin
          SetLength(ResidueRecs,Length(ResidueRecs)+1);
          with ResidueRecs[High(ResidueRecs)] do
            begin
            Residue:=tmpresidues[r];
            ASA:=0;
            ContactSurfaces:=nil;
            ContactResidues:=nil;
            tmpatoms:=Residue.AllAtoms;
            alen:=Length(AtomRecs);
            SetLength(AtomRecs,Length(tmpatoms)+alen);
            for a:=0 to High(tmpatoms) do
              begin
              AtomRecs[a+alen].Atom:=tmpatoms[a];
              AtomRecs[a+alen].ASA:=0;
              AtomRecs[a+alen].ResidueRec:=High(ResidueRecs);
              end;
            end;
          end;
      end;
    end;
end;

procedure CalcResASAs(var Group:TGroup;Rad:TFloat);

var
  f:Integer;

begin
  with Group do
    begin
    CalcASAs(AtomRecs,Rad);
    for f:=0 to High(ResidueRecs) do
      ResidueRecs[f].ASA:=0;
    for f:=0 to High(AtomRecs) do
      ResidueRecs[AtomRecs[f].ResidueRec].ASA:=
        ResidueRecs[AtomRecs[f].ResidueRec].ASA+AtomRecs[f].ASA;
    end;
end;

function CalcHull(Atoms:TAtoms;Rad:TFloat):TCuboid;

var
  f:Integer;
  c1,c2:TCoord;

begin
  if Atoms=nil then
    begin
    c1:=NullVector;
    c2:=NullVector;
    end
  else
    begin
    c1:=Atoms[0].Coords;
    c2:=c1;
    for f:=1 to High(Atoms) do
      begin
      c1:=Min(c1,Atoms[f].Coords);
      c2:=Max(c2,Atoms[f].Coords);
      end;
    for f:=0 to 2 do
      begin
      c1[f]:=c1[f]-Rad;
      c2[f]:=c2[f]+Rad;
      end;
    end;
  Result[0]:=c1;
  Result[1]:=c2;
end;

procedure SetupResiduesContact(var Group:TGroup;Rad:TFloat);

var f:Integer;

begin
  with Group do
    for f:=0 to High(ResidueRecs) do
      with ResidueRecs[f] do
        begin
        ContactHull:=CalcHull(Residue.AllAtoms,Rad);
        IsolatedASA:=Sum(CalcAsas(Residue.AllAtoms,Rad));
        ContactSurfaces:=nil;
        ContactResidues:=nil;
        end;
end;

function Intersect(Cuboid1,Cuboid2:TCuboid):Boolean;

  var c1,c2:TCoord;

  begin
    c1:=Max(Cuboid1[0],Cuboid2[0]);
    c2:=Min(Cuboid1[1],Cuboid2[1]);
    Result:=(c1[0]<=c2[0]) and (c1[1]<=c2[1]) and (c1[2]<=c2[2]);
  end;

procedure CalcContacts(var Group1,Group2:TGroup;Rad:TFloat);

var
  f,r1,r2,len:Integer;
  sr,sr1,sr2:TFloat;
  atoms1,atoms2:TAtoms;

begin
  SetupResiduesContact(Group1,Rad);
  SetupResiduesContact(Group2,Rad);
  for r1:=0 to High(Group1.ResidueRecs) do
    for r2:=0 to High(Group2.ResidueRecs) do
      if Intersect(Group1.ResidueRecs[r1].ContactHull,
                   Group2.ResidueRecs[r2].ContactHull) then
        begin
        atoms1:=Group1.ResidueRecs[r1].Residue.AllAtoms;
        atoms2:=Group2.ResidueRecs[r2].Residue.AllAtoms;
        len:=Length(atoms1);
        SetLength(atoms1,Length(atoms1)+Length(atoms2));
        for f:=0 to High(atoms2) do atoms1[f+len]:=atoms1[f];
        sr:=Sum(CalcASAs(atoms1,Rad));
        sr1:=Group1.ResidueRecs[r1].IsolatedASA;
        sr2:=Group2.ResidueRecs[r2].IsolatedASA;
        //Contact area is half the lost ASA when in contact
        sr:=(sr1+sr2-sr)/2;
        AddToArray(sr,Group1.ResidueRecs[r1].ContactSurfaces);
        AddToArray(r2,Group1.ResidueRecs[r1].ContactResidues);
        AddToArray(sr,Group2.ResidueRecs[r2].ContactSurfaces);
        AddToArray(r1,Group2.ResidueRecs[r2].ContactResidues);
        end;
end;

procedure CalcASA(var Group:TGroup;Rad:TFloat);

var
  f:Integer;

begin
  CalcResASAs(Group,Rad);
  with Group do
    begin
    ASA:=0;
    for f:=0 to High(ResidueRecs) do
      ASA:=ASA+ResidueRecs[f].ASA;
    end;
end;

procedure ReportSurfaces(var Groups:TGroups;Rad:TFloat;Sl:TStringList);

var g,r:Integer;
    index:TIntegers;
    contacts:TFloats;

procedure WriteRec(const Rec:TResidueRec);

begin
  if Rec.ASA>=MinSurf then
    Sl.Add(Rec.Residue.Parent.Name+#9+
          Rec.Residue.Name+#9+
          IntToStr(Rec.Residue.Id)+#9+
          FloatToStr(Rec.ASA));
end;

begin
  for g:=0 to High(Groups) do
    begin
    Sl.Add('');
    Sl.Add('Group: '+Groups[g].Name);
    CalcASA(Groups[g],Rad);
    Sl.Add('Total ASA: '+FloatToStr(Groups[g].ASA));
    if SortResidues then
      begin
      SetLength(contacts,Length(Groups[g].ResidueRecs));
      for r:=0 to High(contacts) do
        contacts[r]:=-Groups[g].ResidueRecs[r].ASA;
      index:=QSAscendingIndex(contacts);
      for r:=0 to High(Groups[g].ResidueRecs) do
        WriteRec(Groups[g].ResidueRecs[index[r]]);
      end
    else
      for r:=0 to High(Groups[g].ResidueRecs) do
        WriteRec(Groups[g].ResidueRecs[r]);
    end;

end;


{ TPdbSurface }

procedure TPdbSurface.DoRun;
var
  pdblayerman:TPDBLayerMan;
  pdblayer:TPDBLayer;
  srec:TSearchRec;

  //parameters
  proberadius:Single;     //water molecule radius
  filemask:string;

  groupids:array of TSimpleStrings;
    //chain identifiers for groups; nil if all one group

  groups:TGroups;

  report:TStringList;

procedure SetParameters;

var
  tmp:TSimpleStrings;
  f:Integer;

begin
  //defaults
  GoldenSphere:=GoldenSpiralPoints(50);//TODO: number of points in argument
  MinHashCell:=2;
  MinSurf:=0;
  SortResidues:=False;
  proberadius:=1.4;
  groupids:=nil;

  if HasOption('w','width') then
    minhashcell:=StringToFloat(GetOptionValue('w','width'));

  if HasOption('m','minsurf') then
    minsurf:=StringToFloat(GetOptionValue('m','minsurf'));

  if HasOption('p','proberadius') then
    minsurf:=StringToFloat(GetOptionValue('p','proberadius'));

  if HasOption('o','ordering') then
    SortResidues:=True;

  if HasOption('x','exclude') then
    ExcludeRes:=SplitString(GetOptionValue('x','exclude'),',');

  if HasOption('u','useonly') then
    begin
      if Uppercase(GetOptionValue('u','useonly'))='AMINOACIDS' then
        begin
        UseOnlyRes:=nil;
        for f:=0 to High(AAData) do
          AddToArray(AAData[f].TLCode,UseOnlyRes);
        end
      else UseOnlyRes:=SplitString(GetOptionValue('u','useonly'),',');
    end;


  if HasOption('g','groups') then
    begin
      tmp:=SplitString(GetOptionValue('g','groups'),';');
      SetLength(groupids,Length(tmp));
      for f:=0 to High(groupids) do
        groupids[f]:=SplitString(tmp[f],',');
    end;
end;


procedure BuildGroups;

var f,g:Integer;

begin
  if groupids=nil then
    begin
    SetLength(groups,1);
    groups[0].Chains:=pdblayer.Molecule.Groups;
    groups[0].Name:=pdblayer.Molecule.Name;
    end
  else
    begin
    SetLength(groups,Length(groupids));
    for f:=0 to High(groupids) do
      begin
      SetLength(groups[f].Chains,Length(groupids[f]));
      groups[f].Name:=pdblayer.Molecule.Name+':'+FlattenStrings(groupids[f],',');
      for g:=0 to High(groupids[f]) do
        groups[f].Chains[g]:=pdblayer.GetChain(groupids[f,g]);
      end;
    end;
  for f:=0 to High(groups) do
    PopulateGroup(groups[f]);
end;

var f:Integer;

begin
  // quick check parameters
  if HasOption('h','help') or (ParamCount<1) then begin //TODO: check parameters
    WriteHelp;
    Terminate;
    Exit;
  end;

  report:=TStringList.Create;

  //configuration
  LoadAtomData; //for VdW radius
  LoadAAData;   //for sequence and listing aminoacids
  pdblayerman:=TPDBLayerMan.Create(Config.MonomersPath);
  SetParameters;
  filemask:=ParamStr(1);

  if FindFirst(filemask,faAnyFile,srec)=0 then
    repeat
    report.Clear;
    pdblayerman.LoadLayer(srec.Name);
    pdblayer:=pdblayerman.LayerByIx(0);
    report.Add('Processing file:'+srec.Name);

    BuildGroups;

    if HasOption('s','surface')
     then ReportSurfaces(groups,proberadius,report);
    pdblayerman.ClearLayers;
    pdblayer:=nil;
    until FindNext(srec)<>0;
  FindClose(srec);
  pdblayerman.Free;
  if HasOption('f','file') then
    report.SaveToFile(GetOptionValue('f','file'))
  else
    for f:=0 to report.Count-1 do
      WriteLn(report.Strings[f]);

  // stop program loop
  report.Free;
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
  writeln('pdbsurface filemask [options]');
  writeln;
  writeln('filemask is the name of pdb file or files (e.g. 1dxg.pdb, *.pdb, ...).');
  writeln('filemask Must be the first argument');
  writeln;
  writeln('Optional arguments:');
  writeln;
  writeln('-g --groups: semi-colon separated list of comma separated');
  writeln('             lists of chain ids. E.g. A,B;C,D');
  writeln;
  writeln('-s --surface: report group ASA total and surface residues');
  writeln;
  writeln('-c --contacts: report inter group residue contacts');
  writeln;
  writeln('-n --neighbours: report neighbour contacts between surface residues');
  writeln;
  writeln('-f --file filename: save report to file (default: write to console)');
  writeln;
  writeln('-w --width value: set the minimum width of hashing grid cells');
  writeln('                 (e.g. -w 5 for a 5A width). Default 2A.');
  writeln;
  writeln('-m --minsurf value: Minimum ASA cutoff for a residue to ');
  writeln('                    count as surface or contact. Default 0A^2.');
  writeln;
  writeln('-o --ordering: Order residues by ASA (most to least)');
  writeln;
  writeln('-p --proberadius value: value of probe molecule radius. Default 1.4A.');
  writeln;
  writeln('-u --useonly list: comma separated list of residue types to use');
  writeln('                   (e.g. ALA,ARG).');
  writeln('                   useonly supersedes --exclude.');
  writeln('                   If followed by "aminoacids" (without quotes) use');
  writeln('                   all 20 aminoacids.');
  writeln;
  writeln('-x --exclude list: comma separated list of residue types to ignore.');


end;

var
  Application: TPdbSurface;

{$R *.res}

begin
  Application:=TPdbSurface.Create(nil);
  Application.Run;
  Application.Free;
end.

