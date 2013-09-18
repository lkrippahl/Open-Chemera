{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain. Enjoy.
********************************************************************************
Author: Ludwig Krippahl
Date: 14.7.2011
Purpose:
  Surface calculation, processing alignments, measuring coevolution
Requirements:
  LCL for some units (e.g. FileUtil)
Revisions:
To do:
  Add hidrogens (in pdb molecule)
  fix residue indexing for sequence (starts at 1, should start at first residue)


*******************************************************************************}
program pdbtools;

{$mode objfpc}{$H+}

uses
  {$IFDEF UNIX}{$IFDEF UseCThreads}
  cthreads,
  {$ENDIF}{$ENDIF}
  interfaces, Classes, SysUtils, basetypes, CustApp, surface, pdbmolecules,
  molecules, geomhash, quicksort, pdbparser, mmcifparser, oclconfiguration,
  alignment, stringutils, progress, pdbsurface, contactprediction;

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

  { TPdbTools }
  TPdbTools = class(TCustomApplication)
  protected
    procedure DoRun; override;
    procedure DoContacDescriptors(Rad,MinSurf:TFloat);

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

  MSAlignment:TMSA;       //for msa mapping
  Organisms:TSimpleStrings;
                          //for sorting the lines in the MSA mapping

function UseRes(Name:string):Boolean;

begin
  Result:=((UseOnlyRes=nil) or (LastIndexOf(Name,UseOnlyRes)>=0)) and
           (LastIndexOf(Name,ExcludeRes)<0);
end;

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
      with ResidueRecs[f] do
        begin
        ASA:=0;
        IsolatedASA:=Sum(CalcASAs(Residue.AllAtoms,Rad));
        end;
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

function CalcResidueContact(Res1,Res2:TResidueRec;Rad:TFloat):TFloat;

var
  sr1,sr2,sr:TFloat;
  atoms1,atoms2:TAtoms;
  len,f:Integer;

begin
  atoms1:=Res1.Residue.AllAtoms;
  atoms2:=Res2.Residue.AllAtoms;
  len:=Length(atoms1);
  SetLength(atoms1,Length(atoms1)+Length(atoms2));
  for f:=0 to High(atoms2) do atoms1[f+len]:=atoms1[f];
  sr:=Sum(CalcASAs(atoms1,Rad));
  sr1:=Res1.IsolatedASA;
  sr2:=Res2.IsolatedASA;
  //Contact area is half the lost ASA when in contact
  Result:=(sr1+sr2-sr)/2;
end;


procedure CalcContacts(var Group1,Group2:TGroup;Rad:TFloat);

var
  r1,r2:Integer;
  sr:TFloat;

begin
  SetupResiduesContact(Group1,Rad);
  SetupResiduesContact(Group2,Rad);
  for r1:=0 to High(Group1.ResidueRecs) do
    for r2:=0 to High(Group2.ResidueRecs) do
      if Intersect(Group1.ResidueRecs[r1].ContactHull,
                   Group2.ResidueRecs[r2].ContactHull) then
        begin
        sr:=CalcResidueContact(Group1.ResidueRecs[r1],
            Group2.ResidueRecs[r2],Rad);
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

function SurfaceResidueIxs(Group:TGroup;MinSurf:TFloat):TIntegers;

var f,curr:Integer;

begin
  SetLength(Result,Length(Group.ResidueRecs));
  curr:=0;
  for f:=0 to High(Group.ResidueRecs) do
      if Group.ResidueRecs[f].ASA>=MinSurf then
        begin
        Result[curr]:=f;
        Inc(curr);
        end;
  SetLength(Result,curr)
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
          FloatToStrF(Rec.ASA,ffFixed,1,1)+#9+
          FloatToStrF(Rec.IsolatedASA,ffFixed,1,1));
end;

begin
  for g:=0 to High(Groups) do
    begin
    Sl.Add('');
    Sl.Add('Group: '+Groups[g].Name);
    CalcASA(Groups[g],Rad);
    Sl.Add('Total ASA: '+FloatToStrF(Groups[g].ASA,ffFixed,1,1));
    Sl.Add('Chain'+#9+'Res'+#9+'ID'+#9+'Surf'+#9+'ASA');
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

procedure CalcNeighbours(Group:TGroup;Rad:TFloat);

var
  r1,r2:Integer;
  sr:TFloat;


begin
  SetupResiduesContact(Group,Rad);
  with Group do
    for r1:=0 to High(ResidueRecs)-1 do
      for r2:=r1+1 to High(ResidueRecs) do
         if Intersect(Group.ResidueRecs[r1].ContactHull,
                      Group.ResidueRecs[r2].ContactHull) then
          begin
          sr:=CalcResidueContact(Group.ResidueRecs[r1],
              Group.ResidueRecs[r2],Rad);
          AddToArray(sr,Group.ResidueRecs[r1].ContactSurfaces);
          AddToArray(r2,Group.ResidueRecs[r1].ContactResidues);
          AddToArray(sr,Group.ResidueRecs[r2].ContactSurfaces);
          AddToArray(r1,Group.ResidueRecs[r2].ContactResidues);
          end;
end;

procedure ReportNeighbours(Groups:TGroups;Rad:TFloat;Sl:TStringList);

var g,f,r1:Integer;
    s:string;

begin
  for g:=0 to High(Groups) do
    begin
    CalcNeighbours(Groups[g],Rad);
    Sl.Add('Neighbours: '+Groups[g].Name);

    for r1:=0 to High(Groups[g].ResidueRecs) do
      with Groups[g].ResidueRecs[r1] do
        begin
        s:=Residue.Name+#9+IntToStr(Residue.Id);
        for f:=0 to High(ContactResidues) do
          if ContactSurfaces[f]>MinSurf then
            s:=s+#9+IntToStr(Groups[g].ResidueRecs[ContactResidues[f]].Residue.Id)
                +#9+FloatToStrF(ContactSurfaces[f],ffFixed,1,1);
        Sl.Add(s);
        end;
    end;
end;

procedure ReportContacts(var Groups:TGroups;Rad:TFloat;Sl:TStringList);

function CalcJointASA(Ix1,Ix2:Integer):TFloat;

var
  f,len:Integer;
  tmpatoms:TAtoms;
begin
  len:=Length(Groups[Ix1].AtomRecs);
  SetLength(tmpatoms,len+Length(Groups[Ix2].AtomRecs));
  for f:=0 to High(Groups[Ix1].AtomRecs) do
    tmpatoms[f]:=Groups[Ix1].AtomRecs[f].Atom;
  for f:=0 to High(Groups[Ix2].AtomRecs) do
    tmpatoms[f+len]:=Groups[Ix2].AtomRecs[f].Atom;
  Result:=Sum(CalcASAs(tmpatoms,Rad));
end;

procedure WriteRecs(const Rec1,Rec2:TResidueRec;Cont:TFloat);

begin
  if Cont>=MinSurf then
    Sl.Add(Rec1.Residue.Parent.Name+#9+
          Rec1.Residue.Name+#9+
          IntToStr(Rec1.Residue.Id)+#9+
          Rec2.Residue.Parent.Name+#9+
          Rec2.Residue.Name+#9+
          IntToStr(Rec2.Residue.Id)+#9+
          FloatToStrF(Cont,ffFixed,1,1));
end;

var g1,g2,r,f:Integer;
    index,resix1,resix2:TIntegers;
    contacts:TFloats;
    jointASA:TFloat;

begin
  //Need group ASAs for interface calculation
  for g1:=0 to High(Groups) do
    if Groups[g1].ASA<0 then //check if already calculated (if -s option)
      CalcASA(Groups[g1],Rad);

  for g1:=0 to High(Groups)-1 do
    for g2:=g1+1 to High(Groups) do
      begin
      Sl.Add('Contacts: '+Groups[g1].Name+' to '+Groups[g2].Name);
      CalcContacts(Groups[g1],Groups[g2],Rad);

      jointASA:=CalcJointAsa(g1,g2);
      jointASA:=(Groups[g1].ASA+Groups[g2].ASA-jointASA)/2;
        //contact surface is half the loss in ASA when joining the groups
      Sl.Add('Contact surface (A^2): '+FloatToStrF(jointASA,ffFixed,1,1));
      resix1:=nil;
      resix2:=nil;
      contacts:=nil;
      with Groups[g1] do
        for r:=0 to High(ResidueRecs) do
          with ResidueRecs[r] do
            for f:=0 to High(ContactResidues) do
              begin
              AddToArray(r,resix1);
              AddToArray(ContactResidues[f],resix2);
              AddToArray(-ContactSurfaces[f],contacts); //negative, for sorting
              end;
      if SortResidues then
        begin
        index:=QSAscendingIndex(contacts);
        for r:=0 to High(index) do
          WriteRecs(Groups[g1].ResidueRecs[resix1[index[r]]],
                    Groups[g2].ResidueRecs[resix2[index[r]]],
                    -contacts[index[r]]);
        end
      else
        for r:=0 to High(resix1) do
          WriteRecs(Groups[g1].ResidueRecs[resix1[r]],Groups[g2].ResidueRecs[resix2[r]],
                    -contacts[r]);
      end;
end;

function ChainSequence(Chain:TMolecule):string;

var
  res:TMolecules;
  r:Integer;
  c,s:string;
  rid:Integer;

begin
  res:=Chain.Groups;
  s:='';
  for r:=0 to High(res) do
    begin
    rid:=res[r].Id;
    while Length(s)<rid-1 do s:=s+'X';
    if Length(s)<rid then
      begin
      c:=AAOneLetterCode(res[r].Name);
      if c<>'' then
        s:=s+c
      else s:=s+'X';
      end;
    end;
  Result:=s;
end;

procedure ExtractSequence(Name:string;Groups:TGroups;Report:TStringList);

function GetChains:TMolecules;

var
  f,g:Integer;

begin
  Result:=nil;
  for f:=0 to High(Groups) do
    with Groups[f] do for g:=0 to High(Chains) do
      begin
      SetLength(Result,Length(Result)+1);
      Result[High(Result)]:=Chains[g];
      end;
end;

var
  chains:TMolecules;
  f:Integer;
  s:string;

begin
  chains:=GetChains;
  for f:=0 to High(chains) do
    begin
    s:=ChainSequence(chains[f]);
    Report.Add('> '+Name+':'+chains[f].Name);
    Report.Add(s);
    end;
end;

function MapOrgsToMSA(MSA:TMSA;Orgs:TSimpleStrings):TIntegers;

var f,g:Integer;

begin
  if Orgs=nil then
    Result:=nil
  else
    begin
    SetLength(Result,Length(Orgs));
    for f:=0 to High(Orgs) do
      begin
      Result[f]:=-1;
      for g:=0 to High(MSA.SequenceIds) do
        if Pos(Orgs[f],MSA.SequenceIds[g])>0 then
          begin
          Result[f]:=g;
          Break;
          end;
      end;
    end;
end;

procedure MatchMSA(Chain:TMolecule;MSA:TMSA;Orgs:TSimpleStrings;Report:TStringList);

var
  ix,ixres:Integer;
  map,osmap:TIntegers;
  seq:string;
  f,g:Integer;
  res:TMolecule;

begin
  if MSA.Alignment=nil then
    begin
    Report.Add('No alignment to process');
    Exit;
    end;
  seq:=ChainSequence(Chain);
  osmap:=MapOrgsToMSA(MSA,Orgs);
  if LastIndexOf(-1,osmap)>=0 then
    begin
    Report.Add('Missing organisms');
    for f:=0 to High(osmap) do
      if osmap[f]<0 then Report.Add(Organisms[f]);
    end;

  ix:=FindMatch(seq,MSA,map);

  if ix>=0 then
    begin
    Report.Add('MSA Match: '+Chain.Name+': matches sequence '+IntToStr(ix+1)+
      ' (in first position)');
    for f:=0 to High(Chain.Groups) do
    if Chain.Groups[f].Id>0 then //TODO:improve this, can be negative, but must change sequence extraction...
      begin
      res:=Chain.Groups[f];
      ixres:=map[res.ID-1];
      seq:=msa.Alignment[ix][ixres];
      if osmap=nil then
        begin
        for g:=0 to High(msa.Alignment) do
          if g<>ix then seq:=seq+msa.Alignment[g][ixres];
        end
      else
        for g:=0 to High(osmap) do
          if osmap[g]>=0 then
            seq:=seq+msa.Alignment[osmap[g]][ixres];
      if UseRes(res.Name) then Report.Add(res.name+#9+IntToStr(res.ID)+#9+seq);
      end;
    end
  else
    Report.Add('MSA Match: '+Chain.Name+' does not match any sequence');
end;



{ TPdbTools }

procedure TPdbTools.DoRun;
var
  pdblayerman:TPDBModelMan;
  pdblayer:TPDBModel;
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
  Organisms:=nil;
  MSAlignment.Alignment:=nil;
  MSAlignment.SequenceIds:=nil;
  if HasOption('width') then
    minhashcell:=StringToFloat(GetOptionValue('width'));

  if HasOption('minsurf') then
    minsurf:=StringToFloat(GetOptionValue('minsurf'));

  if HasOption('probe') then
    minsurf:=StringToFloat(GetOptionValue('probe'));

  if HasOption('order') then
    SortResidues:=True;

  if HasOption('exclude') then
    ExcludeRes:=SplitString(GetOptionValue('exclude'),',');

  if HasOption('useonly') then
    begin
      if Uppercase(GetOptionValue('useonly'))='AMINOACIDS' then
        begin
        UseOnlyRes:=nil;
        for f:=0 to High(AAData) do
          AddToArray(AAData[f].TLCode,UseOnlyRes);
        end
      else UseOnlyRes:=SplitString(GetOptionValue('useonly'),',');
    end;

  if HasOption('msa') then
    MSAlignment:=ReadMSA(GetOptionValue('msa'));

  if HasOption('organisms') then
    Organisms:=ReadAsSimpleStrings(GetOptionValue('organisms'));

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

var
  c:string;
  chain:TMolecule;
  command:string;

  procedure DoPDBStats;

  var f:Integer;

  begin
    filemask:=ParamStr(1);
    report.Clear;
    if FindFirst(filemask,faAnyFile,srec)=0 then
      repeat
        pdblayerman.LoadLayer(ExtractFilePath(filemask)+srec.Name);
        pdblayer:=pdblayerman.LayerByIx(0);
        report.Add('Processing file:'+srec.Name);

        BuildGroups;

        if HasOption('s','surface') then
          ReportSurfaces(groups,proberadius,report);

        if HasOption('c','contacts') then
          ReportContacts(groups,proberadius,report);

        if HasOption('q','sequence') then
          ExtractSequence(pdblayer.Molecule.Name,groups,report);

        if HasOption('m','matchmsa') then
          begin
          c:=GetOptionValue('m','matchmsa');
          chain:=pdblayer.GetChain(c);
          if chain=nil then
            report.Add('Match MSA: chain not found')
          else  MatchMsa(chain,MSAlignment,Organisms, report);
          end;
        if HasOption('n','neighbours') then
          ReportNeighbours(groups,proberadius,report);

        pdblayerman.ClearLayers;
        pdblayer:=nil;
      until FindNext(srec)<>0;
    FindClose(srec);
    pdblayerman.Free;
    if HasOption('file') then
      report.SaveToFile(GetOptionValue('file'))
    else
      for f:=0 to report.Count-1 do
        WriteLn(report.Strings[f]);
  end;

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
  pdblayerman:=TPDBModelMan.Create(Config.MonomersPath);
  SetParameters;
  command:=ParamStr(1);
  if command='pdb' then
    DoPDBStats
  else if command='desc' then
    DoContacDescriptors(proberadius,minsurf);
  // stop program loop
  report.Free;
  Terminate;

end;

procedure TPdbTools.DoContacDescriptors(Rad,MinSurf: TFloat);

var
  descriptorgen:TDescriptorGen;
  paramfile:string;
  cdparams:TCDParameters;

procedure FlushSurfaces(Index:Integer);

var
  surfrep:TSurfaceReport;
  f:Integer;

begin
  surfrep:=descriptorgen.GetSurfaceData(Index);
  with surfrep do
  for f:=0 to High(Residues) do
    begin
    WriteLn(f,':',Residues[f].Name,#9,
      Residues[f].ID,#9,
      IsolatedSurface[f],#9,
      IsolatedSidechainSurface[f],#9,
      SASurface[f],#9,
      SASidechainSurface[f],#9);
    end;
end;

procedure WriteReport;

var
  f:Integer;

begin
  with cdparams do
    begin
    WriteLn(TargetPDB);
    WriteLn(ProbePDB);
    WriteLn(FlattenStrings(TargetIDs,', '));
    WriteLn(FlattenStrings(ProbeIDs,', '));
    for f:=0 to High(TargetIndexes) do
      Write(TargetIndexes[f],';');
    WriteLn;
    for f:=0 to High(ProbeIndexes) do
      Write(ProbeIndexes[f],';');
    end;
  FlushSurfaces(0);
  FlushSurfaces(1);
end;

var
  f,g,c:Integer;


begin
  descriptorgen:=TDescriptorGen.Create;
  paramfile:=ParamStr(2);
  cdparams:=ReadParameters(paramfile);

  descriptorgen.GenerateAllDescriptors(cdparams);

  WriteReport;
  descriptorgen.Free;

end;

constructor TPdbTools.Create(TheOwner: TComponent);
begin
  inherited Create(TheOwner);
  StopOnException:=True;
end;

destructor TPdbTools.Destroy;
begin
  inherited Destroy;
end;

procedure TPdbTools.WriteHelp;
begin
  writeln('For pdb analysis:');
  writeln('pdbtools pdb filemask [options]');
  writeln;
  writeln('filemask is the name of pdb file or files (e.g. 1dxg.pdb, *.pdb, ...).');
  writeln('filemask Must be the first argument after the keyword pdb');
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
  writeln('-q --sequence: extract sequence from all selected chains');
  writeln('              containing aminoacids (X for missing residues).');
  writeln;
  writeln('-m --mapmsa: maps each residue of the identified chain to the');
  writeln('             corresponding MSA column. Requires the chain ID and');
  writeln('             a loaded MSA file. E.g. --msa=msa.fasta -m A');
  writeln;
  writeln('-n --neighbours: for each residue in each group, report all residues');
  writeln('                 in the same group in contact, and contact area');
  writeln;
  writeln('--file filename: save report to file (default: write to console)');
  writeln;
  writeln('--width value: set the minimum width of hashing grid cells');
  writeln('                 (e.g. -width=5 for a 5A width). Default 2A.');
  writeln;
  writeln('--minsurf value: Minimum ASA cutoff for a residue to ');
  writeln('                    count as surface or contact. Default 0A^2.');
  writeln;
  writeln('--order: Order residues by ASA (most to least)');
  writeln;
  writeln('--probe: value of probe molecule radius. Default 1.4A.');
  writeln;
  writeln('--useonly: comma separated list of residue types to use');
  writeln('                   (e.g. ALA,ARG).');
  writeln('                   useonly supersedes --exclude.');
  writeln('                   If followed by "aminoacids" (without quotes) use');
  writeln('                   all 20 aminoacids.');
  writeln;
  writeln('--exclude: comma separated list of residue types to ignore.');
  writeln;
  writeln('--msa filename: loads an msa file (fasta format)');
  writeln;
  writeln('--organisms filename: loads list of organisms to sort the msa mapping');
  writeln;
  writeln('For computing contact descriptors:');
  writeln('pdbtools desc parameters');
  writeln('Parmeters file contains (TODO)');
  writeln('Surface exposure is computed with the whole structure, not just the specified chain');
  writeln('The chain specified is the first chain with that chainID; all other repeats are ignored');
  writeln('The first sequence in the MSA file must correspond to the sequence of the specified chain');
  writeln('Can also use optional parameters, but must come *after* the main arguments');



end;

var
  Application: TPdbTools;

{$R *.res}

begin
  Application:=TPdbTools.Create(nil);
  Application.Run;
  Application.Free;
end.

