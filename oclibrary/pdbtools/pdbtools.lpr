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
  alignment, stringutils, progress, pdbsurface, contactprediction, ebiblastp,
  sequence, fasta,geomutils, base3ddisplay, molfit, molutils;

type

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
    // procedure PDBAnalysis; TODO
    procedure MergeBLASTP(MaxSeqs:Integer=2000);
    // merges a set of xml files with results from BLASTP(EBI)
    procedure SurfaceCharges(Headers:Boolean);
    // computes the surface charges from a pdb file
    // NOTE Assumes charge is in the occupancy field
    procedure BuildComplex;
    // generates complex from unbound and bound structures
    procedure ExportChainSequences;
    procedure Silluette;
    procedure Process;
  public
    constructor Create(TheOwner: TComponent); override;
    destructor Destroy; override;
    procedure WriteHelp; virtual;
  end;

var
  SortResidues:Boolean;   //sort residue contacts and surface, largest first
  MinSurf:Single;         //ASA cutoff for reporting residues at surface



{ TPdbTools }

procedure TPdbTools.DoRun;
var

  report:TStringList;
  command:string;
  proberadius:TFLoat;
    {

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

//        if HasOption('q','sequence') then
//          ExtractSequence(pdblayer.Molecule.Name,groups,report);

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
                                           }
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
//  pdblayerman:=TPDBModelMan.Create(Config.MonomersPath);
//  SetParameters;
  command:=ParamStr(1);
  proberadius:=1.4;
 { if command='pdb' then
    DoPDBStats
  else }if command='desc' then
    DoContacDescriptors(proberadius,minsurf)
  else if command='seq' then
    ExportChainSequences
  else if command='merge' then
    MergeBLASTP
  else if command='scharge' then
    SurfaceCharges(HasOption('headers'))
  else if command='build' then
    BuildComplex
  else if command='silluette' then
    Silluette
  else if command='process' then
    Process;
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
  writeln;
  WriteLn('Res'+#9+'ID'+#9+'Isolated'+#9+'Isolated SC'+#9+'SASurf'+#9+'SAS SC');
  with surfrep do
  for f:=0 to High(Residues) do
    begin
    WriteLn(Residues[f].Name,#9,
      Residues[f].ID,#9,
      FloatToStrF(IsolatedSurface[f],ffFixed,6,1),#9,
      FloatToStrF(IsolatedSidechainSurface[f],ffFixed,6,1),#9,
      FloatToStrF(SASurface[f],ffFixed,6,1),#9,
      FloatToStrF(SASidechainSurface[f],ffFixed,6,1),#9);
    end;
end;

procedure FlushContact(Res1,Res2:TMolecule;Contact:TFloat);

begin
  WriteLn(Res1.Name,#9,Res1.ID,#9,
          Res2.Name,#9,Res2.ID,#9,
          FloatToStrF(Contact,ffFixed,6,1));
end;

procedure FlushInnerContacts(StructureIndex:Integer);

var
  r1,r2:Integer;
  contacts:TMatrix;
  surfrep:TSurfaceReport;
begin
  WriteLn();
  WriteLn('Full, ',StructureIndex+1);
  contacts:=descriptorgen.FullInnerContacts(StructureIndex);
  surfrep:=descriptorgen.GetSurfaceData(StructureIndex);
  for r1:=0 to High(contacts)-1 do
    for r2:=r1+1 to High(contacts) do
      if contacts[r1,r2]>0 then
        FlushContact(surfrep.Residues[r1],surfrep.Residues[r2],contacts[r1,r2]);
  WriteLn('Sidechain, ',StructureIndex+1);
  contacts:=descriptorgen.SidechainInnerContacts(StructureIndex);
  for r1:=0 to High(contacts)-1 do
    for r2:=r1+1 to High(contacts) do
      if contacts[r1,r2]>0 then
        FlushContact(surfrep.Residues[r1],surfrep.Residues[r2],contacts[r1,r2]);

end;


procedure FlushCrossContacts;

var
  rt,rp:Integer;
  contacts:TMatrix;
  surft,surfp:TSurfaceReport;
begin
  WriteLn();
  WriteLn('Full cross contacts');
  contacts:=descriptorgen.FullCrossContacts;
  surft:=descriptorgen.GetSurfaceData(TargetIndex);
  surfp:=descriptorgen.GetSurfaceData(ProbeIndex);
  for rt:=0 to High(contacts) do
    for rp:=0 to High(contacts[0]) do
      if contacts[rt,rp]>0 then
        FlushContact(surft.Residues[rt],surfp.Residues[rp],contacts[rt,rp]);

  WriteLn('Sidechain contacts');
  contacts:=descriptorgen.SidechainCrossContacts;
  for rt:=0 to High(contacts) do
    for rp:=0 to High(contacts[0]) do
      if contacts[rt,rp]>0 then
        FlushContact(surft.Residues[rt],surfp.Residues[rp],contacts[rt,rp]);

end;


procedure FlushNeighbours(Index:Integer);

var
  f,n:Integer;
  s:string;
  res:TMolecule;
  oix1,oix2:Integer;

begin
  WriteLn;
  WriteLn('Neighbours ',Index);
  with descriptorgen.Descriptors do
    begin
    for f:=0 to High(Neighbours[Index]) do
      begin
      s:='';
      oix1:=SelectedIXs[Index,f];
      for n:=0 to High(Neighbours[Index,f]) do
        begin
        res:=Selected[Index,Neighbours[Index,f,n]];
        oix2:=SelectedIXs[Index,Neighbours[Index,f,n]];
        s:=s+#9+res.Name+#9+IntToStr(res.ID)+#9+
          IntToStr(Round(descriptorgen.SidechainInnerContacts(Index)[oix1,oix2]));
        end;
      WriteLn(Selected[Index,f].Name,#9,Selected[Index,f].ID,#9,s)
      end;
    end;
end;

procedure FlushSequences(DBIndex,Index:Integer);

var
  f:Integer;

begin
  WriteLn;
  WriteLn('Sequences ',Index);
  WriteLn('Database:',descriptorgen.Descriptors.MSADBNames[DBIndex]);
  with descriptorgen.Descriptors do
    for f:=0 to High(ResidueMSAs[Index]) do
      WriteLn(Selected[Index,f].Name,#9,Selected[Index,f].ID,#9,ResidueMSAs[DBIndex,Index,f]);
end;


procedure WriteReport;

var
  f:Integer;

begin
  with cdparams do
    begin
    WriteLn(TargetPDBFile);
    WriteLn(ProbePDBFile);
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
  FlushInnerContacts(0);
  FlushInnerContacts(1);
  FlushCrossContacts;
  FlushNeighbours(0);
  FlushNeighbours(1);
  for f:=0 to High(descriptorgen.Descriptors.MSADBNames) do
    begin
    FlushSequences(f,0);
    FlushSequences(f,1);
    end;
end;

var
  f,g,c:Integer;
  //sl:TStringList;


begin
  descriptorgen:=TDescriptorGen.Create(50);
  paramfile:=ParamStr(2);
  cdparams:=ReadParameters(paramfile);

  if HasOption('s','surface') then
    descriptorgen.DistanceDescriptors(cdparams)
  else
    descriptorgen.GenerateAllDescriptors(cdparams);
  descriptorgen.SaveDescriptors(HasOption('headers'), HasOption('c','contacts'));
  if not(HasOption('s','surface') or HasOption('c','contacts')) then
    descriptorgen.SaveRaw();
end;

procedure TPdbTools.MergeBLASTP(MaxSeqs:Integer=2000);

var
  fils:TSimpleStrings;
  seqsset:TOCSequencesSet;
  tmpseqs:TOCSequences;
  f:Integer;
  fastawriter:TFastaWriter;
  tmpfastareader,fastareader:TFastaReader;

begin
  SetLength(fils,ParamCount-1);
  SetLength(seqsset,ParamCount-2);
  fastareader:=TFastaReader.Create(ParamStr(2));
  for f:=3 to ParamCount do
    begin
    tmpseqs:=ReadEBIBlast(ParamStr(f));     // try ebi blastp xml format
    if tmpseqs=nil then                     // try as fasta file
      begin
      tmpfastareader:=TFastaReader.Create(ParamStr(f));
      tmpseqs:=tmpfastareader.AsOCSequences;
      tmpfastareader.Free;
      end;
    seqsset[f-3]:=tmpseqs;
    seqsset[f-3,0].Sequence:=fastareader.SequenceByID(seqsset[f-3,0].ID);
    end;

  seqsset:=FilterByCommonOrganisms(seqsset);

  fastareader.Free;
  fastawriter:=TFastaWriter.Create;

  for f:=3 to ParamCount do
    begin
    fastawriter.clear;
    if Length(seqsset[f-3])>MaxSeqs then
      SetLength(seqsset[f-3],MaxSeqs);
    fastawriter.AppendSeqs(seqsset[f-3],True);
    fastawriter.SaveToFile(ChangeFileExt(ExtractFileName(ParamStr(f)),'.fas'));
    end;

  fastawriter.Free;
end;

procedure TPdbTools.SurfaceCharges(Headers: Boolean);

var
  surman:TPdbSurface;
  pdblayerman:TPDBModelMan;
  pdb:string;
  atoms:TAtoms;
  asas:TFloats;
  totasa,totcharge,surfcharge:TFloat;
  f:Integer;
  dipole,center:TCoord;

begin
  pdblayerman:=TPDBModelMan.Create(Config.MonomersPath);
  pdb:=ParamStr(2);
  pdblayerman.LoadLayer(pdb,pdbOccTemp);
  surman:=TPDBSurface.Create(50);
  atoms:=pdblayerman.LayerByIx(0).Molecule.AllAtoms;
  asas:=surman.CalcASAs(atoms);
  if Headers then
    WriteLn('PDB'+#9+'Total ASA (A^2)'+#9+'Total Charge (e)'+#9+
            'Surface Charge (e)'+#9+'Charge Density (e/A^2)'+#9+'Dipole (D)');
  totasa:=Sum(asas);
  totcharge:=0;
  surfcharge:=0;
  //compute geometri center
  center:=NullVector;
  for f:=0 to High(atoms) do
      center:=Add(center,atoms[f].Coords);
  center:=Scaled(center,1/Length(atoms));

  dipole:=NullVector;
  for f:=0 to High(atoms) do
    begin
    totcharge:=totcharge+atoms[f].Charge;
    dipole:=Add(dipole,Scaled(Subtract(atoms[f].Coords,center),atoms[f].Charge));
    if asas[f]>0 then
      surfcharge:=surfcharge+atoms[f].Charge;
    end;
  WriteLn(pdb+#9+FloatToStrF(totasa,ffFixed,3,1)+
              #9+FloatToStrF(totcharge,ffFixed,3,1)+
              #9+FloatToStrF(surfcharge,ffFixed,3,1)+
              #9+FloatToStrF(surfcharge/totasa,ffExponent,3,1)+
              #9+FloatToStrF(Norm(dipole)*0.20819434,ffExponent,3,1));
  surman.Free;
  pdblayerman.Free;
end;

procedure TPdbTools.BuildComplex;

var
  pdb1,pdb2,complex:string;
  mols1,mols2,molsc,fitmol1,fitmol2,builtc:TMolecule;
  pdblayerman:TPDBModelMan;
  chainlist:TSimpleStrings;
  submat:TSubMatrix;
  fit1,fit2:TFitResult;
  folder:string;
  exclusions:TIntegers;
  minscore:TFloat;

  function GetLayer(PdbChains:string):TMolecule;

  var
    pdb,chains:string;
    ix,f:Integer;
    layer:TPDBModel;

  begin
    ix:=Pos('|',PdbChains);
    if ix>0 then
      begin
      pdb:=Copy(PdbChains,1,ix-1);
      chains:=Copy(PdbChains,ix+1,Length(PdbChains));
      end
    else
      begin
      pdb:=PdbChains;
      chains:=''
      end;
  pdblayerman.ClearLayers;
  pdblayerman.LoadLayer(pdb);
  layer:=pdblayerman.LayerByIx(0);
  layer.DeleteResidues([resNonAA]);
  if chains<>'' then
    begin
    chainlist:=SplitChars(chains);
    Result:=layer.CopyChains(chainlist);
    end
  else
    Result:=TMolecule.CopyFrom(layer.Molecule,nil);

  end;

  function NewMoleculeFromFit(Original,Template:TMolecule;Fit:TFitResult):TMolecule;

  var
    f:Integer;
    chn:TMolecule;

  begin
    Result:=TMolecule.Create('',0,nil);
    for f:=0 to High(Fit.Map) do
      if Fit.Map[f]>=0 then
        begin
        chn:=Original.GetGroup(f);
        chn.Name:=Template.GetGroup(Fit.Map[f]).Name;
        Result.AddGroup(TMolecule.CopyFrom(chn,Result));
        end;
  end;

  procedure WriteReport;

  var
    report:string;
    f,g,count:Integer;

  begin
    report:=ParamStr(2)+#9+ParamStr(3)+#9+ParamStr(4)+#9;
    if (fit1.rmsd>0) and (fit2.rmsd>0) then
      begin
      report:=report+FlattenStrings(fitmol1.ListGroupNames,'')+#9+
                     FlattenStrings(fitmol2.ListGroupNames,'')+#9+
                     FloatToStrF(fit1.Rmsd,ffFixed,3,3)+#9+
                     FloatToStrF(fit2.Rmsd,ffFixed,3,3);
      report:=report+#9;
      for f:=0 to High(fit1.Map) do
        report:=report+IntToStr(fit1.Map[f])+' ';
      report:=report+#9;
      for f:=0 to High(fit2.Map) do
        report:=report+IntToStr(fit2.Map[f])+' ';
      report:=report+IntToStr(CountResidueMatches(fit1))+' ';
      report:=report+IntToStr(CountResidueMatches(fit2))+' ';
      end
    else
      begin
      if fit1.rmsd<0 then
        report:=report+'Failed on '+ParamStr(2)+#9;
      if fit2.rmsd<0 then
        report:=report+'Failed on '+ParamStr(3)+#9;
      end;
    WriteLn(report);
  end;

var f:Integer;

begin
  if ParamCount<8 then
    begin
    WriteLn('Insufficient parameters '+IntToStr(ParamCount));
    WriteLn('(pdb1 pdb2 complex submat out1 out2 outcomplex [minscore]');
    Exit();
    end;
  minscore:=0.7;
  if ParamCount>=9 then
    minscore:=StrToFloat(ParamStr(9));
  pdb1:=ParamStr(2);
  pdb2:=ParamStr(3);
  complex:=ParamStr(4);
  pdblayerman:=TPDBModelMan.Create(Config.MonomersPath);
  mols1:=GetLayer(pdb1);
  mols2:=GetLayer(pdb2);
  molsc:=GetLayer(complex);
  submat:=ReadBLASTMatrix(ParamStr(5));

  fit1:=MagicFit(mols1,molsc,SubMat,nil,minscore);
  exclusions:=nil;
  for f:=0 to High(fit1.Map) do
      if fit1.Map[f]>=0 then AddToArray(fit1.Map[f],exclusions);
  fit2:=MagicFit(mols2,molsc,SubMat,exclusions,minscore);
  if (fit1.Rmsd>0) and (fit2.Rmsd>0) then
    begin
    fitmol1:=NewMoleculeFromFit(mols1,molsc,fit1);
    fitmol2:=NewMoleculeFromFit(mols2,molsc,fit2);
    SaveToPDB(fitmol1,ParamStr(6));
    SaveToPDB(fitmol2,ParamStr(7));
    fitmol1.Transform(fit1.Center,fit1.Rotation,fit1.Translation);
    fitmol2.Transform(fit2.Center,fit2.Rotation,fit2.Translation);
    builtc:=TMolecule.Create('',0,nil);
    for f:=0 to fitmol1.GroupCount-1 do
      builtc.AddGroup(TMolecule.CopyFrom(fitmol1.GetGroup(f),builtc));
    for f:=0 to fitmol2.GroupCount-1 do
      builtc.AddGroup(TMolecule.CopyFrom(fitmol2.GetGroup(f),builtc));
    builtc.RenumberAtoms;
    SaveToPDB(builtc,ParamStr(8));
    end;
  WriteReport;
  pdblayerman.Free;
  builtc.Free;
  fitmol1.Free;
  fitmol2.Free;
end;

procedure TPdbTools.ExportChainSequences;
var
  chainids:TSimpleStrings;
  pdb,pdbname:string;
  f:Integer;
  pdblayerman:TPDBModelMan;
  chains:TMolecules;
begin
  pdblayerman:=TPDBModelMan.Create(Config.MonomersPath);
  SetLength(chainids,ParamCount-2);
  pdb:=ParamStr(2);

  //Delete extensions (.pdb or .pdb.gz, tipically)
  pdbname:=ExtractFileName(pdb);
  while (Length(pdbname)>1) and (Pos('.',pdbname)>0) do
    delete(pdbname,Pos('.',pdbname),Length(pdbname));

  for f:=3 to ParamCount do
    chainids[f-3]:=ParamStr(f);
  pdblayerman.LoadLayer(pdb);
  chains:=pdblayerman.GetChains(0,chainids);
  for f:=0 to High(chains) do
    begin
    WriteLn('>',pdbname+'_'+chains[f].name);
    WriteLn(ChainSequence(chains[f]));
    end;
  pdblayerman.Free;
end;

procedure TPdbTools.Silluette;
//computes the silluette of the pdb file in a number of directions


  function Area(var Coords:TCoords;const Rads:TFloats;MaxRad:TFloat):Integer;

  var
    f:Integer;
    hasher:TGeomHasher;
    minc,maxc:TCoord;
    x,y:TFloat;

  begin
    for f:=0 to High(Coords) do
      Coords[f,2]:=0;
    hasher:=TGeomHasher.Create(Coords,MaxRad,Rads);
    minc:=Min(Coords);
    maxc:=Max(Coords);
    Result:=0;
    x:=minc[0]+0.5-maxrad;
    while x<=maxc[0]+1+maxrad do
      begin
      y:=minc[1]+0.5-maxrad;
      while y<=maxc[1]+1+maxrad do
        begin
        if hasher.IsInnerPoint(Coord(x,y,0)) then
          Inc(Result);
        y:=y+1;
        end;
      x:=x+1;
      end;
    hasher.Free;
  end;

  function UniformSampleAngles(NumAxes:Integer):TQuaternions;

  //Generates rotation quaternions evenly spread around sphere

var
  sphere:TCoords;
  f,g,ix:Integer;
  zaxis:TCoord;
begin
  if NumAxes<=1 then
    begin
    SetLength(Result,1);
    Result[0]:=IdentityQuaternion;
    end
  else
    begin
    zaxis:=Coord(0,0,1);
    //get ZSteps^2 points uniformly distributed around a sphere
    sphere:=GoldenSpiralPoints(NumAxes);
    SetLength(Result,Length(sphere));
    ix:=0;
    zaxis:=Coord(0,0,1);
    for g:=0 to High(sphere) do
      Result[g]:=RotationTo(zaxis,sphere[g]);
    end;
end;

var
  pdb:string;
  angles:Integer;
  added:TFloat;
  coords,rotcoords:TCoords;
  rotations:TQuaternions;
  pdblayerman:TPDBModelMan;
  rads:TFloats;
  maxrad:TFloat;
  areas:TIntegers;
  f:Integer;

begin
  pdblayerman:=TPDBModelMan.Create(Config.MonomersPath);
  pdb:=ParamStr(2);
  pdblayerman.LoadLayer(pdb);
  angles:=StrToInt(ParamStr(3));
  rotations:=UniformSampleAngles(angles);
  added:=StrToFloat(ParamStr(4));
  coords:=ListCoords(pdblayerman.LayerByIx(0).Molecule);
  rads:=ListRadii(pdblayerman.LayerByIx(0).Molecule);
  for f:=0 to High(rads) do rads[f]:=rads[f]+added;
  maxrad:=Max(rads);
  SetLength(areas,length(rotations));
  for f:=0 to High(rotations) do
    begin
    rotcoords:=Rotate(coords,rotations[f]);
    areas[f]:=Area(rotcoords,rads,maxrad);
    end;
  if HasOption('headers') then
    WriteLn('PDB'+#9+'Angles'+#9+'Average'+#9+'Median'+#9+'Max'+#9+'Min');
  WriteLn(pdb+#9+IntToStr(Round(angles))+#9+
                 IntToStr(Round(Average(areas)))+#9+
                 IntToStr(Round(Median(areas)))+#9+
                 IntToStr(Round(Max(areas)))+#9+
                 IntToStr(Round(Min(areas))));
end;

procedure TPdbTools.Process;
{ For processing a pdb file:
  pdbtools process pdb [--writestats] [--renamechains] [--save pdb]
                       [--aaonly] [--headers]');
  writestats: writes pdb stats (after processing)
  renamechains: renames chains consecutively A...Z
  aaonly: deletes all non amino acid residues and empty chains
  headers: writes stats headers
  deletechains Labels: deletes all chains with those Labels, separated by ; (after renaming)');
}

var
  pdb,outfile:string;
  layer:TMolecule;
  pdblayerman:TPDBModelMan;
  mol:TMolecule;

procedure DeleteChains;

var
  f:Integer;
  labels:TSimpleStrings;

begin
  labels:=SplitString(GetOptionValue('deletechains'),';');
  for f:=0 to High(labels) do
    pdblayerman.LayerByIx(0).DeleteChainsByName(labels[f]);
end;

var f:Integer;

begin
  pdblayerman:=TPDBModelMan.Create(Config.MonomersPath);
  pdb:=ParamStr(2);
  mol:=pdblayerman.LoadLayer(pdb);
  if HasOption('aaonly') then
    pdblayerman.LayerByIx(0).DeleteNonAAResidues;
  if HasOption('renamechains') then
    pdblayerman.LayerByIx(0).RenameChains;
  if HasOption('deletechains') then
    DeleteChains;
  if HasOption('headers') then
    WriteLn('PDB'+#9+'Chain'+#9+'Residues');
  if HasOption('writestats') then
    for f:=0 to mol.GroupCount-1 do
      WriteLn(pdb+#9+mol.GetGroup(f).Name+#9+IntToStr(mol.GetGroup(f).GroupCount));
  if HasOption('save') then
    begin
    outfile:=GetOptionValue('save');
    SaveToPdb(mol,outfile);
    end;
  pdblayerman.Free;
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
  {writeln('For pdb analysis:');
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
  writeln;}

  writeln('For computing contact descriptors:');
  writeln('pdbtools desc parameterFile [--headers -c -s]');
  writeln('Parmeters file contains (TODO)');
  writeln('Surface exposure is computed for the specified chains but considering');
  writeln('all the structure, not just the specified chains');
  writeln('Each chain specified is the first chain with that chainID; all other repeats are ignored');
  writeln('The first sequence in the MSA file must correspond to the sequence of the specified chain');
  writeln('Can also use optional parameters, but must come *after* the main arguments:');
  writeln('--headers Include headers line in output');
  writeln('-c --contacts output only actual contacts');
  writeln('-s --surface compute only surface areas');
  writeln();
  writeln('For exporting sequences:');
  writeln('pdbtools seq pdbFile chainIds');
  writeln('E.g. pdbtools seq 1dxg.pdb.gz A B');
  writeln('Writes .fas files in current folder');
  writeln();
  writeln('For merging BLAST results, filtering by organism:');
  writeln('pdbtools merge query res1 res2 ...');
  writeln('E.g. pdbtools merge query res1.xml res2.xml');
  writeln('query is a FASTA file containing the query sequences identified by their ID');
  writeln('Writes .fas files in current folder containing all sequences from organisms');
  writeln('present in all BLAST files, with one sequence for each organism.');
  writeln('Query result files can be either in EBI xml format or Fasta format with organisms');
  writeln('indicated on the ID line with OS="species" or [organism=species]');
  writeln('However, if in fasta format then the first sequence must have the query ID (although the');
  writeln('sequence itself may be empty)');
  writeln();
  writeln('For computing surface charges:');
  writeln('pdbtools scharge [--headers] pdbq_file');
  writeln('NOTE: Assumes that charge is in the occupancy field');
  writeln();
  writeln('For creating a complex from unbound and bound structures:');
  writeln('pdbtools build pdb1[|chains] pdb2[|chains] complex[|chains] submatrix out1 out2 outcomplex [minscore]');
  writeln('Output is out1 out2 outcomplex files with the selected chains and built complex');
  writeln();
  writeln('For computing the silluette of a pdb file:');
  writeln('pdbtools silluette pdb num_angles added_radius [-headers]');
  writeln();
  writeln('For processing a pdb file:');
  writeln('pdbtools process pdb [-writestats -renamechains -save pdb]');
  writeln('                     [-aaonly -headers -deletechain Label]');
  writeln('writestats: writes pdb stats (after processing)');
  writeln('renamechains: renames chains consecutively A...Z');
  writeln('aaonly: deletes all non amino acid residues and empty chains');
  writeln('deletechains Labels: deletes all chains with those Labels, separated by ; (after renaming)');



end;

var
  Application: TPdbTools;

{$R *.res}

begin
  Application:=TPdbTools.Create(nil);
  Application.Run;
  Application.Free;
end.

