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
  sequence, fasta;

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
    procedure ExportChainSequences;

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
    MergeBLASTP;

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
  descriptorgen.SaveDescriptors(HasOption('h','headers'), HasOption('c','contacts'));
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


end;

var
  Application: TPdbTools;

{$R *.res}

begin
  Application:=TPdbTools.Create(nil);
  Application.Run;
  Application.Free;
end.

