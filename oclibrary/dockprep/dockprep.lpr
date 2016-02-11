{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 2014 02 21
Purpose:
  prepares docking orders and structure files, reads dock files


  for docking:
    Path options

    -pdbpath path
      path to pdb files, passed to jobfile loader

    -outputpath
      path for output files (for all reports, export, etc) TODO: implement on all

    Structure file options:

    -t target.pdb file to load (mandatory)
    -p probe.pdb file to load  (mandatory)
    -tout target.pdb file to save (mandatory if processing options)
    -pout probe.pdb fle to save (mandatory if processing options)
    -tchains target chains to select (e.g. ACD)
    -pchains probe chains to select
    -center
      center target and probe molecules
    -randomize
      randomize probe orientation
    -cut SurfaceCutoff
      cut cores on exposed residues with at least SurfaceCutoff exposed surface

    Job file options

    -j jobfile.xml (mandatory)


    -append
      append to jobs file instead of replacing (default)

    -contacts contactConstraintsFile
      file containing a list of residue contacts:
        chainID residueID chainID residueID distance num_models minoverlap
        the first residue is assumed to be target, the second is probe

    -unconstrained nummodels[:minoverlap]
      adds an unconstrained docking set. Nummodels is mandatory, minoverlap
      (after colon) is optional

    -maxcontacts num
      number of contacts to read from contacts file. If absent or 0 reads all

    -baseangles num
      base number of angles for selection. Default 2000

    -shapepoints num
      number of points for defining molecule. Default 20

    -displacement num
      maximum displacement between axes, default 10

    -maxaxes num
      maximum number of axes after selection (overrides displacement)

    -symmetry num
      integer for using symmetry constraints on angles

    -symmetryplane num
      float for symmetry plane constraint, with plane width. Plane is perpendicular to axis

  Modifying job files

    -modify
      this flag is necessary to load and modify a job file
      { TODO : add an identifyer for the job to edit? }

  Scores options (can be added with -modify)


    -rmsd complex.pdb
      computes alfa-carbon rmsd to the complex. Currently only works if target, probe and
      complex are the same pdb file
      The superposition is made on the target molecule, the RMSD is computed for the probe

    -intrmsd complex.pdb
      computes the interface residues RMSD using the alfa carbon of all residues within
      5A (the interface) or the distance specified in intdist

    -intdist distance_in_A

    (both require a substitution matrix and a minimum alignment match)

    -submat filename

    -minalign (value from 0 to 1)


    -delscores
      deletes all additional scores
      TODO: enable specifying score id and deleting only those

    -contactscore contactmatrix
      loads contactmatrix, computes contact score
      outputs to console! TODO: this should be added to xml to compute in bigger... (experimental only)

  Reading results

    list models in each constraintset
    -table ordersfile
    displays on screen, each constraint in one group of columns, then model id, overlap, other scores


    list all models
    -models ordersfile
    displays on screen table with single list of models and all scores (removes duplicates)

    Export pdbs
    -export filename
    -all
       exports all pdb files as ordersfile_run_constreintset_id.pdb
       this option ignores others but uses
       -firstid
       -lastid
       for specifying range

    -outfile filename.pdb
    -constraintset number
    -model  id

    TODO: score filters not implemented
    -minscore scoreix
    -maxscore scoreix
      (if scoreix<0 then it refers to geometry)

    summary:
    -summary [h]
        outputs summary to screen, with headers if h
        -maxmodels NN;NN2;...
          maximum number of models to use per constraint
        -rmsdixs Nprobe;Nint
          indexes of the rmsd scores to be used for probe and interface
        -list
          lists ids of acceptable models (0 based)

    ***TESTING***
    ***This only works in development***
    -test


Requirements:
Revisions:
To do:
*******************************************************************************}

program dockprep;

{$mode objfpc}{$H+}

uses
  {$IFDEF UNIX}{$IFDEF UseCThreads}
  cthreads,
  {$ENDIF}{$ENDIF}
  Classes, SysUtils, docktasks, molecules, pdbmolecules, molfit,
  oclconfiguration, alignment, geomutils, stringutils, quicksort, geomhash,
  basetypes, base3ddisplay, progress, CustApp, rotations,molutils,protinter
  { you can add units after this };

type

  { TDockPrep }

  TDockPrep = class(TCustomApplication)
  protected
    FDockMan:TDockOrdersManager;
    procedure DoRun; override;
  public
    constructor Create(TheOwner: TComponent); override;
    destructor Destroy; override;
    procedure WriteHelp; virtual;
  end;

{ TDockPrep }

procedure RunTests;

var
  modman:TPDBModelMan;
  pdb:TPDBModel;
  axes,coords:TCoords;
  distmatrix:TMatrix;
  sl:TStringList;
  x,y:Integer;
  s:string;
  matrix:TMatrix;
begin
  //distance matrix generation test
  {modman:=TPDBModelMan.Create(Config.MonomersPath);
  WriteLn('CreatedOK');
  modman.LoadLayer('..\testfiles\1RPO-A-prob.pdb');
  pdb:=modman.LayerByIx(0);
  CenterMolecule(pdb.Molecule);
  coords:=ListCoords(pdb.Molecule);
  axes:=BaseSampleAxes(1000);
  matrix:=AxialDisplacementMatrix(axes,coords);
  sl:=TStringList.Create;
  for x:=0 to High(matrix) do
    begin
    s:='';
    for y:=0 to High(matrix[x]) do
      begin
      s:=s+FloatToStrF(matrix[x,y],ffFixed,7,2)+';';
      end;
    SetLength(s,Length(s)-1);
    sl.Add(s);
    end;
  sl.SaveToFile('..\testfiles\Matrix.tsv')}

end;

procedure TDockPrep.DoRun;
var
  f,maxcontacts:Integer;
  targetfile,probefile:string;
  targetout,probeout:string;
  jobfile:string;
  tchains,pchains:string;
  dormsd,center:Boolean;
  rotprobe,appendjobs:Boolean;
  submatfile:string;
  intdist,minalign:TFloat;
  contactfile:string;
  unconstrainedmodels,unconstrainedoverlap:Integer;
  submat:TSubMatrix;
  maxmodels:TIntegers;
  tmps:TSimpleStrings;

  //angle generation
  baseangles,shapepoints,maxaxes:Integer;
  displacement:TFloat;

  //symmetry constraints
  symmetry:Integer;
  symmetryplane:TFloat;

  LigandScore,InterfaceScore:Integer;
  pdbpath,outputpath:string;


  procedure GetOptions;

  var
    tmp:string;
    tmps:TSimpleStrings;
    f:Integer;

  begin
    minalign:=0.9;

    //angle generation
    baseangles:=2000;
    maxaxes:=2000;
    shapepoints:=20;
    displacement:=10;

    //symmetry constraints
    symmetry:=-1;
    symmetryplane:=-1;

    intdist:=5.0;

    pdbpath:=GetOptionValue('pdbpath');
    if pdbpath<>'' then pdbpath:=IncludeTrailingPathDelimiter(pdbpath);
    outputpath:=GetOptionValue('outputpath');
    if outputpath<>'' then outputpath:=IncludeTrailingPathDelimiter(outputpath);
    targetfile:=GetOptionValue('t');
    probefile:=GetOptionValue('p');
    jobfile:=GetOptionValue('j');
    targetout:=GetOptionValue('tout');
    probeout:=GetOptionValue('pout');
    tchains:=GetOptionValue('tchains');
    pchains:=GetOptionValue('pchains');
    center:=HasOption('center');
    rotprobe:=HasOption('randomize');
    submatfile:=GetOptionValue('submat');
    dormsd:=false;
    if submatfile<>'' then
      begin
      submat:=ReadBLASTMatrix(submatfile);
      dormsd:=HasOption('rmsd') or HasOption('intrmsd') or HasOption('fullrmsd');
      end;

    appendjobs:=HasOption('append');
    contactfile:=GetOptionValue('contacts');
    maxcontacts:=0;

    tmp:=GetOptionValue('baseangles');
    if tmp<>'' then baseangles:=StrToInt(tmp);
     tmp:=GetOptionValue('maxaxes');
    if tmp<>'' then maxaxes:=StrToInt(tmp);
    tmp:=GetOptionValue('shapepoints');
    if tmp<>'' then shapepoints:=StrToInt(tmp);
    tmp:=GetOptionValue('displacement');
    if tmp<>'' then displacement:=StrToFloat(tmp);
    tmp:=GetOptionValue('symmetryplane');
    if tmp<>'' then symmetryplane:=StrToFloat(tmp);

    tmp:=GetOptionValue('symmetry');
    if tmp<>'' then symmetry:=StrToInt(tmp);
    tmp:=GetOptionValue('maxcontacts');
    if tmp<>'' then maxcontacts:=StrToInt(tmp);
    tmp:=GetOptionValue('minalign');
    if tmp<>'' then minalign:=StrToFloat(tmp);
    tmp:=GetOptionValue('intdist');
    if tmp<>'' then intdist:=StrToFloat(tmp);


    tmp:=GetOptionValue('maxmodels');
    if tmp<>'' then
        begin
        tmps:=SplitString(tmp,';');
        SetLength(maxmodels,Length(tmps));
        for f:=0 to High(tmps) do
          maxmodels[f]:=StrToInt(tmps[f]);
        end
    else
      begin
      SetLength(maxmodels,1);
      maxmodels[0]:=5000;
      end;
    tmp:=GetOptionValue('rmsdixs');
    if tmp<>'' then
        begin
        tmps:=SplitString(tmp,';');
        LigandScore:=StrToInt(tmps[0]);
        InterfaceScore:=StrToInt(tmps[1]);
        end
    else
      begin
      LigandScore:=-1;
      InterfaceScore:=-1;
      end;

    unconstrainedmodels:=5000;
    //No unconstrained models by default if
    if HasOption('contacts') or HasOption('symmetry') then
      unconstrainedmodels:=-1;
    unconstrainedoverlap:=200;
    tmp:=GetOptionValue('unconstrained');
    if tmp<>'' then
      begin
      tmps:=SplitString(tmp,':');
      unconstrainedmodels:=StrToInt(tmps[0]);
      if Length(tmps)>1 then
        unconstrainedoverlap:=StrToInt(tmps[1]);
      end;
  end;

  procedure Error(Msg:string);

  begin
    WriteLn();
    WriteLn(Msg);
    WriteLn();
    WriteLn(targetfile);
    WriteLn(probefile);
    WriteLn(targetout);
    WriteLn(probeout);
    WriteLn(tchains);
    WriteLn(pchains);
    WriteLn(center);
    WriteLn(rotprobe);;

    Terminate;
  end;

  procedure ShowDockResults(FileName:string);

  var
    table:TSimpleStrings;
    f:Integer;

  begin
    table:=ReadAsTable(FileName);
    for f:=0 to High(table) do
      WriteLn(table[f]);
  end;

  procedure ExportStructures;

  var
    outfile:string;
    constraintset,modelid,minscoreix,maxscoreix:Integer;
    model:TMolecule;
    runs:TDockRuns;

    procedure ExportAllStructures(Run:TDockRun);

    var
      cs,id:Integer;
      modelman:TPDBModelMan;
      probe,target:TMolecule;
      sl:TStringList;
      pdbname:string;
      tmp,jobname:string;
      firstid,lastid:Integer;

    begin
      jobname:=ChangeFileExt(ExtractFileName(jobfile),'');
      tmp:=GetOptionValue('firstid');
      if tmp<>'' then firstid:=StrToInt(tmp)
        else firstid:=0;
      tmp:=GetOptionValue('lastid');
      if tmp<>'' then lastid:=StrToInt(tmp)
        else lastid:=High(Run.ConstraintSets[0].DockModels);

      sl:=TStringList.Create;
      modelman:=TPDBModelMan.Create(Config.MonomersPath);
      target:=modelman.LoadLayer(pdbpath+Run.TargetFile);
      probe:=modelman.LoadLayer(pdbpath+Run.ProbeFile);
      for cs:=0 to High(Run.ConstraintSets) do
        for id:=firstid to lastid do
          begin
          model:=GetModel(target,probe,Run,cs,id);
          pdbname:=jobname+'_'+IntToStr(cs)+'_'+IntToStr(id)+'.pdb';
          SaveToPDB(model,outputpath+pdbname);
          sl.add(pdbname);
          model.Free;
          end;
        modelman.Free;
        sl.SaveToFile(outputpath+ChangeFileExt(jobname,'')+'_models.txt');
        sl.Free;
    end;

  begin
    runs:=LoadOrders(jobfile);
    if HasOption('all') then
      ExportAllStructures(runs[0])
    else
      begin
      outfile:=GetOptionValue('outfile');
      constraintset:=StrToInt(GetOptionValue('constraintset'));
      modelid:=StrToInt(GetOptionValue('model'));
      minscoreix:=StrToInt(GetOptionValue('minscore'));
      maxscoreix:=StrToInt(GetOptionValue('maxscore'));
      if HasOption('model') then
        model:=GetModel(runs[0],constraintset,modelid);
      SaveToPDB(model,outputpath+outfile);
      model.Free;
      end;
    FreeOrders(runs);
  end;

  procedure AddRmsd;

  begin
    if HasOption('rmsd') then
      FDockMan.AddRMSDScore(GetOptionValue('rmsd'),submat,minalign,[rmsdFitTarget,rmsdScoreProbe],'Probe RMSD');
    if HasOption('intrmsd') then
      FDockMan.AddInterfaceRMSDScore(GetOptionValue('intrmsd'),intdist,submat,minalign);
    if HasOption('fullrmsd') then
      FDockMan.AddRMSDScore(GetOptionValue('fullrmsd'),submat,minalign,
        [rmsdFitTarget,rmsdFitProbe,rmsdScoreTarget,rmsdScoreProbe],'Full RMSD');
  end;

  procedure ExportModels;

  var
    infile:string;
    f,g:Integer;
    runs:TDockRuns;
    constixs,modelixs:TIntegers;
    s:string;
    sl:TStringList;
    scoreids:TSimpleStrings;

    function GetDescription(const DockRun:TDockRun;
        ConstIx,ModelIx:Integer;ScoreIds:TSimpleStrings):string;

    var f:Integer;

    begin
      Result:=IntToStr(ConstIx)+#9+IntToStr(ModelIx);
      with DockRun.ConstraintSets[ConstIx] do
        begin
        with DockModels[ModelIx] do
          begin
          Result:=Result+#9+FloatToStr(Rotation[0])
                        +#9+FloatToStr(Rotation[1])
                        +#9+FloatToStr(Rotation[2])
                        +#9+FloatToStr(Rotation[3]);
          Result:=Result+#9+FloatToStr(TransVec[0])
                        +#9+FloatToStr(TransVec[1])
                        +#9+FloatToStr(TransVec[2]);
          Result:=Result+#9+IntToStr(OverlapScore);
          for f:=0 to High(ScoreIds) do
            Result:=Result+#9+FloatToStrF(
              ScoreVal(DockRun.ConstraintSets[ConstIx],
                ModelIx,ScoreIds[f],-1), ffFixed,0,3);
          end;
        end;
    end;

  begin
    sl:=TStringList.Create;
    infile:=GetOptionValue('models');
    runs:=LoadOrders(infile);
    for f:=0 to High(runs) do
      begin
      sl.Add('Run '+IntToStr(f));
      if runs[f].ConstraintSets<>nil then
        begin
        scoreids:=ComputedScores(runs[f].ConstraintSets[0]);
        s:='Const.'+#9+'Model'+#9+'r'+#9+'i'+#9+'j'+#9+'k'+#9+'X'+#9+'Y'+#9+'Z'+#9+'Overlap';
        for g:=0 to High(scoreids) do
          s:=s+#9+scoreids[g];
        SortedModelList(runs[f], constixs, modelixs);
        sl.Add(s);
        for g:=0 to High(constixs) do
          sl.Add(GetDescription(runs[f],constixs[g],modelixs[g],scoreids));
        end;
      end;
    sl.SaveToFile(ChangeFileExt(infile,'')+'_table.tsv');
    FreeOrders(runs);
    sl.Free;
  end;


  procedure EditJob;

  begin
    FDockMan.LoadJobFile;
    if HasOption('delscores') then
      FDockMan.DeleteScores;
    if dormsd then
      begin
      FDockMan.PreparePartners;
      AddRmsd;
      end;
    FDockMan.SaveOrders(False);
  end;

  procedure Summary;

  var
    runs:TDockRuns;
    idlist:string;

    procedure WriteHeaders;

    var f:Integer;
        s:string;

    begin
      s:='';
      for f:=0 to High(maxmodels) do
        s:=s+'High '+IntToStr(maxmodels[f])+#9+'Good '+IntToStr(maxmodels[f])+#9+
          'Acceptable '+IntToStr(maxmodels[f])+#9+'Min RMSD '+IntToStr(maxmodels[f])+#9+
          'Min IntRMSD '+IntToStr(maxmodels[f])+#9;
      //WriteLn(#9+FlattenStrings(ComputedScores(runs[0].ConstraintSets[0]),#9+#9));
      WriteLn('Name'+#9+s+'Constraints TC'+#9+'Digitization TC'+#9+'Domain TC'+#9+'Scoring TC'+#9+
        'Total Rotations'+#9+'Done rotations');
    end;

    function Counts(NumModels:Integer; out AcceptList:string):string;

    var
      c,m,hi,good,accept:Integer;
      lscore,iscore,miniscore,minlscore:TFloat;
      buildlist:Boolean;

    begin
      { TODO : Auto detect probe and interface rmsd scores if LigandScore or InterfaceScore <0 }
      {LigandScore=0;
      InterfaceScore=1;}

      buildlist:=HasOption('list');
      hi:=0;
      good:=0;
      accept:=0;
      minlscore:=100000;
      miniscore:=100000;
      AcceptList:='';
      for c:=0 to High(runs[0].ConstraintSets) do
        for m:=0 to High(runs[0].ConstraintSets[c].ScoreResults[0].ScoreVals) do
          begin
          lscore:=runs[0].ConstraintSets[c].ScoreResults[LigandScore].ScoreVals[m];
          iscore:=runs[0].ConstraintSets[c].ScoreResults[InterfaceScore].ScoreVals[m];
          if lscore<minlscore then minlscore:=lscore;
          if iscore<miniscore then miniscore:=iscore;
          if (lscore<1) or (iscore<1) then Inc(hi)
          else if (lscore<5) or (iscore<2) then Inc(good)
          else if (lscore<10) or (iscore<4) then Inc(accept);
          if buildlist and ((lscore<10) or (iscore<4)) then
            AcceptList:=AcceptList+#9+IntToStr(m);
          if m>=NumModels then Break;
          end;
      Result:=IntToStr(hi)+#9+IntToStr(good)+#9+IntToStr(accept)+
              #9+FloatToStrF(minlscore,ffFixed,2,2)+#9+FloatToStrF(miniscore,ffFixed,2,2);
    end;

  var
    f:Integer;
    rep:string;

  begin
    try
    runs:=LoadOrders(jobfile);
    except
      WriteLn('**'+jobfile+' not found**');
      Halt;
    end;
    if GetOptionValue('summary')='h' then
      WriteHeaders;
    rep:=ExtractFileName(jobfile)+#9;
    for f:=0 to High(maxmodels) do
      rep:=rep+Counts(maxmodels[f],idlist)+#9;
    WriteLn(rep,runs[0].Stats.ConstraintsTickCount,#9,
        runs[0].Stats.DigitizationTickCount,#9,
        runs[0].Stats.DomainTickCount,#9,
        runs[0].Stats.ScoringTickCount, #9,
        Length(runs[0].Rotations),#9,
        runs[0].CompletedRotations,idlist);
    FreeOrders(runs);
  end;

  procedure ContactScore;

  var
    runs:TDockRuns;
    run:TDockRun;
    cs,id:Integer;
    modelman:TPDBModelMan;
    probe,target:TMolecule;
    sl:TStringList;
    pdbname:string;
    matname:string;
    contactmat:TSubMatrix;
    scores:TFLoats;
  begin
    matname:=GetOptionValue('contactscore');
    contactmat:=ReadBLASTMatrix(matname);
    runs:=LoadOrders(jobfile);
    run:=runs[0];
    modelman:=TPDBModelMan.Create(Config.MonomersPath);
    target:=modelman.LoadLayer(pdbpath+run.TargetFile);
    probe:=modelman.LoadLayer(pdbpath+run.ProbeFile);
    for cs:=0 to High(run.ConstraintSets) do
      begin
      scores:=ResidueContactScore(target,probe,run.ConstraintSets[cs].DockModels,
              contactmat,3);
      for id:=0 to High(scores) do
          WriteLn(scores[id]);
      end;
    modelman.Free;
    FreeOrders(runs);
  end;

begin
  Randomize;

  { TODO : Check best way of doing this }
  DecimalSeparator:='.';

  // parse parameters
  if HasOption('h','help') then begin
    WriteHelp;
    Terminate;
    Exit;
  end;

  { add your program here }
  LoadAtomData;
  LoadAAData;

  GetOptions;
  if HasOption('test') then
    RunTests
  else if jobfile='' then
    Error('Job file is mandatory.')
  else
    begin
    FDockMan:=TDockOrdersManager.Create(jobfile);
    if HasOption('contactscore') then
      ContactScore;
    if HasOption('summary') then
      Summary;
    if HasOption('export') then
      ExportStructures;
    if HasOption('models') then
      ExportModels;
    if HasOption('table') then
      ShowDockResults(GetOptionValue('table'));
    if HasOption('modify') then
      EditJob;
    if HasOption('t') and HasOption('p') then
      begin
      if ((tchains<>'') or (pchains<>'') or center or rotprobe) and
        ((targetout='') or (probeout='')) then
      Error('No output files for saving modified structures.')
      else
        begin
        FDockMan:=TDockOrdersManager.Create(jobfile);
        FDockMan.NewDockRun;
        FDockMan.PreparePartners(targetfile,probefile,tchains,pchains,center,rotprobe);
        if targetout<>'' then SaveToPDB(FDockMan.Target,targetout)
          else targetout:=targetfile;
        if probeout<>'' then SaveToPDB(FDockMan.Probe,probeout)
          else probeout:=probefile;
        FDockMan.SetPDBFiles(targetout,probeout);
        FDockMan.BuildShapes;
        FDockMan.BuildRotations(shapepoints, baseangles, displacement,symmetry,maxaxes);
        if dormsd then
          AddRmsd;
        if contactfile<>'' then
          FDockMan.AddResidueContacts(contactfile, maxcontacts);
        if symmetryplane>0 then
          FDockMan.AddSymmetryConstraint(symmetryplane);
        if unconstrainedmodels>0 then
          FDockMan.AddNullConstraintSet(unconstrainedmodels,unconstrainedoverlap);
        FDockMan.SaveOrders(appendjobs);
        end;
      end;
      FDockMan.Free;
    end;

  // stop program loop
  Terminate;
end;

constructor TDockPrep.Create(TheOwner: TComponent);
begin
  inherited Create(TheOwner);
  StopOnException:=True;
end;

destructor TDockPrep.Destroy;
begin
  inherited Destroy;
end;

procedure TDockPrep.WriteHelp;
begin
  { add your help code here }
  writeln('Usage: ',ExeName,' -t "target.pdb/target.pdb.gz" -p "probe.pdb/probe.pdb.gz" -j jobfile');
  WriteLn('Docking Options:');
  WriteLn('-tchains TargetChains');
  WriteLn('-pchains ProbeChains');
  WriteLn('  Specify which chains to use, generating new pdb files (-tout -pout)');
  WriteLn('-tout pdbfile');
  WriteLn('-pout pdbfile');
  WriteLn('  Mandatory if center, chains or randomize');
  WriteLn('-center : centers both target and probe');
  WriteLn('-randomize : randomizes orientation of probe');
  WriteLn('-fixedzrotation Angle');
  WriteLn('  probe will be rotated only this angle (degrees) for each rotation axis');
  WriteLn('-symmetry margin');
  WriteLn('  probe will be restricted to a plane perpendicular to the rotation axis, intersecting the center of the target');
  WriteLn('-contacts ContactConstraintsFile ');
  WriteLn('-minalign minimum alignment match for *rmsd, from 0 to 1');
  WriteLn('-submat substitution_matrix.txt');
  WriteLn('-rmsd template.pdb : scores probe by fitting target');
  WriteLn('-intrmsd template.pdb : scores interface region');
  WriteLn('-fullrmsd template.pdb :scores whole complex');
  WriteLn('-intdist interface distance, A');
  WriteLn('-append');
  WriteLn('-maxcontacts');
  WriteLn('-unconstrained models:minoverlap');
  WriteLn('  default is 5000 models, 200 of min overlap');
  WriteLn('  or no models if -contacts is specified');
  WriteLn('Processing existing files (requires existing -j file):');
  WriteLn('-summary: stats for scores');
  WriteLn('   -maxmodels: ; separated values for maximum number of models considered in summary');
  WriteLn('   -rmsdixs: ; separated values for the index of probe and interface rmsd');

  {  'export') then
    'models') then
    'table') then
    'modify') then

}

end;

var
  Application: TDockPrep;
begin
  Application:=TDockPrep.Create(nil);
  Application.Run;
  Application.Free;
end.

