{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 2014 02 21
Purpose:
  prepares docking orders and structure files, reads dock files


  for docking:
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
      5A (the interface)



    -delscores
      deletes all additional scores
      TODO: enable specifying score id and deleting only those

  Reading results

    list models in each constraintset
    -table filename
    displays on screen, each constraint in one group of columns, then model id, overlap, other scores


    list all models
    -models filename
    displays on screen table with single list of models and all scores (removes duplicates)

    exports pdb
    -export filename
    -outfile filename.pdb
    -constraintset number
    -model  id

    summary:
    -summary [h]
    outputs summary to screen, with headers if h



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
  Classes, SysUtils, docktasks, molecules, pdbmolecules, oclconfiguration, geomutils,
  stringutils, quicksort, basetypes, CustApp
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

procedure TDockPrep.DoRun;
var
  f,maxcontacts:Integer;
  targetfile,probefile:string;
  targetout,probeout:string;
  jobfile:string;
  tchains,pchains:string;
  center:Boolean;
  rotprobe,appendjobs:Boolean;
  rmsdtemplate,intrmsdtemplate:string;
  intdistance:TFloat;
  contactfile:string;
  unconstrainedmodels,unconstrainedoverlap:Integer;


  procedure GetOptions;

  var
    tmp:string;
    tmps:TSimpleStrings;

  begin
    intdistance:=5;
    { TODO : add this to some parameter }

    targetfile:=GetOptionValue('t');
    probefile:=GetOptionValue('p');
    jobfile:=GetOptionValue('j');
    targetout:=GetOptionValue('tout');
    probeout:=GetOptionValue('pout');
    tchains:=GetOptionValue('tchains');
    pchains:=GetOptionValue('pchains');
    center:=HasOption('center');
    rotprobe:=HasOption('randomize');
    rmsdtemplate:=GetOptionValue('rmsd');
    intrmsdtemplate:=GetOptionValue('intrmsd');
    appendjobs:=HasOption('append');
    contactfile:=GetOptionValue('contacts');
    maxcontacts:=0;
    tmp:=GetOptionValue('maxcontacts');
    if tmp<>'' then maxcontacts:=StrToInt(tmp);
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
    infile,outfile:string;
    constraintset,modelid:Integer;
    model:TMolecule;
    runs:TDockRuns;

  begin
    infile:=GetOptionValue('export');
    outfile:=GetOptionValue('outfile');
    constraintset:=StrToInt(GetOptionValue('constraintset'));
    modelid:=StrToInt(GetOptionValue('model'));
    runs:=LoadOrders(infile);
    model:=GetModel(runs[0],constraintset,modelid);
    SaveToPDB(model,outfile);
    FreeOrders(runs);
    model.Free;
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
    if (rmsdtemplate<>'') or (intrmsdtemplate<>'') then
      begin
      FDockMan.PreparePartners;
      if rmsdtemplate<>'' then
        FDockMan.AddRMSDScore(rmsdtemplate);
      if intrmsdtemplate<>'' then
        FDockMan.AddInterfaceRMSDScore(intrmsdtemplate,intdistance);
      end;
    FDockMan.SaveOrders(False);
  end;

  procedure Summary;

  var
    runs:TDockRuns;

    procedure WriteHeaders;

    begin
      //WriteLn(#9+FlattenStrings(ComputedScores(runs[0].ConstraintSets[0]),#9+#9));
      WriteLn('Name'+#9+'High'+#9+'Good'+#9+'Acceptable'+#9);
    end;

    function Counts(MaxModels:Integer):string;


    const
      LigandScore=0;
      InterfaceScore=1;

    var
      c,m,hi,good,accept:Integer;
      lscore,iscore:TFloat;

    begin
      hi:=0;
      good:=0;
      accept:=0;
      for c:=0 to High(runs[0].ConstraintSets) do
        for m:=0 to High(runs[0].ConstraintSets[c].ScoreResults[0].ScoreVals) do
          begin
          lscore:=runs[0].ConstraintSets[c].ScoreResults[LigandScore].ScoreVals[m];
          iscore:=runs[0].ConstraintSets[c].ScoreResults[InterfaceScore].ScoreVals[m];
          if (lscore<1) or (iscore<1) then Inc(hi)
          else if (lscore<5) or (iscore<2) then Inc(good)
          else if (lscore<10) or (iscore<4) then Inc(accept);
          if m>=MaxModels then Break;
          end;
      Result:=IntToStr(hi)+#9+IntToStr(good)+#9+IntToStr(accept);
    end;

  begin
    runs:=LoadOrders(jobfile);
    if GetOptionValue('summary')='h' then
      WriteHeaders;
    WriteLn(jobfile{+#9+Counts(999)},#9,runs[0].Stats.ConstraintsTickCount+
        runs[0].Stats.DigitizationTickCount +
        runs[0].Stats.DomainTickCount +
        runs[0].Stats.ScoringTickCount);
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

  GetOptions;
  if jobfile='' then
    Error('Job file is mandatory.')
  else
    begin
    FDockMan:=TDockOrdersManager.Create(jobfile);
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
        if rmsdtemplate<>'' then
          FDockMan.AddRMSDScore(rmsdtemplate);
        if contactfile<>'' then
          FDockMan.AddResidueContacts(contactfile, maxcontacts);
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
  writeln('Usage: ',ExeName,' -T"target.pdb/target.pdb.gz" -P"probe.pdb/target.pdb.gz" '+
          '--tchainsTargetChains --pchainsProbeChains --contactsContactConstraintsFile --RMSDtemplate.pdb');
  writeln('Usage: ',ExeName,' -h');
end;

var
  Application: TDockPrep;
begin
  Application:=TDockPrep.Create(nil);
  Application.Run;
  Application.Free;
end.

