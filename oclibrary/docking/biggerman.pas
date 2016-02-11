{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 2014 02 24
Purpose:
  Manage dockings
Requirements:
Revisions:
To do:
*******************************************************************************}

unit biggerman;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, bogie, docktasks, molecules,
  geomutils, molutils, oclconfiguration, linegrids,rmsd;


type

  { TBiGGERManager }

  TBiGGERManager=class
    private
      FJobs:TDockRuns;
      FCurrentJobFile:string;

      procedure RunJob(Ix:Integer);
      procedure RunGeometric(Ix:Integer);
      function ComputeRMSD(const ScoreDef:TScoreDef;const DockModels:TDockModels):TFloats;
      function ComputeScore(const ScoreDef:TScoreDef;const DockModels:TDockModels):TFloats;
    public
      Verbosity:Integer;
      constructor Create;
      procedure LoadJobs(JobFile:string);
      procedure SaveCurrent(JobFile:string);
      procedure RunJobs;
      procedure ExportGrids(FileName:string);
      procedure Free;
  end;

implementation

{ TBiGGERManager }

procedure TBiGGERManager.RunJob(Ix: Integer);


function Computed(ScoreId:string;ScoreResults:TScoreResults):Boolean;

var f:Integer;

begin
  Result:=False;
  for f:=0 to High(ScoreResults) do
    if ScoreResults[f].ScoreId=ScoreID then
      begin
      Result:=True;
      Break;
      end;
end;

var score,constraint :Integer;

begin
  with FJobs[Ix] do
    begin
    if CompletedRotations<Length(Rotations) then
      RunGeometric(Ix);
    for score:=0 to High(ScoreDefs) do
      for constraint:=0 to High(ConstraintSets) do
        if not Computed(ScoreDefs[score].Name, ConstraintSets[constraint].ScoreResults) then
          with ConstraintSets[constraint] do
            begin
            WriteLn('Computing ',ScoreDefs[score].Name,' of ',constraint);
            SetLength(ScoreResults,Length(ScoreResults)+1);
            with ScoreResults[High(ScoreResults)] do
              begin
              ScoreId:=ScoreDefs[score].Name;
              ScoreVals:=ComputeScore(ScoreDefs[score],ConstraintSets[constraint].DockModels);
              end;
            end;
    end;
end;

procedure TBiGGERManager.RunGeometric(Ix: Integer);

var
  targetrads,proberads:TFloats;
  targetcoords,probecoords:TCoords;
  initialtime,lasttime,currenttime,lastsave:Int64;
  totalsecs:TFloat;
  dockman:TDockManager;
  s:string;
  f:Integer;
  tick:DWORD;
  rotationcount:Integer;

begin
  with FJobs[ix] do
    begin
    //Load structures and get atoms
    tick:=GetTickCount;
    targetrads:=Target.Rads;
    targetcoords:=Target.Coords;
    proberads:=Probe.Rads;
    probecoords:=Probe.Coords;
    //setup dockmanager with constraints and build all rotations
    dockman:=TDockManager.Create(targetcoords,targetrads,probecoords,proberads,
             Resolution);
    dockman.ImportConstraintSets(FJobs[Ix].ConstraintSets);
    dockman.BuildTargetGrid;
    dockman.BuildRotations(FJobs[ix]);
    dockman.CurrentRotation:=CompletedRotations-1;

    initialtime:=GetTickCount;
    lasttime:=initialtime;
    lastsave:=initialtime;
    rotationcount:=0;
    while not dockman.Dock(1,True) do
      begin
      Inc(rotationcount);
      with Stats do
        begin
        DigitizationTickCount:=DigitizationTickCount+dockman.DigitizationTicks;
        ConstraintsTickCount:=ConstraintsTickCount+dockman.ConstraintTicks;
        DomainTickCount:=DomainTickCount+dockman.DomainTicks;
        ScoringTickCount:=ScoringTickCount+dockman.ScoringTicks;
        //TestedModelsCount,
        //InsertedModelsCount:Integer;
        end;
      currenttime:=GetTickCount;

      if currenttime-lastsave>FJobs[ix].SecondsBetweenSaves*1000 then
        begin
        WriteLn('Autosaving');
        lastsave:=currenttime;
        CompletedRotations:=dockman.CurrentRotation+1;
        dockman.ExportResults(FJobs[Ix].ConstraintSets);
        SaveCurrent(FCurrentJobFile);
        end;

      if currenttime-lasttime>5000 then
          begin
          totalsecs:=(currenttime-initialtime)/1000;
          WriteLn('Current Job: ',FCurrentJobFile,' (',
            IntToStr(ix+1),' of ',IntToStr(Length(FJobs))+')');
          WriteLn('Angle: ',dockman.CurrentRotation,' of ',Length(Rotations),'; Total time:',
            FloatToStrF(totalsecs,ffFixed,0,2),'s Average:',
            FloatToStrF(totalsecs/rotationcount,ffFixed,0,2),
            's Estimated:',
            FloatToStrF(Length(Rotations)*totalsecs/rotationcount,ffFixed,0,2));
          lasttime:=GetTickCount;
          s:='';
          for f:=0 to dockman.ModelGroupsCount-1 do
            s:=s+IntToStr(dockman.Models(f)[0].OverlapScore)+';';
          Writeln('Best models:',s);
          end;
      end;
    CompletedRotations:=dockman.CurrentRotation+1;
    dockman.ExportResults(FJobs[Ix].ConstraintSets);
    Stats.TotalTickCount:=Stats.TotalTickCount+GetTimeInteval(tick);
    end;
end;

function TBiGGERManager.ComputeRMSD(const ScoreDef: TScoreDef;
  const DockModels: TDockModels): TFloats;

var
  rmsdcalc:TRMSDCalculator;
  f,g:Integer;
  tmpcoords:TCoords;
  rdef:TRMSDScoreDef;
  fixed,mobile:TCoords;

  function GetCoords(const Original:TCoords; ProbeStart:Integer;Before,After:Boolean):TCoords;

  begin
    if Before then
      Result:=Copy(Original,0,ProbeStart-1)
    else Result:=nil;
    if After then
      AppendToArray(Copy(Original,ProbeStart,Length(Original)),Result);
  end;

begin
  WriteLn(ScoreDef.Name);
  rmsdcalc:=TRMSDCalculator.Create;
  rdef:=TRMSDScoreDef(ScoreDef);
  tmpcoords:=Copy(rdef.Predicted,0,Length(rdef.Predicted));
  SetLength(Result,Length(DockModels));
  WriteLn(Length(DockModels));
  for f:=0 to High(DockModels) do
    begin
    for g:=rdef.ProbeStart to High(tmpcoords) do
      tmpcoords[g]:=Add(Rotate(rdef.Predicted[g],DockModels[f].Rotation),DockModels[f].TransVec);
    rmsdcalc.Initialize;
    fixed:=GetCoords(rdef.Actual,Rdef.ProbeStart,rdef.FitTarget,rdef.FitProbe);
    mobile:=GetCoords(tmpcoords,Rdef.ProbeStart,rdef.FitTarget,rdef.FitProbe);
    rmsdcalc.SetCoords(fixed,mobile);
    rmsdcalc.Minimise(100000,0.001);
    fixed:=GetCoords(rdef.Actual,Rdef.ProbeStart,rdef.ScoreTarget,rdef.ScoreProbe);
    mobile:=GetCoords(tmpcoords,Rdef.ProbeStart,rdef.ScoreTarget,rdef.ScoreProbe);
    mobile:=rmsdcalc.TransformedCoords(mobile);
    Result[f]:=StaticRMSD(fixed,mobile);
    end;
  rmsdcalc.Free;
end;

function TBiGGERManager.ComputeScore(const ScoreDef: TScoreDef;
  const DockModels: TDockModels): TFloats;
begin
  case ScoreDef.ScoreType of
    ScoreRMSD:Result:=ComputeRMSD(ScoreDef,DockModels);
  end;

end;

constructor TBiGGERManager.Create;
begin
  inherited;
  Verbosity:=10;
end;

procedure TBiGGERManager.LoadJobs(JobFile: string);
begin
  FJobs:=docktasks.LoadOrders(JobFile);
  FCurrentJobFile:=JobFile;
end;

procedure TBiGGERManager.SaveCurrent(JobFile: string);
begin
  docktasks.SaveOrders(FJobs,JobFile,False);
end;

procedure TBiGGERManager.RunJobs;

var f:Integer;

begin
  for f:=0 to High(FJobs) do
    RunJob(f);
end;

procedure TBiGGERManager.ExportGrids(FileName: string);
var
  targetrads,proberads:TFloats;
  targetcoords,probecoords:TCoords;
  dockman:TDockManager;
  ix:Integer;
  sl:TStringList;

  procedure DescribeGrid(DockGrid:TDockingGrid);

  var x,y,z:Integer;
  begin
    sl.Add('::Core');
    with DockGrid.Core do
      for x:=0 to High(Grid) do
        for y:=0 to High(Grid[x]) do
          for z:=0 to High(Grid[x,y]) do
            Sl.Add(IntToStr(x)+','+IntToStr(y)+','+IntToStr(Grid[x,y,z,0])+':'+IntToStr(Grid[x,y,z,1]));
    sl.Add('::Surface');
    with DockGrid.Surf do
      for x:=0 to High(Grid) do
        for y:=0 to High(Grid[x]) do
          for z:=0 to High(Grid[x,y]) do
            Sl.Add(IntToStr(x)+','+IntToStr(y)+','+IntToStr(Grid[x,y,z,0])+':'+IntToStr(Grid[x,y,z,1]));

  end;

begin
  sl:=TStringList.Create;
  for ix:=0 to High(FJobs) do with FJobs[ix] do
    begin
    targetrads:=Target.Rads;
    targetcoords:=Target.Coords;
    proberads:=Probe.Rads;
    probecoords:=Probe.Coords;
    //setup dockmanager with constraints
    dockman:=TDockManager.Create(targetcoords,targetrads,probecoords,proberads,
             Resolution);
    dockman.ImportConstraintSets(FJobs[Ix].ConstraintSets);
    dockman.BuildTargetGrid;
    dockman.BuildProbeGrid(IdentityQuaternion);
    sl.Add(IntToStr(ix));
    sl.Add(':'+TargetFile);
    DescribeGrid(dockman.TargetGrid);
    sl.Add(':'+ProbeFile);
    DescribeGrid(dockman.ProbeGrid);
    dockman.Free;
    end;
  sl.SaveToFile(FileName);
  sl.Free;
end;

procedure TBiGGERManager.Free;
begin
  if Self<>nil then
    begin
    FreeOrders(FJobs);
    inherited;
    end;
end;

end.

