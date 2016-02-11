{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 31.3.2011
Purpose:
  Geometric docking module
  Name from original algorithm (Boolean Geometric Interaction Evaluation)
Revisions:
TO DO:
  split up
*******************************************************************************}
unit bogie;

{$mode objfpc}{$H+}

interface

uses SysUtils,basetypes,linegrids,dockdomains,surface,geomutils,docktasks;

type

   { TModelManager }

  TModelManager=class
    private
    FModels:TDockModels;
    FMinOverlap:Integer;
    FLastModel:Integer;
    public
    CurrentRotation:TQuaternion;
    GridScale:TFloat;
    ProbeGridTransVec:TCoord;
    property Models:TDockModels read FModels;
    property MinimumOverlap:Integer read FMinOverlap;
    constructor Create(NumModels,MinOverlap:Integer;DockModels:TDockModels);
    function AddModel(Score,X,Y,Z:Integer):Integer;
    procedure ImportModels(const CSet:TConstraintSet);
    procedure ExportModels(var CSet:TConstraintSet);

  end;

  TModelManagers=array of TModelManager;


  { TDockManager }

  TDockManager=class
     private
     //These are copies so they cannot be changed outside
     FTargetCoords,FProbeCoords:TCoords;
     FTargetRads,FProbeRads:TFloats;
     FResolution:TFloat;
     //add coord and radius for soft docking

     FRotations:TQuaternions;
     FAxes:TCoords;
     FCurrentRotation:Integer;

     FTargetGrid,FProbeGrid:TDockingGrid;

     FConstraintSets:TConstraintSets;
     FModelManagers:TModelManagers;

     procedure ClearModelManagers;

     public
     //Statistics
     DigitizationTicks,ConstraintTicks,DomainTicks,ScoringTicks:Integer;

     property TargetGrid:TDockingGrid read FTargetGrid;
     property ProbeGrid:TDockingGrid read FProbeGrid;
     property CurrentRotation:Integer read FCurrentRotation write FCurrentRotation;
     function TotalRotations:Integer;
     function Models(DockIx:Integer):TDockModels;
     function ModelGroupsCount:Integer;
     constructor Create(TargCoords:TCoords;TargRads:TFloats;
                        ProbCoords:TCoords;ProbRads:TFloats;
                        AResolution:TFLoat);
     procedure Free;
     procedure ImportConstraintSets(var CSets:TConstraintSets);
     procedure ExportResults(var CSets:TConstraintSets);
     procedure BuildTargetGrid;
     procedure BuildProbeGrid(Rotation:TQuaternion);
     procedure BuildRotations(Job:TDockRun);
     function Dock(NumSteps:Integer=1;CalcStats:Boolean=False):Boolean;
      //Returns true if finished



  end;

implementation

{ TModelManager }

constructor TModelManager.Create(NumModels, MinOverlap: Integer;DockModels:TDockModels);

var f:Integer;

begin
  inherited Create;
  SetLength(FModels,NumModels);
  for f:=0 to High(FModels) do
    if f<Length(DockModels) then
      FModels[f]:=DockModels[f]
    else
      begin
      FModels[f].OverlapScore:=-1;
      FModels[f].TransVec:=NullVector;
      FModels[f].Rotation:=IdentityQuaternion;
      end;
  FLastModel:=High(DockModels);
  if FLastModel>=0 then
      FMinOverlap:=Min(FModels[FLastModel].OverlapScore,MinOverlap)
  else FMinOverlap:=MinOverlap;
end;

function TModelManager.AddModel(Score, X, Y, Z: Integer): Integer;

var
  ix1,ix2:Integer;

begin
  if Score>FMinOverlap then
    { TODO : Interaction tests here }
    begin
    ix1:=0;
    while FModels[ix1].OverlapScore>=Score do
      ix1:=ix1+1;
      //This must stop if Score>FMinOverlap
    if FLastModel<High(FModels) then
      Inc(FLastModel);
    ix2:=FLastModel;
    while ix2>ix1 do
      begin
      FModels[ix2]:=FModels[ix2-1];
      Dec(ix2);
      end;
    FModels[ix1].OverlapScore:=Score;
    FModels[ix1].TransVec:=Add(ProbeGridTransvec,
                              Coord(X*GridScale,Y*GridScale,Z*GridScale));
    FModels[ix1].Rotation:=CurrentRotation;
    if FLastModel=High(FModels) then
      FMinOverlap:=FModels[FLastModel].OverlapScore;
    Result:=FMinOverlap;

    end;
end;

procedure TModelManager.ImportModels(const CSet: TConstraintSet);

var f:Integer;

begin
  for f:=0 to High(CSet.DockModels) do
    begin
    if f>High(FModels) then Break;
    FModels[f]:=CSet.DockModels[f];
    end;
   FMinOverlap:=Max(FMinOverlap,FModels[High(FModels)].OverlapScore);
end;

procedure TModelManager.ExportModels(var CSet: TConstraintSet);
begin
  CSet.DockModels:=Copy(FModels,0,Length(FModels));
end;


{ TDockManager }

procedure TDockManager.ClearModelManagers;

var f:Integer;

begin
  for f:=0 to High(FModelManagers) do
    FModelManagers[f].Free;
  FModelManagers:=nil;
end;

function TDockManager.TotalRotations: Integer;
begin
  Result:=Length(FRotations);
end;

function TDockManager.Models(DockIx: Integer): TDockModels;
begin
  Result:=FModelManagers[DockIx].Models;
end;

function TDockManager.ModelGroupsCount: Integer;
begin
  Result:=Length(FModelManagers);
end;

constructor TDockManager.Create(TargCoords:TCoords;TargRads:TFloats;
                        ProbCoords:TCoords;ProbRads:TFloats;
                        AResolution:TFLoat);
begin
  inherited Create;

  Assert((Length(TargCoords)=Length(TargRads)) and (Length(ProbCoords)=Length(ProbRads)),
        'Mismatching coordinate and radius arrays');

  FTargetCoords:=Copy(TargCoords,0,Length(TargCoords));
  FProbeCoords:=Copy(ProbCoords,0,Length(ProbCoords));
  FTargetRads:=Copy(TargRads,0,Length(TargRads));
  FProbeRads:=Copy(ProbRads,0,Length(ProbRads));

  FResolution:=AResolution;

  FTargetGrid:=TDockingGrid.Create(FResolution);
  FProbeGrid:=TDockingGrid.Create(FResolution);
  FModelManagers:=nil;
end;

procedure TDockManager.Free;
begin
  FTargetGrid.Free;
  FProbeGrid.Free;
  ClearModelManagers;
  inherited Free;
end;

procedure TDockManager.ImportConstraintSets(var CSets: TConstraintSets);

var f:Integer;

begin
  FConstraintSets:=CSets;
  ClearModelManagers;
  SetLength(FModelManagers,Length(FConstraintSets));
  for f:=0 to High(FModelManagers) do
    begin
    FModelManagers[f]:=TModelManager.Create(FConstraintSets[f].NumModels,
                                            FConstraintSets[f].MinOverlap,
                                            FConstraintSets[f].DockModels);
    FModelManagers[f].GridScale:=FResolution;
    FModelManagers[f].ImportModels(CSets[f]);
    end;
end;

procedure TDockManager.ExportResults(var CSets: TConstraintSets);

var f:Integer;

begin
  for f:=0 to High(FModelManagers) do
    FModelManagers[f].ExportModels(CSets[f]);
end;

procedure TDockManager.BuildTargetGrid;
begin
  FTargetGrid.BuildFromSpheres(FTargetCoords,FTargetRads);
end;

procedure TDockManager.BuildProbeGrid(Rotation: TQuaternion);

var
  tmpcoords:TCoords;

begin
  tmpcoords:=Rotate(FProbeCoords,Rotation);
  FProbeGrid.BuildFromSpheres(tmpcoords,FProbeRads);
end;

procedure TDockManager.BuildRotations(Job: TDockRun);

var f:Integer;

begin
  FRotations:=Job.Rotations;
  SetLength(FAxes,Length(FRotations));
  for f:=0 to High(Job.AxisIndexes) do
    FAxes[f]:=Job.Axes[Job.AxisIndexes[f]];
end;

function TDockManager.Dock(NumSteps: Integer; CalcStats: Boolean): Boolean;

var
  dockdomain:TDockDomain;
  f:Integer;
  ticktime:DWORD;

begin
  while (FCurrentRotation<High(FRotations)) and (NumSteps>0) do
    begin
    Inc(FCurrentRotation);
    Dec(NumSteps);

    if CalcStats then ticktime:=GetTickCount;

    BuildProbeGrid(FRotations[FCurrentRotation]);

    if CalcStats then
      begin
      DigitizationTicks:=GetTimeInteval(ticktime);
      ConstraintTicks:=0;
      DomainTicks:=0;
      ScoringTicks:=0;
      end;

    for f:=0 to High(FConstraintSets) do
      begin
      FModelManagers[f].CurrentRotation:=FRotations[FCurrentRotation];
      dockdomain:=TDockDomain.Create(FTargetGrid,FProbeGrid);
      dockdomain.ConstraintManager.ImportConstraints(FConstraintSets[f].Constraints,
        FRotations[FCurrentRotation],FAxes[FCurrentRotation]);

      if CalcStats then
        ConstraintTicks:=ConstraintTicks+GetTimeInteval(ticktime);

      dockdomain.AssignModelManager(@FModelManagers[f].AddModel);
      dockdomain.MinimumOverlap:=FModelManagers[f].MinimumOverlap;
      dockdomain.RemoveCores:=True;
      dockdomain.BuildInitialDomain;

      if CalcStats then
        DomainTicks:=DomainTicks+GetTimeInteval(ticktime);

      FModelManagers[f].ProbeGridTransVec:=Subtract(FProbeGrid.TransVec,FTargetGrid.TransVec);
      dockdomain.Score;
      if CalcStats then
        ScoringTicks:=ScoringTicks+GetTimeInteval(ticktime);

      dockdomain.Free;
      end;
    end;
  Result:=FCurrentRotation=High(FRotations);
end;


end.

