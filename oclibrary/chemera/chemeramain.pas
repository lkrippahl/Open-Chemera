{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 9.1.2011
Purpose:
  This is the main form for Chemera.
Requirements:
Revisions:
To do:
*******************************************************************************}

unit chemeramain;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, FileUtil, Forms, Controls, Graphics, Dialogs, Menus,
  pdbmolecules, oglform, displayobjects, LCLProc, oclconfiguration, molecules,
  molutils, geomutils, basetypes, linegrids, base3ddisplay, dockdomains, bogie;

type

  { TCmMainForm }

  TCmMainForm = class(TForm)
    MainMenu: TMainMenu;
    FileMn: TMenuItem;
    MenuItem1: TMenuItem;
    MenuItem2: TMenuItem;
    MenuItem3: TMenuItem;
    MenuItem4: TMenuItem;
    MenuItem5: TMenuItem;
    MenuItem6: TMenuItem;
    MenuItem7: TMenuItem;
    OpenDlg: TOpenDialog;
    OpenMn: TMenuItem;
    procedure FormActivate(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure MenuItem2Click(Sender: TObject);
    procedure MenuItem3Click(Sender: TObject);
    procedure MenuItem4Click(Sender: TObject);
    procedure MenuItem5Click(Sender: TObject);
    procedure MenuItem6Click(Sender: TObject);
    procedure MenuItem7Click(Sender: TObject);
    procedure OpenMnClick(Sender: TObject);
  private
    { private declarations }
    FInitializing:Boolean;
    FDisplay:TOpenGlForm;
    FDispMan:TDisplayManager;
    FMolecules:TPdbModelMan;

    //Test functions
    procedure TestGrids;
    procedure TestRotationSampling;
    procedure TestBogie;
    procedure OptimizeAddedRad;
    procedure TestConstraints;
    procedure TestAngles;

    procedure InitChemera;
  public
    { public declarations }
  end; 

var
  CmMainForm: TCmMainForm;

implementation

{$R *.lfm}

{ TCmMainForm }

procedure TCmMainForm.FormCreate(Sender: TObject);
begin
  FInitializing:=True;
end;

procedure TCmMainForm.MenuItem2Click(Sender: TObject);
begin
  TestGrids;
end;

procedure TCmMainForm.MenuItem3Click(Sender: TObject);
begin
  TestRotationSampling;
end;

procedure TCmMainForm.MenuItem4Click(Sender: TObject);
begin
  TestBogie;
end;

procedure TCmMainForm.MenuItem5Click(Sender: TObject);
begin
  OptimizeAddedRad;
end;

procedure TCmMainForm.MenuItem6Click(Sender: TObject);
begin
  TestConstraints;
end;

procedure TCmMainForm.MenuItem7Click(Sender: TObject);
begin
  TestAngles;
end;

procedure TCmMainForm.OpenMnClick(Sender: TObject);

var mol:TMolecule;

begin
  OpenDlg.Filter:='';//Pdb file|*.pdb;GZipped PDB file|*.gz';
  if OpenDlg.Execute then
    begin
    mol:=FMolecules.LoadLayer(OpenDlg.FileName);
    mol.Transform(Simmetric(FindCenter(mol)));
    FDispMan.Attach(mol);
    FDispMan.Render;
    end;
end;


procedure TCmMainForm.TestGrids;

var
  targetrads,proberads:TFloats;
  targetcoords,probecoords:TCoords;
  targetgrid,probegrid:TDockingGrid;
  domain:TDockDomain;
  tick1,tick2:DWORD;
  models:TModelManager;
  f:Integer;
  target,probe:TMolecule;
  MaxIters:Integer;

begin   {
  //target:=FMolecules.LoadLayer('C:\My Documents\programs\fpc\oclibrary\testfiles\DX-A.pdb');
  //probe:=FMolecules.LoadLayer('C:\My Documents\programs\fpc\oclibrary\testfiles\DX-B.pdb');
  target:=FMolecules.LoadLayer('C:\My Documents\programs\fpc\oclibrary\testfiles\1JMJ-A.pdb');
  probe:=FMolecules.LoadLayer('C:\My Documents\programs\fpc\oclibrary\testfiles\2CN0_HL.pdb');


  target.Transform(Simmetric(FindCenter(target)));
  probe.Transform(Simmetric(FindCenter(probe)));

  targetrads:=Add(ListRadii(target),1.4);
  targetcoords:=ListCoords(target);

  proberads:=Add(ListRadii(probe),1.4);
  probecoords:=ListCoords(probe);


  models:=TModelManager.Create(100,300);
  models.GridScale:=1;
  Writeln('Workin...');
  targetgrid:=TDockingGrid.Create(1);
  targetgrid.BuildFromSpheres(targetcoords,targetrads);


  MaxIters:=1;

  for f:=1 to MaxIters do
    begin
    tick1:=GetTickCount;
    probegrid:=TDockingGrid.Create(1);
    probegrid.BuildFromSpheres(probecoords,proberads);


    domain:=TDockDomain.Create(targetgrid,probegrid,0);
    domain.MinimumOverlap:=models.MinimumOverlap;
    domain.AssignModelManager(@models.AddModel);
    domain.RemoveCores:=True;
    domain.BuildInitialDomain;
    domain.Score;
      writeLn(models.Models[0].Score,' (',models.Models[0].TransVec[0],',',
              models.Models[0].TransVec[1],',',models.Models[0].TransVec[2],')');
    if f<MaxIters then
      begin
      domain.free;
      probegrid.free;
      end;
    end;
  tick2:=GetTickCount;

  domain.CalcDomainStats;

  WriteLn(FloatToStrF((tick2-tick1)/1000,ffFixed,4,3));
  WriteLn(domain.Size,' cells');
  DebugLn(IntToStr(Length(targetcoords)),' atoms');
  {for f:=0 to High(models.Models) do
    with models.Models[f] do
      writeLn(Score,' (',TransVec[0],',',TransVec[1],',',TransVec[2],')');
  }

  FDispMan.Attach(target);
  FDispMan.Attach(domain.Grid,domain.TranslateToTarget,targetgrid.Resolution,RGBAColor(0.7,0,0,0.5));
  FDispMan.Render;
  FDisplay.Refresh; }
end;


procedure TCmMainForm.TestRotationSampling;

var
  cs1,cs2,cs3:TCoords;
  rotations:TQuaternions;
  mindists:TFloats;
  rms:TFloat;
  f,g:Integer;
  sl:TStringList;
  mol:TMolecule;
begin
  SetLength(cs1,3);
  cs1[0]:=Coord(1,0,0);
  cs1[1]:=Coord(0,1,0);
  cs1[2]:=Coord(0,0,1);
  {cs1[3]:=Coord(-1,0,0);
  cs1[4]:=Coord(0,-1,0);
  cs1[5]:=Coord(0,0,-1);}

  mol:=FMolecules.LoadLayer('C:\My Documents\programs\fpc\oclibrary\testfiles\DX-A.pdb');
  mol.Transform(Simmetric(FindCenter(mol)));
  FDispMan.Attach(mol);
  cs1:=ListCoords(mol);

  rotations:=UniformSampleAngles(24);
  //rotations:=EulerSampleZYX(24);
  MinDists:=FilledFloats(Length(rotations),1000000);

  for f:=0 to High(Rotations)-1 do
    begin
    cs2:=Rotate(cs1,rotations[f]);
    {mol.Transform(rotations[f]);
    FDispMan.Render;
    FDisplay.Refresh;
    mol.Transform(Conjugated(rotations[f]));
    Application.ProcessMessages;}
    WriteLn(f,' of ',High(Rotations));
    for g:=f+1 to High(Rotations) do
      begin
      cs3:=Rotate(cs1,rotations[g]);
      rms:=StaticRMSD(cs2,cs3);
      if rms<MinDists[f] then
        MinDists[f]:=rms;
      if rms<MinDists[g] then
        MinDists[g]:=rms;              {mol.Transform(rotations[f]);
    FDispMan.Render;
    FDisplay.Refresh;
    mol.Transform(Conjugated(rotations[f]));
    Application.ProcessMessages;}
      end;
    end;
  sl:=TStringList.Create;
  for f:=0 to High(mindists) do
    sl.Add(FloatToStr(mindists[f]));
  sl.SaveToFile('C:\My Documents\programs\fpc\oclibrary\testfiles\anglesample.txt');
  sl.Free;
end;

procedure TCmMainForm.TestBogie;
var
  targetrads,proberads:TFloats;
  targetcoords,probecoords:TCoords;
  targetgrid,probegrid:TDockingGrid;
  domain:TDockDomain;
  tick1,tick2:DWORD;
  models:TModelManager;
  f:Integer;
  target,probe:TMolecule;
  dockman:TDockManager;
  rot:TQuaternion;
  tran:TCoord;
  continue:Boolean;

begin{
  //target:=FMolecules.LoadLayer('C:\My Documents\programs\fpc\oclibrary\testfiles\DX-A.pdb');
  //probe:=FMolecules.LoadLayer('C:\My Documents\programs\fpc\oclibrary\testfiles\DX-B.pdb');
  target:=FMolecules.LoadLayer('C:\My Documents\programs\fpc\oclibrary\testfiles\1JMJ-A.pdb');
  probe:=FMolecules.LoadLayer('C:\My Documents\programs\fpc\oclibrary\testfiles\2CN0_HL.pdb');

  target.Transform(Simmetric(FindCenter(target)));
  probe.Transform(Simmetric(FindCenter(probe)));
  FDispMan.Attach(target);
  FDispMan.Attach(probe);

  targetrads:=Add(ListRadii(target),1.35);
  targetcoords:=ListCoords(target);

  proberads:=Add(ListRadii(probe),1.35);
  probecoords:=ListCoords(probe);

  writeln(length(probecoords),' ',length(targetcoords));

  dockman:=TDockManager.Create(targetcoords,TargetRads,ProbeCoords,ProbeRads,1,5000,600);
  dockman.BuildTargetGrid;
  dockman.BuildRotations(24);
  continue:=True;
  tick1:=GetTickCount;
  while continue do
    begin
    continue:=not dockman.Dock(5);
    rot:=dockman.Models[0].ProbeRotation;
    tran:=dockman.Models[0].TransVec;
    //probe.Transform(rot);
    //probe.Transform(tran);
    tick2:=GetTickCount-tick1;
    WriteLn(dockman.CurrentRotation,' : ',
            FloatToStrF(tick2/1000,ffFixed,4,3),'s (',
            FloatToStrF(tick2/60000*dockman.TotalRotations/dockman.CurrentRotation,ffFixed,4,2),
            'm) Score:',dockman.Models[0].Score{,
            ' RMSD: ',FloatToStrF(RMSD(probecoords,ListCoords(probe)),ffFixed,4,3)});
    {FDispMan.Render;
    FDisplay.Refresh;
    probe.Transform(Simmetric(tran));
    probe.Transform(Conjugated(rot));}
    Application.ProcessMessages;
    end;}
end;

procedure TCmMainForm.OptimizeAddedRad;
var
  targetrads,proberads:TFloats;
  tmp,targetcoords,probecoords:TCoords;
  targetgrid,probegrid:TDockingGrid;
  dockdomain:TDockDomain;
  tick1,tick2:DWORD;
  models:TModelManager;
  f:Integer;
  target,probe:TMolecule;
  dockman:TDockManager;
  rot:TQuaternion;
  tran:TCoord;
  rad,dist:TFloat;
  sl:TStringList;

begin{
  sl:=TStringList.Create;

  target:=FMolecules.LoadLayer('C:\My Documents\programs\fpc\oclibrary\testfiles\DX-A.pdb');
  probe:=FMolecules.LoadLayer('C:\My Documents\programs\fpc\oclibrary\testfiles\DX-B.pdb');
  //target.Transform(Simmetric(FindCenter(target)));
  //probe.Transform(Simmetric(FindCenter(probe)));
  FDispMan.Attach(target);
  FDispMan.Attach(probe);

  targetcoords:=ListCoords(target);
  probecoords:=ListCoords(probe);
  writeln(length(probecoords),' ',length(targetcoords));

  for f:=0 to 50 do
    begin
    rad:=f/20;
    //rad:=1.3;
    targetrads:=Add(ListRadii(target),rad);
    proberads:=Add(ListRadii(probe),rad);

    dockman:=TDockManager.Create(targetcoords,TargetRads,ProbeCoords,ProbeRads,1,100,300);
    dockman.BuildTargetGrid;
    dockman.BuildEulerRotations(2);
//    dockman.BuildNoRotations;
    dockman.Dock(1);
    rot:=dockman.Models[0].ProbeRotation;
    tran:=dockman.Models[0].TransVec;
    {
    FDispMan.Attach(dockman.TargetGrid.Core.Grid,Simmetric(dockman.TargetGrid.TransVec),dockman.TargetGrid.Resolution,RGBAColor(0.7,0,0,1));
    FDispMan.Attach(dockman.ProbeGrid.Core.Grid,Simmetric(Add(tran,dockman.TargetGrid.TransVec)),dockman.TargetGrid.Resolution,RGBAColor(0,1,0,1));

    probe.Transform(rot);
    probe.Transform(Subtract(dockman.ProbeGrid.TransVec,Add(tran,dockman.TargetGrid.TransVec)));
    }
    tmp:=ListCoords(probe);


    tmp:=Rotate(tmp,rot);
    tmp:=Add(tmp,tran);;
    dist:=RMSD(tmp,probecoords);
    Application.ProcessMessages;
    WriteLn(rad,' ',dist,' ',dockman.Models[0].Score);
    sl.Add(FloatToStr(rad)+#9+FloatToStr(dist)+#9+IntToStr(dockman.Models[0].Score));

    {dockdomain:=TDockDomain.Create(dockman.TargetGrid,dockman.ProbeGrid);
    dockdomain.MinimumOverlap:=300;
    dockdomain.RemoveCores:=True;
    dockdomain.BuildInitialDomain;
    FDispMan.Attach(dockdomain.Grid,dockdomain.TranslateToTarget,dockman.TargetGrid.Resolution,RGBAColor(0.7,0,0,0.5));


    FDispMan.Render;
    FDisplay.Refresh;}


    dockman.Free;
    end;

  sl.SaveToFile('C:\My Documents\programs\fpc\oclibrary\testfiles\rmsds.txt');
  sl.Free;}
end;

procedure TCmMainForm.TestConstraints;

var
  targetrads,proberads:TFloats;
  targetcoords,probecoords:TCoords;
  targetgrid,probegrid:TDockingGrid;
  domain:TDockDomain;
  tick1,tick2:DWORD;
  models:TModelManager;
  f:Integer;
  target,probe:TMolecule;
  MaxIters:Integer;
  tcpoints,pcpoints:TCoords;
  rot:TQuaternion;
  tran:TCoord;
begin
  {target:=FMolecules.LoadLayer('C:\My Documents\programs\fpc\oclibrary\testfiles\DX-A.pdb');
  probe:=FMolecules.LoadLayer('C:\My Documents\programs\fpc\oclibrary\testfiles\DX-B.pdb');}
  target:=FMolecules.LoadLayer('C:\My Documents\programs\fpc\oclibrary\testfiles\t1\7cei-A.pdb');
  probe:=FMolecules.LoadLayer('C:\My Documents\programs\fpc\oclibrary\testfiles\t1\7cei-B.pdb');

  rot:=Quaternion(-0.080,-0.437,-0.441,-0.780);
  probe.Transform(rot);


  {target:=FMolecules.LoadLayer('C:\My Documents\programs\fpc\oclibrary\testfiles\1JMJ-A.pdb');
  probe:=FMolecules.LoadLayer('C:\My Documents\programs\fpc\oclibrary\testfiles\2CN0_HL.pdb');}

  {target.Transform(Simmetric(FindCenter(target)));
  probe.Transform(Simmetric(FindCenter(probe)));}

  targetrads:=Add(ListRadii(target),1.35);
  targetcoords:=ListCoords(target);

  proberads:=Add(ListRadii(probe),1.35);
  probecoords:=ListCoords(probe);

  tcpoints:=ListCoords(target.GetGroup(0).GetGroupById(55));
  pcpoints:=ListCoords(probe.GetGroup(0).GetGroupById(538));
  {for f:=0 to High(pcpoints) do
    writeln(pcpoints[f,0],'+',pcpoints[f,1],'+',pcpoints[f,2]);
  }
  models:=TModelManager.Create(10,0,nil);

  models.GridScale:=1;
  Writeln('Workin...');
  targetgrid:=TDockingGrid.Create(1);
  targetgrid.BuildFromSpheres(targetcoords,targetrads);


  MaxIters:=1;
  tick1:=GetTickCount;

  for f:=1 to MaxIters do
    begin
    probegrid:=TDockingGrid.Create(1);
    probegrid.BuildFromSpheres(probecoords,proberads);
    models.ProbeGridTransVec:=Subtract(ProbeGrid.TransVec,TargetGrid.TransVec);



    tick1:=GetTickCount;
    domain:=TDockDomain.Create(targetgrid,probegrid,0);
    //domain.ConstraintManager.LinearDistanceConstraint(tcpoints,pcpoints,10);
    domain.ConstraintManager.EuclideanDistanceConstraint(tcpoints,pcpoints,10);

    domain.MinimumOverlap:=models.MinimumOverlap;
    domain.AssignModelManager(@models.AddModel);
    domain.RemoveCores:=True;

    domain.BuildInitialDomain;
    //domain.Score;
    if f<MaxIters then
      begin
      domain.free;
      probegrid.free;
      end;
    end;
  tick2:=GetTickCount;
  WriteLn('Time: ',FloatToStrF((tick2-tick1)/1000,ffFixed,4,3));


  writeLn(models.Models[0].OverlapScore,' (',models.Models[0].TransVec[0],',',
         models.Models[0].TransVec[1],',',models.Models[0].TransVec[2],')');

  domain.CalcDomainStats;

  WriteLn(domain.Size,' cells');
  DebugLn(IntToStr(Length(targetcoords)),' atoms');
  for f:=0 to High(models.Models) do
    with models.Models[f] do
      writeLn(OverlapScore,' (',TransVec[0],',',TransVec[1],',',TransVec[2],')');


  FDispMan.Attach(target.GetGroup(0).GetGroupById(55));
  FDispMan.Attach(domain.Grid,domain.TranslateToTarget,targetgrid.Resolution,RGBAColor(0.7,0,0,0.5));

  rot:=models.Models[0].Rotation;
  tran:=models.Models[0].TransVec;
  //tran:=Coord(3.732,26.433,10.850);

  probe.Transform(tran);
  SaveToPDB(target,'C:\My Documents\programs\fpc\oclibrary\testfiles\t1\7cei-AM0.pdb');
  SaveToPDB(probe,'C:\My Documents\programs\fpc\oclibrary\testfiles\t1\7cei-BM0.pdb');


  //FDispMan.Attach(probe{.GetGroup(0).GetGroupById(528)});


  FDispMan.Render;
  FDisplay.Refresh;
end;

procedure TCmMainForm.TestAngles;

var
  targetcoords,probecoords:TCoords;
  f,r,iters:Integer;
  target:TMolecule;
  MaxIters:Integer;
  tmp:TFloat;
  minrmsds:TFloats;
  rot:TQuaternion;
  rots:TQuaternions;
  tran:TCoord;
  ZStep:Integer;

begin
  target:=FMolecules.LoadLayer('C:\My Documents\programs\fpc\oclibrary\testfiles\DX-A.pdb');
  target.Transform(Simmetric(FindCenter(target)));

  targetcoords:=ListCoords(target);


  ZStep:=24;
  rots:=UniformSampleAngles(ZStep);

  MaxIters:=100;

  minrmsds:=FilledFloats(MaxIters,1000000);
  for f:=0 to High(minrmsds) do
    begin
    rot:=Quaternion(Random-0.5,Random-0.5,Random-0.5,Random-0.5);
    Normalize(rot);
    probecoords:=Rotate(ListCoords(target),rot);
    for r:=0 to High(rots) do
      begin
      tmp:=StaticRMSD(Rotate(probecoords,rots[r]),targetcoords);
      if tmp<minrmsds[f] then minrmsds[f]:=tmp;
      end;
    WriteLn(minrmsds[f]);
    end;
  WriteLn('Max: ',Max(minrmsds));
  target.Free;

end;

procedure TCmMainForm.InitChemera;

procedure TestIntersect;

var
  l1,l2: TGridLine;
  f:Integer;

begin
  SetLength(l1,3);
  l1[0,0]:=4;
  l1[0,1]:=8;
  l1[1,0]:=10;
  l1[1,1]:=12;
  l1[2,0]:=15;
  l1[2,1]:=18;

  SetLength(l2,2);
  l2[0,0]:=4;
  l2[0,1]:=25;
  l2[1,0]:=20;
  l2[1,1]:=21;

  l1:=Intersect(l1,l2);
  for f:=0 to High(l1) do
    writeln(l1[f,0],' ',l1[f,1]);

end;

begin
  FInitializing:=False;

  //Load configuration data
  LoadAtomData;
  LoadAAData;
  FDisplay:=TOpenGLForm.Create(Self);
  FDisplay.Show;
  FDispMan:=TDisplayManager.Create(FDisplay);

  //DEBUG
  //FDispMan.Test;

  FMolecules:=TPdbModelMan.Create(Config.MonomersPath);

  //TestIntersect;

end;

procedure TCmMainForm.FormActivate(Sender: TObject);
begin
  if FInitializing then InitChemera;
end;

end.

