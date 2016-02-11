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
  molutils, geomutils, basetypes, linegrids, base3ddisplay, dockdomains,
  bogie, povray,selections, surface, rotations;

type

  { TCmMainForm }

  TCmMainForm = class(TForm)
    MainMenu: TMainMenu;
    FileMn: TMenuItem;
    MenuItem1: TMenuItem;
    MenuItem10: TMenuItem;
    MenuItem11: TMenuItem;
    MenuItem12: TMenuItem;
    MenuItem13: TMenuItem;
    MenuItem2: TMenuItem;
    MenuItem3: TMenuItem;
    MenuItem4: TMenuItem;
    MenuItem5: TMenuItem;
    MenuItem6: TMenuItem;
    MenuItem7: TMenuItem;
    MenuItem8: TMenuItem;
    MenuItem9: TMenuItem;
    SavePOVmi: TMenuItem;
    OpenDlg: TOpenDialog;
    OpenMn: TMenuItem;
    procedure FormActivate(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure FormDropFiles(Sender: TObject; const FileNames: array of String);
    procedure MenuItem10Click(Sender: TObject);
    procedure MenuItem11Click(Sender: TObject);
    procedure MenuItem12Click(Sender: TObject);
    procedure MenuItem13Click(Sender: TObject);
    procedure MenuItem2Click(Sender: TObject);
    procedure MenuItem3Click(Sender: TObject);
    procedure MenuItem4Click(Sender: TObject);
    procedure MenuItem5Click(Sender: TObject);
    procedure MenuItem6Click(Sender: TObject);
    procedure MenuItem7Click(Sender: TObject);
    procedure MenuItem8Click(Sender: TObject);
    procedure SavePOVmiClick(Sender: TObject);
    procedure OpenMnClick(Sender: TObject);
  private
    { private declarations }
    FInitializing:Boolean;
    FDisplay:TOpenGlForm;
    FDispMan:TDisplayManager;
    FMolecules:TPdbModelMan;

    //Test functions
    procedure TestShapePoints;
    procedure TestGrids;
    procedure TestRotationSampling;
    procedure TestBogie;
    procedure OptimizeAddedRad;
    procedure TestConstraints;
    procedure TestSymmetryConstraints;
    procedure TestAngles;
    procedure TestShowRotation;

    procedure LoadFile(FileName:string);
    procedure InitChemera;
  public
    { public declarations }
  end; 

var
  CmMainForm: TCmMainForm;

implementation

{$R *.lfm}

//For testing only; remove
function OldUniformSampleAngles(ZSteps:Integer;Out Axes:TCoords):TQuaternions;

//Rotates around Z axis ZSteps, for each rotates pole to one point in sphere

var
  sphere:TCoords;
  f,g,ix:Integer;
  rot:TFloat;
  zaxis:TCoord;
begin
  if ZSteps=1 then
    begin
    SetLength(Result,1);
    Result[0]:=IdentityQuaternion;
    end
  else
    begin
    zaxis:=Coord(0,0,1);
    //get ZSteps^2 points uniformly distributed around a sphere
    sphere:=GoldenSpiralPoints(ZSteps*ZSteps div 2);
    Axes:=sphere;
    SetLength(Result,Length(sphere)*ZSteps);
    ix:=0;
    for f:=0 to ZSteps-1 do
      begin
      rot:=2*PI/ZSteps*f;
      for g:=0 to High(sphere) do
        begin
        //rotate around z axis by rot then rotate to point z axis to sphere point
        Result[ix]:=Multiply(RotationTo(zaxis,sphere[g]),RotationQuaternion(sphere[g],rot));
        Inc(ix);
        end;
      end;
    end;
end;


{ TCmMainForm }

procedure TCmMainForm.FormCreate(Sender: TObject);
begin
  FInitializing:=True;
end;

procedure TCmMainForm.FormDropFiles(Sender: TObject;
  const FileNames: array of String);

var f:Integer;

begin
  for f:=0 to High(FileNames) do
    LoadFile(FileNames[f]);
end;

procedure TCmMainForm.MenuItem10Click(Sender: TObject);

var f:Integer;
    atomsets:TAtomSettingsArray;

begin
  atomsets:=FDispMan.Settings.MolSettings[0].AtomSettings;
  for f:=0 to High(atomsets) do
    atomsets[f].IsVisible:=AtomIsAABackbone(atomsets[f].Atom);
  FDispMan.Render;
end;

procedure TCmMainForm.MenuItem11Click(Sender: TObject);
begin
  TestShowRotation;
end;

procedure TCmMainForm.MenuItem12Click(Sender: TObject);
begin
  TestSymmetryConstraints;
end;

procedure TCmMainForm.MenuItem13Click(Sender: TObject);
begin
  TestShapePoints;
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

procedure TCmMainForm.MenuItem8Click(Sender: TObject);

var
  Bmp:TBitMap;
  f:Integer;
  s:string;

begin
  for f:=1 to 72 do
    begin
    FDisplay.RotateMatrix(5,0);
    Bmp:=FDisplay.GetImage;
    s:=IntToStr(f);
    while Length(s)<4 do s:='0'+s;
    Bmp.SaveToFile('G:\test'+s+'.bmp');
    Bmp.Free;
    end;
end;

procedure TCmMainForm.SavePOVmiClick(Sender: TObject);

var
  c,r,a:Integer;
  colfrac,x,y,z:TFloat;
  pve:TPovRayExport;
  layer:TPDBModel;
  chain,res:TMolecule;
  atom:TAtom;
  mat:TGluMatrix;
begin
  mat:=FDisplay.ModelViewMat;
  pve:=TPovRayExport.Create();
  FDisplay.CameraPos(x,y,z);
  pve.AddCodeLine('#declare View_Transf = transform { matrix <');
  pve.AddCodeLine(Format('%f, %f, %f,',[mat[0],mat[1],mat[2]]));
  pve.AddCodeLine(Format('%f, %f, %f,',[mat[4],mat[5],mat[6]]));
  pve.AddCodeLine(Format('%f, %f, %f,',[mat[8],mat[9],mat[10]]));
  pve.AddCodeLine(Format('0, 0, %f>',[-FDisplay.ZDist]));
  pve.AddCodeLine('inverse}');

  pve.AddCodeLine(Format('#declare camera_pos = <%f,%f,%f>;',[x,y,z]));
  pve.AddCodeLine('#include "scene.pov"');
  //pve.AddCodeLine('camera {look_at <0,0,0> location camera_pos right -x*image_width/image_height sky <0,0,1>}');
  layer:=FMolecules.LayerByIx(0);
  for c:=0 to layer.ChainCount-1 do
    begin
    chain:=layer.GetChain(c);
    colfrac:=0.1+c/(layer.ChainCount+3);
    pve.AddCodeLine(Format('#declare BallPigment = pigment {color <%f, %f, %f>};',
                    [colfrac,colfrac,colfrac]));
    for r:=0 to chain.GroupCount-1 do
        begin
        res:=chain.GetGroup(r);
        for a:=0 to res.AtomCount-1 do
          begin
          atom:=res.GetAtom(a);
          pve.AddCodeLine(Format('sphere{ <%f, %f, %f>, %f texture { BallFinish BallPigment }}',
                    [atom.Coords[0],atom.Coords[1],atom.Coords[2],atom.Radius]));
          end;
        end;
    end;

  pve.SaveToFile('D:\test.pov');
end;

procedure TCmMainForm.OpenMnClick(Sender: TObject);

begin
  OpenDlg.Filter:='';//Pdb file|*.pdb;GZipped PDB file|*.gz';
  if OpenDlg.Execute then
    LoadFile(OpenDlg.FileName);
end;

procedure TCmMainForm.TestShapePoints;
//generates a set of points describing the shape of the protein

var
  mol:TMolecule;
  coords,selcoords,axes:TCoords;
  selatoms:TIntegers;
  atoms:TAtoms;
  f:Integer;
  rotman:TRotationManager;
  rads:TFloats;

begin
  //mol := FMolecules.LoadLayer('C:\My Documents\Research\papers\2015-WCB\data\cleandimers\1OAC-A-prob.pdb');
  mol := FMolecules.LoadLayer('C:\My Documents\Research\papers\2015-WCB\data\cleantrimers\1AA0.pdb');
  //mol := FMolecules.LoadLayer('C:\My Documents\Research\papers\2015-WCB\data\cleantrimers\1CA4.pdb');
  CenterMolecule(mol);
  coords:=mol.AllCoords;
  selatoms:=SpacedPoints(coords,40);
  atoms := mol.AllAtoms;
  writeln(length(atoms));
  for f:=0 to High(atoms) do atoms[f].Radius:=0.1;
  for f:=0 to high(selatoms) do atoms[selatoms[f]].Radius:=1;
  SetLength(selcoords,Length(selatoms));
  for f:=0 to High(selatoms) do selcoords[f]:=coords[selatoms[f]];
  FDispMan.Attach(mol);
  rotman:=TRotationManager.Create(selcoords);
  //rotman.InitializeAxes(2000);
  //rotman.ComputeAxisDistances;
  //rotman.SelectAxes(6);
  //rotman.SelectAxes(6);
  //rotman.SelectAxesByMaxmin(8,5000);
  //rotman.ForceSingleRotation(-1);
  //rotman.GenerateQuaternions(16,8);
  //rotman.GenerateQuaternions(rotman.AxialDistance,rotman.AxialDistance/2);

  //rotman.MaxDistQuaternions(3000,10);
  rotman.MaxDistQuaternions(3000,5,2);

  axes := Multiply(rotman.GetSelectedAxes,100);
  rads := FilledFloats(Length(axes),1);
  for f:=0 to High(axes) do
    //rads[f]:=rotman.MinDist(f);
    rads[f]:=rotman.RotationCount(f);
  WriteLn(Round(Sum(rads)));
  rads:=Multiply(rads,1/Max(rads));
  FDispMan.Attach(axes,rads,RGBAColor(0.2,0,0,1),RGBAColor(1,0,0,1),RGBAColor(0.2,0.2,0.2,1));

  rads:=FilledFloats(length(rotman.Rotations),0.2);
  SetLength(axes,length(rads));
  for f:=0 to High(axes) do
    axes[f]:=Rotate(selcoords[0],rotman.Rotations[f]);

  FDispMan.Attach(axes,rads,RGBAColor(0,0,0.2,1),RGBAColor(0,0,1,1),RGBAColor(0.2,0.2,0.2,1));
  FDispMan.Render;
  FDisplay.Refresh;
  rotman.Free;
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

{Tests all rotations of a sphere of 10 points against a set of random rotations
 Returns the maximum of the minimum distance found}




var
  cs1,cs2,cs3,tmp:TCoords;
  tc:TCoord;
  rms:TFloat;
  rotations:TQuaternions;
  rot:TQuaternion;
  maxrms,minrms:TFloat;
  f,g:Integer;
  MaxIters:Integer;

begin
  cs1:=Multiply(GoldenSpiralPoints(10),10);
  //rotations:=UniformSampleAngles(24,tmp);
  rotations:=OldUniformSampleAngles(24,tmp);

  MaxIters:=100000;
  maxrms:=0;
  for f:=1 to MaxIters do
    begin
    tc:=Coord(Random, Random, Random);

    //uniform sampling of quaternions
    //from http://planning.cs.uiuc.edu/node198.html
    rot:=Quaternion(Sqrt(1-tc[0])*Sin(2*Pi*tc[1]),
                    Sqrt(1-tc[0])*Cos(2*Pi*tc[1]),
                    Sqrt(tc[0])*Sin(2*Pi*tc[2]),
                    Sqrt(tc[0])*Cos(2*Pi*tc[2]));
    Normalize(rot);
    cs2:=Rotate(cs1,rot);
    minrms:=10000000;
    for g:=0 to High(rotations) do
      begin
      rms:=StaticRMSD(Rotate(cs1,rotations[g]),cs2);
      if rms<minrms then minrms:=rms;
      end;
    if f mod 100 = 0 then WriteLn(f,' ',MaxIters,' ',minrms,'-',maxrms);
    if minrms>maxrms then
      maxrms:=minrms;
    end;
  WriteLn('Max: ',maxrms);
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

procedure TCmMainForm.TestSymmetryConstraints;

var
  targetrads,proberads:TFloats;
  axis:TCoord;
  targetcoords,probecoords:TCoords;
  targetgrid,probegrid:TDockingGrid;
  domain:TDockDomain;
  tick1,tick2:DWORD;
  models:TModelManager;
  f:Integer;
  target,probe:TMolecule;
  MaxIters:Integer;
  tran:TCoord;
  rot:TQuaternion;
begin
  {target:=FMolecules.LoadLayer('C:\My Documents\programs\fpc\oclibrary\testfiles\DX-A.pdb');
  probe:=FMolecules.LoadLayer('C:\My Documents\programs\fpc\oclibrary\testfiles\DX-B.pdb');}
  target:=FMolecules.LoadLayer('C:\My Documents\programs\fpc\oclibrary\testfiles\1RPO-A-targ.pdb');
  probe:=FMolecules.LoadLayer('C:\My Documents\programs\fpc\oclibrary\testfiles\1RPO-A-prob.pdb');

  axis:=Coord(0,0,1);
  Normalize(axis);

  targetrads:=Add(ListRadii(target),1.35);
  targetcoords:=ListCoords(target);

  proberads:=Add(ListRadii(probe),1.35);
  probecoords:=Rotate(ListCoords(probe),RotationQuaternion(axis,Pi));

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
    domain.ConstraintManager.PlaneConstraint(Coord(0,0,0),axis,3);

    domain.MinimumOverlap:=models.MinimumOverlap;
    domain.AssignModelManager(@models.AddModel);
    domain.RemoveCores:=True;

    domain.BuildInitialDomain;
    domain.Score;
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


  FDispMan.Attach(target);
  FDispMan.Attach(domain.Grid,domain.TranslateToTarget,targetgrid.Resolution,RGBAColor(0.7,0,0,0.5));
  //FDispMan.Attach(probegrid.Surf.Grid,Subtract(models.Models[0].TransVec,ProbeGrid.TransVec),probegrid.Resolution,RGBAColor(0.2,0.5,0.5,0.5));
  rot:=RotationQuaternion(axis,Pi);
  tran:=models.Models[0].TransVec;
  writeln(tran[0],' ',tran[1],' ',tran[2]);
  probe.Transform(rot);
  probe.Transform(tran);
  SaveToPDB(target,'C:\My Documents\programs\fpc\oclibrary\testfiles\t1\7cei-AM0.pdb');
  SaveToPDB(probe,'C:\My Documents\programs\fpc\oclibrary\testfiles\t1\7cei-BM0.pdb');


  FDispMan.Attach(probe);


  FDispMan.Render;
  FDisplay.Refresh;
end;

procedure TCmMainForm.TestAngles;

var
  targetcoords,probecoords,target,tmpcoords:TCoords;
  f,r,iters:Integer;
  //target:TMolecule;
  MaxIters:Integer;
  tmp:TFloat;
  minrmsds,rmsds:TFloats;
  rot:TQuaternion;
  rots:TQuaternions;
  tran:TCoord;
  ZStep:Integer;
  sl:TStringList;
  minixs,hist:TIntegers;
  s:string;

begin
  //target:=FMolecules.LoadLayer('C:\My Documents\programs\fpc\oclibrary\testfiles\DX-A.pdb');
  //target.Transform(Simmetric(FindCenter(target)));
  target:=Multiply(GoldenSpiralPoints(10),10);
  ZStep:=48;
  rots:=UniformSampleAngles(ZStep,tmpcoords);
  //rots:=OldUniformSampleAngles(ZStep,tmpcoords);
  minrmsds:=FilledFloats(Length(rots),1000000);
  minixs:=FilledInts(Length(rots),-1);
  sl:=TStringList.Create;
  for f:=0 to High(rots)-1 do
    begin
    targetcoords:=Rotate(target,rots[f]);
    hist:=FilledInts(40,0);
    for r:=f+1 to High(rots) do
      begin
      probecoords:=Rotate(target,rots[r]);
      tmp:=StaticRMSD(probecoords,targetcoords);
      if round(tmp*10)<=High(hist) then
        Inc(hist[round(tmp*10)]);
      if minrmsds[f]>tmp then
        begin
        minrmsds[f]:=tmp;
        minixs[f]:=r;
        end;
      if minrmsds[r]>tmp then
        begin
        minrmsds[r]:=tmp;
        minixs[r]:=f;
        end;
      end;
    s:='';
    for r:=0 to High(hist) do
      s:=s+IntToStr(hist[r])+#9;
    sl.Add(s);
    if f mod 100 = 0 then
      Writeln(f,' of ',High(rots));
    end;
  WriteLn('Max: ',Max(minrmsds));
  //target.Free;
  sl.SaveToFile('C:\My Documents\programs\fpc\oclibrary\testfiles\anglehistogram.txt');
  sl.Clear;
  for f:=0 to High(minrmsds) do
    sl.add(FloatToStr(minrmsds[f])+#9+IntToStr(minixs[f]));
  sl.SaveToFile('C:\My Documents\programs\fpc\oclibrary\testfiles\anglepairmindist.txt');
  sl.free;

end;

procedure TCmMainForm.TestShowRotation;
begin

end;

procedure TCmMainForm.LoadFile(FileName: string);

var mol:TMolecule;

begin
    mol:=FMolecules.LoadLayer(FileName);
    FMolecules.LayerByIx(0).DeleteWater;
    mol.Transform(Simmetric(FindCenter(mol)));
    FDispMan.Attach(mol);
    FDispMan.Render;
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
  FDisplay.OnDropFiles:=@FormDropFiles;
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

