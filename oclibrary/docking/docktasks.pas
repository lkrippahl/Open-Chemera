{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 2014 02 18
Purpose:
  Handles jobs, results and input/output of orders and result files for docking
  runs (bigger and dockprep)
Requirements:
Revisions:
To do:
    TDockOrdersManager should handle different docking runs per file, but currently
    the way it deals with probe and target structures doesn«t allow this in
    general.

*******************************************************************************}
unit docktasks;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes,geomutils, molecules, pdbmolecules, molutils,
  stringutils, oclconfiguration, quicksort, laz2_DOM, laz2_XMLRead, laz2_XMLWrite;

const
  //current orders version
  OrdersVersion=1;

  //constraint types
  CNull=-1;           //for errors in reading
  CUnconstrained=0;
  CLinearDistance=1;
  CEuclideanDistance=2;
  CContactCount=3;

  //Score types
  ScoreNone=0;
  ScoreRMSD=1;
  ScoreContacts=1;

type
  TDockModel=record
    TransVec:TCoord;
    Rotation:TQuaternion;
    OverlapScore:Integer;
  end;

  TDockModels=array of TDockModel;
    //the ID of each model is its index

  TScoreResult=record
    ScoreID:string;
    ScoreVals:TFloats;
      //same index as corresponding models
      //the id for each is the index (0, ... )
    ComputedCount:Integer;
  end;

  TScoreResults=array of TScoreResult;

  TConstraintDef=record
    Name:string;
    ProbePoints,TargetPoints:TCoords;
    Distance:TFloat;
    case ConstraintType:Integer of
      //CLinearDistance:
      //CEuclideanDistance:
      CContactCount: ( MinContacts:Integer;  MaxContacts:Integer);
  end;

  TConstraintDefs=array of TConstraintDef;

  TConstraintSet=record
    Name:string;
    Constraints:TConstraintDefs;
      //Constraints may be empty for unconstrained docking
    NumModels,MinOverlap:Integer;
    DockModels:TDockModels;
    ScoreResults:TScoreResults;
  end;


    //these constraints are added in series to the same constraint manager for
    //conjunction of constraints

  TConstraintSets=array of TConstraintSet;
    //each element is an independent set of constraints, leading to an independent
    //docking run in paralel (to use the same grids and rotations)

  //Score definitions

  { TScoreDef }

  TScoreDef=class
    private
      FName:string;
    public
      property Name:string read FName;
      constructor Create(AName:string);virtual;
      function ScoreType:Integer;virtual;
  end;

  { TRMSDScoreDef }

  TRMSDScoreDef=class(TScoreDef)
    private
      FPredictedCoords,FActualCoords:TCoords;
      FFirstProbeCoord:Integer;
    public
      FitTarget,FitProbe,ScoreTarget,ScoreProbe:Boolean;
      constructor Create(AName:string);override;
      property Predicted:TCoords read FPredictedCoords;
      property Actual:TCoords read FActualCoords;
      property ProbeStart:Integer read FFirstProbeCoord;
      procedure SetCoords(Preds,Actuals:TCoords;FirstProbe:Integer);
      function ScoreType:Integer;override;
  end;

  { TContactScoreDef }

  TContactScoreDef=class(TScoreDef)
    private
      FTargetGroups,FProbeGroups:TCoordGroups;
      FTargetIndexes,FProbeIndexes:TIntegers;
      FContactScores:TFloats;
    public
      constructor Create(AName:string);override;
      property TargetGroups:TCoordGroups read FTargetGroups;
      property ProbeGroups:TCoordGroups read FProbeGroups;
      property TargetIndexes:TIntegers read FTargetIndexes;
      property ProbeIndexes:TIntegers read FProbeIndexes;
      property ContactScores:TFloats read FContactScores;
      procedure AddTargetGroup(Group:TCoords);
      procedure AddProbeGroup(Group:TCoords);
      procedure AddContact(TargetIx,ProbeIx:Integer;Score:TFloat);
      function ScoreType:Integer;override;
  end;


  TScoreDefs=array of TScoreDef;

  TDockStats=record
    DigitizationTickCount,
    ConstraintsTickCount,
    DomainTickCount,
    ScoringTickCount,
    TotalTickCount,
    TestedModelsCount,
    InsertedModelsCount:Integer;
  end;

  TDockRun=record
    TargetFile,ProbeFile:string;
    Resolution,AddedRadius:TFloat;
    NumAxisSteps,
    SecondsBetweenSaves:Integer;
    CompletedAngles,TotalAngles:Integer;
    ConstraintSets:TConstraintSets;
    ScoreDefs:TScoreDefs;
    Stats:TDockStats;
  end;

  TDockRuns=array of TDockRun;

  { TDockOrdersManager }

  TDockOrdersManager=class
    private
      FTarget,FProbe:TMolecule;
      FModelMan:TPDBModelMan;
      FDockRuns:TDockRuns;
      FCurrentDock:Integer;
      FJobFile:string;
      function GetLayer(LayerName:string):TPDBModel;
    public
      property Target:TMolecule read FTarget;
      property Probe:TMolecule read FProbe;
      constructor Create(AJobFile:string);
      procedure Free;
      procedure NewDockRun;
      procedure PreparePartners(TargetFile,ProbeFile,TargetChains,ProbeChains:string;
                                Center:Boolean;RandomizeProbe:Boolean);overload;
      procedure PreparePartners(JobIndex:Integer=0);overload;
        //prepares the docking partners for that job
      procedure SetPDBFiles(TargetFile,ProbeFile:string);
      procedure AddResidueContacts(ContactFile:string;MaxContacts:Integer=0);
        //text file containing one contact per line
        //each line has target chain space target res ID space probe chain space probe res ID space distance
        //each line is an independent distance constraint between the atoms in those residues

      procedure AddRMSDScore(TemplateFile:string);
      { TODO : only works if same as target and probe files
        Must implement chain mapping with sequence comparisons to generalize
        Currently is fixed on fitting target and measuing probe}

      procedure AddInterfaceRMSDScore(TemplateFile:string;Dist:TFloat);
      { TODO : only works if same as target and probe files
        Must implement chain mapping with sequence comparisons to generalize }

      procedure AddNullConstraintSet(NumModels,MinOverlap:Integer);
        //for unconstrained docking
      procedure LoadJobFile;
      procedure SaveOrders(Append:Boolean);
      procedure DeleteScores(Id: string='');


  end;

  procedure SaveOrders(const Orders:TDockRuns;FileName:string;Append:Boolean=False);
  function LoadOrders(FileName:string):TDockRuns;
  function ReadAsTable(FileName:string):TSimpleStrings;
  procedure ZeroStats(var Stats:TDockStats);
  function GetModel(const Run:TDockRun;ConstSetIx,ModelId:Integer):TMolecule;
  procedure FreeOrders(var Orders:TDockRuns);
  procedure SortedModelList(const DockRun:TDockRun; out ConstSetIxs,ModelIxs:TIntegers);
  function CountModels(const DockRun:TDockRun):Integer;
  function ComputedScores(const ConstSet:TConstraintSet):TSimpleStrings;
  function ScoreVal(const ConstSet:TConstraintSet;ModelIx:Integer;
                              ScoreId:string;Default:TFloat=0):TFloat;

implementation

procedure SaveOrders(const Orders:TDockRuns;FileName:string;Append:Boolean=False);

var
 doc: TXMLDocument;
 root, node, run: TDOMNode;
 DockRun:TDockRun;


 procedure AddTextNode(Name,Val:string;Parent:TDOMNode);

 var tmp,txt:TDOMNode;

 begin
    tmp:=Doc.CreateElement(Name);
    txt:=Doc.CreateTextNode(Val);
    tmp.AppendChild(txt);
    Parent.AppendChild(tmp);
 end;

 function AddDockParameters(Parent:TDOMNode):TDOMNode;
 begin
    Result:=Doc.CreateElement('DockParameters');
    AddTextNode('TargetFile',DockRun.TargetFile,Result);
    AddTextNode('ProbeFile',DockRun.ProbeFile,Result);
    AddTextNode('Resolution',FloatToStrF(DockRun.Resolution,ffFixed,0,2),Result);
    AddTextNode('AddedRadius',FloatToStrF(DockRun.AddedRadius,ffFixed,0,2),Result);
    AddTextNode('NumAxisSteps',IntToStr(DockRun.NumAxisSteps),Result);
    AddTextNode('SecondsBetweenSaves',IntToStr(DockRun.SecondsBetweenSaves),Result);
    AddTextNode('CompletedRotations',IntToStr(DockRun.CompletedAngles),Result);
    AddTextNode('TotalRotations',IntToStr(DockRun.TotalAngles),Result);
    Parent.AppendChild(Result);
 end;

 function CreatePointSet(Points:TCoords;PointSetName:string):TDOMNode;

 var
   point:TDOMNode;
   f:Integer;

 begin
    Result:=Doc.CreateElement(PointSetName);
    for f:=0 to High(Points) do
      begin
      point:=Doc.CreateElement('Point');
      TDOMElement(point).SetAttribute('X',FloatToStrF(Points[f,0],ffFixed,0,3));
      TDOMElement(point).SetAttribute('Y',FloatToStrF(Points[f,1],ffFixed,0,3));
      TDOMElement(point).SetAttribute('Z',FloatToStrF(Points[f,2],ffFixed,0,3));
      Result.AppendChild(point);
      end;

 end;

 procedure AddContactScore(Parent:TDOMNode;Score:TContactScoreDef);

 var
   tmp,contact:TDOMNode;
   f:Integer;

 begin
    tmp:=Doc.CreateElement('ContactScore');
    AddTextNode('ScoreID',Score.Name,tmp);
    for f:=0 to High(Score.TargetGroups) do
      tmp.AppendChild(CreatePointSet(Score.TargetGroups[f],'TargetGroup'));
    for f:=0 to High(Score.ProbeGroups) do
      tmp.AppendChild(CreatePointSet(Score.ProbeGroups[f],'ProbeGroup'));
    for f:=0 to High(Score.TargetIndexes) do
      begin
      contact:=Doc.CreateElement('Contact');
      TDOMElement(contact).SetAttribute('TargetIx',IntToStr(Score.TargetIndexes[f]));
      TDOMElement(contact).SetAttribute('ProbeIx',IntToStr(Score.ProbeIndexes[f]));
      TDOMElement(contact).SetAttribute('ScoreVal',FloatToStrF(Score.ContactScores[f],ffFixed,0,3));
      tmp.AppendChild(contact);
      end;
    Parent.AppendChild(tmp);
  end;

 function BooleanToStr(Val:Boolean):string;

 begin
    if Val then Result:='1'
    else Result:='0';
 end;

 procedure AddRMSDScore(Parent:TDOMNode;Score:TRMSDScoreDef);

 var
   tmp,pointset:TDOMNode;
   f:Integer;

 begin
    tmp:=Doc.CreateElement('RMSDScore');
    AddTextNode('ScoreID',Score.Name,tmp);
    tmp.AppendChild(CreatePointSet(Score.Predicted,'Predicted'));
    tmp.AppendChild(CreatePointSet(Score.Actual,'Template'));
    TDOMElement(tmp).SetAttribute('ScoreTarget',BooleanToStr(Score.ScoreTarget));
    TDOMElement(tmp).SetAttribute('ScoreProbe',BooleanToStr(Score.ScoreProbe));
    TDOMElement(tmp).SetAttribute('FitTarget',BooleanToStr(Score.FitTarget));
    TDOMElement(tmp).SetAttribute('FitProbe',BooleanToStr(Score.FitProbe));
    AddTextNode('FirstProbePoint',IntToStr(Score.ProbeStart),tmp);
    Parent.AppendChild(tmp);
  end;

  procedure AddCoordinate(Parent:TDOMNode;CoordName:string;Coor:TCoord);

  var tmp:TDOMElement;
  begin
    tmp:=Doc.CreateElement(CoordName);
    tmp.SetAttribute('X',FloatToStrF(Coor[0],ffFixed,0,3));
    tmp.SetAttribute('Y',FloatToStrF(Coor[1],ffFixed,0,3));
    tmp.SetAttribute('Z',FloatToStrF(Coor[2],ffFixed,0,3));
    Parent.AppendChild(tmp);
  end;

  procedure AddQuaternion(Parent:TDOMNode;QuatName:string;Quat:TQuaternion);

  var tmp:TDOMElement;
  begin
    tmp:=Doc.CreateElement(QuatName);
    tmp.SetAttribute('r',FloatToStrF(Quat[0],ffFixed,0,3));
    tmp.SetAttribute('i',FloatToStrF(Quat[1],ffFixed,0,3));
    tmp.SetAttribute('j',FloatToStrF(Quat[2],ffFixed,0,3));
    tmp.SetAttribute('k',FloatToStrF(Quat[3],ffFixed,0,3));
    Parent.AppendChild(tmp);
  end;


  procedure AddConstraintSet(Parent:TDOMNode;ConstSet:TConstraintSet);

  var
    setnode,constnode,modelset,model,scoreset,score:TDOMNode;
    pointset,pair,point:TDOMNode;
    f,g:Integer;

  begin
    setnode:=Doc.CreateElement('JointConstraintSet');
    TDOMElement(setnode).SetAttribute('Name',ConstSet.Name);
    TDOMElement(setnode).SetAttribute('NumModels',IntToStr(ConstSet.NumModels));
    TDOMElement(setnode).SetAttribute('MinOverlap',IntToStr(ConstSet.MinOverlap));
    for f:=0 to High(ConstSet.Constraints) do
      case ConstSet.Constraints[f].ConstraintType of
        CEuclideanDistance:
          begin
          constnode:=Doc.CreateElement('DistanceConstraint');
          TDOMElement(constnode).SetAttribute('Name',ConstSet.Constraints[f].Name);
          constnode.AppendChild(CreatePointSet(ConstSet.Constraints[f].TargetPoints,'TargetPoints'));
          constnode.AppendChild(CreatePointSet(ConstSet.Constraints[f].ProbePoints,'ProbePoints'));
          AddTextNode('Distance',FloatToStrF(ConstSet.Constraints[f].Distance,ffFixed,0,3),constnode);
          setnode.AppendChild(constnode);
          end;
      end;
    if ConstSet.DockModels<>nil then
      begin
      modelset:=Doc.CreateElement('ComputedModels');
      for f:=0 to High(ConstSet.DockModels) do
        begin
        model:=Doc.CreateElement('Model');
        AddTextNode('Contact',IntToStr(ConstSet.DockModels[f].OverlapScore),model);
        AddCoordinate(model,'Translation',ConstSet.DockModels[f].TransVec);
        AddQuaternion(model,'Rotation',ConstSet.DockModels[f].Rotation);
        modelset.AppendChild(model);
        end;
      setnode.AppendChild(modelset);
      end;

    for f:=0 to High(ConstSet.ScoreResults) do
      with ConstSet.ScoreResults[f] do
        begin
        scoreset:=Doc.CreateElement('ComputedScore');
        TDOMElement(scoreset).SetAttribute('Name',ScoreID);
        for g:=0 to High(ScoreVals) do
          begin
          score:=Doc.CreateElement('ModelScore');
          TDOMElement(score).SetAttribute('Value',FloatToStrF(ScoreVals[g],ffFixed,0,3));
          scoreset.AppendChild(score);
          end;
        setnode.AppendChild(scoreset);
        end;
    Parent.AppendChild(setnode);
  end;


  procedure AddDockStats(Parent:TDOMNode;const Stats:TDockStats);

  var statnode:TDOMNode;

  begin
    statnode:=Doc.CreateElement('Stats');
    with Stats do
      begin
      AddTextNode('DigitizationTickCount',IntToStr(DigitizationTickCount),statnode);
      AddTextNode('ConstraintsTickCount',IntToStr(ConstraintsTickCount),statnode);
      AddTextNode('DomainTickCount',IntToStr(DomainTickCount),statnode);
      AddTextNode('ScoringTickCount',IntToStr(ScoringTickCount),statnode);
      AddTextNode('TotalTickCount',IntToStr(TotalTickCount),statnode);
      AddTextNode('TestedModelsCount',IntToStr(TestedModelsCount),statnode);
      AddTextNode('InsertedModelsCount',IntToStr(InsertedModelsCount),statnode);
      end;
    Parent.AppendChild(statnode);
  end;


var
  f,g:Integer;

begin
  try
    if Append and FileExists(FileName) then
      begin
      ReadXMLFile(doc, FileName);
      root:= Doc.DocumentElement;
      end
    else
      begin
      doc:=TXMLDocument.Create;
      root:=Doc.CreateElement('Orders');
      Doc.Appendchild(root);
      end;
    TDOMElement(root).SetAttribute('Version',IntToStr(OrdersVersion));

    for g:=0 to High(Orders) do
      begin
      DockRun:=Orders[g];
      run:=Doc.CreateElement('DockRun');
      root.AppendChild(run);
      AddDockStats(run,DockRun.Stats);
      AddDockParameters(run);

      for f:=0 to High(DockRun.ConstraintSets) do
        AddConstraintSet(run,DockRun.ConstraintSets[f]);

      for f:=0 to High(DockRun.ScoreDefs) do
        case DockRun.ScoreDefs[f].ScoreType of
          ScoreRMSD:AddRMSDScore(run,TRMSDScoreDef(DockRun.ScoreDefs[f]));
        end;
      end;

    WriteXMLFile(Doc, FileName);
  finally
    Doc.Free;
  end;
end;

function LoadOrders(FileName: string): TDockRuns;

var
  orders:TDOMNode;
  doc:TXMLDocument;

  function ReadTextNode(Name:string;Parent:TDOMNode):string;

  begin
    //WriteLn(Name,':',Parent.NodeName);
    Result:=Parent.FindNode(Name).FirstChild.NodeValue;
    //WriteLn(Name,':',Result);
  end;

  function ReadBooleanAttribute(Node:TDOMNode;Name:string):Boolean;

  var s:string;

  begin
    s:=UpperCase(TDOMElement(Node).GetAttribute(Name));
    Result:=(s='TRUE') or (s='1');
  end;



 function ReadPointSet(PointSet:TDOMNode):TCoords;overload;

 var
   point:TDOMNode;
   f:Integer;

 begin
    SetLength(Result,PointSet.ChildNodes.Count);
    for f:=0 to High(Result) do
      begin
      point:=pointset.ChildNodes.Item[f];
      Result[f,0]:=StrToFloat(TDOMElement(point).GetAttribute('X'));
      Result[f,1]:=StrToFloat(TDOMElement(point).GetAttribute('Y'));
      Result[f,2]:=StrToFloat(TDOMElement(point).GetAttribute('Z'));
      end;
 end;


 function ReadPointSet(Node:TDOMNode;PointSetName:string):TCoords;overload;

 var
   pointset:TDOMNode;

 begin
    pointset:=Node.FindNode(PointSetName);
    Result:=ReadPointSet(pointset);
 end;

  procedure ReadDockParameters(DockIx:Integer;Node:TDomNode);
  begin
    with Result[DockIx] do
      begin
      TargetFile:=ReadTextNode('TargetFile',Node);
      ProbeFile:=ReadTextNode('ProbeFile',Node);
      Resolution:=StrToFLoat(ReadTextNode('Resolution',Node));
      AddedRadius:=StrToFLoat(ReadTextNode('AddedRadius',Node));
      NumAxisSteps:=StrToInt(ReadTextNode('NumAxisSteps',Node));
      SecondsBetweenSaves:=StrToInt(ReadTextNode('SecondsBetweenSaves',Node));
      CompletedAngles:=StrToInt(ReadTextNode('CompletedRotations',Node));
      TotalAngles:=StrToInt(ReadTextNode('TotalRotations',Node));
      end
  end;

  function ReadCoordinate(Parent:TDOMNode;CoordName:string):TCoord;

  var tmp:TDOMNode;

  begin
    tmp:=Parent.FindNode(CoordName);
    Result[0]:=StrToFloat(TDOMElement(tmp).GetAttribute('X'));
    Result[1]:=StrToFloat(TDOMElement(tmp).GetAttribute('Y'));
    Result[2]:=StrToFloat(TDOMElement(tmp).GetAttribute('Z'));
  end;

  function ReadQuaternion(Parent:TDOMNode;QuatName:string):TQuaternion;

  var tmp:TDOMNode;

  begin
    tmp:=Parent.FindNode(QuatName);
    Result[0]:=StrToFloat(TDOMElement(tmp).GetAttribute('r'));
    Result[1]:=StrToFloat(TDOMElement(tmp).GetAttribute('i'));
    Result[2]:=StrToFloat(TDOMElement(tmp).GetAttribute('j'));
    Result[3]:=StrToFloat(TDOMElement(tmp).GetAttribute('k'));
  end;



  procedure AddConstraintSet(DockIx:Integer;Node:TDomNode);

  var
    f,ix,g:Integer;
    child,model:TDOMNode;
    nodename:string;

  begin
    SetLength(Result[DockIx].ConstraintSets,Length(Result[DockIx].ConstraintSets)+1);
    with Result[DockIx].ConstraintSets[High(Result[DockIx].ConstraintSets)] do
      begin
      Constraints:=nil;
      DockModels:=nil;
      ScoreResults:=nil;
      Name:=TDOMElement(Node).GetAttribute('Name');
      NumModels:=StrToInt(TDOMElement(Node).GetAttribute('NumModels'));
      MinOverlap:=StrToInt(TDOMElement(Node).GetAttribute('MinOverlap'));
      for f:=0 to Node.ChildNodes.Count-1 do
        begin
        child:=Node.ChildNodes.Item[f];
        nodename:=UpperCase(child.NodeName);
        if nodename='DISTANCECONSTRAINT' then
          begin
          ix:=Length(Constraints);
          SetLength(Constraints,ix+1);
          Constraints[ix].ConstraintType:=CEuclideanDistance;
          Constraints[ix].Name:=TDOMElement(child).GetAttribute('Name');
          Constraints[ix].TargetPoints:=ReadPointSet(child,'TargetPoints');
          Constraints[ix].ProbePoints:=ReadPointSet(child,'ProbePoints');
          Constraints[ix].Distance:=StrToFloat(ReadTextNode('Distance',child));
          end
        else if nodename='COMPUTEDMODELS' then
          begin
          SetLength(DockModels,child.ChildNodes.Count);
          for ix:=0 to High(DockModels) do
             begin
             model:=child.ChildNodes.Item[ix];
             DockModels[ix].OverlapScore:=StrToInt(ReadTextNode('Contact',model));
             DockModels[ix].TransVec:=ReadCoordinate(model,'Translation');
             DockModels[ix].Rotation:=ReadQuaternion(model,'Rotation');
             end;
          end
        else if nodename='COMPUTEDSCORE' then
          begin
          SetLength(ScoreResults,Length(ScoreResults)+1);
          ix:=High(ScoreResults);
          ScoreResults[ix].ScoreId:=TDOMElement(child).GetAttribute('Name');
          SetLength(ScoreResults[ix].ScoreVals,child.ChildNodes.Count);
          for g:=0 to High(ScoreResults[ix].ScoreVals) do
             ScoreResults[ix].ScoreVals[g]:=
              StrToFLoat(TDOMElement(child.ChildNodes.Item[g]).GetAttribute('Value'));
          end;
        end;
      end;
  end;

 procedure AddContactScore(JobIx:Integer; Node:TDOMNode);

 var
  tmp,child:TDOMNode;
  childname:string;
  f,ix:Integer;
  score:TContactScoreDef;

 begin
   with Result[JobIx] do
     begin
     ix:=Length(ScoreDefs);
     SetLength(ScoreDefs,ix+1);
     score:=TContactScoreDef.Create(ReadTextNode('ScoreID',Node));
     ScoreDefs[ix]:=score;
     for f:=0 to node.ChildNodes.count-1 do
        begin
        child:=node.ChildNodes.Item[f];
        childname:=UpperCase(child.NodeName);
        if childname='TARGETGROUP' then
          score.AddTargetGroup(ReadPointSet(child))
        else if childname='PROBEGROUP' then
          score.AddProbeGroup(ReadPointSet(child))
        else if childname='CONTACT' then
          score.AddContact(StrToInt(TDOMElement(child).GetAttribute('TargetIx')),
                           StrToInt(TDOMElement(child).GetAttribute('ProbeIx')),
                           StrToFloat(TDOMElement(child).GetAttribute('ScoreVal')));
        end;
     end;
  end;

  procedure AddRMSDScore(JobIx:Integer;Node:TDOMNode);

  var
   tmp,pointset:TDOMNode;
   f,ix:Integer;

  begin
    with Result[JobIx] do
      begin
      ix:=Length(ScoreDefs);
      SetLength(ScoreDefs,ix+1);
      ScoreDefs[ix]:=TRMSDScoreDef.Create(ReadTextNode('ScoreID',Node));
      TRMSDScoreDef(ScoreDefs[ix]).SetCoords(ReadPointSet(Node,'Predicted'),
                                             ReadPointSet(Node,'Template'),
                                StrToInt(ReadTextNode('FirstProbePoint',Node)));
      with TDOMElement(Node) do
        begin
        if HasAttribute('FitTarget') then
          TRMSDScoreDef(ScoreDefs[ix]).FitTarget:=ReadBooleanAttribute(Node,'FitTarget');
        if HasAttribute('FitProbe') then
          TRMSDScoreDef(ScoreDefs[ix]).FitProbe:=ReadBooleanAttribute(Node,'FitProbe');
        if HasAttribute('ScoreTarget') then
          TRMSDScoreDef(ScoreDefs[ix]).ScoreTarget:=ReadBooleanAttribute(Node,'ScoreTarget');
        if HasAttribute('ScoreProbe') then
          TRMSDScoreDef(ScoreDefs[ix]).ScoreProbe:=ReadBooleanAttribute(Node,'ScoreProbe');
        end;
      end;
  end;

  procedure ReadDockStats(Parent:TDOMNode;var Stats:TDockStats);

  var statnode:TDOMNode;

  begin
    statnode:=Parent.FindNode('Stats');
    if statnode=nil then
      ZeroStats(Stats)
    else with Stats do
      begin
      DigitizationTickCount:=StrToInt(ReadTextNode('DigitizationTickCount',statnode));
      ConstraintsTickCount:=StrToInt(ReadTextNode('ConstraintsTickCount',statnode));
      DomainTickCount:=StrToInt(ReadTextNode('DomainTickCount',statnode));
      ScoringTickCount:=StrToInt(ReadTextNode('ScoringTickCount',statnode));
      TotalTickCount:=StrToInt(ReadTextNode('TotalTickCount',statnode));
      TestedModelsCount:=StrToInt(ReadTextNode('TestedModelsCount',statnode));
      InsertedModelsCount:=StrToInt(ReadTextNode('InsertedModelsCount',statnode));
      end;
    Parent.AppendChild(statnode);
  end;



var
  f,g:Integer;
  node,currdock:TDOMNode;
  nodename:string;
begin
  ReadXMLFile(doc, FileName);
  orders:=doc.DocumentElement;
  if (orders<>nil) and (UpperCase(orders.NodeName)='ORDERS') then
    begin
    SetLength(Result,orders.ChildNodes.Count);
    for f:=0 to orders.ChildNodes.Count-1 do
      begin
      Result[f].ConstraintSets:=nil;
      Result[f].ScoreDefs:=nil;
      currdock:=orders.ChildNodes.Item[f];
      ReadDockStats(currdock,Result[f].Stats);
      for g:=0 to currdock.ChildNodes.count-1 do
        begin
        node:=currdock.ChildNodes.Item[g];
        nodename:=UpperCase(node.NodeName);
        if nodename='DOCKPARAMETERS' then
          ReadDockParameters(f,node)
        else if nodename='JOINTCONSTRAINTSET' then
          AddConstraintSet(f,node)
        else if nodename='RMSDSCORE' then
          AddRMSDScore(f,node);
        end;
      end;
    end;
  doc.Free;
end;

function ReadAsTable(FileName: string): TSimpleStrings;

var
  jobs:TDockRuns;
  f,g,h,i,nummodels:Integer;
  s,s2:string;

begin
  jobs:=LoadOrders(FileName);
  Result:=nil;
  AddToArray(FileName,Result);
  for f:=0 to High(jobs) do
    with jobs[f] do
      begin
      AddToArray('',Result);
      AddToArray(TargetFile+' to '+ProbeFile,Result);
      nummodels:=0;

      //headers
      s:='';
      s2:='';
      for g:=0 to High(ConstraintSets) do
        begin
        s:=s+ConstraintSets[g].Name+#9+#9;
        s2:=s2+'ID'+#9+'Overlap'+#9;
        nummodels:=Max(nummodels,Length(ConstraintSets[g].DockModels));
        for h:=0 to High(ConstraintSets[g].ScoreResults) do
          begin
          s:=s+#9;
          s2:=s2+ConstraintSets[g].ScoreResults[h].ScoreId+#9;
          end;
        end;
      AddToArray(s,Result);
      AddToArray(s2,Result);

      //values
      for g:=0 to nummodels-1 do
        begin
        s:='';
        for h:=0 to High(ConstraintSets) do
          begin
          if g<=High(ConstraintSets[h].DockModels) then
            s:=s+IntToStr(g)+#9+IntToStr(ConstraintSets[h].DockModels[g].OverlapScore)+#9
          else s:=s+#9+#9 ;
          for i:=0 to High(ConstraintSets[h].ScoreResults) do
            if g<=High(ConstraintSets[h].ScoreResults[i].ScoreVals) then
              s:=s+FloatToStrF(ConstraintSets[h].ScoreResults[i].ScoreVals[g],ffFixed,0,3)+#9
            else s:=s+#9
          end;
        AddToArray(s,Result);
        end;
      end;
end;

procedure ZeroStats(var Stats: TDockStats);
begin
  with Stats do
    begin
    DigitizationTickCount:=0;
    ConstraintsTickCount:=0;
    DomainTickCount:=0;
    ScoringTickCount:=0;
    TotalTickCount:=0;
    TestedModelsCount:=0;
    InsertedModelsCount:=0;
    end;
end;

function GetModel(const Run: TDockRun; ConstSetIx, ModelId: Integer): TMolecule;

var
  modelman:TPDBModelMan;
  f:Integer;
  rot:TQuaternion;
  transvec:TCoord;
  probe,target,chain:TMolecule;

begin
  modelman:=TPDBModelMan.Create(Config.MonomersPath);

  target:=modelman.LoadLayer(Run.TargetFile);
  probe:=modelman.LoadLayer(Run.ProbeFile);
  Result:=TMolecule.CopyFrom(target,nil);
  rot:=Run.ConstraintSets[ConstSetIx].DockModels[ModelId].Rotation;
  transvec:=Run.ConstraintSets[ConstSetIx].DockModels[ModelId].TransVec;

  for f:=0 to High(probe.Groups) do
    begin
    chain:=TMolecule.CopyFrom(probe.Groups[f],nil);
    chain.Transform(rot);
    chain.Transform(transvec);
    Result.AddGroup(chain);
    end;
  modelman.Free;
end;

procedure FreeOrders(var Orders: TDockRuns);

var f,g:Integer;

begin
  for f:=0 to High(Orders) do
    for g:=0 to High(Orders[f].ScoreDefs) do
      Orders[f].ScoreDefs[g].Free;
  Orders:=nil;
end;

procedure SortedModelList(const DockRun: TDockRun; out ConstSetIxs,
  ModelIxs: TIntegers);


var
  c,m,f,count,checkfrom,lastoverlap,overlap,ix:Integer;
  sortedixs,mixs,cixs:TIntegers;
  overlaps:TFloats;
  isrepeat:Boolean;

  function IsSame(Cix1,Mix1,Cix2,Mix2:Integer):Boolean;

  begin
    Result:=Distance(DockRun.ConstraintSets[Cix1].DockModels[Mix1].Rotation,
                     DockRun.ConstraintSets[Cix2].DockModels[Mix2].Rotation)<=TINY;
    if Result then
      Result:=Distance(DockRun.ConstraintSets[Cix1].DockModels[Mix1].TransVec,
                       DockRun.ConstraintSets[Cix2].DockModels[Mix2].TransVec)<=TINY;


  end;

begin
  count:=CountModels(DockRun);
  SetLength(ConstSetIxs,count);
  SetLength(ModelIxs,count);
  SetLength(cixs,count);
  SetLength(mixs,count);
  SetLength(overlaps,count);
  count:=0;

  //flatten list and sort by geometric overlap
  with DockRun do
    for c:=0 to High(ConstraintSets) do
      with ConstraintSets[c] do
        for m:=0 to High(DockModels) do
          begin
          cixs[count]:=c;
          mixs[count]:=m;
          overlaps[count]:=-DockModels[m].OverlapScore;
          Inc(count);
          end;
  sortedixs:=QSAscendingIndex(overlaps);

  //run through sorted list and check those with same overlap score
  //to skip repeats
  checkfrom:=Length(sortedixs);
  count:=0;
  lastoverlap:=-1;
  for m:=0 to High(sortedixs) do
    begin
    isrepeat:=False;
    ix:=sortedixs[m];
    overlap:=DockRun.ConstraintSets[cixs[ix]].DockModels[mixs[ix]].OverlapScore;
    //check all with same overlap value to test if repeated
    if overlap=lastoverlap then
      for f:=checkfrom to count-1 do
        if IsSame(cixs[ix],mixs[ix],ConstSetIxs[f],ModelIxs[f]) then
          begin
          isrepeat:=True;
          Break;
          end;
    if not isrepeat then
      begin
      if overlap <> lastoverlap then
        begin
        checkfrom:=count;
        lastoverlap:=overlap
        end;
      ConstSetIxs[count]:=cixs[ix];
      ModelIxs[count]:=mixs[ix];
      Inc(count);
      end;
    end;
  SetLength(ConstSetIxs,count);
  SetLength(ModelIxs,count);
end;

function CountModels(const DockRun:TDockRun): Integer;

var
  f:Integer;

begin
  Result:=0;
  for f:=0 to High(DockRun.ConstraintSets) do
    Result:=Result+Length(DockRun.ConstraintSets[f].DockModels);
end;

function ComputedScores(const ConstSet: TConstraintSet): TSimpleStrings;

var
  f:Integer;

begin
  SetLength(Result,Length(ConstSet.ScoreResults));
  for f:=0 to High(Result) do
    Result[f]:=ConstSet.ScoreResults[f].ScoreId;
end;

function ScoreVal(const ConstSet: TConstraintSet; ModelIx: Integer;
  ScoreId: string; Default: TFloat): TFloat;

var
  scoreix:Integer;

begin
  scoreix:=High(ConstSet.ScoreResults);
  while (scoreix>=0) and (ConstSet.ScoreResults[scoreix].ScoreId<>ScoreId) do
    Dec(scoreix);
  if scoreix>=0 then
    Result:=ConstSet.ScoreResults[scoreix].ScoreVals[ModelIx]
  else Result:=Default;
end;

{ TContactScoreDef }

constructor TContactScoreDef.Create(AName: string);
begin
  inherited Create(AName);
end;

procedure TContactScoreDef.AddTargetGroup(Group: TCoords);
begin
  SetLength(FTargetGroups,Length(FTargetGroups)+1);
  FTargetGroups[High(FTargetGroups)]:=Copy(Group,0,Length(Group));
end;

procedure TContactScoreDef.AddProbeGroup(Group: TCoords);
begin
  SetLength(FProbeGroups,Length(FProbeGroups)+1);
  FProbeGroups[High(FProbeGroups)]:=Copy(Group,0,Length(Group));
end;

procedure TContactScoreDef.AddContact(TargetIx, ProbeIx: Integer; Score: TFloat);
begin
  AddToArray(TargetIx,FTargetIndexes);
  AddToArray(ProbeIx,FProbeIndexes);
  AddToArray(Score,FContactScores);
end;

function TContactScoreDef.ScoreType: Integer;
begin
  Result:=inherited ScoreType;
end;

{ TRMSDScoreDef }

constructor TRMSDScoreDef.Create(AName: string);
begin
  inherited Create(AName);
  ScoreTarget:=True;
  ScoreProbe:=True;
  FitTarget:=True;
  FitProbe:=True;
end;

procedure TRMSDScoreDef.SetCoords(Preds, Actuals: TCoords;
  FirstProbe: Integer);
begin
  FPredictedCoords:=Copy(Preds,0,Length(Preds));
  FActualCoords:=Copy(Actuals,0,Length(Actuals));
  FFirstProbeCoord:=FirstProbe;
end;

function TRMSDScoreDef.ScoreType: Integer;
begin
  Result:=ScoreRMSD;
end;

{ TScoreDef }

constructor TScoreDef.Create(AName: string);
begin
  inherited Create;
  FName:=AName;
end;

function TScoreDef.ScoreType: Integer;
begin
  Result:=ScoreNone;
end;


{ TDockOrdersManager }

function TDockOrdersManager.GetLayer(LayerName: string): TPDBModel;

begin
  Result:=FModelMan.LayerByFilename(LayerName);
  if Result=nil then
    begin
    FModelMan.LoadLayer(LayerName);
    Result:=FModelMan.LayerByFilename(LayerName);
    end;
end;

constructor TDockOrdersManager.Create(AJobFile:string);
begin
  inherited Create;
  FModelMan:=TPDBModelMan.Create(Config.MonomersPath);
  FJobFile:=AJobFile;
end;

procedure TDockOrdersManager.Free;
begin
  if Self<>nil then
    begin
    FreeOrders(FDockRuns);
    inherited;
    end;
end;

procedure TDockOrdersManager.NewDockRun;
begin
  SetLength(FDockRuns,Length(FDockRuns)+1);
  FCurrentDock:=High(FDockRuns);
  with FDockRuns[FCurrentDock] do
    begin
    Resolution:=1;
    AddedRadius:=1.35;
    NumAxisSteps:=24;
    SecondsBetweenSaves:=1000;
    CompletedAngles:=0;
    ConstraintSets:=nil;
    ScoreDefs:=nil;
    CompletedAngles:=0;
    TotalAngles:=-1;
    ZeroStats(Stats);
    end;
end;

procedure TDockOrdersManager.PreparePartners(TargetFile, ProbeFile, TargetChains,
  ProbeChains: string; Center: Boolean; RandomizeProbe: Boolean);

var
  tchains,pchains:TSimpleStrings;
  tlayer,player:TPDBModel;
  quat:TQuaternion;

  f:Integer;

begin
  FTarget.Free;
  FProbe.Free;
  FTarget:=nil;
  FProbe:=nil;
  tlayer:=GetLayer(TargetFile);
  player:=GetLayer(ProbeFile);


  //if no chains selected, then whole protein
  if TargetChains='' then
    tchains:=tlayer.ListChains
  else tchains:=SplitChars(TargetChains);
  if ProbeChains='' then
    pchains:=player.ListChains
  else pchains:=SplitChars(ProbeChains);

  FTarget:=tlayer.CopyChains(tchains);
  FProbe:=player.CopyChains(pchains);

  if Center then
    begin
    CenterMolecule(FTarget);
    CenterMolecule(FProbe);
    end;

 if RandomizeProbe then
    begin
    quat:=Quaternion(Random-0.5,Random-0.5,Random-0.5,Random-0.5);
    Normalize(quat);
    FProbe.Transform(quat);
    end;
end;

procedure TDockOrdersManager.PreparePartners(JobIndex:Integer=0);

begin
 FCurrentDock:=JobIndex;
 with FDockRuns[FCurrentDock] do
   PreparePartners(TargetFile, ProbeFile, '', '', False, False);
end;

procedure TDockOrdersManager.SetPDBFiles(TargetFile, ProbeFile: string);

begin
  FDockRuns[FCurrentDock].TargetFile:=TargetFile;
  FDockRuns[FCurrentDock].ProbeFile:=ProbeFile;
end;

procedure TDockOrdersManager.AddResidueContacts(ContactFile: string;
  MaxContacts: Integer);

var
  sl:TStringList;
  contact:TSimpleStrings;
  f:Integer;


begin
  sl:=TStringList.Create;
  sl.LoadFromFile(ContactFile);
  if MaxContacts=0 then MaxContacts:=sl.Count
  else MaxContacts:=Min(MaxContacts,sl.Count);
  for f:=0 to MaxContacts-1 do
    begin
    contact:=SplitString(sl.Strings[f],' ');
    with FDockRuns[FCurrentDock] do
      begin
      SetLength(ConstraintSets,Length(ConstraintSets)+1);
      ConstraintSets[High(ConstraintSets)].NumModels:=StrToInt(contact[5]);
      ConstraintSets[High(ConstraintSets)].MinOverlap:=StrToInt(contact[6]);
      ConstraintSets[High(ConstraintSets)].Name:=sl.Strings[f];
      SetLength(ConstraintSets[High(ConstraintSets)].Constraints,1);
      with ConstraintSets[High(ConstraintSets)].Constraints[0] do
        begin
        Name:='Euclidean '+contact[0]+contact[1]+':'+contact[2]+contact[3];
        TargetPoints:=ListCoords(FTarget.GetGroup(contact[0]).GetGroupById(StrToInt(contact[1])));
        ProbePoints:=ListCoords(FProbe.GetGroup(contact[2]).GetGroupById(StrToInt(contact[3])));
        Distance:=StrToFloat(contact[4]);
        ConstraintType:=CEuclideanDistance;
        end;
      end;
    end;
  sl.Free;
end;

procedure TDockOrdersManager.AddRMSDScore(TemplateFile: string);
      { TODO : only works if same as target and probe files
        Must implement chain mapping with sequence comparisons to generalize }

var
  layer:TPDBModel;
  realcoords,predcoords,tmpcoords:TCoords;
  firstprobe,f:Integer;
  rmsd:TRMSDScoreDef;

begin
  predcoords:=nil;
  realcoords:=nil;
  layer:=GetLayer(TemplateFile);
  //build target coords
  for f:=0 to High(FTarget.Groups) do
    begin
    AppendToArray(ListCoords(FTarget.Groups[f],'CA'),predcoords);
    AppendToArray(ListCoords(layer.Molecule.GetGroup(FTarget.Groups[f].Name),'CA'),
                    realcoords);
    end;
  firstprobe:=Length(predcoords);

  //build probe coords
  for f:=0 to High(FProbe.Groups) do
    begin
    AppendToArray(ListCoords(FProbe.Groups[f],'CA'),predcoords);
    AppendToArray(ListCoords(layer.Molecule.GetGroup(FProbe.Groups[f].Name),'CA'),
                    realcoords);
    end;
  WriteLn(length(predcoords),' ', Length(realcoords));
  with FDockRuns[FCurrentDock] do
    begin
    rmsd:=TRMSDScoreDef.Create('RMSD to '+ExtractFileName(TemplateFile));
    rmsd.SetCoords(predcoords,realcoords,firstprobe);
    rmsd.FitTarget:=True;
    rmsd.FitProbe:=False;
    rmsd.ScoreTarget:=False;
    rmsd.ScoreProbe:=True;
    SetLength(ScoreDefs,Length(ScoreDefs)+1);
    ScoreDefs[High(ScoreDefs)]:=rmsd;
    end;
end;

procedure TDockOrdersManager.AddInterfaceRMSDScore(TemplateFile: string;Dist:TFloat);
var
  layer:TPDBModel;
  realcoords,predcoords,tmpcoords:TCoords;
  firstprobe,f:Integer;
  templatet,templatep:TMolecules;
  interfacet,interfacep:TMolecules;
  rmsd:TRMSDScoreDef;

begin
  predcoords:=nil;
  realcoords:=nil;
  templatet:=nil;
  templatep:=nil;
  layer:=GetLayer(TemplateFile);
  //build template interface groups
  for f:=0 to High(FTarget.Groups) do
    AppendGroupsToArray(templatet,layer.GetChain(FTarget.Groups[f].Name).AllTerminalGroups);
  for f:=0 to High(FProbe.Groups) do
    AppendGroupsToArray(templatep,layer.GetChain(FProbe.Groups[f].Name).AllTerminalGroups);
  GroupsInContact(templatet, templatep, Dist, interfacet, interfacep);

  for f:=0 to High(interfacet) do
    begin
    AppendToArray(ListCoords(interfacet[f],'CA'),realcoords);
    AppendToArray(ListCoords(GetResidue(FTarget,interfacet[f].Parent.Name,interfacet[f].Id),'CA'),predcoords);
    end;
  firstprobe:=Length(predcoords);
  for f:=0 to High(interfacep) do
    begin
    AppendToArray(ListCoords(interfacep[f],'CA'),realcoords);
    AppendToArray(ListCoords(GetResidue(FProbe,interfacep[f].Parent.Name,interfacep[f].Id),'CA'),predcoords);
    end;
  with FDockRuns[FCurrentDock] do
    begin
    rmsd:=TRMSDScoreDef.Create('Interface RMSD to '+
      ExtractFileName(TemplateFile)+' at '+FloatToStrF(Dist,ffFixed,0,3)+'A');
    rmsd.SetCoords(predcoords,realcoords,firstprobe);
    rmsd.FitTarget:=True;
    rmsd.FitProbe:=True;
    rmsd.ScoreTarget:=True;
    rmsd.ScoreProbe:=True;
    SetLength(ScoreDefs,Length(ScoreDefs)+1);
    ScoreDefs[High(ScoreDefs)]:=rmsd;
    end;
end;

procedure TDockOrdersManager.AddNullConstraintSet(NumModels,MinOverlap:Integer);

begin
  with FDockRuns[FCurrentDock] do
    begin
    SetLength(ConstraintSets,Length(ConstraintSets)+1);
    ConstraintSets[High(ConstraintSets)].NumModels:=NumModels;
    ConstraintSets[High(ConstraintSets)].MinOverlap:=MinOverlap;
    ConstraintSets[High(ConstraintSets)].Name:='Unconstrained';
    ConstraintSets[High(ConstraintSets)].Constraints:=nil;
    end;
end;

procedure TDockOrdersManager.LoadJobFile;
begin
  FDockRuns:=LoadOrders(FJobFile);
end;

procedure TDockOrdersManager.SaveOrders(Append: Boolean);

begin
    docktasks.SaveOrders(FDockRuns,FJobFile,Append);
end;

procedure TDockOrdersManager.DeleteScores(Id: string);

{ TODO : Allow deletion of only those scores matching ID }

var f,g:Integer;

begin
  for f:=0 to High(FDockRuns) do
    with FDockRuns[f] do
      begin
      for g:=0 to High(ScoreDefs) do
        ScoreDefs[g]:=nil;
      ScoreDefs:=nil;
      for g:=0 to High(ConstraintSets) do
        ConstraintSets[g].ScoreResults:=nil;
      end;
end;

end.

