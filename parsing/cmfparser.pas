unit cmfparser;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, filebuffer, basetypes, molecules,protein;

type

 TCMFResidue=record
    Name:String;
    Atoms:array of string;
    end;

  TCMFFile=class
    private
    FStructure:TProteinSet;
    FSelectionForm:Pointer;
    FTempCP1,FTempCP2,FTempCP3:Pointer;
    FResTable:array of TCMFResidue;
    FList:TOCIntegers;
    FFileName:string;
    FNames:TOCStrings;
    FBuffer:TFileBuffer;
    FDockOffset:Integer;
    FFileVersion:Integer;
    FDockDataVersion:Integer;
    procedure ReadAtomIdCallback(Atom,Residue,Chain,Protein:Pointer; Tag:Integer);
    procedure WriteAtomIdCallBack(Atom,Residue,Chain,Protein:Pointer; Tag:Integer);
    procedure AddName(s:string);
    procedure AddInteger(I:Integer);
    function GetStructure: TProteinSet;
    procedure SetStructure(const Value: TProteinSet);
    procedure ClearResidueTable;
    procedure ClearLists;
    procedure BuildBasics;
    procedure WriteBasicsCallback(Atom,Residue,Chain,Protein:Pointer;Tag:Integer);
    procedure WriteCoordsCallBack(Atom,Residue,Chain,Protein:Pointer;Tag:Integer);
    procedure ReadCoordsCallBack(Atom,Residue,Chain,Protein:Pointer;Tag:Integer);
    function AddRes(ResName:string):Integer;
    function AddAtom(ResName,AtomName:string; var Rindex:Integer):Integer;
    procedure InitializeWriting;
    procedure WriteDockingInfo;
    procedure WriteDock(Dock:TMPDock);
    procedure WriteSolutions(DockData:TMDDockData);
    procedure WriteClusters(DockData:TMDDockData);
    procedure WriteValues(DockData:TMDDockData);
    procedure ReadDockingInfo;
    procedure ReadDock(Dock:TMPDock);
    procedure ReadSolutions(DockData:TMDDockData);
    procedure ReadClusters(DockData:TMDDockData);
    procedure ReadValues(DockData:TMDDockData);
    procedure ReadDockParameters(DockData:TMDDockData);
    procedure WriteResTable;
    procedure WriteList;
    procedure WriteNames;
    procedure WriteCoords;
    function CheckFile:Boolean;
    procedure ReadResTable;
    procedure ReadLists;
    procedure ReadCoords;
    procedure ReadNames;
    {version 102}
    procedure ReadProteinInfo;
    procedure WriteProteinInfo;
    {version 103}
    procedure WriteProteinSettings(Protein: TMPProtein;SaveDs:Boolean);
    procedure ReadProteinSettings(Protein: TMPProtein;LoadDs:Boolean);
    procedure ReadSettings;
    procedure WriteSettings(SaveDisplaySettings:Boolean);
    procedure ReadStoredSelections;
    procedure WriteStoredSelections;
    procedure ReadAll;
    procedure WriteAll(SaveDisplaySettings:Boolean);
    procedure WriteObjects;
    procedure ReadObjects;
    public
    property Structure:TMPStructures read GetStructure write SetStructure;
    constructor Create(AStructure:TMPStructures; AFileName:string; ASelForm:Pointer);
    procedure Free;
    procedure Load;
    procedure Save(SaveDisplaySettings:Boolean);
    procedure LoadFromBuffer(Buffer:TCMFBuffer);
    procedure SaveToBuffer(Buffer:TCMFBuffer;SaveDisplaySettings:Boolean);

    end;

procedure CMFSaveStructure(FileName:string;Structure:TMPStructures;
  SaveSettings:Boolean;SelectionForm:Pointer=nil);
procedure CMFLoadStructure(FileName:string;Structure:TMPStructures;
  SelectionForm:Pointer=nil);

implementation

uses SysUtils,MolBufferUtils,CMolMain,CMSelection;

const
  VersionString='CMOLFILE';
  {Version 100 Replaced 6 - Feb - 2001
    Modification: atom User Field in Read/Write Coords}

  {Version 101 Replaced 7-Feb- 2001
    Modif: Added
    TMPProteinInfo=record
    Header,Title,Compound,Source,
    KeyWords,Experimental,Author,Journal,Remarks:TMTStrings;
    UserComments:TMTStrings;
    end;
  to TMPProtein in ReadProteinInfo and WriteProteinInfo;
  Removed TMPDock.Parameters, moved parameter information to
  Dock.Dockdata.Info, and changed  ReadDockParameters in ReadDock
  }

  {Version 102 Replaced 27 Aug 2001
    Added:
    procedure WriteProteinSettings;
    procedure ReadProteinSettings;
    procedure ReadSettings;
    procedure WriteSettings;
    to store connections, display settings
    Both called from Load and Save

    TCMFFile.ReadAtomIdCallback
    TCMFFile.WriteAtomIdCallBack
    to write/retrieve ID information
    }

  {Version 103 replaced 5 Nov 2001
    Corrected XRotation on Mol3DCalc
    Changed ReadSolutions}

  {Version 104 replaced 15 jun 2005. added WriteObjects
  to write protein 3d and calc objects
  Added FIDObjects}

  VersionNumber=105;
  FIDResData=0;
  FIDList=1;
  FIDProtChainNames=2;
  FIDAtomCoords=3;
  FIDDockingInfo=4;
  FIDProteinInfo=5;
  FIDSettings=6;
  FIDStoredSelections=7;
  FIDObjects=8;

procedure CMFSaveStructure(FileName:string;Structure:TMPStructures;SaveSettings:Boolean;
  SelectionForm:Pointer=nil);

var CMFile:TCMFFile;

begin
  CMFile:=TCMFFile.Create(Structure,FileName,SelectionForm);
  CMFile.Save(SaveSettings);
  CMFile.Free;
end;

procedure CMFLoadStructure(FileName:string;Structure:TMPStructures;
  SelectionForm:Pointer);

var CMFile:TCMFFile;

begin
  CMFile:=TCMFFile.Create(Structure,FileName,SelectionForm);
  CMFile.Load;
  CMFile.Free;
end;


{ TCMFFile }

function TCMFFile.AddAtom(ResName, AtomName: string; var RIndex:Integer): Integer;

var I:Integer;

begin
  Result:=0;
  I:=AddRes(ResName);
  RIndex:=i;
  with FResTable[i] do
    begin
    while (Result<Length(Atoms)) and (AtomName<>Atoms[Result]) do
      Inc(Result);
    if Result>High(Atoms) then
      begin
      SetLength(Atoms,Length(Atoms)+1);
      Atoms[Result]:=AtomName;
      end;
    end;
end;

function TCMFFile.AddRes(ResName: string): Integer;
begin
  Result:=0;
  while (Result<Length(FResTable)) and (ResName<>FResTable[Result].Name) do
    Inc(Result);
  if Result>High(FResTable) then
    begin
    SetLength(FResTable,Length(FResTable)+1);
    Pointer(FResTable[Result].Atoms):=nil;
    FResTable[Result].Name:=ResName;
    end
end;

procedure TCMFFile.BuildBasics;
begin
  ClearLists;
  if (FStructure<>nil) then
    begin
    AddInteger(FStructure.Count+FStructure.DockProteinsCount);
    FStructure.EnumerateDockingAtoms(WriteBasicsCallback,0);
    FStructure.EnumerateAtoms(WriteBasicsCallback,0);
    end;
end;

procedure TCMFFile.WriteBasicsCallback(Atom, Residue, Chain,
  Protein: Pointer; Tag: Integer);

var ai,ri:Integer;

begin
  ai:=AddAtom(TMPResidue(Residue).Name,TMPAtom(Atom).Name,ri);
  if Protein<>FTempCP1 then
    begin
    AddInteger(TMPProtein(Protein).ChainNumber);
    AddName(TMPProtein(Protein).FileName);
    FTempCP1:=Protein;
    end;
  if Chain<>FTempCP2 then
    begin
    AddInteger(TMPChain(Chain).NumResidues);
    AddName(TMPChain(Chain).Name);
    FTempCP2:=Chain;
    end;
  if Residue<>FTempCP3 then
    begin
    AddInteger(ri);
    AddInteger(TMPResidue(Residue).Id);
    AddInteger(TMPResidue(Residue).NumAtoms);
    FTempCP3:=Residue;
    end;
  AddInteger(ai);
end;

procedure TCMFFile.ClearLists;
begin
  FTempCP1:=nil;
  FTempCP2:=nil;
  FTempCP3:=nil;
  FList:=nil;
  FNames:=nil;
  ClearResidueTable;
end;

procedure TCMFFile.ClearResidueTable;

var f:Integer;

begin
  if FResTable<>nil then
    for f:=0 to High(FResTable) do
      FResTable[f].Atoms:=nil;
  FResTable:=nil;
end;

constructor TCMFFile.Create(AStructure: TMPStructures; AFileName:string;ASelForm:Pointer);
begin
  inherited Create;
  FSelectionForm:=ASelForm;
  FStructure:=AStructure;
  FBuffer:=TCMFBuffer.Create;
  FFileName:=AFileName;
end;

procedure TCMFFile.Free;
begin
  if Self<>nil then
    begin
    ClearLists;
    ClearResidueTable;
    FBuffer.Free;
    inherited Free;
    end;
end;

function TCMFFile.GetStructure: TMPStructures;
begin
  Result:=FStructure;
end;

procedure TCMFFile.InitializeWriting;
begin
  FBuffer.ClearBuffer;
  FBuffer.WriteString(VersionString);
  FBuffer.WriteInteger(VersionNumber);
end;

procedure TCMFFile.Save(SaveDisplaySettings: Boolean);

begin
  if (FStructure<>nil) and (FFileName<>'') then
    begin
    WriteAll(SaveDisplaySettings);
    FBuffer.FlushToFile(FFileName);
    FStructure.Changed:=False;
    end;
end;

procedure TCMFFile.SetStructure(const Value: TMPStructures);
begin
  FStructure:=Value;
end;

procedure TCMFFile.WriteResTable;

var Res,At:Integer;

begin
  FBuffer.WriteInteger(FIDResData);
  FBuffer.WriteInteger(Length(FResTable));
  if FResTable<>nil then
    for Res:=0 to High(FResTable) do
      begin
      FBuffer.WriteString(FResTable[Res].Name);
      FBuffer.WriteInteger(Length(FResTable[Res].Atoms));
      if FResTable[Res].Atoms<>nil then
        for At:=0 to High(FResTAble[Res].Atoms) do
          FBuffer.WriteString(FResTable[Res].Atoms[At]);
      end;
end;

procedure TCMFFile.AddInteger(I: Integer);
begin
  SetLength(FList,Length(FList)+1);
  FList[High(FList)]:=I;
end;

procedure TCMFFile.AddName(s: string);
begin
  SetLength(FNames,Length(FNames)+1);
  FNames[High(FNames)]:=s;
end;


procedure TCMFFile.WriteList;

var f:Integer;

begin
  FBuffer.WriteInteger(FIDList);
  FBuffer.WriteInteger(Length(FList));
  if FList<>nil then
    for f:=0 to High(FList) do
      FBuffer.WriteInteger(FList[f]);
end;

procedure TCMFFile.WriteNames;

var f:Integer;

begin
  FBuffer.WriteInteger(FIDProtChainNames);
  FBuffer.WriteInteger(Length(FNames));
  if FNames<>nil then
    for f:=0 to High(FNames) do
      FBuffer.WriteString(FNames[f]);
end;


procedure TCMFFile.WriteCoordsCallBack(Atom, Residue, Chain,
  Protein: Pointer; Tag: Integer);
begin
  with TMPAtom(Atom) do
    begin
    FBuffer.WriteSingle(Coords[1]);
    FBuffer.WriteSingle(Coords[2]);
    FBuffer.WriteSingle(Coords[3]);
    {Version 101}
    FBuffer.WriteSingle(User);
    end;
end;

procedure TCMFFile.WriteCoords;
begin
  FBuffer.WriteInteger(FIDAtomCoords);
  FStructure.EnumerateDockingAtoms(WriteCoordsCallBack,0);
  FStructure.EnumerateAtoms(WriteCoordsCallBack,0);
end;

procedure TCMFFile.Load;

begin
  if (FStructure<>nil) and (FFileName<>'') then
    begin
    FBuffer.LoadFromFile(FFileName);
    ReadAll;
    end;
end;

function TCMFFile.CheckFile: Boolean;

var
  Ver:string;

begin
  FBuffer.ReadString(Ver);
  FBuffer.ReadInteger(FFileVersion);
  Result:=(Ver=VersionString) and (FFileVersion<=VersionNumber);
end;

procedure TCMFFile.ReadCoords;
begin
  if FStructure<>nil then
    begin
    FStructure.EnumerateDockingAtoms(ReadCoordsCallback,0);
    FStructure.EnumerateAtoms(ReadCoordsCallback,0);
    end;
end;

procedure TCMFFile.ReadLists;

var
  p,np,c,nc,r,nr,a,na,rt,rid,at:Integer;
  Prot:TMPProtein;
  Chain:TMPChain;
  Res:TMPResidue;
  Atom:TMPAtom;
  AtomNumber:Integer;

begin
  FBuffer.ReadInteger(np);
  if np>0 then
  begin
  FBuffer.ReadInteger(np);
  if np>0 then for p:=1 to np do
    begin
    Prot:=FStructure.AddMissingProtein;
    AtomNumber:=0;
    FBuffer.ReadInteger(nc);
    if nc>0 then for c:=1 to nc do
      begin
      Chain:=Prot.AddChain;
      FBuffer.ReadInteger(nr);
      if nr>0 then for r:=1 to nr do
        begin
        FBuffer.ReadInteger(rt);
        FBuffer.ReadInteger(rid);
        Res:=Chain.AddResidue(FResTable[rt].Name,rid);
        FBuffer.ReadInteger(na);
        if na>0 then for a:=1 to na do
          begin
          Inc(AtomNumber);
          FBuffer.ReadInteger(at);
          Atom:=Res.AddAtom(FResTable[rt].Atoms[at]);
          Atom.Id:=AtomNumber;
          end;
        end;
      Chain.FixTemplates;
      Chain.BuildConnections;
      end;
    end;
    end;
end;

procedure TCMFFile.ReadNames;

var
  p,np,c,nc:Integer;
  Prot:TMPProtein;
  Chain:TMPChain;
  s:string;

begin
  FBuffer.ReadInteger(np);
  if np>0 then
  begin
  np:=FStructure.DockProteinsCount+FStructure.Count-1;
  if np>=0 then for p:=0 to np do
    begin
    Prot:=FStructure.DockProteinByIndex(p);
    FBuffer.ReadString(s);
    Prot.FileName:=s;
    nc:=Prot.ChainNumber-1;
    if nc>=0 then for c:=0 to nc do
      begin
      Chain:=Prot.ChainByIndex(c);
      FBuffer.ReadString(s);
      Chain.Name:=S;
      end;
    end;
    end;
end;


procedure TCMFFile.ReadResTable;

var
  nr,r,na,a:Integer;
  s:string;

begin
  FBuffer.ReadInteger(nr);
  if nr>0 then
    begin
    SetLength(FResTable,nr);
    for r:=0 to nr-1 do
      begin
      FBuffer.ReadString(s);
      FResTable[r].Name:=s;
      Pointer(FResTable[r].Atoms):=nil;
      FBuffer.ReadInteger(na);
      SetLength(FResTable[r].Atoms,na);
      if na>0 then for a:=0 to na-1 do
        begin
        FBuffer.ReadString(s);
        FResTable[r].Atoms[a]:=s;
        end;
      end;
    end;
end;

procedure TCMFFile.ReadCoordsCallBack(Atom, Residue, Chain,
  Protein: Pointer; Tag: Integer);

var
  s:Single;
  c:TMTCoords;

begin
  FBuffer.ReadSingle(s);
  c[1]:=s;
  FBuffer.ReadSingle(s);
  c[2]:=s;
  FBuffer.ReadSingle(s);
  c[3]:=s;
  TMPAtom(Atom).Coords:=c;
  {Version 101}
  if FFileVersion>100 then
    begin
    FBuffer.ReadSingle(s);
    TMPAtom(Atom).User:=s;
    end;
end;

procedure TCMFFile.WriteDockingInfo;

var
  Count,f:Integer;
  Dock:TMPDock;

begin
  Count:=FStructure.DockCount;
  if Count>0 then
    begin
    FBuffer.WriteInteger(FIDDockingInfo);
    FBuffer.WriteInteger(Count);
    for f:=0 to Count-1 do
      begin
      Dock:=FStructure.GetDock(f);
      WriteDock(Dock);
      end;
    end;
end;

procedure TCMFFile.WriteSolutions(DockData: TMDDockData);

var
  f,Count,t:Integer;
  c:TMTCoords;

begin
  Count:=DockData.SolutionCount;
  FBuffer.WriteInteger(Count);
  if Count>0 then
    for f:=0 to Count-1 do
      begin
      c:=DockData.SolutionPosition(f);
      for t:=1 to 3 do
        FBuffer.WriteSingle(c[t]);
      c:=DockData.SolutionOrientation(f);
      for t:=1 to 3 do
        FBuffer.WriteSingle(c[t]);
      end;
end;

procedure TCMFFile.WriteDock(Dock: TMPDock);

procedure WriteFlags;

var
  Flags:TMTIntegers;
  f:Integer;

begin
  Flags:=nil;
  Flags:=Dock.DockData.ExportFlags;
  FBuffer.WriteInteger(Length(Flags));
  for f:=0 to High(Flags) do
    FBuffer.WriteInteger(Flags[f]);
  Flags:=nil;
end;

procedure WriteDockInfo;

begin
  with Dock.DockData.Info do
    begin
    WriteIntegers(FBuffer,TotalSteps);

    WriteStrings(FBuffer,UserComments);
    FBuffer.WriteSingle(Resolution);
    FBuffer.WriteSingle(AddedRadius);
    FBuffer.WriteInteger(MinOverlap);
    FBuffer.WriteInteger(MaxStored);
    FBuffer.WriteInteger(Smoothing);
    FBuffer.WriteSingle(AngularStep);
    MolBufferUtils.WriteCoords(FBuffer,FirstAngle);
    MolBufferUtils.WriteCoords(FBuffer,LastAngle);
    WriteIntegers(FBuffer,FiltersUsed);
    WriteSingles(FBuffer,FilterCutoffValues);

    WriteStrings(FBuffer,Machines);
    WriteSingles(FBuffer,ComputationTime);
    WriteSingles(FBuffer,Received);
    WriteSingles(FBuffer,Complete);
    WriteSingles(FBuffer,Delivered);
    WriteIntegers(FBuffer,StepCounts);
    FBuffer.WriteInteger(ProbeCuts);
    FBuffer.WriteInteger(TargetCuts);
    FBuffer.WriteInteger(MinExposure);
    {Version 102 of DockData}
    FBuffer.WriteInteger(Origin);
    {Version 103 of DockData}
    FBuffer.WriteInteger(Symmetry);
    FBuffer.WriteInteger(RotationType);
    FBuffer.WriteInteger(MonomerCount);
    end
end;

begin
  FBuffer.WriteString(Dock.Name);
  FBuffer.WriteInteger(Dock.CopyTarget);
  FBuffer.WriteInteger(Dock.CopyProbe);
  FBuffer.WriteInteger(Dock.DockData.Version);
  WriteFlags;
  WriteSolutions(Dock.DockData);
  WriteClusters(Dock.DockData);
  WriteValues(Dock.DockData);
  {Version 101 of DockData}
  WriteDockInfo;
end;

procedure TCMFFile.WriteClusters(DockData: TMDDockData);

var
  f,Count,g:Integer;
  List:TMTIntegers;

begin
  Pointer(List):=nil;
  Count:=DockData.ClusterCount;
  FBuffer.WriteInteger(Count);
  if Count>0 then
    for f:=0 to Count-1 do
      begin
      List:=DockData.GetCluster(f);
      FBuffer.WriteInteger(Length(List));
      if List<>nil then
        for g:=0 to High(List) do
          FBuffer.WriteInteger(List[g]);
      List:=nil;
      end;
end;

procedure TCMFFile.WriteValues(DockData: TMDDockData);

var
  f,Count:Integer;
  Data:TMDBasicData;

begin
  Count:=DockData.DataCount;
  FBuffer.WriteInteger(Count);
  if Count>0 then
    for f:=0 to Count-1 do
      begin
      Data:=DockData.GetDataByIndex(f);
      FBuffer.WriteInteger(Data.DataType);
      Data.WriteToBuffer(FBuffer);
      end;
end;

procedure TCMFFile.ReadClusters(DockData: TMDDockData);

var
  f,Count,g:Integer;
  List:TMTIntegers;

begin
  Pointer(List):=nil;
  FBuffer.ReadInteger(Count);
  DockData.ClearClusters;
  if Count>0 then
    for f:=0 to Count-1 do
      begin
      List:=DockData.GetCluster(f);
      FBuffer.ReadInteger(g);
      if g>0 then
        begin
        SetLength(List,g);
        for g:=0 to High(List) do
          FBuffer.ReadInteger(List[g]);
        DockData.AppendCluster(List);
        end;
      List:=nil;
      end;
end;

procedure TCMFFile.ReadDock(Dock: TMPDock);

var
  s:string;
  i:Integer;

procedure ReadFlags;

var
  Flags:TMTIntegers;
  f:Integer;

begin
  Pointer(Flags):=nil;
  FBuffer.ReadInteger(f);
  SetLength(Flags,f);
  for f:=0 to High(Flags) do
    FBuffer.ReadInteger(Flags[f]);
  Dock.DockData.ImportFlags(Flags);
  Flags:=nil;
end;

procedure ReadDockInfo;

var Info:TMDDockDataInfo;

begin
  Info:=Dock.DockData.Info;
  with Info do
    begin
    TotalSteps:=ReadIntegers(FBuffer);
    UserComments:=ReadStrings(FBuffer);

    Resolution:=FBuffer.GetSingle;
    AddedRadius:=FBuffer.GetSingle;
    MinOverlap:=FBuffer.GetInteger;
    MaxStored:=FBuffer.GetInteger;
    Smoothing:=FBuffer.GetInteger;
    AngularStep:=FBuffer.GetSingle;
    FirstAngle:=MolBufferUtils.ReadCoords(FBuffer);
    LastAngle:=MolBufferUtils.ReadCoords(FBuffer);
    FiltersUsed:=ReadIntegers(FBuffer);
    FilterCutoffValues:=ReadSingles(FBuffer);

    Machines:=ReadStrings(FBuffer);
    ComputationTime:=ReadSingles(FBuffer);
    Received:=ReadSingles(FBuffer);
    Complete:=ReadSingles(FBuffer);
    Delivered:=ReadSingles(FBuffer);
    StepCounts:=ReadIntegers(FBuffer);
    FBuffer.ReadInteger(ProbeCuts);
    FBuffer.ReadInteger(TargetCuts);
    FBuffer.ReadInteger(MinExposure);
    {Version 102 of DockData}
    if FDockDataVersion>=102 then
      FBuffer.ReadInteger(Origin);
    {Version 103 of DockData}
    if FDockDataVersion>=103 then
      begin
      FBuffer.ReadInteger(Symmetry);
      FBuffer.ReadInteger(RotationType);
      FBuffer.ReadInteger(MonomerCount);
      end;
    end;
  Dock.DockData.Info:=Info;
end;

begin
  FBuffer.ReadString(s);
  Dock.Name:=s;
  FBuffer.ReadInteger(i);
  if i>=0 then i:=i+FDockOffset;
  Dock.CopyTarget:=i;
  FBuffer.ReadInteger(i);
  if i>=0 then i:=i+FDockOffset;
  Dock.CopyProbe:=i;

  {Version 102}
  if FFileVersion<=101 then ReadDockParameters(Dock.DockData);

  FBuffer.ReadInteger(i);
  FDockDataVersion:=i;
  ReadFlags;
  ReadSolutions(Dock.DockData);
  ReadClusters(Dock.DockData);
  ReadValues(Dock.DockData);

  {Version 101 of DockData}
  if FDockDataVersion>=101 then
    ReadDockInfo;
  Dock.DockData.MaxSolutions:=Dock.DockData.Info.MaxStored;
end;

procedure TCMFFile.ReadDockingInfo;
var
  Count,f:Integer;
  Dock:TMPDock;

begin
  FBuffer.ReadInteger(Count);
  if Count>0 then
    begin
    FDockOffset:=FStructure.DockCount;
    for f:=0 to Count-1 do
      begin
      Dock:=FStructure.AddDock;
      ReadDock(Dock);
      end;
    end;
end;

procedure TCMFFile.ReadSolutions(DockData: TMDDockData);

var
  f,Count,t:Integer;
  c1,c2:TMTCoords;
  s:Single;

begin
  FBuffer.ReadInteger(Count);
  DockData.ClearSolutions;
  if Count>0 then
    begin
    DockData.SetSolutionLength(Count);
    for f:=0 to Count-1 do
      begin
      for t:=1 to 3 do
        begin
        FBuffer.ReadSingle(s);
        c1[t]:=s;
        end;
      for t:=1 to 3 do
        begin
        FBuffer.ReadSingle(s);
        c2[t]:=s;
        end;
      {Version 104, correcting XRotation errors}
      if (FFileVersion<104) and
        (DockData.Info.Origin=MDOriginBigger)  then
          c2[1]:=-c2[1];
      DockData.SetSolution(f,c1,c2);
      end;
    end;
end;

procedure TCMFFile.ReadValues(DockData: TMDDockData);

var
  DType,f,Count:Integer;
  Data:TMDBasicData;

begin
  FBuffer.ReadInteger(Count);
  if Count>0 then
    for f:=0 to Count-1 do
      begin
      FBuffer.ReadInteger(DType);
      Data:=NewDataByType(DType);
      if Data<>nil then
        begin
        Data.ReadFromBuffer(FBuffer);
        DockData.AddData(Data);
        end;
      end;
end;


procedure TCMFFile.ReadDockParameters(DockData:TMDDockData);

var Info:TMDDockDataInfo;

begin
  {Outdated by version 102, 7-feb-2001
  Retained only for compatibility with CMF 100, 101}
  Info:=DockData.Info;
  with FBuffer do
  with Info do
    begin
    GetInteger{Version};
    Resolution:=GetSingle;
    AddedRadius:=GetSingle;
    AngularStep:=GetInteger/180*PI{was in degrees, now in Rad};
    MinOverlap:=GetInteger;
    MaxStored:=GetInteger;
    TotalSteps:=MolBufferUtils.ReadIntegers(FBuffer);
    UserComments:=MolBufferUtils.ReadStrings(FBuffer);
  end;
  DockData.Info:=Info;
end;

procedure TCMFFile.WriteProteinInfo;

var
  Protein:TMPProtein;
  p,np:Integer;


begin
  FBuffer.WriteInteger(FIDProteinInfo);
  np:=FStructure.DockProteinsCount+FStructure.Count-1;
  if np>=0 then for p:=0 to np do
    begin
    Protein:=FStructure.DockProteinByIndex(p);
    with Protein.ProteinInfo do
      begin
      WriteStrings(FBuffer,Header);
      WriteStrings(FBuffer,Title);
      WriteStrings(FBuffer,Compound);
      WriteStrings(FBuffer,Source);
      WriteStrings(FBuffer,KeyWords);
      WriteStrings(FBuffer,Experimental);
      WriteStrings(FBuffer,Author);
      WriteStrings(FBuffer,Journal);
      WriteStrings(FBuffer,Remarks);
      WriteStrings(FBuffer,UserComments);
      end;
    end;
end;

procedure TCMFFile.ReadProteinInfo;

var
  Info:TMPProteinInfo;
  Protein:TMPProtein;
  p,np:Integer;

begin
  np:=FStructure.DockProteinsCount+FStructure.Count-1;
  if np>=0 then for p:=0 to np do
    begin
    Protein:=FStructure.DockProteinByIndex(p);
    with Info do
      begin
      Header:=ReadStrings(FBuffer);
      Title:=ReadStrings(FBuffer);
      Compound:=ReadStrings(FBuffer);
      Source:=ReadStrings(FBuffer);
      KeyWords:=ReadStrings(FBuffer);
      Experimental:=ReadStrings(FBuffer);
      Author:=ReadStrings(FBuffer);
      Journal:=ReadStrings(FBuffer);
      Remarks:=ReadStrings(FBuffer);
      UserComments:=ReadStrings(FBuffer);
      end;
    Protein.ProteinInfo:=Info;
    end;
end;

procedure TCMFFile.ReadProteinSettings(Protein: TMPProtein;LoadDs:Boolean);

begin
  if LoadDs then
      Protein.ReadSettings(FBuffer);
  Protein.ReadConnections(FBuffer);
end;

procedure TCMFFile.WriteProteinSettings(Protein: TMPProtein;SaveDS:Boolean);
begin
  if SaveDs then
    Protein.WriteSettings(FBuffer);
  Protein.WriteConnections(FBuffer);
end;

procedure TCMFFile.ReadSettings;
var
  Protein:TMPProtein;
  p,np:Integer;
  LoadDS:Boolean;
begin
  FBuffer.ReadBoolean(LoadDS);
  np:=FStructure.DockProteinsCount+FStructure.Count-1;
  if np>=0 then for p:=0 to np do
    begin
    Protein:=FStructure.DockProteinByIndex(p);
    Protein.EnumerateAtoms(ReadAtomIdCallback,0);
    ReadProteinSettings(Protein,LoadDS);
    end;
  FStructure.FixDockProteins;
  if LoadDS then MainForm.ReadCMFSettings(FBuffer,FFileVersion);
end;

procedure TCMFFile.WriteSettings(SaveDisplaySettings:Boolean);

var
  Protein:TMPProtein;
  p,np:Integer;

begin
  FBuffer.WriteInteger(FIDSettings);
  FBuffer.WriteBoolean(SaveDisplaySettings);
  np:=FStructure.DockProteinsCount+FStructure.Count-1;
  if np>=0 then for p:=0 to np do
    begin
    Protein:=FStructure.DockProteinByIndex(p);
    Protein.EnumerateAtoms(WriteAtomIdCallback,0);
    WriteProteinSettings(Protein,SaveDisplaySettings);
    end;
  if SaveDisplaySettings then MainForm.WriteCMFSettings(FBuffer);
end;

procedure TCMFFile.ReadAtomIdCallback(Atom, Residue, Chain,
  Protein: Pointer; Tag: Integer);

var f:Integer;

begin
  FBuffer.ReadInteger(f);
  TMPAtom(Atom).Id:=f;
  if TMPProtein(Protein).HighestId<f then TMPProtein(Protein).HighestId:=f;
end;

procedure TCMFFile.WriteAtomIdCallBack(Atom, Residue, Chain,
  Protein: Pointer; Tag: Integer);
begin
  FBuffer.WriteInteger(TMPAtom(Atom).Id);
end;

procedure TCMFFile.ReadStoredSelections;
begin
  if FSelectionForm<>nil then
    TSelectForm(FSelectionForm).ReadStoredSelections(FBuffer)
  else
    begin
    FSelectionForm:=TSelectForm.Create(nil,nil);
    TSelectForm(FSelectionForm).ReadStoredSelections(FBuffer);
    TSelectForm(FSelectionForm).Free;
    FSelectionForm:=nil;
    end
end;

procedure TCMFFile.WriteStoredSelections;
begin
  if FSelectionForm<>nil then
    begin
    FBuffer.WriteInteger(FIDStoredSelections);
    TSelectForm(FSelectionForm).WriteStoredSelections(FBuffer);
    end;
end;

procedure TCMFFile.LoadFromBuffer(Buffer: TCMFBuffer);
var Tmp:TCMFBuffer;

begin
  Tmp:=FBuffer;
  FBuffer:=Buffer;
  ReadAll;
  FBuffer:=Tmp;
end;

procedure TCMFFile.SaveToBuffer(Buffer: TCMFBuffer;SaveDisplaySettings:Boolean);

var Tmp:TCMFBuffer;

begin
  Tmp:=FBuffer;
  FBuffer:=Buffer;
  WriteAll(SaveDisplaySettings);
  FBuffer:=Tmp;
end;

procedure TCMFFile.ReadAll;

var Field:Integer;

begin
  ClearLists;
  FDockOffset:=0;
  if CheckFile then
    begin
    while not FBuffer.EndOfBuffer do
      begin
      FBuffer.ReadInteger(Field);
      case Field of
        FIDResData:ReadResTable;
        FIDList:ReadLists;
        FIDProtChainNames:ReadNames;
        FIDAtomCoords:ReadCoords;
        FIDDockingInfo:ReadDockingInfo;
        {version 102}
        FIDProteinInfo:ReadProteinInfo;
        {version 103}
        //After V105 objects are incompatible, and ReadSettings had to be changed
        // not to read objects, so older files are truncated here
        FIDSettings:if FFileVersion>=105 then ReadSettings else Break;
        FIDStoredSelections:if FFileVersion>=105 then ReadStoredSelections else Break;
        {Version 105}
        FIDObjects:ReadObjects;
        end;
      end;
    end;
end;

procedure TCMFFile.WriteAll(SaveDisplaySettings:Boolean);
begin
  if FStructure<>nil then
    begin
    BuildBasics;
    InitializeWriting;
    WriteDockingInfo;
    WriteResTable;
    WriteList;
    WriteNames;
    WriteCoords;
    {Version 102};
    WriteProteinInfo;
    {version 103}
    WriteSettings(SaveDisplaySettings);
    WriteStoredSelections;
    {version 105}
    WriteObjects;
    end;
end;

procedure TCMFFile.WriteObjects;

var
  p,f:Integer;
  Prot:TMPProtein;

begin
  FBuffer.WriteInteger(FIDObjects);
  for p:=0 to FStructure.Count-1 do
    FStructure.GetProtein(p).WriteObjects(FBuffer);
end;

procedure TCMFFile.ReadObjects;
var
  p,f:Integer;
  Prot:TMPProtein;

begin
  for p:=0 to FStructure.Count-1 do
    FStructure.GetProtein(p).ReadObjects(FBuffer);
end;

end.
