{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 21.12.2013
Purpose:
  Classify potential contacts
Requirements:
Revisions:
To do: Still in experimental stage.
*******************************************************************************}

program cclassif;

{$mode objfpc}{$H+}

uses
  {$IFDEF UNIX}{$IFDEF UseCThreads}
  cthreads,
  {$ENDIF}{$ENDIF}
  Classes, SysUtils, CustApp, lrndata, basetypes, filebuffer, stringutils,
  quicksort, lrnhistogram, lrnconfig, lrnnaivebayes, lrnneural, geomutils;

const
  //Feature modes
  FMSelect=1;
  FMAggregate=2;

type

  { TContactsApp }

  TContactsApp = class(TCustomApplication)
  protected
    procedure DoRun; override;
    procedure ClassifyData(DM:TDataManager;ClassCol:string);
  public
    constructor Create(TheOwner: TComponent); override;
    destructor Destroy; override;
    procedure WriteHelp; virtual;
    procedure TestScript;
    procedure CompileBinaryData(FoldFileName,SourceFolder,DestFolder,ClassCol:string;TrainData:string='');
    procedure CreateHistograms(DataFile,HistFolder:string;NumFolds:Integer);
    procedure DisplayHistogramStats(HistFolder:string;NumFolds:Integer);
    procedure NBFeatureSearchBeam(DataFile,HistFolder:string;
                NumFolds,BeamWidth:Integer;Mode,MaxMAP,HistBins,NumCycles:Integer;
                ReportFile:string;
                Disjoint:Boolean=False;Exclusions:TSimpleStrings=nil);
    procedure NBEvaluate(HistFolder,DataFile,FeaturesFile,ReportFile:string;
                ContactsFolder:string='');
      //if contactsfolder<>'' then exports sorted contacts for each complex to
      //that folder, with the original file name

    procedure SelectAggregates(InFile,OutFile:string);
      //file formats: tab separated, first column is score (preserved), following column names
      //can contain repeated tabs (empty strings are ignored)
      //outputs all those that do not have repeated features
    procedure Aggregate(DataFile,AggregateFile,OutFile:string;DropColumns:Boolean);
      //Adds to the data file new columns with the average of those indicated in the aggregatefile
      //Aggregate file format: tab separated, first column is score (preserved), following column names

  end;

{ TContactsApp }

procedure TContactsApp.DoRun;
var
  ErrorMsg: String;
  MaxMAP:Integer;
begin
  // quick check parameters
  ErrorMsg:=CheckOptions('h','help');
  if ErrorMsg<>'' then begin
    ShowException(Exception.Create(ErrorMsg));
    Terminate;
    Exit;
  end;

  // parse parameters
  if HasOption('h','help') then begin
    WriteHelp;
    Terminate;
    Exit;
  end;

  { add your program here }
   LrnVerbose:=True;

   //Actual procedure

  {Compile binary data}
  {CompileBinaryData('C:\My Documents\Research\coevol\SVM-experiment\G-trainingfolds.txt',
                      'H:\Coevol\data','H:\Coevol\bindata\FullASA','Full Contact ASA');

  CreateHistograms('H:\Coevol\bindata\FullASA\data.bin','H:\Coevol\bindata\FullASA',5);}


  {Feature selection}
  {procedure NBFeatureSearchBeam(
    DataFile: binary file with training data
    NumFolds: number of cross validation folds
    BeamWidth: number of combinations kept at each search iteration
    Mode: select or aggregate. Aggregate is not very good
    MaxMAP: ignored. Used to be the maximum number of true contacts considered for MAP but
          now is equal to the number of true contacts;
    HistBins: bins for report histogram
    NumCycles: number of search cycles
    ReportFile:string;
    Disjoint:Boolean=False: true to forbid repeated occurences of same feature in retained lists
    Exclusions:TSimpleStrings=nil: features to excluse from search}

  {NBFeatureSearchBeam('H:\Coevol\bindata\FullASA5\data.bin',
    'H:\Coevol\bindata\FullASA5',5,20,FMSelect,100,30,45,
    'H:\Coevol\MAPAllBasicRPrecisionPen5A.txt',False);}



  {Compile test data}

  {CompileBinaryData('C:\My Documents\Research\coevol\SVM-experiment\G-testset.txt',
                      'H:\Coevol\data','H:\Coevol\bindata\FullASA5A64_test','Full Contact ASA',
                      'H:\Coevol\bindata\FullASA5\data.bin');}

  {create histograms from original data but in single fold}
  {CreateHistograms('H:\Coevol\bindata\FullASA5A64\data.bin','H:\Coevol\bindata\FullASA5A64_test',1);}

  {Test with test data}
  {NBEvaluate('H:\Coevol\bindata\FullASA5A64_test',
              'H:\Coevol\bindata\FullASA5A64_test\data.bin',
              'H:\Coevol\bindata\FullASA5A64_test\SelectedFeatures.txt',
              'H:\Coevol\bindata\FullASA5A64_test\Test.txt',
              'H:\Coevol\bindata\fullASA5A64_test\contacts');}
  {NBEvaluate('H:\Coevol\bindata\FullASA5_test',
              'H:\Coevol\bindata\FullASA5_test\data.bin',
              'H:\Coevol\bindata\FullASA5_test\SelectedFeatures.txt',
              'H:\Coevol\bindata\FullASA5_test\Test.txt',
              'H:\Coevol\bindata\fullASA5_test\contacts');}
  {last argument is to folder where to save contact files to use in docking; ranked by classification}

  //DisplayHistogramStats('H:\Coevol\5A_38cutoff\bindata\FullASA5',5);

  {***************************
  ****** Unbound tests *******
  ****************************}


  {CompileBinaryData('C:\My Documents\Research\coevol\SVM-experiment\Unbound\G-testset.txt',
                      'H:\_Research\Coevol\data',
                      'H:\_Research\Coevol\FullASA5_unbound_test',
                      'Full Contact ASA',
                      'H:\_Research\Coevol\5A_38cutoff\bindata\FullASA5\data.bin');}

  //Copy selected features and 0.hist from FullASA5_test

  NBEvaluate('H:\_Research\Coevol\FullASA5_unbound_test',
              'H:\_Research\Coevol\FullASA5_unbound_test\data.bin',
              'H:\_Research\Coevol\FullASA5_unbound_test\SelectedFeatures.txt',
              'H:\_Research\Coevol\FullASA5_unbound_test\Test.txt',
              'H:\_Research\Coevol\FullASA5_unbound_test\contacts');


  // stop program loop


  Terminate;
end;

procedure TContactsApp.ClassifyData(DM: TDataManager; ClassCol: string);

//Drops all other classification data

const ColsToDrop:array [0..5] of string = ('Sidechain Contact ASA',
      'Full Contact ASA','Sidechain Contact ASA (all to one)',
      'Full Contact ASA (all to one)','Sidechain Contact ASA (all to all)',
      'Full Contact ASA (all to all)');


var
  colix,f,g:Integer;
  table:TLearnData;
  dropcols:TSimpleStrings;


begin
  colix:=DM.Tables[0].ColumnIndex(ClassCol);
  for f:=0 to High(DM.Tables) do
    begin
    table:=DM.Tables[f];
    SetLength(table.DataClasses,table.Rows);
    for g:=0 to High(table.DataClasses) do
      begin
      if TFloatDataColumn(table.Columns[colix]).Vals[g]>0 then
        table.DataClasses[g]:=1
      else
        table.DataClasses[g]:=0;
      end;
    WriteLn(table.Name,': ',CountInArray(1,Table.DataClasses));
    end;
  dropcols:=nil;
  for f:=0 to High(ColsToDrop) do
    AddToArray(ColsToDrop[f],dropcols);
  DM.DropColumns(dropcols);
end;

constructor TContactsApp.Create(TheOwner: TComponent);
begin
  inherited Create(TheOwner);
  StopOnException:=True;
end;

destructor TContactsApp.Destroy;
begin
  inherited Destroy;
end;

procedure TContactsApp.WriteHelp;
begin
  { add your help code here }
  writeln('Usage: ',ExeName,' -h');
end;

procedure TContactsApp.TestScript;

var
  files:TSimpleStrings;
  dataman:TDataManager;
  sl:TStringList;
  selcols:TSimplesTrings;
  seltables:TIntegers;
  cont,nocont:TMatrix;

var f:Integer;

begin
  files:=nil;
  AddToArray('H:\Coevol\data\1ay7.dat',files);
  AddToArray('H:\Coevol\data\1bvn.dat',files);
  dataman:=TDataManager.Create();
  dataman.LoadTables(files,'Target ID-Chain-Res Name-Res ID:Probe ID-Chain-Res Name-Res ID');
  WriteLn('ToClassify');
  ClassifyData(dataman,'Sidechain Contact ASA');
  selcols:=nil;
  seltables:=nil;
  AddToArray('SCOTCH Score(50full) (all to one)',selcols);
  AddToArray(0,seltables);
  AddToArray(1,seltables);
  WriteLn('ToScale');
  dataman.ScaleToMinusPlusOne;
  WriteLn('Scaled');
  cont:=dataman.SelectData(selcols,seltables,1);
  nocont:=dataman.SelectData(selcols,seltables,0);
  {for f:=0 to High(cont) do
    WriteLn('c:'+#9+FloatToStr(cont[f,0]));
  WriteLn(IntToStr(Length(cont)));
  WriteLn(IntToStr(Length(nocont)));}
end;

procedure TContactsApp.CompileBinaryData(FoldFileName, SourceFolder,
  DestFolder,ClassCol: string;TrainData:string);

var
  sl:TStringList;
  fname:string;
  f,ix:Integer;
  filenames:TSimpleStrings;
  tmp:TSimpleStrings;
  data:TDataManager;
  subs,divs:TFloats;


begin
  if TrainData<>'' then
    //get scaling vectors. Assumes all tables in training data have the same
    begin
    data:=TDataManager.Create;
    data.LoadFromFile(TrainData);
    WriteLn('Loaded ',Length(data.Tables),' tables of training data');
    subs:=Copy(data.Tables[0].Subtract,0,Length(data.Tables[0].Subtract));
    divs:=Copy(data.Tables[0].Divide,0,Length(data.Tables[0].Divide));
    data.Free;
    end;
  sl:=TStringList.Create;
  sl.LoadFromFile(FoldFileName);
  filenames:=nil;
  for f:=0 to sl.Count-1 do
    begin
    tmp:=SplitString(sl.Strings[f],#9);
    if Length(tmp)>=1 then
      AddToArray(IncludeTrailingPathDelimiter(SourceFolder)+tmp[0]+'.dat',filenames);
    end;
  sl.Free;

  data:=TDataManager.Create;
  data.LoadTables(filenames,'Target ID-Chain-Res Name-Res ID:Probe ID-Chain-Res Name-Res ID');
  ClassifyData(data,ClassCol);
  if TrainData='' then
    data.ScaleToMinusPlusOne
  else data.ScaleToVectors(subs,divs);
  data.SaveToFile(IncludeTrailingPathDelimiter(DestFolder)+'data.bin');
  data.Free;
end;

procedure TContactsApp.CreateHistograms(DataFile,HistFolder:string;NumFolds:Integer);
{ TODO : Only works for 2 classes, 0 and 1
creates folds in order of data files }

  procedure SaveHistSets(DataMan:TDataManager;DataCols:TSimpleStrings;
                        TableIxs:TIntegers;FileName:string);

  var
    datamat:TMatrix;
    f,dc:Integer;
    histset:THistMatrix;

  begin
    SetLength(histset,2,Length(DataCols));
    for dc:=0 to 1 do
      begin
      datamat:=DataMan.SelectData(DataCols,TableIxs,dc);
      WriteLn(Length(datamat[0]));
      for f:=0 to High(DataCols) do
        begin
        histset[dc,f]:=TFloatHistogram.Create(datamat[f]);
        TFloatHistogram(histset[dc,f]).MinVal:=-1.1;
        TFloatHistogram(histset[dc,f]).MaxVal:=1.1;
        histset[dc,f].CalcHistogram(-1,1);
        histset[dc,f].Normalize;
        histset[dc,f].ConvertToLog;
        WriteLn(dc,': Computed ',f,' of ',Length(DataCols),' with ',histset[dc,f].BinCount);
        end;
      end;
    SaveHistMatrix(histset,FileName);


    WriteLn('Saved to ',FileName);
  end;


var
  data:TDataManager;
  trainixs,testixs:TIntegers;
  f,g,foldsize:Integer;
  trainstr,teststr,fname:string;
  colnames:TSimpleStrings;


begin
  data:=TDataManager.Create;
  WriteLn('Loading');
  data.LoadFromFile(DataFile);
  WriteLn('Loaded');
  colnames:=data.Tables[0].GetColumnNames;
  for f:=0 to High(colnames) do WriteLn(colnames[f]);
  foldsize:=Round(Length(data.Tables)/NumFolds);
  for f:=0 to NumFolds-1 do
    begin
    trainixs:=nil;
    testixs:=nil;
    trainstr:='';
    teststr:='';

    for g:=0 to High(data.Tables) do
      if (Trunc(g/foldsize)=f) and (NumFolds>1) then
        begin
        AddToArray(g,testixs);
        teststr:=teststr+' '+IntToStr(g);
        end
      else
        begin
        AddToArray(g,trainixs);
        trainstr:=trainstr+' '+IntToStr(g);
        end;

    WriteLn('Fold '+IntToStr(f));
    WriteLn('Train:'+trainstr);
    WriteLn('Test:'+teststr);

    fname:='';
    if NumFolds=1 then
      fname:='0'
    else
      for g:=0 to NumFolds-1 do
        if g<>f then fname:=fname+IntToStr(g);

    SaveHistSets(data,colnames,trainixs,
                 IncludeTrailingPathDelimiter(HistFolder)+fname+'.hist');
    end;
end;

procedure TContactsApp.DisplayHistogramStats(HistFolder: string;
  NumFolds: Integer);

var
  f,g,h:Integer;
  fname:string;
  hists:THistMatrix;
begin
  for f:=0 to NumFolds-1 do
      begin
      fname:='';
      for g:=0 to NumFolds-1 do
        if g<>f then fname:=fname+IntToStr(g);
      LoadHistMatrix(hists,
        IncludeTrailingPathDelimiter(HistFolder)+fname+'.hist');
      WriteLn(hists[0,0].BinCount,' ',hists[1,0].BinCount);
      FreeHistMatrix(hists);
      end;
end;


procedure TContactsApp.NBFeatureSearchBeam(DataFile,HistFolder:string;
                              NumFolds,BeamWidth:Integer;
                              Mode,MaxMAP,HistBins,NumCycles:Integer;
                              ReportFile:string;
                              Disjoint:Boolean;Exclusions:TSimpleStrings);


//ATENTION: MAXMap is ignored, equals sum of true classes

{ TODO : Only works for 2 classes, 0 and 1
creates folds in order of data files }
{beam search, adapted to combinations:
  first, scan all features isolated, pick N best, keep features and score
  add N best to Selected
  for each step:
    empty Tried list
    for each S in Selected
      for each F in all features
        if sorted([S,F]) not in Tried
          add [S,F] to Tried, score
    Pick N best in Tried
    Run till end, keeping each set of Best and respective scores
}

  type
    TSelectedList=record
      List:array of TIntegers;
      Best:array of TIntegers;
      Hists:array of TIntegers;
      Scores:TFloats;
    end;

    TSelectedLists=array of TSelectedList;


  function GetFold(NumTables,Fold:Integer):TIntegers;

  var
    f:Integer;
    frac:TFloat;

  begin
    Result:=nil;
    frac:=NumTables/NumFolds;
    for f:=0 to NumTables-1 do
      if Trunc(f/frac)=Fold then
        AddToArray(f,Result);
  end;

var HistSets:array of THistMatrix;

  procedure LoadHistSets;

  var
    f,g:Integer;
    fname:string;

  begin
    SetLength(HistSets,NumFolds);
    for f:=0 to NumFolds-1 do
      begin
      fname:='';
      for g:=0 to NumFolds-1 do
        if g<>f then fname:=fname+IntToStr(g);
      LoadHistMatrix(HistSets[f],
        IncludeTrailingPathDelimiter(HistFolder)+fname+'.hist');
      end;
  end;

  function MeanAveragePrecision(Data:TDataManager;Columns:TSimpleStrings;
                                Feats:TIntegers;out Best,Hist:TIntegers):TFloat;
  var
    histmat:THistMatrix;


    procedure SetHists(TestTables:TIntegers;FoldIx:Integer;ColNames:TSimpleStrings);

      var
        classix,featix:Integer;
        f,tableix:Integer;
        tmpmat:TMatrix;
        averages:TFloats;
        trainset:TIntegers;

      begin
        if Mode=FMSelect then
          begin
          SetLength(histmat,2,Length(Feats));
          for classix:=0 to 1 do
            for featix:=0 to High(Feats) do
              histmat[classix,featix]:=HistSets[FoldIx,classix,Feats[featix]]
          end
        else
          begin
          SetLength(histmat,2,1);
          SetLength(trainset,Length(Data.Tables)-Length(TestTables));
          f:=0;
          for tableix:=0 to High(Data.Tables) do
            if not IsInArray(tableix,TestTables) then
              begin
              trainset[f]:=tableix;
              Inc(f);
              end;
          for classix:=0 to 1 do
            begin
            tmpmat:=Data.SelectData(ColNames,trainset,classix);
            averages:=FilledFloats(Length(tmpmat[0]),0);
            for featix:=0 to High(tmpmat) do
              averages:=Sum(averages,tmpmat[featix]);
            averages:=Multiply(averages,1/Length(tmpmat));
            histmat[classix,0]:=TFloatHistogram.Create(averages);
            TFloatHistogram(histmat[classix,0]).MaxVal:=1.1;
            TFloatHistogram(histmat[classix,0]).MinVal:=-1.1;
            histmat[classix,0].CalcHistogram(-1,1);
            histmat[classix,0].ConvertToLog;
            end;
          end;
      end;

  var
    f,count,countrel,foldix,classix,featix,tableix:Integer;
    sortedix,tables,singletable:TIntegers;
    scoremat:TMatrix;
    tmpvals:TFLoats;
    trueclasses,nbclasses,sortresults:TIntegers;
    averageprec,tmp:TFloat;
    scorediff:TFloats;
    colnames:TSimpleStrings;
    datmat:TMatrix;



  begin
    Best:=nil;
    Hist:=FilledInts(HistBins,0);
    SetLength(singletable,1);
    Result:=0;
    count:=0;
    SetLength(colnames,Length(Feats));
    for f:=0 to High(Feats) do colnames[f]:=Columns[Feats[f]];
    for foldix:=0 to NumFolds-1 do
      begin
      tables:=GetFold(Length(Data.Tables),foldix);
      SetHists(tables,foldix,colnames);
      //classify each complex
      for tableix:=0 to High(tables) do
        begin
        //get data and true classification
        singletable[0]:=tables[tableix];
        datmat:=Data.SelectData(colnames,singletable);
        if (Mode=FMAggregate) and (Length(datmat)>1) then
          begin
          for f:=1 to High(datmat) do
            datmat[0]:=Sum(datmat[0],datmat[f]);
          datmat[0]:=Multiply(datmat[0],1/Length(datmat));
          SetLength(datmat,1);
          end;
        trueclasses:=Data.Tables[tables[tableix]].DataClasses;
        //Naive Bayes
        nbclasses:=NBClassify(histmat,datmat,scoremat);
        //Compute MAP
        SetLength(scorediff,Length(nbclasses));
        for f:=0 to High(scorediff) do
          scorediff[f]:=scoremat[0,f]-scoremat[1,f];//lower is more towards class 1, for ascending sort
        sortedix:=QSAscendingIndex(scorediff);
        averageprec:=0;
        countrel:=0;
        Inc(count);
        {for f:=0 to High(trueclasses) do
          if trueclasses[sortedix[f]]>0 then
            begin
            if countrel=0 then AddToArray(f,Best);
            Inc(countrel);
            Inc(Hist[Trunc(f/Length(trueclasses)*HistBins)]);
            averageprec:=averageprec+countrel/(f+1);
            if countrel>=MaxMAP then Break;
            end;
        if countrel>0 then
          Result:=Result+averageprec/countrel;}
        {R Precision}
        MaxMAP:=Sum(trueclasses);
        for f:=0 to MaxMap-1 do
          if trueclasses[sortedix[f]]>0 then
            begin
            if countrel=0 then AddToArray(f,Best);
            Inc(countrel);
            Inc(Hist[Trunc(f/Length(trueclasses)*HistBins)]);
            end;
         if countrel=0 then
           for f:=MaxMap to High(sortedix) do
             if trueclasses[sortedix[f]]>0 then
               begin
               AddToArray(f,Best);
               countrel:=MaxMAP-f;
               break
               end;
         Result:=Result+countrel/MaxMAP;
        end;
      if Mode=FMAggregate then
        for f:=0 to 1 do histmat[f,0].Free;
      end;
    Result:=Result/count;
  end;

  function IsInSelection(List:TIntegers;const SelectedList:TSelectedList;LastIx:Integer):Boolean;

  var f,g:Integer;

  begin
    Result:=False;
    for f:=0 to LastIx do
     if IsEqual(List,SelectedList.List[f]) then
        begin
        Result:=True;
        Break;
        end;
  end;

  function HasElementsInSelection(List:TIntegers;const SelectedList:TSelectedList;LastIx:Integer):Boolean;

  var f,g,h:Integer;

  begin
    Result:=False;
    for f:=0 to LastIx do
     begin
     for g:=0 to High(List) do
      if IsInArray(List[g],SelectedList.List[f]) then
        begin
        Result:=True;
        Break;
        end;
     if Result then Break;
     end;
  end;


  procedure AddLists(Data:TDataManager;Columns:TSimpleStrings;Exclusion:TIntegers;
                     Base:TIntegers;var SelectedList:TSelectedList);

  var
    unsortedlist,sortedlist:TIntegers;
    f,col,rowix:Integer;

  begin
    SetLength(unsortedlist,Length(Base)+1);
    for f:=0 to High(Base) do unsortedlist[f]:=Base[f];
    SetLength(sortedlist,Length(unsortedlist));
    rowix:=Length(SelectedList.List);
    SetLength(SelectedList.List,rowix+Length(Columns));
    SetLength(SelectedList.Scores,rowix+Length(Columns));
    SetLength(SelectedList.Best,rowix+Length(Columns));
    SetLength(SelectedList.Hists,rowix+Length(Columns));
    for col:=0 to High(Columns) do
      begin
      if not IsInArray(col,Base) and not IsInArray(col,Exclusion) then
        begin
        unsortedlist[High(unsortedlist)]:=col;
        sortedlist:=QSSorted(unsortedlist);
        if not IsInSelection(sortedlist,SelectedList,rowix-1) then
          begin
          SelectedList.List[rowix]:=sortedlist;
          SelectedList.Scores[rowix]:=-MeanAveragePrecision(Data,Columns,sortedlist,
                                        SelectedList.Best[rowix],SelectedList.Hists[rowix]);
          inc(rowix);
          end;
        end;
      end;

    SetLength(SelectedList.List,rowix);
    SetLength(SelectedList.Scores,rowix);
    SetLength(SelectedList.Best,rowix);
    SetLength(SelectedList.Hists,rowix);

  end;


  function ExpandSelection(Data:TDataManager;Columns:TSimpleStrings;
                           SelectedList:TSelectedList;BeamWidth:Integer;
                           ExclusionList:TIntegers):TSelectedList;

  var
    tempselection:TSelectedList;
    ixs:TIntegers;
    f,currentix:Integer;
    actualwidth:Integer;

  begin
    tempselection.List:=nil;
    tempselection.Scores:=nil;

    for f:=0 to High(SelectedList.List) do
      begin
      writeln('add:',f);
      AddLists(Data,Columns,ExclusionList,SelectedList.List[f],tempselection);
      WriteLn(f,' ',Length(tempselection.List),' total');
      end;
    actualwidth:=Min(BeamWidth,Length(tempselection.List));
    ixs:=QSAscendingIndex(tempselection.Scores);
    SetLength(Result.List,actualwidth);
    SetLength(Result.Scores,actualwidth);
    SetLength(Result.Best,actualwidth);
    SetLength(Result.Hists,actualwidth);
    if not Disjoint then
      for f:=0 to actualwidth-1 do
        begin
        Result.Scores[f]:=-tempselection.Scores[ixs[f]];
        Result.List[f]:=tempselection.List[ixs[f]];
        Result.Best[f]:=tempselection.Best[ixs[f]];
        Result.Hists[f]:=tempselection.Hists[ixs[f]];
        end
    else
      begin
      currentix:=-1;
      for f:=0 to High(ixs) do
        begin
        if not HasElementsInSelection(tempselection.List[ixs[f]],
                  Result,currentix) then
          begin
          Inc(currentix);
          Result.Scores[currentix]:=-tempselection.Scores[ixs[f]];
          Result.List[currentix]:=tempselection.List[ixs[f]];
          Result.Best[currentix]:=tempselection.Best[ixs[f]];
          Result.Hists[currentix]:=tempselection.Hists[ixs[f]];
          if (currentix)>=actualwidth-1 then
            Break;
          end;
        end;

      actualwidth:=currentix+1;
      SetLength(Result.List,actualwidth);
      SetLength(Result.Scores,actualwidth);
      SetLength(Result.Best,actualwidth);
      SetLength(Result.Hists,actualwidth);
      end;

    WriteLn;
    WriteLn('Actual width:',actualwidth,' from ',Length(ixs));
    WriteLn('Top Score: ', Result.Scores[0]);
    for f:=0 to High(Result.List[0]) do
      WriteLn(Columns[Result.List[0,f]]);
    WriteLn;
  end;


var
  data:TDataManager;
  selectedlists:TSelectedLists;
  allcols:TSimpleStrings;
  feat:Integer;
  sl:TStringList;
  s:string;
  exclusionlist:TIntegers;

  procedure WriteReport;

  var
    f,g,h:Integer;
    fname:string;


  begin
    sl:=TStringList.Create;
    for f:=1 to feat do
    for h:=0 to High(SelectedLists[f].Scores) do
      begin
      s:=IntToStr(f)+#9+FloatToStr(selectedlists[f].Scores[h])+#9;
      for g:=0 to High(allcols) do
        begin
        if g<Length(selectedlists[f].List[h]) then
          s:=s+allcols[selectedlists[f].List[h,g]];
        s:=s+#9;
        end;
      for g:=0 to High(selectedlists[f].Best[h]) do
        s:=s+#9+IntToStr(selectedlists[f].Best[h,g]);
      for g:=0 to 10 do s:=s+#9;
      for g:=0 to High(selectedlists[f].Hists[h]) do
        s:=s+#9+IntToStr(selectedlists[f].Hists[h,g]);
      sl.Add(s);
      end;
    sl.SaveToFile(ReportFile);
    sl.Free;
  end;

var f,g:Integer;

begin
  data:=TDataManager.Create;
  WriteLn('Loading');
  data.LoadFromFile(DataFile);
  LoadHistSets;
  WriteLn('Loaded');
  allcols:=data.Tables[0].GetColumnNames;
  SetLength(selectedlists,NumCycles+1);


  //Initialize first to empty
  SetLength(selectedlists[0].List,1);
  SetLength(selectedlists[0].Scores,1);
  selectedlists[0].List[0]:=nil;
  selectedlists[0].Scores[0]:=0;
  exclusionlist:=data.Tables[0].GetColumnIndexes(Exclusions);

  for feat:=1 to NumCycles do
    begin
    writeln(feat);
    selectedlists[feat]:=ExpandSelection(data,allcols,selectedlists[feat-1],BeamWidth,exclusionlist);
    WriteReport;
    //Exclusion prevents the combination of previously selected features
    {if Disjoint then
      begin
      for f:=0 to High(selectedlists[feat].List) do
        for g:=0 to High(selectedlists[feat].List[f]) do
         AddToArray(selectedlists[feat].List[f,g],exclusionlist);
      end;}
    end;
  for f:=0 to High(HistSets[0]) do
    for g:=0 to High(HistSets[0,f]) do
      begin
      HistSets[0,f,g].Free;
      HistSets[1,f,g].Free;
      end;
  data.Free;


end;

procedure TContactsApp.NBEvaluate(HistFolder, DataFile, FeaturesFile,
  ReportFile: string; ContactsFolder: string);

var HistSets:array of THistMatrix;

  procedure LoadHistSets;

  var
    fname:string;

  begin
    SetLength(HistSets,1);
    LoadHistMatrix(HistSets[0],
        IncludeTrailingPathDelimiter(HistFolder)+'0.hist');
  end;

  function RPrecision(Data:TDataManager;ColNames:TSimpleStrings;
                                out Best:TIntegers; WriteContacts:Boolean):TFloat;
  var
    histmat:THistMatrix;


    procedure SetHists(Feats:TIntegers);

      var
        classix,featix:Integer;

      begin
        SetLength(histmat,2,Length(ColNames));
        for classix:=0 to 1 do
          for featix:=0 to High(Feats) do
            histmat[classix,featix]:=HistSets[0,classix,Feats[featix]]
      end;

    procedure ExportContacts(TableIx:Integer;SortedIxs:TIntegers;Scores:TFloats);

    var
      f:Integer;
      sl:TStringList;
      s:string;
      labels:TSimpleStrings;
      classes:TIntegers;
      seps,sepst,sepsp:TSimpleStrings;

    begin
      labels:=Data.Tables[TableIx].GetLabels;
      classes:=Data.Tables[tableix].DataClasses;
      sl:=TStringList.Create;
      for f:=0 to High(SortedIxs) do
       begin
       seps:=SplitString(labels[SortedIxs[f]],':');
       sepst:=SplitString(seps[0],'-');
       sepsp:=SplitString(seps[1],'-');
       s:=sepst[High(sepst)-2]+' '+sepst[High(sepst)]+' '+
          sepsp[High(sepsp)-2]+' '+sepsp[High(sepsp)]+' '+
          ' 5 200 200';
       sl.Add(s+' '+labels[SortedIxs[f]]+' '+FloatToStrF(Scores[SortedIxs[f]],ffFixed,0,3)
              +' '+IntToStr(classes[SortedIxs[f]]));

       end;

      sl.SaveToFile(IncludeTrailingPathDelimiter(ContactsFolder)+
        ExtractFileName(Data.Tables[TableIx].Name));
      sl.Free;
    end;

  var
    f,count,countrel,r,classix,featix,tableix:Integer;
    sortedix,tables,singletable:TIntegers;
    scoremat:TMatrix;
    tmpvals:TFLoats;
    trueclasses,nbclasses,sortresults:TIntegers;
    averageprec,tmp:TFloat;
    scorediff:TFloats;
    datmat:TMatrix;



  begin
    Best:=nil;
    SetLength(singletable,1);
    Result:=0;
    count:=0;
    SetHists(data.Tables[0].GetColumnIndexes(ColNames));
    //classify each complex
    for tableix:=0 to High(data.Tables) do
      begin
      //get data and true classification
      singletable[0]:=tableix;
      datmat:=Data.SelectData(ColNames,singletable);
      trueclasses:=Data.Tables[tableix].DataClasses;
      //Naive Bayes
      nbclasses:=NBClassify(histmat,datmat,scoremat);

      SetLength(scorediff,Length(nbclasses));
      for f:=0 to High(scorediff) do
        scorediff[f]:=scoremat[0,f]-scoremat[1,f];//lower is more towards class 1, for ascending sort
      sortedix:=QSAscendingIndex(scorediff);

      averageprec:=0;
      countrel:=0;
      Inc(count);

      {R Precision}
      r:=Sum(trueclasses);
      for f:=0 to r-1 do
        if trueclasses[sortedix[f]]>0 then
          begin
          if countrel=0 then AddToArray(f,Best);
          Inc(countrel);
          end;
       if countrel=0 then
         for f:=r to High(sortedix) do
           if trueclasses[sortedix[f]]>0 then
             begin
             AddToArray(f,Best);
             countrel:=r-f;
             break
             end;
       Result:=Result+countrel/r;

       if WriteContacts then
        ExportContacts(tableix,sortedix,scorediff);


      end;
    Result:=Result/count;
  end;



var
  data:TDataManager;
  cols:array of TSimpleStrings;
  sl:TStringList;
  s:string;
  scores:TFloats;
  bests:array of TIntegers;
  count,f,colcount,g:Integer;

begin
  data:=TDataManager.Create;
  WriteLn('Loading');
  data.LoadFromFile(DataFile);
  colcount:=Length(data.Tables[0].Columns);
  LoadHistSets;
  WriteLn('Loaded');
  sl:=TStringList.Create;
  sl.LoadFromFile(FeaturesFile);
  count:=0;
  for f:=0 to sl.Count-1 do
    if sl.Strings[f]<>'' then Inc(count);
  SetLength(cols,count);
  SetLength(scores,count);
  SetLength(bests,count);
  for f:=0 to High(cols) do
    begin
    cols[f]:=nil;
    scores[f]:=0;
    bests[f]:=nil;
    end;
  count:=0;
  for f:=0 to sl.Count-1 do
    if sl.Strings[f]<>'' then
      begin
      cols[count]:=SplitString(TrimmedBlanks(sl.Strings[f]),#9);
      Inc(count);
      end;
  for f:=0 to High(cols[0]) do WriteLn(cols[0,f]);
  sl.Clear;

  //Add headers
  s:='';
  for g:=0 to colcount do s:=s+#9;
  for g:=0 to High(data.Tables) do
    s:=s+#9+ExtractFileName(data.Tables[g].Name);
  sl.Add(s);

  //Add total count
  s:='';
  for g:=0 to colcount do s:=s+#9;
  for g:=0 to High(data.Tables) do
    s:=s+#9+IntToStr(Length(data.Tables[g].DataClasses));
  sl.Add(s);

  //Add true counts
  s:='';
  for g:=0 to colcount do s:=s+#9;
  for g:=0 to High(data.Tables) do
    s:=s+#9+IntToStr(Sum(data.Tables[g].DataClasses));
  sl.Add(s);

  for f:=0 to High(cols) do
    begin
    writeln(f);
    scores[f]:=RPrecision(data,cols[f],bests[f],
                (ContactsFolder<>'') and (f=0));
    writeln(scores[f]);
    s:=IntToStr(Length(cols[f]))+#9+FloatToStr(scores[f]);
    for g:=0 to colcount-1 do
      if g<=High(cols[f]) then s:=s+#9+cols[f,g]
      else s:=s+#9;
    for g:=0 to High(bests[f]) do s:=s+#9+IntToStr(bests[f,g]);
    sl.Add(s);
    end;
  sl.SaveToFile(ReportFile);
  sl.Free;
  data.Free;
  for f:=0 to High(HistSets) do
    FreeHistMatrix(HistSets[f]);
end;

procedure TContactsApp.SelectAggregates(InFile, OutFile: string);

var
  aggregates:array of TSimpleStrings;

  function IsInAggregate(ss:TSimpleStrings):Boolean;

  var
    f,g:Integer;

  begin
    Result:=False;
    for f:=0 to High(aggregates) do
      begin
      for g:=0 to High(ss) do
        if IsInArray(ss[g],aggregates[f]) then
          begin
          Result:=True;
          Break;
          end;
      if Result then Break;
      end;
  end;

var
  sl1,sl2:TStringList;
  f,g:Integer;
  s:string;
  spl,spl2:TSimpleStrings;


begin
  aggregates:=nil;
  sl1:=TStringList.Create;
  sl2:=TStringList.Create;
  sl1.LoadFromFile(InFile);
  for f:=0 to sl1.Count-1 do
    begin
    spl:=SplitString(sl1.Strings[f],#9);
    spl2:=nil;
    for g:=1 to High(spl) do
      if spl[g]<>'' then
        AddToArray(spl[g],spl2);
    if not IsInAggregate(spl2) then
      begin
      SetLength(aggregates,Length(aggregates)+1);
      aggregates[High(aggregates)]:=Copy(spl2,0,Length(spl2));
      if Length(spl2)>0 then
        begin
        s:=spl[0];
        for g:=0 to High(spl2) do
          s:=s+#9+spl2[g];
        sl2.Add(s);
        end;
      end;
    end;
  sl2.SaveToFile(OutFile);
  sl1.Free;
  sl2.Free;
end;

procedure TContactsApp.Aggregate(DataFile, AggregateFile, OutFile: string;
                         DropColumns:Boolean);

var
  data:TDataManager;
  sl:TStringList;
  f,tab,col,g:Integer;
  spl,toaggregate,todrop:TSimpleStrings;
  tmp,averages:TFloats;
  s:string;
  newcol:TFloatDataColumn;

begin
  todrop:=nil;
  WriteLn('Loading');
  data:=TDataManager.Create;
  data.LoadFromFile(DataFile);
  sl:=TStringList.Create;
  sl.LoadFromFile(AggregateFile);
  for f:=0 to sl.Count-1 do
    begin
    WriteLn('Adding ',f+1,' of ',sl.Count);
    spl:=SplitString(sl.Strings[f],#9);
    toaggregate:=nil;
    for g:=1 to High(spl) do
      AddToArray(spl[g],toaggregate);
    if Length(toaggregate)>1 then
      for tab:=0 to High(data.Tables) do
        begin
        tmp:=data.Tables[tab].GetColumn(toaggregate[0]);
        s:=toaggregate[0];
        if tab=0 then AddToArray(toaggregate[0],todrop);
        averages:=Copy(tmp,0,Length(tmp));
        for col:=1 to High(toaggregate) do
          begin
          tmp:=data.Tables[tab].GetColumn(toaggregate[col]);
          averages:=Sum(averages,tmp);
          s:=s+'+'+toaggregate[col];
          if tab=0 then AddToArray(toaggregate[col],todrop);
          end;
        averages:=Multiply(averages,1/Length(toaggregate));
        newcol:=TFloatDataColumn.Create;
        newcol.Vals:=Copy(averages,0,Length(averages));
        newcol.Name:=s;
        data.Tables[tab].AddColumn(newcol,0,1);
        end;
    end;
  if DropColumns then
    begin
    data.DropColumns(todrop);
    Writeln('Drop ',Length(todrop),' remaining ',Length(data.Tables[0].Columns));
    end;
  Writeln('Saving');
  data.SaveToFile(OutFile);
  data.Free;
  sl.Free;
end;


var
  Application: TContactsApp;
begin
  Application:=TContactsApp.Create(nil);
  Application.Title:='Contacts Classifier';
  Application.Run;
  Application.Free;
end.

