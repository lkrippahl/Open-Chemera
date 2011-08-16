{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 23.6.2011
Purpose:
  Console application for parsing MSA files
Requirements:
Revisions:
To do:
  Work in progress...
*******************************************************************************}

program msaparser;

{$mode objfpc}{$H+}

uses
  {$IFDEF UNIX}{$IFDEF UseCThreads}
  cthreads,
  {$ENDIF}{$ENDIF}
  Classes, SysUtils, msfparser, clustalparser, alignment, CustApp, debugutils,
  msablocksparser, msaconstraints, process, basetypes, stringutils;


type

  { TMSAParser }

  TMSAParser = class(TCustomApplication)
  protected
    procedure DoRun; override;
  public
    constructor Create(TheOwner: TComponent); override;
    destructor Destroy; override;
    procedure WriteHelp; virtual;
  end;

{ TMSAParser }

procedure TMSAParser.DoRun;

function BaliScores(Reference,Test:string;var Report:string):TFloats;

const
  BUFFER=2048;

var
  process:TProcess;
  sl:TStringList;
  m:TMemoryStream;
  f,bytesread,count:Integer;
  formatsettings:TFormatSettings;
  tmp:TFLoat;

begin
  sl:=TStringList.Create;
  m:=TMemoryStream.Create;
  process:=TProcess.Create(nil);
  formatsettings.DecimalSeparator:='.';

  process.CommandLine:='bali_score.exe '+Reference+' '+Test;
  process.Options:=process.Options+[poUsePipes];
  process.Execute;
  bytesread:=0;
  while process.Running do
    begin
    m.SetSize(bytesread+BUFFER);
    count:=process.Output.Read((m.Memory+bytesread)^,BUFFER);
    Inc(bytesread,count);
    m.SetSize(bytesread);
    Sleep(100);
    end;

  m.SetSize(bytesread+BUFFER);
  count:=process.Output.Read((m.Memory+bytesread)^,BUFFER);
  Inc(bytesread,count);
  m.SetSize(bytesread);

  sl.LoadFromStream(m);
  Result:=nil;
  //for f:=0 to sl.Count-1 do WriteLn(sl.Strings[f]);
  for f:=4 to sl.Count-1 do
    begin
    if Pos('SP',sl.Strings[f])>0 then
      begin
      //WriteLn('SP');
      tmp:=StrToFloat(Copy(sl.Strings[f],Pos('=',sl.Strings[f])+2,Length(sl.Strings[f])),formatsettings);
      Report:=Report+#9+FloatToStr(tmp);
      Break;
      end;
    AddToArray(StrToFloat(sl.Strings[f],formatsettings),Result);
    //Writeln(StrToFloat(sl.Strings[f]));
    end;
  process.Free;
  sl.Free;
  m.Free;
end;

procedure RunOne(BaseName:string);

var
  scorealn,scorecp:TFloats;
  blocks:TBlockSet;
  msa:TMSA;
  totscorecp,totscorealn:Single;
  f,count:Integer;
  reportstr:string;

begin
  //fix msf header
  if not FileExists(BaseName+'-head.msf') then
    begin
    msa:=ReadMSF(BaseName+'.msf');
    SaveMSF(BaseName+'-head.msf',MSA);
    end;
  //create msf from the clustal
  if not FileExists(BaseName+'-aln.msf') then
    begin
    msa:=ReadClustal(BaseName+'.aln');
    SaveMSF(BaseName+'-aln.msf',msa);
    end;
  //create cp fixed .msf from clustal
  if not FileExists(BaseName+'-aln-cp.msf') then
    begin
    msa:=ReadClustal(BaseName+'.aln');
    blocks:=LoadBlocks(BaseName+'.align');
    MergeBlocks(msa,blocks);
    SaveMSF(BaseName+'-aln-cp.msf',msa);
    end;
  Sleep(100);  //is this necessary? TO DO: test it
  reportstr:=BaseName;

  //score relative to bali msf
  scorealn:=BaliScores(BaseName+'-head.msf',BaseName+'-aln.msf',reportstr);
  scorecp:=BaliScores(BaseName+'-head.msf',BaseName+'-aln-cp.msf',reportstr);
  totscorecp:=0;
  totscorealn:=0;
  count:=0;
  for f:=0 to High(scorealn) do
    if scorealn[f]<>scorecp[f] then
      begin
      inc(count);
      totscorecp:=totscorecp+scorecp[f];
      totscorealn:=totscorealn+scorealn[f];
      //WriteLn(IntToStr(f)+':'+FloatToStr(scorealn[f])+'->'+FloatToStr(scorecp[f]));
      end;
  reportstr:=reportstr+#9+IntToStr(count)+#9+FloatToStr(totscorecp/totscorealn);
  WriteLn(reportstr);
end;

function ImprovementSum(FileName:string):Single;

var
  sl:TStringList;
  f:Integer;
  s:string;
  formatsettings:TFormatSettings;

begin
  formatsettings.DecimalSeparator:='.';
  sl:=TStringList.Create;
  sl.LoadFromFile(FileName);
  Result:=0;
  for f:=0 to sl.Count-1 do
    begin
    s:=sl.Strings[f];
    if Pos('# Total Scores (per sequence pair)    :',s)>0 then
      begin
      Delete(s,1,Pos(':',s));
      //subtract original score
      Result:=Result-StrToFloat(Deblank(Copy(s,1,Pos('|',s)-1)),formatsettings);
      Delete(s,1,Pos('|',s));
      //add corrected score
      Result:=Result+StrToFloat(Deblank(Copy(s,1,Pos('|',s)-1)),formatsettings);
      end;
    end;
  sl.Free;
end;

procedure CopyBest(FileName:string);

function TotPBScore(Sl:TStringList):Integer;

var
  f:Integer;
  s:string;

begin
  Result:=0;
  for f:=0 to Sl.Count-1 do
    begin
    s:=Sl.Strings[f];
    if Pos('best=',s)>0 then
      begin
      Delete(s,1,Pos('init:',s)+4);
      Result:=Result-StrToInt(Copy(s,1,Pos(';',s)-1));
      Delete(s,1,Pos('best=',s)+4);
      Result:=Result+StrToInt(Copy(s,1,Pos(';',s)-1));
      end;
    end;
end;


var
  sl,slkeep:TStringList;
  best:Integer;


begin
  sl:=TStringList.Create;
  slkeep:=TStringList.Create;
  slkeep.LoadFromFile('test1\'+FileName);
  sl.LoadFromFile('test2\'+FileName);
  if TotPBScore(sl)>TotPBScore(slkeep) then
    slkeep.Assign(sl);
  sl.LoadFromFile('test3\'+FileName);
  if TotPBScore(sl)>TotPBScore(slkeep) then
    slkeep.Assign(sl);
  Delete(FileName,1,2);
  slkeep.SaveToFile(FileName);
  WriteLn(Copy(FileName,1,Pos('.',FileName)-1),#9,TotPBScore(slkeep));
  sl.Free;
  slkeep.Free;
end;

var
  infile, outfile, blockfile,command:string;
  f:Integer;
  msa,msacp:TMSA;
  cpmsa:TCPMSA;
  blocks:TBlockSet;
  srec:TSearchRec;
  s:string;

 begin
  // Check arguments
  if HasOption('h','help') or (ParamCount<1) then begin
    WriteHelp;
    Terminate;
    Exit;
  end;

  command:=Uppercase(ParamStr(1));
  if command='C' then
    begin
    infile:=ParamStr(2);
    outfile:=ParamStr(3);
    msa:=ReadClustal(infile);
    SaveMSF(outfile,msa);
    end
  else if command='H' then
    begin
    infile:=ParamStr(2);
    outfile:=ParamStr(3);
    msa:=ReadMSF(infile);
    SaveMSF(outfile,msa);
    end
  else if command='M' then
    begin
    infile:=ParamStr(2);
    blockfile:=ParamStr(3);
    outfile:=ParamStr(4);
    msa:=ReadMSF(infile);
    blocks:=LoadBlocks(blockfile);
    MergeBlocks(MSA,blocks);
    SaveMSF(outfile,msa);
    WriteLn('Merge OK');
    end
  else if command='CP' then
    begin
    cpmsa:=TCPMSA.Create;
    infile:=ParamStr(2);
    outfile:=ParamStr(3);
    msa:=ReadMSF(infile);
    cpmsa.AssignMSA(MSA);
    if ParamStr(4)<>'' then
      cpmsa.MinShiftable:=StrToInt(ParamStr(4));
    cpmsa.GenerateBlocks(outfile+'.blocks');
    SaveMSAToFile(MSA,outfile+'.msa');
    cpmsa.Free;
    end
  else if command='LIST' then
    begin
    cpmsa:=TCPMSA.Create;
    if FindFirst('*.blocks',faAnyFile,srec)=0 then
      repeat
      blocks:=LoadBlocks(srec.Name);
      for f:=0 to High(blocks.Blocks) do
        begin
        if f=0 then s:='| tee '
          else s:='| tee -a ';
        WriteLn('./main '+srec.name+' '+ChangeFileExt(srec.Name,'.msa')+
          ' Gonnet.const '+IntToStr(f)+' 0 - 500 '+s+ChangeFileExt(srec.Name,'.align'))
        end;
      until FindNext(srec)<>0;
    cpmsa.Free;
    end
  else if command='CPALL' then
    begin
    cpmsa:=TCPMSA.Create;
    if FindFirst('*.aln',faAnyFile,srec)=0 then
      repeat
       WriteLn(srec.Name);
      infile:=srec.Name;
      outfile:=ChangeFileExt(srec.Name,'');
      msa:=ReadClustal(infile);
      cpmsa.AssignMSA(MSA);
      cpmsa.GenerateBlocks(outfile+'.blocks');
      SaveMSAToFile(MSA,outfile+'.msa');
      until FindNext(srec)<>0;
    cpmsa.Free;
    end
  else if command='BALI' then
    begin
    if FindFirst('*.align',faAnyFile,srec)=0 then
      repeat
      if FileExists(ChangeFileExt(srec.Name,'.msf')) and
          FileExists(ChangeFileExt(srec.Name,'.aln')) then
        RunOne(ChangeFileExt(srec.Name,''));
      until FindNext(srec)<>0;
      FindClose(srec);
    end
  else if command='IMPROV' then
  //only for Marco, # commented blocks file
    begin
    if FindFirst('*.align',faAnyFile,srec)=0 then
      repeat
      WriteLn(FloatToStr(ImprovementSum(srec.Name)));
      until FindNext(srec)<>0;
      FindClose(srec);
    end
  else if command='PB' then
  //only for PB align files, to pick best
    begin
    if FindFirst('test1\*.align',faAnyFile,srec)=0 then
      repeat
      CopyBest(srec.Name);
      until FindNext(srec)<>0;
      FindClose(srec);
    end;

  // stop program loop
  Terminate;
end;

constructor TMSAParser.Create(TheOwner: TComponent);
begin
  inherited Create(TheOwner);
  StopOnException:=True;
end;

destructor TMSAParser.Destroy;
begin
  inherited Destroy;
end;

procedure TMSAParser.WriteHelp;
begin
  { add your help code here }
  WriteLn('Usage: ',ExtractFileName(ExeName),'command params');
  WriteLn('Command: c (convert MSA files)');
  WriteLn('Params: infile.aln outfile.msf ');
  WriteLn('Command: h (add header to msf file)');
  WriteLn('Params: infile.msf outfile.msf');
  WriteLn('Command: m (merge blocks into a .msf file)');
  WriteLn('Params: infile.msf blockfile outfile.msf');
  WriteLn('Command: cp (create blocks and alignment file)');
  WriteLn('Params: infile.msf outfile minshiftable');
  WriteLn('(outfile has no extension, will be used for .blocs and .msa;');
  WriteLn('minshiftable is the minimum number of lines that can be shifted,');
  WriteLn('default is 20% of number of sequences)');

end;

var
  Application: TMSAParser;

{$R *.res}

begin
  Application:=TMSAParser.Create(nil);
  Application.Run;
  Application.Free;
end.

