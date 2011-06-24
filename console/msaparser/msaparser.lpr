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
  msablocksparser, msaconstraints
  { you can add units after this };

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
var
  infile,outfile, blockfile,command:string;
  MSA:TMSA;
  cpmsa:TCPMSA;

begin
  // Check arguments
  if HasOption('h','help') or (ParamCount<1) then begin
    WriteHelp;
    Terminate;
    Exit;
  end;

  command:=Uppercase(ParamStr(1));
  WriteLn(command);
  if command='C' then
    begin
    infile:=ParamStr(2);
    outfile:=ParamStr(3);
    MSA:=ReadClustal(infile);
    SaveMSF(outfile,MSA);
    end
  else if command='H' then
    begin
    infile:=ParamStr(2);
    outfile:=ParamStr(3);
    MSA:=ReadMSF(infile);
    SaveMSF(outfile,MSA);
    end
  else if command='M' then
    begin
    infile:=ParamStr(2);
    blockfile:=ParamStr(3);
    outfile:=ParamStr(4);
    MSA:=ReadMSF(infile);
    MergeBlocks(MSA,blockfile);
    SaveMSF(outfile,MSA);
    WriteLn('Merge OK');
    end
  else if command='CP' then
    begin
    cpmsa:=TCPMSA.Create;
    infile:=ParamStr(2);
    outfile:=ParamStr(3);
    MSA:=ReadMSF(infile);
    cpmsa.AssignMSA(MSA);
    if ParamStr(4)<>'' then
      cpmsa.MinShiftable:=StrToInt(ParamStr(4));
    cpmsa.GenerateBlocks(outfile+'.blocks');
    SaveMSAToFile(MSA,outfile+'.msa');
    cpmsa.Free;
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

