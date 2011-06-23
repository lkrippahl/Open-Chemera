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
  msablocksparser
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
  else if command='M' then
    begin
    infile:=ParamStr(2);
    blockfile:=ParamStr(3);
    outfile:=ParamStr(4);
    MSA:=ReadClustal(infile);
    MergeBlocks(MSA,blockfile);
    SaveMSF(outfile,MSA);
    WriteLn('Merge OK');

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
  WriteLn('Params: infile.aln outfile.msf (MSA formats are determined by the extension)');
  WriteLn('Command: m (merge blocks into a .msf file)');
  WriteLn('Params: infile.aln blockfile outfile.msf');

end;

var
  Application: TMSAParser;

{$R *.res}

begin
  Application:=TMSAParser.Create(nil);
  Application.Run;
  Application.Free;
end.

