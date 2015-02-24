{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 1.7.2011
Purpose:
  BiGGER main program file

  -j ordersfile.xml
     mandatory, job file

  -e exportfile.txt
     optional: exports grids and terminates execution.

Requirements:
Revisions:
To do:
*******************************************************************************}

program bigger;

{$mode objfpc}{$H+}

uses
  {$IFDEF UNIX}{$IFDEF UseCThreads}
  cthreads,
  {$ENDIF}{$ENDIF}
  Classes, SysUtils, geomhash, CustApp, linegrids, dockdomains, geomutils,
  basetypes, dockconstraints, bogie, surface, docktasks, biggerman,
  oclconfiguration, alignment, rmsd, cgradient, base3ddisplay, progress;

type

  { TBigger }

  TBigger = class(TCustomApplication)
  protected
    procedure DoRun; override;
  public
    constructor Create(TheOwner: TComponent); override;
    destructor Destroy; override;
    procedure WriteHelp; virtual;
  end;

{ TBigger }

procedure TBigger.DoRun;
var
  biggerman:TBiGGERManager;
  jobfile:string;

begin
  DecimalSeparator:='.';
  { TODO : Find a better solution for the decimal separator thing }

  // parse parameters
  if HasOption('h','help') then begin
    WriteHelp;
    Terminate;
    Exit;
  end;

  { add your program here }

  jobfile:=GetOptionValue('j');
  if jobfile<>'' then
    begin
    biggerman:=TBiGGERManager.Create;
    biggerman.LoadJobs(jobfile);
    if HasOption('e') then
      biggerman.ExportGrids(GetOptionValue('e'))
    else
      begin
      biggerman.RunJobs;
      biggerman.SaveCurrent(jobfile);
      biggerman.Free;
      end;
    end;

  // stop program loop
  Terminate;
end;

constructor TBigger.Create(TheOwner: TComponent);
begin
  inherited Create(TheOwner);
  StopOnException:=True;
end;

destructor TBigger.Destroy;
begin
  inherited Destroy;
end;

procedure TBigger.WriteHelp;
begin
  { add your help code here }
  writeln('Usage: ',ExeName,' -h');
end;

var
  Application: TBigger;

{$R *.res}

begin
  Application:=TBigger.Create(nil);
  Application.Run;
  Application.Free;
end.
