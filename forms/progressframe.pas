{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 9.1.2011
Purpose:
  Displays progress reports

Requirements:
Revisions:
To do:
  Currently only working on first task on the list. Change ProgressCallback
*******************************************************************************}

unit progressframe;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, FileUtil, LResources, Forms,progress, ComCtrls, StdCtrls;

type

  { TProgressFrame }

  TProgressFrame = class(TFrame)
    CancelBt: TButton;
    TaskTimeLbl: TLabel;
    TaskTitleLbl: TLabel;
    TaskPBar: TProgressBar;
  private
    { private declarations }
    FRefreshWait:Real; //time to wait for next refresh, in days
    FLastRefresh:Real;
    FLastTaskCount:Integer;
    procedure ProgressCallBack;
    procedure SetHook;
  public
    { public declarations }
    constructor Create(AOwner: TComponent); override;
    procedure Free;
  end; 

implementation

{ TProgressFrame }

procedure TProgressFrame.ProgressCallBack;

var
  n:Real;
  tasks:TTasks;

begin
  n:=Now;
  tasks:=ListTasks;
  if (Length(tasks)<>FLastTaskCount) or (n-FLastRefresh>FRefreshWait) then
    begin
    FLastRefresh:=n;
    FLastTaskCount:=Length(tasks);
    tasks:=ListTasks;
    //TODO: add support for multiple tasks
    if tasks<>nil then
      begin
      TaskTitleLbl.Caption:=tasks[0].Title;
      TaskPBar.Position:=Round(tasks[0].Progress*TaskPBar.Max);
      CancelBt.Enabled:=tasks[0].CanCancel;
      end
    else
      begin
      TaskTitleLbl.Caption:='No tasks running';
      TaskPBar.Position:=0;
      CancelBt.Enabled:=False;
      end;
    Application.ProcessMessages;
    end;
end;

constructor TProgressFrame.Create(AOwner: TComponent);
begin
  inherited Create(AOwner);
  SetHook;
  FRefreshWait:=1/24/3600/2; // half a second;
  FLastRefresh:=0;
  FLastTaskCount:=-1;
end;

procedure TProgressFrame.Free;
begin
  if Self<>nil then
    RemoveObjectHook(@ProgressCallback);
end;

procedure TProgressFrame.SetHook;
begin
  AddObjectHook(@ProgressCallBack);
end;

initialization
  {$I progressframe.lrs}

end.

