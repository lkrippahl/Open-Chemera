{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 9.1.2011
Purpose:
  Reports progress of running tasks. Should be used by all units running
  lengthy processes to help display progress and cancel tasks.

  All memory management of taks is run in this unit. Units outside should not use
  constructors from these classes, nor free the objects.

  All event processing, such as Application.Processmessages, must be called from
  the registered hooks. They are not called from this unit.

Requirements:
Revisions:
To do:
  SinglesToMSA is ignoring the TrimToQuery parameter. Rewrite
  SumOfPairScores is not calculating gap scores
*******************************************************************************}

unit progress;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils; 

type

  { TRunningTask }
  TProgressHook_Object=procedure of object;
  TOHooks=array of TProgressHook_Object;
  TProgressHook=procedure;
  THooks=array of TProgressHook;

  TRunningTask=class
  private
    FCanCancel:Boolean;
    FTitle:string;
    FStartTime:TDateTime;
    FProgress:Real; // from 0 to 1
  public
    AbortRequested:Boolean;
    property CanCancel:Boolean read FCanCancel;
    property Progress:Real read FProgress;
    property Title:string read FTitle;
    constructor Create(ACancel:Boolean;ATitle:string);
    procedure Step(s:Real);
      // Step calls ll hooks
  end;
  TTasks=array of TRunningTask;

  procedure AddObjectHook(hook:TProgressHook_Object);
    // for hooks that are methods of objects
  procedure AddHook(hook:TProgressHook);
    // for other hooks
  function NewTask(Cancancel:Boolean;Title:string):TRunningTask;
  procedure FreeTask(task:TRunningTask);
  procedure CallHooks;
  procedure RemoveObjectHook(hook:TProgressHook_Object);
  procedure RemoveHook(hook:TProgressHook);
  function ListTasks:TTasks;

implementation

var
  RunningTasks:TTasks;
  Hooks:THooks;
  ObjectHooks:TOHooks;

function TaskIndex(Task:TRunningTask):Integer;

begin
  Result:=High(RunningTasks);
  while (Result>=0) and (RunningTasks[Result]<>Task) do
    Dec(Result);
end;

procedure AddObjectHook(Hook: TProgressHook_Object);
begin
  SetLength(ObjectHooks,Length(ObjectHooks)+1);
  ObjectHooks[High(ObjectHooks)]:=Hook;
end;

procedure AddHook(Hook: TProgressHook);
begin
  SetLength(Hooks,Length(ObjectHooks)+1);
  Hooks[High(Hooks)]:=Hook;
end;

function NewTask(CanCancel: Boolean; Title: string): TRunningTask;
begin
  Result:=TRunningTask.Create(CanCancel,Title);
  SetLength(RunningTasks,Length(RunningTasks)+1);
  RunningTasks[High(RunningTasks)]:=Result;
  CallHooks;
end;

procedure FreeTask(task: TRunningTask);

var ix:Integer;

begin
  ix:=TaskIndex(task);
  if ix>=0 then
    begin
    RunningTasks[ix].Free;
    RunningTasks[ix]:=RunningTasks[High(RunningTasks)];
    SetLength(RunningTasks,Length(RunningTasks)-1);
    end;
  CallHooks;
end;

procedure CallHooks;

var f:Integer;

begin
  for f:=0 to High(Hooks) do
    Hooks[f];
  for f:=0 to High(ObjectHooks) do
    ObjectHooks[f];
end;

procedure RemoveObjectHook(hook: TProgressHook_Object);

var f:Integer;

begin
  for f:=0 to High(ObjectHooks) do
    if ObjectHooks[f]=hook then
      begin
      ObjectHooks[f]:=ObjectHooks[High(ObjectHooks)];
      SetLength(ObjectHooks,Length(ObjectHooks)-1);
      Break;
      end;
end;

procedure RemoveHook(hook: TProgressHook);

var f:Integer;

begin
  for f:=0 to High(Hooks) do
    if Hooks[f]=hook then
      begin
      Hooks[f]:=Hooks[High(Hooks)];
      SetLength(Hooks,Length(Hooks)-1);
      Break;
      end;
end;

function ListTasks: TTasks;
begin
  Result:=RunningTasks;
end;

{ TRunningTask }

constructor TRunningTask.Create(ACancel: Boolean; ATitle: string);
begin
  inherited Create;
  FCanCancel:=ACancel;
  FTitle:=ATitle;
  FStartTime:=Now;
  FProgress:=0;
  AbortRequested:=False;
end;

procedure TRunningTask.Step(s: Real);
begin
  FProgress:=FProgress+s;
  CallHooks;
end;

end.

