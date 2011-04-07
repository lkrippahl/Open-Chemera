unit progress;

{
  to report progress of running tasks. Should be used by all units running
  lengthy processes to help display progress and cancel

  All memory management of taks is run in this unit. Units outside should not use
  constructors from these classes, nor free the objects.

  All event processing, such as Application.Processmessages, must be called from
  the registered hooks. Not called from this unit
}


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
    FAborted:Boolean;
    FTitle:string;
    FStartTime:TDateTime;
    FProgress:Real; // from 0 to 1
    {Not sure if abort hooks are useful...
    FAbortHook:TProgressHook;
    FAbortOHook:TProgressHook_Object;}
  public
    property CanCancel:Boolean read FCanCancel;
    property Progress:Real read FProgress;
    property Title:string read FTitle;
    constructor Create(acancel:Boolean;atitle:string);
    procedure Abort;
    procedure Step(s:Real);
      // StepBy also calls CallHooks to call all hooks
  end;
  TTasks=array of TRunningTask;

  procedure AddObjectHook(hook:TProgressHook_Object);
    // for hooks inside objects
  procedure AddHook(hook:TProgressHook);
    // for hooks outside objects
  function NewTask(cancancel:Boolean;title:string):TRunningTask;
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

function TaskIndex(task:TRunningTask):Integer;

begin
  Result:=High(RunningTasks);
  while (Result>=0) and (RunningTasks[Result]<>task) do
    Dec(Result);
end;

procedure AddObjectHook(hook: TProgressHook_Object);
begin
  SetLength(ObjectHooks,Length(ObjectHooks)+1);
  ObjectHooks[High(ObjectHooks)]:=hook;
end;

procedure AddHook(hook: TProgressHook);
begin
  SetLength(Hooks,Length(ObjectHooks)+1);
  Hooks[High(Hooks)]:=hook;
end;

function NewTask(cancancel: Boolean; title: string): TRunningTask;
begin
  Result:=TRunningTask.Create(cancancel,title);
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

constructor TRunningTask.Create(acancel: Boolean; atitle: string);
begin
  inherited Create;
  FCanCancel:=acancel;
  FTitle:=atitle;
  FStartTime:=Now;
  FProgress:=0;
  FAborted:=False;
end;

procedure TRunningTask.Abort;
begin
  FAborted:=True;
end;

procedure TRunningTask.Step(s: Real);
begin
  FProgress:=FProgress+s;
  CallHooks;
end;

end.

