{*******************************************************************************
This source file is public domain. It is provided as is, with no explicit
or implied warranties regarding anything whatsoever.
********************************************************************************
Author: Ludwig Krippahl
Date: 18-8-2012
Purpose:
  Simplify communication between forms and other unit in one app. Create
  notificaton lists for calling back forms to update displays, etc.
Requirements:
Revisions:
To do: Comments
*******************************************************************************}

unit notifications;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes;

type
  TNotificationCallback=procedure (Msg:TIntegers) of object;
  TNotifications=array of TNotificationCallback;
  TNotifyLists=array of TNotifications;

var
  NotifyLists:array of TNotifications;

function NewList:Integer;
procedure Notify(ListIx:Integer;Msg:TIntegers);
procedure AddNotification(ListIx:Integer;NCB:TNotificationCallback);

implementation

function NewList: Integer;
begin
  SetLength(NotifyLists,Length(NotifyLists)+1);
  NotifyLists[High(NotifyLists)]:=nil;
  Result:=High(NotifyLists);
end;

procedure Notify(ListIx: Integer; Msg: TIntegers);

var f:Integer;

begin
  for f:=0 to High(NotifyLists[ListIx]) do
      NotifyLists[ListIx][f](Msg);
end;

procedure AddNotification(ListIx: Integer; NCB: TNotificationCallback);
begin
  SetLength(NotifyLists[ListIx],Length(NotifyLists[ListIx])+1);
  NotifyLists[ListIx][High(NotifyLists[ListIx])]:=NCB;
end;

end.

