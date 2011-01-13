unit ogremanager;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, dynlibs,Dialogs,thread,ctypes;
type


  TRun=procedure;stdcall;
  TAddSphere=procedure(ID,mat,x,y,z,rad:ctypes.cint32);stdcall;

procedure Setup;
procedure Cleanup;

var
  OgreMan:pointer;
  NewOgre:TNewOgre;
  Run:TRun;
  DllSetup:TRun;
  AddSphere:TAddSphere;
  DllH:Integer;

implementation


procedure Setup;


begin
  OgreMan:=nil;
  DllH:= LoadLibrary('Ogre3Ddll');
  Run:=TRun(GetProcAddress(DllH, 'Run'));
  DllSetup:=TRun(GetProcAddress(DllH, 'DllSetup'));
  AddSphere:=TAddSphere(GetProcAddress(DllH, 'AddSphere'));
  //ShowMessage('setup OK');
  DllSetup();

end;

procedure Cleanup;

begin
  //OgreMan.Free;

end;

end.

