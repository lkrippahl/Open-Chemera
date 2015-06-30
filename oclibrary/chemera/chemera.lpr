program chemera;

{$mode objfpc}{$H+}

uses
  {$IFDEF UNIX}{$IFDEF UseCThreads}
  cthreads,
  {$ENDIF}{$ENDIF}
  Interfaces, // this includes the LCL widgetset
  Forms, chemeramain, displayobjects, molecules, pdbmolecules, basetypes,
  selections, base3ddisplay, displaysettings, oglform, pdbparser,
  lazopenglcontext, oclconfiguration, geomhash, geomutils, moleculetree,
  notifications, molutils, surface, linegrids, povray, dockdomains,
  dockconstraints, alignment, formdocker;

{$R *.res}

begin
  Application.Initialize;
  Application.CreateForm(TCmMainForm, CmMainForm);
  Application.CreateForm(TMolTreeForm, MolTreeForm);
  Application.Run;
end.

