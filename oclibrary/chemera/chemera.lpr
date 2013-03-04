program chemera;

{$mode objfpc}{$H+}

uses
  {$IFDEF UNIX}{$IFDEF UseCThreads}
  cthreads,
  {$ENDIF}{$ENDIF}
  Interfaces, // this includes the LCL widgetset
  Forms, chemeramain, displayobjects, molecules, pdbmolecules, basetypes, selections, base3ddisplay,
displaysettings, oglform, pdbparser, lazopenglcontext,
oclconfiguration, geomhash, moleculetree, notifications;

{$R *.res}

begin
  Application.Initialize;
  Application.CreateForm(TCmMainForm, CmMainForm);
  Application.CreateForm(TMolTreeForm, MolTreeForm);
  Application.Run;
end.

