program chemera;

{$mode objfpc}{$H+}

uses
  {$IFDEF UNIX}{$IFDEF UseCThreads}
  cthreads,
  {$ENDIF}{$ENDIF}
  Interfaces, // this includes the LCL widgetset
  Forms, chemeramain, displayobjects, molecules, pdbmolecules, basetypes, selections, glutwindow, base3ddisplay,
displaysettings;

{$R *.res}

begin
  Application.Initialize;
  Application.CreateForm(TCmMainForm, CmMainForm);
  Application.CreateForm(TDisplaySettingsForm, DisplaySettingsForm);
  Application.Run;
end.

