{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 14-8-2012
Purpose:
  Form to select and browse molecules in a pdblayer
Requirements:
Revisions:
To do:
*******************************************************************************}

unit moleculetree;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, FileUtil, LResources, Forms, Controls, Graphics, Dialogs,
  ComCtrls,pdbmolecules;

type

  { TMolTreeForm }

  TMolTreeForm = class(TForm)
    TreeView1: TTreeView;
  private
    { private declarations }
    FLayerMan:TPdbModelMan;
  public
    { public declarations }
    constructor Create(TheOwner: TComponent;ALayerMan:TPDBModelMan);reintroduce;
    procedure Rebuild;
  end;

var
  MolTreeForm: TMolTreeForm;

implementation

{ TMolTreeForm }

constructor TMolTreeForm.Create(TheOwner: TComponent;ALayerMan: TPDBModelMan);
begin
  inherited Create(TheOwner);
end;

procedure TMolTreeForm.Rebuild;
begin

end;

initialization
  {$I moleculetree.lrs}

end.

