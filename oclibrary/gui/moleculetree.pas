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
    FLayerMan:TPdbLayerMan;
  public
    { public declarations }
    constructor Create(ALayerMan:TPDBLayerMan);
    procedure Rebuild;
  end;

var
  MolTreeForm: TMolTreeForm;

implementation

{ TMolTreeForm }

constructor TMolTreeForm.Create(ALayerMan: TPDBLayerMan);
begin

end;

procedure TMolTreeForm.Rebuild;
begin

end;

initialization
  {$I moleculetree.lrs}

end.

