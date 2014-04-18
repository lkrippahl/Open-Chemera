{*******************************************************************************
This source file is public domain. It is provided as is, with no explicit
or implied warranties regarding anything whatsoever.
********************************************************************************
Author: Ludwig Krippahl
Date: 2014 03 07
Purpose:
  Class to position forms relative to one another
  ----- header ----
  left|center|right
  ------footer ----
Requirements:
Revisions:
To do: Flexible docking, for each form specify all the neighbours
*******************************************************************************}
unit formdocker;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Forms;

const
  formHeader=0;
  formLeft=1;
  formCenter=2;
  formRight=3;
  formFooter=4;

type
  TFormPlacement=record
    FixedWidth,FixedHeight:Boolean;
    MinWidth,MinHeight,MaxWidth,MaxHeight:Integer;
    Form:TForm;
    end;

  TFormManager=class

  end;

implementation

end.

