{*******************************************************************************
This source file is public domain. It is provided as is, with no explicit
or implied warranties regarding anything whatsoever.
********************************************************************************
Author: Ludwig Krippahl
Date: 27.12.2013
Purpose:
  Configuration variables for learning library
Requirements:
Revisions:
To do:
*******************************************************************************}

unit lrnconfig;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils;

procedure LrnReport(s:string);

var
  LrnVerbose:Boolean=False;

implementation

procedure LrnReport(s: string);
begin
  if LrnVerbose then WriteLn(s);
end;

end.

