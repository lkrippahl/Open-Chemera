{*******************************************************************************
This source file is public domain. It is provided as is, with no explicit
or implied warranties regarding anything whatsoever.
********************************************************************************
Author: Ludwig Krippahl
Date: 25 01 2014
Purpose:
  Exporting scenes to POV Ray
Revisions:
*******************************************************************************}

unit povray;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes;

type

  { TPovRayExport }

  TPovRayExport=class
    private
      FCode:TStringList;
    public
      constructor Create;
      procedure Free;virtual;
  end;


implementation

{ TPovRayExport }

constructor TPovRayExport.Create;
begin
  inherited;
  FCode:=TStringList.Create;
end;

procedure TPovRayExport.Free;
begin
  FCode.Free;
  FCode:=nil;
  inherited Free;
end;

end.

