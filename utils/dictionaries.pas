{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 20.4.2011
Purpose:
  Translations for names, such as 3 letter code AAs to 1 letter code.
  Residue properties and classification (polar, charged, etc)


Requirements:
  None, for now, but see TODO

Revisions:
  These data should not be hardcoded. It would be better to store it in data
  files
*******************************************************************************}

unit dictionaries;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils; 

function OneLetterCode(ResName:string):Char;

implementation

function OneLetterCode(ResName:string):Char;

const
  ThreeLetterCodes:array [0..19] of string=('ALA', 'ARG', 'ASN', 'ASP', 'CYS','GLU',
    'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
    'THR', 'TRP', 'TYR', 'VAL');
  OneLetterCodes:array [0..19] of Char=('A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G',
    'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V');

var f:Integer;

begin
  Result:='?';
  f:=High(ThreeLetterCodes);
  ResName:=UpperCase(ResName);
  while (f>=0) and (ResName<>ThreeLetterCodes[f]) do
    Dec(f);
  if f>=0 then Result:=OneLetterCodes[f];
end;

end.

