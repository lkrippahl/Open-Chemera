{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain. Enjoy.
********************************************************************************
Author: Ludwig Krippahl
Date: 9.1.2011
Purpose:
  Functions for loading pdb files into a molecule
Requirements:
Revisions:
To do:
  Everything
*******************************************************************************}
unit pdbloader;


{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, molecules, pdbparser;

function PDB2Molecule(Pdb:TPdbReader):TMolecule;

implementation

function PDB2Molecule(Pdb:TPdbReader):TMolecule;

begin
  Result:=TMolecule.Create;
end;

end.

