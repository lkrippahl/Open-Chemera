{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 9.1.2011
Purpose:
  Classes for atoms and groups of atoms (residues) or groups of grups (chains).
Requirements:
Revisions:
To do:
  RemoveTaggedBonds is only for bonds with atoms that will be deleted
  Do another one for removing bonds alone. Needs another callback, and bonds need
  tags too
  functions needing callbacks: clearbonds, clearallbonds
*******************************************************************************}


unit molecules;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes;

const
  //bond types
  SingleBond=0;
  DoubleBond=1;
  HBond=2;
  Aromatic=3;

  //special tag
  TagNone=0;
  TagToDelete=1;


type
  TOnDeleteCallback=procedure(Tag:Integer) of object;
  TAbstractMolecule=class
    {placeholder class so that atoms molecules can be nested hierarchies}

    protected

    public
    end;

  { TAtom }

  TAtom=class
  protected
    //Data for each atom in molecule. Bond data is in parent group
    FAtomicNumber:Integer;
    FName:string;
    FParent:TAbstractMolecule;

    //Atom properties
    FCoord:TCoord;
    FRadius:TFloat;
    FCharge:TFloat;
    FMass:TFloat;
    FID:Integer;

  public

    //temporary tag for selection, deletion, etc.
    //Should not be relied on to store persistent information. It is meant
    //to be used and discarded in the same function (such as deleting atoms)
    Tag:Integer;
    property Name:string read FName write FName;
    property Coords:TCoord read FCoord write FCoord;
    property Charge:TFloat read FCharge write FCharge;
    property ID:Integer read FID write FID;       //should be unique for each molecule
    property Parent:TAbstractMolecule read FParent;
    constructor Create(AName:string; AID:Integer;AParent:TAbstractMolecule);
  end;

  TAtoms=array of TAtom;

  TAtomBond=record
    Atom1,Atom2:TAtom;
    BondType:Integer;
    Tag:Integer;   //tag is temporary, not suitable for persistent data
  end;
  TAtomBonds=array of TAtomBond;


  { TMolecule }

  TMolecule=class(TAbstractMolecule)
    {Class for nested hierarchy of molecular fragments.
    It can have both atoms (e.g a residue) or groups (e.g. a chain)
    though, chemically, it doesn't make sense to have both}

  protected
    FAtoms:TAtoms;              //Atoms that do not belong in groups

    FGroups:array of TMolecule; //Groups such as monomers. If molecule has
                                //groups it should not have atoms.
                                //Will not assume that it only
                                //can have atoms or groups, but
                                //that makes more sense for
                                //modelling actual molecules

    FBondsTable:TAtomBonds;     //Best if there is only one bonds table
                                //at the top level of the molecule hierarchy, for performance.
                                //But this is not mandatory

    FParent:TMolecule;          // if a monomer in a larger complex. Can make tres of arbitrary depth
    FType:string;               //not fixed. For proteins should be residue, chain, complex
    FName:string;               //identifier, managed externally
    FID:Integer;                //identifier, managed externally

    function AtomIndex(AAtom:TAtom):Integer;      //index in FAtoms
    function GroupIndex(Group:TMolecule):Integer; //Index in FGroups
    procedure AddGroup(Group:TMolecule);          //adds an existing group
    procedure AddAtom(Atom:TAtom);                //adds an existing atom
  public
    property MolType:string read FType write FType;
    property Name:string read FName write FName;
    property ID:Integer read FID write FID;
    constructor Create(AName:string; AID:Integer; AParent:TMolecule);
    function NewGroup(GName:string; GID:Integer):TMolecule;// creates and adds group
    function NewAtom(AName:string; AID:integer):TAtom; //creates and adds atom;
    function AllAtoms:TAtoms;
    function AllBonds:TAtomBonds;
    function GroupCount:Integer;
    function GetGroup(GroupIx:Integer):TMolecule;
    function GetAtom(AtomIx:Integer):TAtom;
    procedure TagAllAtoms(Tag:Integer);
    procedure TagAllBonds(Tag:Integer);
    procedure TagAtoms(Atoms:TAtoms;Tag:Integer);
    procedure TagBonds(Bonds:TAtomBonds;Tag:Integer);
    procedure DeleteTaggedAtoms(Tag:Integer;OnDelete:TOnDeleteCallback=nil);
    procedure RemoveTaggedAtomBonds(Tag:Integer;OnDelete:TOnDeleteCallback=nil);
    procedure RemoveTaggedBonds(Tag:Integer;OnDelete:TOnDeleteCallback=nil);
    procedure DeleteAtoms(Atoms:TAtoms;OnDelete:TOnDeleteCallback=nil);
    function AtomById(AId:Integer):TAtom;    //looks also recursively in groups

    //These procedures clear all groups or atoms, respectively
    //They are meant to be used on empty molecules
    //since they do not callback to ondelete events
    procedure CreateEmptyGroups(Count:Integer);
    procedure CreateEmptyAtoms(Count:Integer);

    procedure ClearAtoms(OnDelete:TOnDeleteCallback=nil);
    procedure ClearGroups(OnDelete:TOnDeleteCallback=nil);
    procedure ClearBonds(OnDelete:TOnDeleteCallback=nil);   //TO DO: Needs implementing callback

    //This is a recursive version of ClearBonds to clear all bonds
    //from offspring too.
    procedure ClearAllBonds(OnDelete:TOnDeleteCallback=nil);//TO DO: Needs implementing callback

    //Clears everything
    procedure CleanUp(OnDelete:TOnDeleteCallback=nil);

    procedure Free;

  end;
  TMolecules=array of TMolecule;


implementation

{ TAtom }

constructor TAtom.Create(AName: string;AID:Integer;AParent:TAbstractMolecule);
begin
  inherited Create;
  FName:=AName;
  FParent:=AParent;
  FID:=AID;
end;

{ TMolecule }

function TMolecule.AtomIndex(AAtom: TAtom): Integer;
begin
  Result:=High(FAtoms);
  while (Result>=0) and (AAtom<>FAtoms[Result]) do
    Dec(Result);
end;

function TMolecule.GroupIndex(Group: TMolecule): Integer;
begin
  Result:=High(FGroups);
  while (Result>=0) and (FGroups[Result]<>Group) do
    Dec(Result);
end;

constructor TMolecule.Create(AName: string; AID: Integer; AParent: TMolecule);
begin
  inherited Create;
  FName:=AName;
  FID:=AID;
  FParent:=AParent;
  FGroups:=nil;
  FBondsTable:=nil;
  FAtoms:=nil;
end;

procedure TMolecule.DeleteTaggedAtoms(Tag:Integer;OnDelete:TOnDeleteCallback=nil);
//Recursive through groups
//Calls removeTaggedBonds first

var
  tmp:TAtoms;
  f,i:Integer;

begin
  RemoveTaggedAtomBonds(Tag); //Remove bonds before deleting atoms
  SetLength(tmp,Length(FAtoms)); //New atoms array
  i:=0;

  //if atom is not to delete, copy into new array
  for f:=0 to High(FAtoms) do
    if FAtoms[f].Tag<>Tag then
      begin
      tmp[i]:=FAtoms[f];
      inc(i);
      end
    else FAtoms[f].Free;

  //Resize array, store and delete atoms in offpring
  SetLength(tmp,i);
  FAtoms:=tmp;
  for f:=0 to High(FGroups) do
    FGroups[f].DeleteTaggedAtoms(Tag);
end;

procedure TMolecule.RemoveTaggedAtomBonds(Tag: Integer;OnDelete:TOnDeleteCallback=nil);

//Recursive through groups
//Removes bonds with at least one tagged atom.
//This procedure is meant for cleaning bonds with atoms about to be deleted

var
  tmp:TAtomBonds;
  f,i:Integer;

begin
  //callback to warn about deletions
  if Assigned(OnDelete) then OnDelete(Tag);

  SetLength(tmp,Length(FBondsTable)); //new bonds table
  i:=0;

  //if not to be deleted, copy to new bonds table
  for f:=0 to High(FBondsTable) do
    if (FBondsTable[f].Atom1.Tag<>Tag) and
       (FBondsTable[f].Atom2.Tag<>Tag) then
        begin
        tmp[i]:=FBondsTable[f];
        Inc(i);
        end;
  //resize and store
  SetLength(tmp,i);
  FBondsTable:=tmp;
  for f:=0 to High(FGroups) do
      //No more callbacks, as one should be enough to warn about all deletions
      //involving tagged atoms
      FGroups[f].RemoveTaggedAtomBonds(Tag);
end;

procedure TMolecule.RemoveTaggedBonds(Tag: Integer;
    OnDelete: TOnDeleteCallback);

//Recursive through groups
//Removes bonds that are tagged

var
  tmp:TAtomBonds;
  f,i:Integer;

begin
  //callback to warn about deletions
  if Assigned(OnDelete) then OnDelete(Tag);

  SetLength(tmp,Length(FBondsTable)); //new bonds table
  i:=0;

  //if not to be deleted, copy to new bonds table
  for f:=0 to High(FBondsTable) do
    if FBondsTable[f].Tag<>Tag then
        begin
        tmp[i]:=FBondsTable[f];
        Inc(i);
        end;

  //resize and store
  SetLength(tmp,i);
  FBondsTable:=tmp;
  for f:=0 to High(FGroups) do
      //No more callbacks, as one should be enough to warn about all deletions
      //of tagged bonds
      RemoveTaggedBonds(Tag);
end;

function TMolecule.NewGroup(GName: string; GID: Integer): TMolecule;
begin
  Result:=TMolecule.Create(GName,GID,Self);
  AddGroup(Result);
end;

procedure TMolecule.AddGroup(Group: TMolecule);
begin
  Group.FParent:=Self;
  SetLength(FGroups,Length(FGroups)+1);
  FGroups[High(FGroups)]:=Group;
end;

function TMolecule.NewAtom(AName: string; AID: integer): TAtom;
begin
  Result:=TAtom.Create(AName,AID,Self);
  AddAtom(Result);
end;

function TMolecule.AllAtoms: TAtoms;

//flattens hierarchy into a single array of atoms
//returns a copy array but each element is the atom object

var
  f,g,i:Integer;
  tmp:TAtoms;

begin
  Result:=Copy(FAtoms,0,Length(FAtoms));

  //Add to atoms at this level all the offspring atoms
  for f:=0 to High(FGroups) do
    begin
    tmp:=FGroups[f].AllAtoms;
    i:=Length(Result);
    SetLength(Result,i+Length(tmp));
    for g:=0 to High(tmp) do
      Result[g+i]:=tmp[g];
    end;
end;

function TMolecule.AllBonds: TAtomBonds;

//flattens hierarchy into a single array of atom bonds
//returns a copy array, so changes in elements do not affect the bonds tables
//in the molecule (but changes in the atoms will, as these are objects)

var
  f,g,i:Integer;
  tmp:TAtomBonds;

begin
  Result:=Copy(FBondsTable,0,Length(FBondsTable));

  //Add to atoms at this level all the offspring atoms
  for f:=0 to High(FGroups) do
    begin
    tmp:=FGroups[f].AllBonds;
    i:=Length(Result);
    SetLength(Result,i+Length(tmp));
    for g:=0 to High(tmp) do
      Result[g+i]:=tmp[g];
    end;
end;

function TMolecule.GroupCount: Integer;
begin
  Result:=Length(FGroups);
end;


function TMolecule.GetGroup(GroupIx:Integer): TMolecule;
begin
  Assert((GroupIx<Length(FGroups)) and (GroupIx>=0),'Invalid group index');
  Result:=FGroups[GroupIx];
end;

function TMolecule.GetAtom(AtomIx: Integer): TAtom;
begin
  Assert((AtomIx<Length(FGroups)) and (AtomIx>=0),'Invalid atom index');
  Result:=FAtoms[AtomIx];
end;

procedure TMolecule.TagAllAtoms(Tag: Integer);

//Tag all atoms, including in offspring, recusively

var f:Integer;

begin
  for f:=0 to High(FAtoms) do FAtoms[f].Tag:=Tag;
  for f:=0 to High(FGroups) do FGroups[f].TagAllAtoms(Tag);
end;

procedure TMolecule.TagAllBonds(Tag: Integer);

//Tag all bonds, including in offspring, recusively

var f:Integer;

begin
  for f:=0 to High(FBondsTable) do FBondsTable[f].Tag:=Tag;
  for f:=0 to High(FGroups) do FGroups[f].TagAllBonds(Tag);
end;


procedure TMolecule.TagAtoms(Atoms: TAtoms; Tag: Integer);

//tag a list of atoms

var f:Integer;

begin
  for f:=0 to High(Atoms) do Atoms[f].Tag:=Tag;
end;

procedure TMolecule.TagBonds(Bonds: TAtomBonds; Tag: Integer);

var f:Integer;

begin
  for f:=0 to High(FBondsTable) do
      FBondsTable[f].Tag:=Tag;
end;

procedure TMolecule.DeleteAtoms(Atoms: TAtoms;OnDelete:TOnDeleteCallback=nil);

//Deletes those atoms from the Atoms set that belong to this molecule
//or to offspring groups. Note that Atoms array becomes unuseable,
//and that atoms outside this molecule will be tagged for deletion but not deleted.

begin
  //reset tags first so that none get deleted that shouldn't
  TagAllAtoms(TagNone);
  TagAtoms(Atoms,TagToDelete);
  DeleteTaggedAtoms(TagToDelete,OnDelete);
end;

function TMolecule.AtomById(AId: Integer): TAtom;
//WARNING: if IDs are repeated, will return the first atom.
//It's up to caller to guarantee no repeated ids.

var f:Integer;

begin
  Result:=nil;
  for f:=0 to High(FAtoms) do
    if FAtoms[f].Id=AId then
      begin
      Result:=FAtoms[f];
      if Result<>nil then Break;
      end;
  if Result=nil then
      for f:=0 to High(FGroups) do
        begin
        Result:=FGroups[f].AtomById(Id);
        if Result<> nil then Break;
        end;
end;

procedure TMolecule.CreateEmptyGroups(Count: Integer);
//Sets molecule with an array of empty groups
//Clears all existing groups first.

var f:Integer;

begin
  ClearGroups;
  SetLength(FGroups,Count);
  for f:=0 to High(FGroups) do
    FGroups[f]:=TMolecule.Create('',f,Self);
end;

procedure TMolecule.CreateEmptyAtoms(Count: Integer);
//Sets molecule with an array of empty atoms
//Clears all existing atoms first.

var f:Integer;

begin
  ClearAtoms;
  SetLength(FAtoms,Count);
  for f:=0 to High(FAtoms) do
    FAtoms[f]:=TAtom.Create('',f,Self);
end;

procedure TMolecule.ClearAtoms(OnDelete:TOnDeleteCallback=nil);

//Clears atoms at this level. Not recursive, does not delete atoms in offespring.

begin
  DeleteAtoms(FAtoms,OnDelete);
  FAtoms:=nil;
end;

procedure TMolecule.ClearGroups(OnDelete:TOnDeleteCallback=nil);

var f:Integer;

begin
  //first, tag all atoms in all offspring groups
  for f:=0 to High(FGroups) do
    FGroups[f].TagAllAtoms(TagToDelete);

  //Warn deletion before deleting anything
  if Assigned(OnDelete) then OnDelete(TagToDelete);

  ClearAllBonds;

  //finally, free each group (atoms are freed along with the group
  for f:=0 to High(FGroups) do
    FGroups[f].Free;
  FGroups:=nil;
end;

procedure TMolecule.ClearBonds(OnDelete:TOnDeleteCallback);
//TO DO: use callback to warn of bond deletion

begin
  //Warn of deletion if callback provided, by tagging all bonds
  if Assigned(OnDelete) then
     begin
     TagBonds(FBondsTable,TagToDelete);
     OnDelete(TagToDelete);
     end;
  FBondsTable:=nil;
end;

procedure TMolecule.ClearAllBonds(OnDelete:TOnDeleteCallback);

//Clears all bonds in this level and, recursively, in all descendants
//useful for cleaning up everything, before deleting all atoms
//TO DO: callback to warn of bond deletion

var f:Integer;

begin
  //if needs to warn about deletion
  if Assigned(OnDelete) then
    begin
    //Tag all bonds first
    TagAllBonds(TagToDelete);
    for f:=0 to High(FGroups) do
      FGroups[f].TagAllBonds(TagToDelete);
    OnDelete(TagToDelete);
    end;
  //Clear bonds with no more callbacks, even in offspring
  ClearBonds;
  for f:=0 to High(FGroups) do
      FGroups[f].ClearAllBonds;
end;

procedure TMolecule.CleanUp(OnDelete:TOnDeleteCallback=nil);
begin
  //Warn first of all deletions if callback provided
  if Assigned(OnDelete) then
      begin
      TagAllAtoms(TagToDelete);
      OnDelete(TagToDelete);
      end;
  ClearAllBonds(OnDelete); //clear bonds first so no time spent checking them for deleted atoms
  ClearAtoms(OnDelete);
  ClearGroups(OnDelete);
end;

procedure TMolecule.Free;

//when freeing

begin
  if Self<>nil then
    begin
    CleanUp;
    inherited;
    end;
end;

procedure TMolecule.AddAtom(Atom: TAtom);
begin
  Atom.FParent:=Self;
  SetLength(FAtoms,Length(FAtoms)+1);
  FAtoms[High(FAtoms)]:=Atom;
end;

end.

