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
  30-4-2012 added Parent property to TMolecule
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
  Classes, SysUtils, basetypes, geomutils;

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
  TMolecule=class;

  { TAtom }

  TAtom=class
  protected
    //Data for each atom in molecule. Bond data is in parent group
    FAtomicNumber:Integer;
    FName:string;
    FParent:TMolecule;

    //Atom properties
    FCoord:TCoord;
    FRadius:TFloat;
    FCharge:TFloat;
    FMass:TFloat;       { TODO : Should this be here? }
    FID:Integer;

  public

    //temporary tag for selection, deletion, etc.
    //Should not be relied on to store persistent information. It is meant
    //to be used and discarded in the same function (such as deleting atoms)
    Tag:Integer;
    property AtomicNumber:Integer read FAtomicNumber write FAtomicNumber;
      //AtomicNumber<=0 means undetermined
    property Name:string read FName write FName;
    property Coords:TCoord read FCoord write FCoord;
    property Charge:TFloat read FCharge write FCharge;
    property Radius:TFloat read FRadius write FRadius;
    property ID:Integer read FID write FID;       //should be unique for each molecule
    property Parent:TMolecule read FParent;
    constructor Create(AName:string; AID:Integer;AParent:TMolecule);
    constructor CopyFrom(AAtom:TAtom;NewParent:TMolecule);
  end;

  TAtoms=array of TAtom;

  TAtomBond=record
    Atom1,Atom2:TAtom;
    BondType:Integer;
    Tag:Integer;   //tag is temporary, not suitable for persistent data
  end;
  TAtomBonds=array of TAtomBond;

  TMolecules=array of TMolecule;

  { TMolecule }

  TMolecule=class
    {Class for nested hierarchy of molecular fragments.
    It can have both atoms (e.g a residue) or groups (e.g. a chain)
    though, chemically, it doesn't make sense to have both}

  protected
    FAtoms:TAtoms;              //Atoms that do not belong in groups

    FGroups:TMolecules;         //Groups such as monomers. If molecule has
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
    procedure AddAtom(Atom:TAtom);                //adds an existing atom
  public
    property GroupAtoms:TAtoms read FAtoms;
    property MolType:string read FType write FType;
    property Name:string read FName write FName;
    property ID:Integer read FID write FID;
    property Parent:TMolecule read FParent;
    property Groups:TMolecules read FGroups;
    constructor Create(AName:string; AID:Integer; AParent:TMolecule);
    constructor CopyFrom(AMolecule:TMolecule;NewParent:TMolecule);
    procedure AddGroup(Group:TMolecule);//adds an existing group
    procedure RenumberAtoms(Start:Integer=1);
    function NewGroup(GName:string; GID:Integer):TMolecule; // creates and adds group
    function NewAtom(AName:string; AID:integer):TAtom;      //creates and adds atom;
    function AllAtoms:TAtoms;
    function AllCoords:TCoords;
    function AllTerminalGroups:TMolecules;
    function AllBonds:TAtomBonds;
    function GroupCount:Integer;
    function AtomCount:Integer;
    function BondCount:Integer;
    function GetGroup(GroupIx:Integer):TMolecule;overload;
    function GetGroup(GroupName:string):TMolecule;overload;
    function GetGroupById(GroupId:Integer):TMolecule;
    function GetAtom(AtomIx:Integer):TAtom;overload;
    function GetAtom(AtomName:string):TAtom;overload;

    function ListGroupNames:TSimpleStrings;

    procedure TagAllAtoms(Tag:Integer);
    procedure TagAllBonds(Tag:Integer);
    procedure TagAtoms(Atoms:TAtoms;Tag:Integer);
    procedure TagBonds(Bonds:TAtomBonds;Tag:Integer);
    procedure DeleteTaggedAtoms(Tag:Integer;OnDelete:TOnDeleteCallback=nil);
    procedure RemoveTaggedAtomBonds(Tag:Integer;OnDelete:TOnDeleteCallback=nil);
    procedure RemoveTaggedBonds(Tag:Integer;OnDelete:TOnDeleteCallback=nil);
    procedure DeleteAtoms(Atoms:TAtoms;OnDelete:TOnDeleteCallback=nil);
    procedure DeleteEmptyGroups;
    function AtomById(AId:Integer):TAtom;    //looks also recursively in groups

    procedure Transform(TranslationVec:TCoord);overload;
    procedure Transform(RotationMat:TRotMatrix);overload;
    procedure Transform(Quat:TQuaternion);overload;
    procedure Transform(TranslationVec:TCoord;RotationMat:TRotMatrix);overload;
    procedure Transform(Center:TCoord;Rotation:TRotMatrix;Translation:TCoord);overload;
    //Subtracts center, rotates then translates

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

  function AppendAtomsToArray(const Original,ToAppend:TAtoms):TAtoms;
  procedure AppendGroupToArray(const Group:TMolecule;var Groups:TMolecules);
  procedure AppendGroupsToArray(var Original:TMolecules;const ToAppend:TMolecules);

implementation

function AppendAtomsToArray(const Original,ToAppend:TAtoms):TAtoms;

var f:Integer;
    lorig,lapp:Integer;
begin
  lorig:=Length(Original);
  lapp:=Length(ToAppend);
  SetLength(Result,lorig+lapp);
  for f:=0 to lorig-1 do
      Result[f]:=Original[f];
  for f:=0 to lapp-1 do
    Result[f+lorig]:=ToAppend[f];
end;

procedure AppendGroupToArray(const Group: TMolecule; var Groups: TMolecules);
begin
  SetLength(Groups,Length(Groups)+1);
  Groups[High(Groups)]:=Group;
end;

procedure AppendGroupsToArray(var Original: TMolecules; const ToAppend: TMolecules);

var f:Integer;
    lorig,lapp:Integer;
begin
  lorig:=Length(Original);
  lapp:=Length(ToAppend);
  SetLength(Original,lorig+lapp);
  for f:=0 to lapp-1 do
    Original[f+lorig]:=ToAppend[f];
end;


{ TAtom }

constructor TAtom.Create(AName: string;AID:Integer;AParent:TMolecule);
begin
  inherited Create;
  FName:=AName;
  FParent:=AParent;
  FID:=AID;
  FAtomicNumber:=-1; //undetermined
end;

constructor TAtom.CopyFrom(AAtom: TAtom; NewParent: TMolecule);
begin
  inherited Create;
  FName:=AAtom.Name;
  FParent:=NewParent;
  FID:=AAtom.FID;
  FAtomicNumber:=AAtom.FAtomicNumber;
  FCoord:=AAtom.FCoord;
  FRadius:=AAtom.FRadius;
  FCharge:=AAtom.FCharge;
  FMass:=AAtom.FMass;
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

constructor TMolecule.CopyFrom(AMolecule: TMolecule;NewParent:TMolecule);

var
  f,bondtot:Integer;
  at1,at2:TAtom;

begin
  inherited Create;
  FName:=AMolecule.FName;
  FID:=AMolecule.FID;
  FParent:=NewParent;
  SetLength(FGroups,Length(AMolecule.FGroups));
  for f:=0 to High(FGroups) do
    FGroups[f]:=TMolecule.CopyFrom(AMolecule.FGroups[f],Self);

  SetLength(FAtoms,Length(AMolecule.FAtoms));
  for f:=0 to High(FAtoms) do
    FAtoms[f]:=TAtom.CopyFrom(AMolecule.FAtoms[f],Self);

  //copy only bonds between atoms in the new molecule
  //this means some bonds may be lost when copying a subset of chains...
  SetLength(FBondsTable,Length(AMolecule.FBondsTable));
  bondtot:=0;
  for f:=0 to High(AMolecule.FBondsTable) do
    begin
    at1:=AtomById(AMolecule.FBondsTable[f].Atom1.ID);
    at2:=AtomById(AMolecule.FBondsTable[f].Atom2.ID);
    if (at1<>nil) and (at2<>nil) then
      begin
      FBondsTable[bondtot].Atom1:=at1;
      FBondsTable[bondtot].Atom2:=at2;
      FBondsTable[bondtot].BondType:=AMolecule.FBondsTable[f].BondType;
      Inc(bondtot);
      end;
    end;
  SetLength(FBondsTable,bondtot);
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

procedure TMolecule.RenumberAtoms(Start: Integer);

var
  atomlist:TAtoms;
  f:Integer;

begin
  atomlist:=AllAtoms;
  for f:=0 to High(atomlist) do
    atomlist[f].ID:=Start+f;
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

function TMolecule.AllCoords: TCoords;
//returns a copy of all coords in all atoms

var
  f,g,i:Integer;
  tmp:TCoords;

begin
  SetLength(Result,Length(FAtoms));
  for f:=0 to High(FAtoms) do Result[f]:=FAtoms[f].Coords;

  //Add to this level all the offspring levels
  for f:=0 to High(FGroups) do
    begin
    tmp:=FGroups[f].AllCoords;
    i:=Length(Result);
    SetLength(Result,i+Length(tmp));
    for g:=0 to High(tmp) do
      Result[g+i]:=tmp[g];
    end;
end;

function TMolecule.AllTerminalGroups: TMolecules;

var
  f:Integer;

begin
  Result:=nil;
  for f:=0 to High(FGroups) do
    if FGroups[f].Groups=nil then
      AppendGroupToArray(FGroups[f],Result)
    else
      AppendGroupsToArray(Result,FGroups[f].AllTerminalGroups);
end;

function TMolecule.AllBonds: TAtomBonds;

//flattens hierarchy into a single array of atom bonds
//returns a copy array, so changes in elements do not affect the bonds tables
//in the molecule (but changes in the atoms will, as these are objects)

var
  f,i,j,tot:Integer;
  tmp:TAtomBonds;

begin
  tot:=BondCount;
  SetLength(Result,tot);
  f:=0;
  while f<=High(FBondsTable) do
    begin
    Result[f]:=FBondsTable[f];
    Inc(f);
    end;

  //Add to atoms at this level all the offspring atoms
  for i:=0 to High(FGroups) do
    begin
    tmp:=FGroups[i].AllBonds;
    j:=0;
    while j<=High(tmp) do
      begin
      Result[f]:=tmp[j];
      Inc(j);
      Inc(f);
      end;
    end;
end;

function TMolecule.GroupCount: Integer;
begin
  Result:=Length(FGroups);
end;

function TMolecule.AtomCount: Integer;
//total number of atoms in molecule and groups

var
  f:Integer;

begin
  Result:=Length(FAtoms);
  for f:=0 to High(FGroups) do
      Result:=Result+FGroups[f].AtomCount;
end;

function TMolecule.BondCount: Integer;

var
  f:Integer;

begin
  Result:=Length(FBondsTable);
  for f:=0 to High(FGroups) do
      Result:=Result+FGroups[f].BondCount;
end;


function TMolecule.GetGroup(GroupIx:Integer): TMolecule;
begin
  Assert((GroupIx<Length(FGroups)) and (GroupIx>=0),'Invalid group index');
  Result:=FGroups[GroupIx];
end;

function TMolecule.GetGroup(GroupName: string): TMolecule;

var f:Integer;

begin
  Result:=nil;
  for f:=0 to High(FGroups) do
  if FGroups[f].Name=GroupName then
    begin
    Result:=FGroups[f];
    Break;
    end;
end;

function TMolecule.GetGroupById(GroupId: Integer): TMolecule;

var f:Integer;

begin
  Result:=nil;
  for f:=0 to High(FGroups) do
  if FGroups[f].ID=GroupID then
    begin
    Result:=FGroups[f];
    Break;
    end;
end;

function TMolecule.GetAtom(AtomIx: Integer): TAtom;
begin
  Assert((AtomIx<Length(FGroups)) and (AtomIx>=0),'Invalid atom index');
  Result:=FAtoms[AtomIx];
end;

function TMolecule.GetAtom(AtomName: string): TAtom;
var f:Integer;

begin
  Result:=nil;
  for f:=0 to High(FAtoms) do
    if FAtoms[f].Name=AtomName then
      begin
      Result:=FAtoms[f];
      Break;
      end;
end;

function TMolecule.ListGroupNames: TSimpleStrings;

var f:Integer;

begin
  SetLength(Result,Length(FGroups));
  for f:=0 to High(FGroups) do
    Result[f]:=FGroups[f].Name;
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
  for f:=0 to High(Bonds) do
      Bonds[f].Tag:=Tag;
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

procedure TMolecule.DeleteEmptyGroups;

var
  tmp:TMolecules;
  f:Integer;

begin
  tmp:=nil;
  for f:=0 to High(FGroups) do
    if FGroups[f].AtomCount>0 then
      begin
      FGroups[f].DeleteEmptyGroups;
      SetLength(tmp,Length(tmp)+1);
      tmp[High(tmp)]:=FGroups[f];
      end
    else FGroups[f].Free;
  FGroups:=tmp;
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

procedure TMolecule.Transform(TranslationVec: TCoord);

var
  ats:TAtoms;
  f:Integer;

begin
  ats:=AllAtoms;
  for f:=0 to High(ats) do
    ats[f].Coords:=Add(ats[f].Coords,TranslationVec);
end;

procedure TMolecule.Transform(RotationMat: TRotMatrix);

var
  ats:TAtoms;
  f:Integer;

begin
  ats:=AllAtoms;
  for f:=0 to High(ats) do
    ats[f].Coords:=Rotate(ats[f].Coords,RotationMat);
end;

procedure TMolecule.Transform(Quat: TQuaternion);

var
  ats:TAtoms;
  f:Integer;

begin
  ats:=AllAtoms;
  for f:=0 to High(ats) do
    ats[f].Coords:=Rotate(ats[f].Coords,Quat);
end;


procedure TMolecule.Transform(TranslationVec: TCoord; RotationMat: TRotMatrix);

var
  ats:TAtoms;
  f:Integer;

begin
  ats:=AllAtoms;
  for f:=0 to High(ats) do
    begin
    ats[f].Coords:=Add(ats[f].Coords,TranslationVec);
    ats[f].Coords:=Rotate(ats[f].Coords,RotationMat);
    end;
end;

procedure TMolecule.Transform(Center: TCoord; Rotation: TRotMatrix;
  Translation: TCoord);

var
  ats:TAtoms;
  f:Integer;

begin
  ats:=AllAtoms;
  for f:=0 to High(ats) do
    begin
    ats[f].Coords:=Subtract(ats[f].Coords,Center);
    ats[f].Coords:=Rotate(ats[f].Coords,Rotation);
    ats[f].Coords:=Add(ats[f].Coords,Translation);
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

