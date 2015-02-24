{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 4.12.2014
Purpose:
  Functions for aligning molecules by sequence and structure.
Requirements:
Revisions:
To do:
*******************************************************************************}
unit molfit;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils,basetypes, geomutils, alignment, rmsd, pdbmolecules, molecules;

type


  TFitResult=record
    Rotation:TRotMatrix;
    Translation,Center:TCoord;
    Rmsd:TFloat;
    AlignmentScores:TFloats;
    Map:TIntegers;
      //which chain of Probe corresponds to which of probe
      //nil if fit failed
    ResidueMaps:array of TIntegers;
      //one map per probe chain indicating corresponding residues in target
  end;

  TFitCoordsResult=record
    ProbeCoords,TargetCoords:TCoords;
    FitResult:TFitResult;
  end;

function AlignSequences(Mol1,Mol2:TMolecule;SubMat:TSubMatrix;MissingMarker:string='X'):TPairwiseAlignment;

function MagicFit(const Probe,Target:TMolecule;
                  const SubMat:TSubMatrix;
                  TargetExclusion:TIntegers=nil;
                  MinSequenceMatch:TFloat=0.7;
                  MinAlfaCarbons:Integer=20):TFitResult;
{Computes the fits from Probe to Target.
 First matches the probe chains to the avaliable (not in TargetExclusion) target
 chains by sequence (best first). Then computes RMSD of matched chains and returns the fit}

function MagicFitAlphaCoords(const Probe,Target:TMolecule;
                             const SubMat:TSubMatrix;
                             TargetExclusion:TIntegers=nil;
                             MinSequenceMatch:TFloat=0.7):TFitCoordsResult;
{Fills ProbeCoords and TargetCoords with the alpha carbon coordinates of the
matching residues after a MagicFit with the values given}
procedure GetMatchingResidues(const Probe,Target:TMolecule;
                             const FitResult:TFitResult;
                             out ProbeRes,TargetRes:TMolecules);
{returns all residues from the two molecules that were matched in the fit result}

function CountResidueMatches(const Fit:TFitResult):Integer;


implementation

function AlignSequences(Mol1,Mol2:TMolecule;SubMat:TSubMatrix;MissingMarker:string='X'):TPairwiseAlignment;

var
  seq1,seq2:string;
  seqixs1,seqixs2:TIntegers;

begin
  seq1:=ChainSequence(Mol1,MissingMarker);
  seq2:=ChainSequence(Mol2,MissingMarker);
  Result:=NeedlemanWunschAlign(seq1,seq2,SubMat);
end;

function MagicFit(const Probe,Target:TMolecule;
                  const SubMat:TSubMatrix;
                  TargetExclusion:TIntegers;
                  MinSequenceMatch:TFloat;
                  MinAlfaCarbons:Integer=20):TFitResult;

var
  alignments:array of array of TPairwiseAlignment;
  //Target chains in first dimension
  //Probe chains in second dimension
  //All alignments are Probe to Target
  map:TIntegers; //target ix for each probe chain
  bestmap:TIntegers;

  procedure BuildAlignmentMatrix;

  var
    f,cp,ct:Integer;
    chp,cht:TMolecule;
    score:TFloat;

  begin
    SetLength(alignments,Target.GroupCount,Probe.GroupCount);
    for ct:=0 to Target.GroupCount-1 do
      if IsInArray(ct,TargetExclusion) then
        for f:=0 to High(alignments[ct]) do
          begin
          alignments[ct,f].Score:=-1;
          alignments[ct,f].Map:=nil;
          end
      else
        begin
        cht:=Target.GetGroup(ct);
        for cp:=0 to Probe.GroupCount-1 do
          begin
          chp:=Probe.GetGroup(cp);
          alignments[ct,cp]:=AlignSequences(chp,cht,SubMat,'');
          //convert score to percent match. Empty missing marker to ignore missing residues
          score:=0;
          if alignments[ct,cp].Map<>nil then
            begin
            for f:=0 to High(alignments[ct,cp].Map) do
              begin
              if alignments[ct,cp].Map[f]>=0 then
                score:=score+1;
              end;
            //alignments[ct,cp].Score:=score/Min(cht.GroupCount,chp.GroupCount);
            alignments[ct,cp].Score:=score/(0.5*(cht.GroupCount+chp.GroupCount));
            end;
          end;
        end;
  end;

  procedure TryFit;

  var
    tar,prb:TCoords;
    totlen,currpos,f,r:Integer;
    tarchain,prbchain:TMolecule;
    tca,pca:TAtom;
    rmsdcalc:TRMSDCalculator;

    tmp:string;

  begin
    totlen:=0;
    for f:=0 to High(map) do
      if map[f]>=0 then
        totlen:=totlen+Length(alignments[map[f],f].Map);
    SetLength(tar,totlen);
    SetLength(prb,totlen);
    currpos:=0;
    for f:=0 to High(map) do
      if map[f]>=0 then
        begin
        tarchain:=Target.GetGroup(map[f]);
        prbchain:=Probe.GetGroup(f);
        for r:=0 to High(alignments[map[f],f].Map) do
          if alignments[map[f],f].Map[r]>=0 then
            begin
            tca:=tarchain.GetGroup(alignments[map[f],f].Map[r]).GetAtom('CA');
            pca:=prbchain.GetGroup(r).GetAtom('CA');
            if (tca<>nil) and (pca<>nil) then
              begin
              tar[currpos]:=tca.Coords;
              prb[currpos]:=pca.Coords;
              Inc(currpos);
              end;
            end;
        end;
    SetLength(tar,currpos);
    SetLength(prb,currpos);


    //check if at least MinAlfaCarbons atoms were used
    if currpos<MinAlfaCarbons then
      Exit;

    rmsdcalc:=TRMSDCalculator.Create;
    rmsdcalc.AddCoordinates(tar,prb);
    rmsdcalc.Minimise(1000,0.01);

    tmp:='';
    for f:=0 to High(map) do
      tmp:=tmp+IntToStr(map[f])+' ';

    if (Result.Rmsd<0) or (Result.Rmsd>rmsdcalc.Rmsd) then
      begin
      Result.Rmsd:=rmsdcalc.Rmsd;
      Result.Center:=rmsdcalc.CenterVec;
      Result.Rotation:=rmsdcalc.PlacedRotation;
      Result.Translation:=rmsdcalc.PlacedTranslation;

      for f:=0 to High(map) do
        begin
        if map[f]>=0 then
          begin
          Result.AlignmentScores[f]:=alignments[map[f],f].Score;
          Result.ResidueMaps[f]:=alignments[map[f],f].Map;
          end
        else
          begin
          Result.AlignmentScores[f]:=-1;
          Result.ResidueMaps[f]:=nil;
          end;
        Result.Map[f]:=map[f];
        end;
      end;
    rmsdcalc.Free;
  end;

  procedure HighestScore(out Tix,Pix:Integer; out Score:TFLoat);

  var f,g:Integer;

  begin
    Tix:=0;
    Pix:=0;
    Score:=alignments[0,0].Score;
    for f:=0 to High(alignments) do
      for g:=0 to High(alignments[0]) do
        if alignments[f,g].Score>Score then
          begin
          Tix:=f;
          Pix:=g;
          Score:=alignments[f,g].Score;
          end;
  end;

  procedure WipeRowCol(Tix,Pix:Integer);

  var f:Integer;

  begin
    for f:=0 to High(alignments) do
      alignments[f,Pix].Score:=-1;
    for f:=0 to High(alignments[0]) do
      alignments[Tix,f].Score:=-1;
  end;

  procedure FindBestMap(Level:Integer);

  var
    tix:Integer;
    score:TFloat;
    foundone:Boolean;

  begin
    if Level>High(map) then
      TryFit
    else
      begin
      foundone:=False;
      for tix:=0 to High(alignments) do
        begin
        if alignments[tix,Level].Score>=MinSequenceMatch then
          begin
          if (Level=0) or ((Level>0) and (not IsInArray(tix,Copy(map,0,Level)))) then
            foundone:=True;
          map[Level]:=tix;
          FindBestMap(Level+1);
          end;
        end;
      if not foundone then
        begin
        map[Level]:=-1;
        FindBestMap(Level+1);
        end;
      end;
  end;

begin
  Result.Rmsd:=-1;
  BuildAlignmentMatrix;
  map:=FilledInts(Probe.GroupCount,-1);
  SetLength(Result.AlignmentScores,Probe.GroupCount);
  SetLength(Result.Map,Probe.GroupCount);
  SetLength(Result.ResidueMaps,Probe.GroupCount);
  FindBestMap(0);
end;

function MagicFitAlphaCoords(const Probe,Target:TMolecule;
                             const SubMat:TSubMatrix;
                             TargetExclusion:TIntegers=nil;
                             MinSequenceMatch:TFloat=0.7):TFitCoordsResult;

var
  f,g,count,currix:Integer;
  pchain,tchain,pres,tres:TMolecule;
  tca,pca:TAtom;

begin
  Result.FitResult:=MagicFit(Probe, Target, SubMat, TargetExclusion, MinSequenceMatch);
  with Result do
    begin
    count:=0;
    for f:=0 to High(FitResult.ResidueMaps) do
      begin
      count:=count+Length(FitResult.ResidueMaps[f]);
      end;
    SetLength(ProbeCoords,count);
    SetLength(TargetCoords,count);
    currix:=0;
    for f:=0 to High(FitResult.Map) do
      begin
      if FitResult.Map[f]>=0 then
        begin
        pchain:=Probe.GetGroup(f);
        tchain:=Target.GetGroup(FitResult.Map[f]);
        if (pchain<>nil) and (tchain<>nil) then  //this should not be necessary...?
          for g:=0 to High(FitResult.ResidueMaps[f]) do
            if FitResult.ResidueMaps[f,g]>=0 then
              begin
              pres:=pchain.GetGroup(g);
              tres:=tchain.GetGroup(FitResult.ResidueMaps[f,g]);
              if (pres<>nil) and (tres<>nil) then
                begin
                pca:=pres.GetAtom('CA');
                tca:=tres.GetAtom('CA');
                if (pca<>nil) and (tca<>nil) then
                  begin
                  ProbeCoords[currix]:=pca.Coords;
                  TargetCoords[currix]:=tca.Coords;
                  Inc(currix);
                  end;
                end;
              end;
        end;
      end;
    SetLength(ProbeCoords,currix);
    SetLength(TargetCoords,currix);
    end;
end;

procedure GetMatchingResidues(const Probe, Target: TMolecule;
  const FitResult: TFitResult; out ProbeRes, TargetRes: TMolecules);

var
  f,g,count,ix:Integer;
  targc,probec:TMolecule;

begin
  count:=0;
  for f:=0 to High(FitResult.ResidueMaps) do
    count:=count+Length(FitResult.ResidueMaps[f]);
  SetLength(ProbeRes,count);
  SetLength(TargetRes,count);
  ix:=0;
  for f:=0 to High(FitResult.ResidueMaps) do
    if FitResult.Map[f]>=0 then
      begin
      targc:=Target.GetGroup(FitResult.Map[f]);
      probec:=Probe.GetGroup(f);
      for g:=0 to High(FitResult.ResidueMaps[f]) do
        if FitResult.ResidueMaps[f,g]>=0 then
          begin
          ProbeRes[ix]:=probec.GetGroup(g);
          TargetRes[ix]:=targc.GetGroup(FitResult.ResidueMaps[f,g]);
          Inc(ix);
          end;
      end;
  SetLength(ProbeRes,ix);
  SetLength(TargetRes,ix);

end;

function CountResidueMatches(const Fit: TFitResult): Integer;

var f,g:Integer;

begin
  Result:=0;
  for f:=0 to High(Fit.ResidueMaps) do
    for g:=0 to High(Fit.ResidueMaps[f]) do
      if Fit.ResidueMaps[f,g]>=0 then Inc(Result);
end;

end.

