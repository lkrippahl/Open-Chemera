{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 25 01 2014
Purpose:
  Domains for geometric docking
Requirements:
Revisions:
To do:

*******************************************************************************}
unit dockdomains;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, linegrids, geomutils, dockconstraints;

type
  TAddModel=function (Score,X,Y,Z:Integer):Integer of object;
  { TDockDomain }

  TDockDomain=class
    private
      FDomainGrid:TDomainGrid;
      FCoreContacts:TLineArray;
        //X oriented array of gridlines where core contacts can occur in the domain
      FTarget,FProbe:TDockingGrid;
        //Target and Probe grids are NOT managed by this class

      //for multigrid speedup
      FCoarseDomain:TDockDomain;
      FCoarseTarget,FCoarseProbe:TDockingGrid;
        //Coarse grids and domain are created and managed by this class (frees
        //coarse grids after freeing coarse domain)
      FConstraintManager:TDockConstraintManager;

      // for InPlaceIntersect
      FIntersect:TGridLine;
      FIntersectHigh:Integer;

      // For scoring and overlaping without creating integer arrays
      FScores:TIntegers;
      FOverlapixs:TIntegers;

      FAddModel:TAddModel;

      procedure InPlaceIntersect(const Line1, Line2: TGridLine;const Displace:Integer);
        //uses FIntersect and FIntersectHigh to avoid creating new lines
        //this removes about 10% in domain generation time (7.55 to 6.7 s tested on 6k atoms)

      procedure SetDomainYLines;
      procedure SetDomainZLines;
      procedure SetZExtremes(DomainX,DomainY:Integer);
      procedure RemoveZCores(DomainX, DomainY: Integer);
      procedure ScoreOverlap(DomainX, DomainY: Integer);
    public
      RemoveCores:Boolean;
      MinimumOverlap:Integer;
      CountXY,CountXYCore,CountX:Integer;
      property ConstraintManager:TDockConstraintManager read FConstraintManager;
      property Grid:TGridPlane read FDomainGrid.Shape.Grid;
      property Size:Integer read FDomainGrid.Shape.TotalCount;
      function TranslateToTarget:TCoord;
      constructor Create(ATarget,AProbe:TDockingGrid;MaxGrids:Integer=0);
      procedure AssignModelManager(AddModel:TAddModel);
      procedure BuildInitialDomain({Constraints:TDockConstraints=nil});
      procedure CalcDomainStats;
      procedure Score;
      procedure Free;

  end;

implementation

procedure OverlapRegion(const ProbeIx1,ProbeIx2,TargetIx1,TargetIx2,ProbeOffset,DomainEnd:Integer;
              out DomainIx1,DomainIx2:Integer);

begin
  DomainIx1:=TargetIx1-(ProbeIx2+ProbeOffset);
  DomainIx2:=DomainIx1+(ProbeIx2-ProbeIx1)+(TargetIx2-TargetIx1);
  if DomainIx1<0 then DomainIx1:=0;
  if DomainIx2>DomainEnd then DomainIx2:=DomainEnd;
end;


procedure GetIndexes(const ProbeEnd,TargetEnd,ProbeOffset,Coord:Integer;
    out ProbeIx1,ProbeIx2,TargetIx1,TargetIx2:Integer);
   //Gets first and last index of target and probe grids for this domain coord

begin
  ProbeIx1:=-ProbeOffset-Coord;
  ProbeIx2:=ProbeEnd;
  TargetIx1:=0;
  TargetIx2:=TargetIx1+ProbeIx2-ProbeIx1;

  if ProbeIx1<0 then
    begin
    TargetIx1:=-ProbeIx1;
    ProbeIx1:=0;
    end;

  If TargetIx2>TargetEnd then
    begin
    ProbeIx2:=ProbeIx2-TargetIx2+TargetEnd;
    TargetIx2:=TargetEnd;
    end;

end;

{ TDockDomain }

procedure TDockDomain.InPlaceIntersect(const Line1, Line2: TGridLine;const Displace:Integer);

var
  ni1,ni2:Integer;
  i1,i2,ll1,ll2:Integer;
  top,bot:Integer;
begin
  FIntersectHigh:=-1;
  if (Line1<>nil) and (Line2<>nil) then
    begin
    ll1:=Length(Line1)+Length(Line2);
    if Length(FIntersect)<ll1 then
      SetLength(FIntersect,ll1);
    i1:=0;
    i2:=0;
    ll1:=Length(Line1);
    ll2:=Length(Line2);
    while (i1<ll1) and (i2<ll2) do
      begin
      top:=Min(Line1[i1,1]+Displace,Line2[i2,1]);
      bot:=Max(Line1[i1,0]+Displace,Line2[i2,0]);
      if top>=bot then
        begin
        Inc(FIntersectHigh);
        FIntersect[FIntersectHigh,0]:=bot;
        FIntersect[FIntersectHigh,1]:=top;
        end;
      if Line1[i1,1]+Displace>=Line2[i2,1] then ni2:=i2+1 else ni2:=i2;
      if Line1[i1,1]+Displace<=Line2[i2,1] then ni1:=i1+1 else ni1:=i1;
      i1:=ni1;
      i2:=ni2;
      end;
    end;
end;

procedure TDockDomain.SetDomainYLines;

var
  x,xx:Integer;
  xdomain:TGridLine;
  tmpmin,tmpmax,mi,ma,p1,p2,t1,t2:Integer;
  coremin,coremax:Integer;
  limit1,limit2:Integer;
  pix1,pix2,tix1,tix2:Integer;

begin
  with FDomainGrid do
    begin
    SetLength(Shape.NonEmpty,Block.DomainEndX+1);
    SetLength(FCoreContacts,Block.DomainEndX+1);
    for x:=0 to High(Shape.NonEmpty) do
      begin
      Shape.NonEmpty[x]:=nil;
      FCoreContacts[x]:=nil;
      end;
    xdomain:=FConstraintManager.XDomain;
    for xx:=0 to High(xdomain) do
      for x:=xdomain[xx,0] to xdomain[xx,1] do
      begin
      Inc(CountX);
      GetLineExtremes(FConstraintManager.YDomainAtX(x),limit1,limit2);
      if limit1>limit2 then
        begin
        Shape.NonEmpty[x]:=nil;
        FCoreContacts[x]:=nil;
        end
      else
        begin
        GetIndexes(Block.ProbeEndX,Block.TargetEndX,Block.XOffset,x,pix1,pix2,tix1,tix2);
        tmpmin:=limit2+1;
        tmpmax:=-1;
        coremin:=limit2+1;
        coremax:=-1;
        while tix1<=tix2 do
          begin
          //Surface overlap Y range
          if (FProbe.Surf.NonEmpty[pix1]<>nil) and (FTarget.Surf.NonEmpty[tix1]<>nil) then
            begin
            GetLineExtremes(FProbe.Surf.NonEmpty[pix1],p1,p2);
            GetLineExtremes(FTarget.Surf.NonEmpty[tix1],t1,t2);
            OverlapRegion(p1,p2,t1,t2,Block.YOffset,limit2,mi,ma);
            tmpmin:=Min(mi,tmpmin);
            tmpmax:=Max(ma,tmpmax);
            end;
          //Core overlap Y range
          if (FProbe.Core.NonEmpty[pix1]<>nil) and (FTarget.Core.NonEmpty[tix1]<>nil) then
            begin
            GetLineExtremes(FProbe.Core.NonEmpty[pix1],p1,p2);
            GetLineExtremes(FTarget.Core.NonEmpty[tix1],t1,t2);
            OverlapRegion(p1,p2,t1,t2,Block.YOffset,limit2,mi,ma);
            coremin:=Min(mi,coremin);
            coremax:=Max(ma,coremax);
            end;
          Inc(tix1);
          Inc(pix1);
          end;
        tmpmin:=Max(limit1,tmpmin);
        tmpmax:=Min(limit2,tmpmax);
        coremin:=Max(limit1,coremin);
        coremax:=Min(limit2,coremax);

        if tmpmax<tmpmin then
          Shape.NonEmpty[x]:=nil
        else
          Shape.NonEmpty[x]:=NewLine(tmpmin,tmpmax);

        if coremax<coremin then
          FCoreContacts[x]:=nil
        else
          FCoreContacts[x]:=NewLine(coremin,coremax);
        end;
      end;
    end;
end;

procedure TDockDomain.SetDomainZLines;

var
  x,y,yy:Integer;

begin
  //large buffer for intersections
  SetLength(FIntersect,100);
  FIntersectHigh:=-1;
  with FDomainGrid do
    begin
    //clear grid and set to nil to distinguish lines where core was processed
    SetLength(Shape.Grid,Block.DomainEndX+1,Block.DomainEndY+1);
    {for x:=0 to High(Shape.Grid) do
      for y:=0 to High(Shape.Grid[0]) do
        Shape.Grid[x,y]:=nil;}


    //Set lines with surface contacts
    FOverlapIxs:=FilledInts(2*FDomainGrid.Block.ProbeEndZ+FDomainGrid.Block.TargetEndZ,0);
      //FOverlapIxs used with accumulator, must span all range of surface, not actual domain
    for x:=0 to High(Shape.NonEmpty) do
      for y:=0 to High(Shape.NonEmpty[x]) do
        for yy:=Shape.NonEmpty[x,y,0] to Shape.NonEmpty[x,y,1] do
          SetZExtremes(x,yy);

    //Set lines with core overlaps
    if RemoveCores then
      begin
      FOverlapIxs:=FilledInts(FDomainGrid.Block.DomainEndZ+1,0);
      for x:=0 to High(FCoreContacts) do
        for y:=0 to High(FCoreContacts[x]) do
          for yy:=FCoreContacts[x,y,0] to FCoreContacts[x,y,1] do
            if Shape.Grid[x,yy]<>nil then
              RemoveZCores(x,yy);
      end;

    end;
end;

procedure TDockDomain.SetZExtremes(DomainX, DomainY: Integer);

var
  tmpmin,tmpmax,mi,ma,p1,p2,t1,t2:Integer;
  z,y,yy,limit1,limit2:Integer;
  pix1,pix2,tix1,tix2:Integer;
  displace:Integer;
  overlap,acc,dif:Integer;
  st,en,wid1,wid2:Integer;

begin
  Inc(CountXY);
  with FDomainGrid do
    begin
    for z:=0 to High(FOverlapIxs) do
      FOverlapIxs[z]:=0;
    GetLineExtremes(FConstraintManager.ZDomainAtXY(DomainX,DomainY),limit1,limit2);
    GetIndexes(Block.ProbeEndX,Block.TargetEndX,Block.XOffset,DomainX,pix1,pix2,tix1,tix2);
    tmpmin:=limit2+1;
    tmpmax:=-1;
    displace:=DomainY+Block.YOffset;
    overlap:=0;
    while tix1<=tix2 do
      begin
      InPlaceIntersect(FProbe.Surf.NonEmpty[pix1],FTarget.Surf.NonEmpty[tix1],displace);
      for y:=0 to FIntersectHigh do
        for yy:=FIntersect[y,0] to FIntersect[y,1] do
          begin
          overlap:=overlap+Min(Fprobe.Surf.CellCounts[pix1,yy-displace],
                               FTarget.Surf.CellCounts[tix1,yy]);
          GetLineExtremes(FProbe.Surf.Grid[pix1,yy-displace],p1,p2);
          GetLineExtremes(FTarget.Surf.Grid[tix1,yy],t1,t2);

          p1:=p1-Block.ProbeEndZ;
          p2:=p2-Block.ProbeEndZ;
          st:=t1-p2;
          en:=t2-p1;
          wid1:=p2-p1;
          wid2:=t2-t1;
          if wid1>wid2 then
            begin
            z:=wid1;
            wid1:=wid2;
            wid2:=z;
            end;
          wid1:=st+wid1;
          wid2:=st+wid2;
          Inc(FOverlapIxs[st]);
          Dec(FOverlapIxs[wid1+1]);
          Dec(FOverlapIxs[wid2+1]);
          Inc(FOverlapIxs[en+1]);
          end;
      Inc(tix1);
      Inc(pix1);
      end;
    if (overlap<MinimumOverlap) then
      Shape.Grid[DomainX,DomainY]:=nil
    else
      begin
      acc:=0;
      dif:=0;
      z:=-1;
      while (z<High(FOverlapIxs)) and (acc<=MinimumOverlap) do
        begin
        Inc(z);
        dif:=dif+FOverlapIxs[z];
        acc:=acc+dif;
        end;
      tmpmin:=z;
      z:=Length(FOverlapIxs);
      dif:=0;
      acc:=0;
      while (z>0) and (acc<=MinimumOverlap) and (z>tmpmin) do
        begin
        Dec(z);
        dif:=dif+FOverlapIxs[z];
        acc:=acc+dif;
        end;
      tmpmax:=z;
      //convert from grid to domain
      tmpmin:=tmpmin+Block.ProbeEndZ+Block.ZOffset;
      tmpmax:=tmpmax+Block.ProbeEndZ+Block.ZOffset;
      tmpmin:=Max(tmpmin,limit1);
      tmpmax:=Min(tmpmax,limit2);
      if tmpmin<=tmpmax then
        Shape.Grid[DomainX,DomainY]:=NewLine(tmpmin,tmpmax)
      else
        Shape.Grid[DomainX,DomainY]:=nil;
      end;
  end;
end;




procedure TDockDomain.RemoveZCores(DomainX, DomainY: Integer);
//assumes a domain is already created

var
  y,yy,acc,limit2:Integer;
  pix1,pix2,tix1,tix2:Integer;
  displace:Integer;
  zline:TGridLine;
  corefound:Boolean;
  pix,tix:Integer;
  pline,tline:TGridLine;
  zmi,zma:Integer;


begin
  Inc(CountXYCore);
  corefound:=False;
  with FDomainGrid do
    begin
    GetIndexes(Block.ProbeEndX,Block.TargetEndX,Block.XOffset,DomainX,pix1,pix2,tix1,tix2);
    displace:=DomainY+Block.YOffset;
    while tix1<=tix2 do
      begin
      if (FProbe.Core.NonEmpty[pix1]<>nil) and (FTarget.Core.NonEmpty[tix1]<>nil) then
        begin
        InPlaceIntersect(FProbe.Core.NonEmpty[pix1],FTarget.Core.NonEmpty[tix1],displace);
        for y:=0 to FIntersectHigh do
          for yy:=FIntersect[y,0] to FIntersect[y,1] do
            begin
            corefound:=True;
            pline:=FProbe.Core.Grid[pix1,yy-displace];
            tline:=FTarget.Core.Grid[tix1,yy];
            for pix:=0 to High(pline) do
              for tix:=0 to High(tline) do
                begin
                zmi:=tline[tix,0]-(pline[pix,1]+FDomainGrid.Block.ZOffset);
                zma:=zmi+pline[pix,1]-pline[pix,0]+tline[tix,1]-tline[tix,0];
                Dec(FOverlapIxs[zmi]);
                Inc(FOverlapIxs[zma+1]);
                end;
            end;
        end;
      Inc(tix1);
      Inc(pix1);
      end;
    if corefound then
      begin
      acc:=1;
      for yy:=0 to High(FOverlapIxs) do
        begin
        acc:=acc+FOverlapIxs[yy];
        FOverlapIxs[yy]:=acc;
        end;
      Shape.Grid[DomainX,DomainY]:=Intersect(Shape.Grid[DomainX,DomainY],IntegersToLine(FOverlapIxs));
      for yy:=0 to High(FOverlapIxs) do
        FOverlapIxs[yy]:=0;
      end;
  end;
end;


procedure TDockDomain.ScoreOverlap(DomainX, DomainY: Integer);


var
  z,zz,y,yy:Integer;
  pix1,pix2,tix1,tix2:Integer;
  zline:TGridLine;
  hizline:Integer;
  displace:Integer;

  procedure CalcOverlaps;

  var
    pix,tix,z,zz:Integer;
    pline,tline:TGridLine;
    pz1,pz2,tz1,tz2:Integer;
    ad,wid1,wid2,st,en,hi:Integer;

  begin
    pline:=FProbe.Surf.Grid[pix1,yy-displace];
    tline:=FTarget.Surf.Grid[tix1,yy];
    hi:=High(FScores);
    for pix:=0 to High(pline) do
      begin
      pz1:=pline[pix,0]+FDomainGrid.Block.ZOffset;
      pz2:=pline[pix,1]+FDomainGrid.Block.ZOffset;
      for tix:=0 to High(tline) do
        begin
        tz1:=tline[tix,0];
        tz2:=tline[tix,1];
        st:=tz1-pz2;
        en:=tz2-pz1;
        wid1:=-1;
        for z:=0 to hizline do
          if en<zline[z,0] then
              Break
          else if ((st<=zline[z,1]) and (en>=zline[z,0])) then
            begin

            if wid1<0 then
              begin
              wid1:=pz2-pz1;
              wid2:=tz2-tz1;
              if wid1>wid2 then
                begin
                ad:=wid1;
                wid1:=wid2;
                wid2:=ad;
                end;
              wid1:=st+wid1;
              wid2:=st+wid2;
              end;

            for zz:=zline[z,0] to zline[z,1] do
               if zz>=st then
                 begin
                 if zz<=wid1 then Inc(FScores[zz],zz-st+1)
                 else if zz<=wid2 then Inc(FScores[zz],wid1-st+1)
                 else if zz<=en then Inc(FSCores[zz],en-zz+1)
                 else Break;
                 end;
            end;
        end;
      end;
  end;


begin
  with FDomainGrid do
    begin
    zline:=Shape.Grid[DomainX,DomainY];
    hizline:=High(zline);
    //Reset relevant scores (kept at FScores to avoid creating array at each iteration)
    for z:=0 to High(zline) do
      for zz:=zline[z,0] to zline[z,1] do
        FScores[zz]:=0;

    GetIndexes(Block.ProbeEndX,Block.TargetEndX,Block.XOffset,DomainX,pix1,pix2,tix1,tix2);
    displace:=DomainY+Block.YOffset;
    while tix1<=tix2 do
      begin
      InPlaceIntersect(FProbe.Surf.NonEmpty[pix1],FTarget.Surf.NonEmpty[tix1],displace);
      for y:=0 to FIntersectHigh do
        for yy:=FIntersect[y,0] to FIntersect[y,1] do
          CalcOverlaps;
      Inc(tix1);
      Inc(pix1);
      end;
    for z:=0 to High(zline) do
      for zz:=zline[z,0] to zline[z,1] do
        if FScores[zz]>MinimumOverlap then
          MinimumOverlap:=FAddModel(FScores[zz],DomainX+Block.XOffset,
                                   DomainY+Block.YOffset,zz+Block.ZOffset);
  end;
end;

function TDockDomain.TranslateToTarget: TCoord;

var c:TCoord;

begin
  c[0]:=(FDomainGrid.Block.XOffset+FDomainGrid.Block.ProbeEndX/2)*FProbe.Resolution;
  c[1]:=(FDomainGrid.Block.YOffset+FDomainGrid.Block.ProbeEndY/2)*FProbe.Resolution;
  c[2]:=(FDomainGrid.Block.ZOffset+FDomainGrid.Block.ProbeEndZ/2)*FProbe.Resolution;
  Result:=Subtract(c,FTarget.TransVec);
end;

constructor TDockDomain.Create(ATarget, AProbe: TDockingGrid; MaxGrids: Integer);
begin
  inherited Create;
  FTarget:=ATarget;
  FProbe:=AProbe;
  RemoveCores:=True;
  MinimumOverlap:=0;
  FConstraintManager:=TDockConstraintManager.Create(FProbe,FTarget);
  FAddModel:=nil;
end;

procedure TDockDomain.AssignModelManager(AddModel: TAddModel);
begin
  FAddModel:=AddModel;
end;

procedure TDockDomain.BuildInitialDomain;
begin
  FDomainGrid.Block:=FConstraintManager.DomainBlock;
  FDomainGrid.Shape.Grid:=nil;
  FDomainGrid.Shape.NonEmpty:=nil;
  SetDomainYLines;
  SetDomainZLines;
end;

procedure TDockDomain.CalcDomainStats;
begin
  ComputeShapeStats(FDomainGrid.Shape);
end;

procedure TDockDomain.Score;

var x,y,yy:Integer;

begin
  if Assigned(FAddModel) then
    begin
    //large buffer for intersections
    SetLength(FIntersect,100);
    FIntersectHigh:=-1;
    FScores:=FilledInts(FDomainGrid.Block.DomainEndZ+1,0);
    with FDomainGrid do
      for x:=0 to High(Shape.NonEmpty) do
        for y:=0 to High(Shape.NonEmpty[x]) do
          for yy:=Shape.NonEmpty[x,y,0] to Shape.NonEmpty[x,y,1] do
            if Shape.Grid[x,yy]<>nil then
              ScoreOverlap(x,yy);
    end;
end;

procedure TDockDomain.Free;
begin
  FConstraintManager.Free;
  Inherited Free;
end;

end.

