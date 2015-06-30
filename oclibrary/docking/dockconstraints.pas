{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain (see README.TXT).
********************************************************************************
Author: Ludwig Krippahl
Date: 25 01 2014
Purpose:
  Special constraint propagation for constrained docking
Requirements:
Revisions:
To do:

*******************************************************************************}

unit dockconstraints;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, linegrids, geomutils, docktasks;

type



  { TDockConstraintManager }

  TDockConstraintManager=class
  private
    FProbeTransVec,FTargetTransVec:TCoord;

    FTargetResolution,FProbeResolution:TFloat;
      //Resolutions should be identical, but best not hard code that

    //Base domain block
    FDomainBlock:TDomainBlock;
    FXDomain:TGridLine;
    FYDomainAtX:TLineArray;
    FZDomainAtXY:TGridPlane;

    procedure ComputeBaseBlock(Probe,Target:TDockingGrid);

    { TODO : Implement constraints... }



  public
    property DomainBlock:TDomainBlock read FDomainBlock;
    constructor Create(Probe,Target:TDockingGrid);
    property XDomain:TGridLine read FXDomain;
    function YDomainAtX(DomainX:Integer):TGridLine;
    function ZDomainAtXY(DomainX,DomainY:Integer):TGridLine;
    procedure LinearDistanceConstraint(const ATargetPoints,AProbePoints:TCoords;Dist:TFloat);
    procedure EuclideanDistanceConstraint(const ATargetPoints,AProbePoints:TCoords;Dist:TFloat);
    procedure PlaneConstraint(const Point,Normal:TCoord;const Margin:TFloat);
    procedure ImportConstraints(Constraints:TConstraintDefs;ProbeRot:TQuaternion; ProbeAxis:TCoord);

  end;

implementation

{ TDockConstraintManager }

procedure TDockConstraintManager.ComputeBaseBlock(Probe, Target: TDockingGrid);

var
  x,y:Integer;

begin
  FProbeTransVec:=Probe.TransVec;
  FTargetTransVec:=Target.TransVec;
  FTargetResolution:=Target.Resolution;
  FProbeResolution:=Probe.Resolution;
  with FDomainBlock do
    begin
      ProbeEndX:=High(Probe.Surf.Grid);
      ProbeEndY:=High(Probe.Surf.Grid[0]);
      ProbeEndZ:=Probe.Surf.ZMax;
      TargetEndX:=High(Target.Surf.Grid);
      TargetEndY:=High(Target.Surf.Grid[0]);
      TargetEndZ:=Target.Surf.ZMax;
      XOffset:=-ProbeEndX;
      YOffset:=-ProbeEndY;
      ZOffset:=-ProbeEndZ;
      DomainEndX:=TargetEndX-XOffset;
      DomainEndY:=TargetEndY-YOffset;
      DomainEndZ:=TargetEndZ-ZOffset;

      FXDomain:=NewLine(0,DomainEndX);
      SetLength(FYDomainAtX,DomainEndX+1);
      SetLength(FZDomainAtXY,DomainEndX+1,DomainEndY+1);
      for x:=0 to High(FYDomainAtX) do
        begin
        FYDomainAtX[x]:=NewLine(0,DomainEndY);
        for y:=0 to High(FZDomainAtXY[x]) do
          FZDomainAtXY[x,y]:=NewLine(0,DomainEndZ);
        end;
    end;

end;

constructor TDockConstraintManager.Create(Probe, Target: TDockingGrid);
begin
  inherited Create;
  ComputeBaseBlock(Probe,Target);
end;

function TDockConstraintManager.YDomainAtX(DomainX: Integer): TGridLine;
begin
  Result:=FYDomainAtX[DomainX];
end;

function TDockConstraintManager.ZDomainAtXY(DomainX, DomainY: Integer
  ): TGridLine;
begin
  Result:=FZDomainAtXY[DomainX,DomainY];
end;

procedure TDockConstraintManager.LinearDistanceConstraint(const ATargetPoints,
  AProbePoints: TCoords; Dist: TFloat);

var
  targetpoints,probepoints:TCoords;
  mit,mat,mip,map:TCoord;
  lx1,lx2,ly1,ly2,lz1,lz2:Integer;
  f,oldix,xx,x,y,yy:Integer;
  d:TFloat;
  tres,pres:TFloat;
begin
  targetpoints:=Add(FTargetTransvec,ATargetPoints);
  probepoints:=Add(FProbeTransvec,AProbePoints);
  tres:=1/FTargetResolution;
  pres:=1/FProbeResolution;

  mit:=Multiply(Min(targetpoints),tres);
  mat:=Multiply(Max(targetpoints),tres);
  mip:=Multiply(Min(probepoints),pres);
  map:=Multiply(Max(probepoints),pres);
  d:=Dist/FTargetResolution;

  with FDomainBlock do
    begin
    lx1:=Round(mit[0]-map[0]-XOffset-0.5-d);
    lx2:=Round(mat[0]-mip[0]-XOffset+0.5+d);
    ly1:=Round(mit[1]-map[1]-YOffset-0.5-d);
    ly2:=Round(mat[1]-mip[1]-YOffset+0.5+d);
    lz1:=Round(mit[2]-map[2]-ZOffset-0.5-d);
    lz2:=Round(mat[2]-mip[2]-ZOffset+0.5+d);
    end;
  Intersect(FXDomain,lx1,lx2);

  //clear FYDomainAtX outside FXDomain
  oldix:=0;
  for x:=0 to High(FXDomain) do
    begin
    for xx:=oldix to FXDomain[x,0]-1 do
      FYDomainAtX[xx]:=nil;
    oldix:=FXDomain[x,1]+1;
    end;
  for xx:=oldix to High(FYDomainAtX) do
    FYDomainAtX[xx]:=nil;


  //Set FYDomainAtX
  for x:=0 to High(FXDomain) do
    for xx:=FXDomain[x,0] to FXDomain[x,1] do
      begin
      Intersect(FYDomainAtX[xx],ly1,ly2);
      //clear XDomainAtXY outside FYDomainAtX
      oldix:=0;
      for y:=0 to High(FYDomainAtX[xx]) do
        begin
        for yy:=oldix to FYDomainAtX[xx,y,0]-1 do
          FZDomainAtXY[xx,yy]:=nil;
        oldix:=FYDomainAtX[xx,y,1]+1;
        end;
      for yy:=oldix to High(FZDomainAtXY[xx]) do
        FZDomainAtXY[xx,yy]:=nil;

      //set FZDomainAtXY
      for y:=0 to High(FYDomainAtX[xx]) do
        for yy:=FYDomainAtX[xx,y,0] to FYDomainAtX[xx,y,1] do
          Intersect(FZDomainAtXY[xx,yy],lz1,lz2);

      end;

end;

procedure TDockConstraintManager.EuclideanDistanceConstraint(
  const ATargetPoints, AProbePoints: TCoords; Dist: TFloat);

var
  targetpoints,probepoints:TCoords;
  mit,mat,mip,map:TCoord;
  lx1,lx2,ly1,ly2,lz1,lz2:Integer;
  f,oldix,xx,x,y,yy:Integer;
  griddist:TFloat;
  distsquare:TFloat;
  tres,pres:TFloat;

  procedure SetYLimits(XCoord:Integer);

  var
    f,g:Integer;
    dif,dy:TFloat;

  begin
    ly1:=FDomainBlock.DomainEndY;
    ly2:=0;

    for f:=0 to High(targetpoints) do
      for g:=0 to High(probepoints) do
        begin
        dif:=Abs(probepoints[g,0]+XCoord-targetpoints[f,0]);
        if dif<=griddist then
          begin
          dy:=Sqrt(distsquare-Sqr(dif));
          dif:=targetpoints[f,1]-probepoints[g,1];
          ly1:=Min(ly1,Round(dif-0.5-dy));
          ly2:=Max(ly2,Round(dif+0.5+dy));
          end;
     end;
  end;

  procedure SetZLimits(XCoord,YCoord:Integer);
  var
    f,g:Integer;
    dif,dz:TFloat;

  begin
    lz1:=FDomainBlock.DomainEndZ;
    lz2:=0;
    for f:=0 to High(targetpoints) do
      for g:=0 to High(probepoints) do
        begin
        dif:=Abs(probepoints[g,0]+XCoord-targetpoints[f,0]);
        if dif<=griddist then
          begin
          dif:=Sqr(dif)-Sqr(probepoints[g,1]+YCoord-targetpoints[f,1]);
          if dif<=distsquare then
            begin
            dz:=Sqrt(distsquare-dif);
            dif:=targetpoints[f,2]-probepoints[g,2];
            lz1:=Min(lz1,Round(dif-0.5-dz));
            lz2:=Max(lz2,Round(dif+0.5+dz));
            end;
          end;
       end;
  end;


begin
  tres:=1/FTargetResolution;
  pres:=1/FProbeResolution;

  targetpoints:=Multiply(Add(FTargetTransvec,ATargetPoints),tres);
  with FDomainBlock do
    probepoints:=Add(Coord(XOffset,YOffset,ZOffset),Multiply(Add(FProbeTransvec,AProbePoints),pres));

  mit:=Min(targetpoints);
  mat:=Max(targetpoints);
  mip:=Min(probepoints);
  map:=Max(probepoints);
  griddist:=Dist/FTargetResolution;
  distsquare:=Sqr(griddist);

  lx1:=Round(mit[0]-map[0]-0.5-griddist);
  lx2:=Round(mat[0]-mip[0]+0.5+griddist);
  Intersect(FXDomain,lx1,lx2);

  //clear FYDomainAtX outside FXDomain
  oldix:=0;
  for x:=0 to High(FXDomain) do
    begin
    for xx:=oldix to FXDomain[x,0]-1 do
      FYDomainAtX[xx]:=nil;
    oldix:=FXDomain[x,1]+1;
    end;
  for xx:=oldix to High(FYDomainAtX) do
    FYDomainAtX[xx]:=nil;


  //Set FYDomainAtX
  for x:=0 to High(FXDomain) do
    for xx:=FXDomain[x,0] to FXDomain[x,1] do
      begin
      SetYLimits(xx);
      Intersect(FYDomainAtX[xx],ly1,ly2);

      //clear XDomainAtXY outside FYDomainAtX
      oldix:=0;
      for y:=0 to High(FYDomainAtX[xx]) do
        begin
        for yy:=oldix to FYDomainAtX[xx,y,0]-1 do
          FZDomainAtXY[xx,yy]:=nil;
        oldix:=FYDomainAtX[xx,y,1]+1;
        end;
      for yy:=oldix to High(FZDomainAtXY[xx]) do
        FZDomainAtXY[xx,yy]:=nil;


      //set FZDomainAtXY
      for y:=0 to High(FYDomainAtX[xx]) do
        for yy:=FYDomainAtX[xx,y,0] to FYDomainAtX[xx,y,1] do
          begin
          SetZLimits(xx,yy);
          Intersect(FZDomainAtXY[xx,yy],lz1,lz2);
          end;

      end;

end;

procedure TDockConstraintManager.PlaneConstraint(const Point, Normal: TCoord;
  const Margin: TFloat);
//NOTE: Normal vector must be normalized.

var
  planepoint,linetoplane:TCoord;
  xflags,yflags,zflags:TIntegers;

  griddist:TFloat;
  dist:TFloat;
  tres,pres:TFloat;

  procedure SetZDomain(xx,yy:Integer);
  //computes distance from point to plane for each point.
  //this seems necessary to give the proper margin...

  var zz,z:Integer;

  begin
    linetoplane:=Coord(planepoint[0]-xx,planepoint[1]-yy,planepoint[2]);
    for z:=0 to High(FZDomainAtXY[xx,yy]) do
      for zz:=FZDomainAtXY[xx,yy,z,0] to FZDomainAtXY[xx,yy,z,1] do
        begin
        linetoplane[2]:=planepoint[2]-zz;
        dist:=Abs(DotProduct(linetoplane,Normal));
        if dist<=Margin then
          zflags[zz]:=1;
        end;
    FZDomainAtXY[xx,yy]:=IntegersToLine(zflags);
    //cleanup zflags
    SetIntegersLine(FZDomainAtXY[xx,yy],zflags,-1);
  end;

  procedure SetYDomain(xx:Integer);

  var y,yy:Integer;

  begin
    for y:=0 to High(FYDomainAtX[xx]) do
      for yy:=FYDomainAtX[xx,y,0] to FYDomainAtX[xx,y,1] do
        begin
        SetZDomain(xx,yy);
        if FZDomainAtXY[xx,yy]<>nil then
          yflags[yy]:=1;
        end;
    FYDomainAtX[xx]:=IntegersToLine(yflags);
    //cleanup yflags
    SetIntegersLine(FYDomainAtX[xx],yflags,-1);
  end;

var xx,x:Integer;

begin

  //scaling factors for the shapes
  tres:=1/FTargetResolution;
  pres:=1/FProbeResolution;

  //plane and line variables
  planepoint:=Subtract(Add(Point,Multiply(FTargetTransvec,tres)),Multiply(FProbeTransVec,pres));
  planepoint:=Subtract(planepoint,Coord(FDomainBlock.XOffset,FDomainBlock.YOffset,FDomainBlock.ZOffset));

  xflags:=FilledInts(FDomainBlock.DomainEndX+1,-1);
  yflags:=FilledInts(FDomainBlock.DomainEndY+1,-1);
  zflags:=FilledInts(FDomainBlock.DomainEndZ+1,-1);

  for x:=0 to High(FXDomain) do
    for xx:=FXDomain[x,0] to FXDomain[x,1] do
      begin
      SetYDomain(xx);
      if FYDomainAtX[xx]<>nil then
        xflags[xx]:=1;
      end;
  FXDomain:=IntegersToLine(xflags);

end;

procedure TDockConstraintManager.ImportConstraints(
  Constraints: TConstraintDefs; ProbeRot: TQuaternion; ProbeAxis:TCoord);

var
  f:Integer;

begin
  for f:=0 to High(Constraints) do
    case Constraints[f].ConstraintType of
      CEuclideanDistance:
        EuclideanDistanceConstraint(Constraints[f].TargetPoints,
                Rotate(Constraints[f].ProbePoints,ProbeRot),Constraints[f].Distance);
      CNormalPlane:PlaneConstraint(Coord(0,0,0),ProbeAxis,Constraints[f].Distance);
    end;
end;


end.

