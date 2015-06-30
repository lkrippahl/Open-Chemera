{*******************************************************************************
This source file is public domain. It is provided as is, with no explicit
or implied warranties regarding anything whatsoever.
********************************************************************************
Author: Ludwig Krippahl
Date: 2014 02 24
Purpose:
  Conjugated gradient annd linear minimization
Revisions:
TODO:
  This was adapted from old Chemera code on Delphi. Probably needs a revision
  and some rewriting
*******************************************************************************}

unit cgradient;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, basetypes, geomutils;

type
  TCGReport=procedure(Iteration:Integer;Value:TFloat) of object;
  TCGGetValueAndDerivative=function(var Point:TFloats;var Deriv:TFloats):TFloat of object;
  TCGGetValue=function(const Point:TFloats):TFloat of object;

  TCGMinimiser=class
  private
    FGetValue:TCGGetValue;
    FGetValueAndDerivative:TCGGetValueAndDerivative;
    FReport:TCGReport;
    FLineCounts:Integer;
    FMins,FMaxs:TFloats;
    FPrecision:TFloat;
  public
    property Report:TCGReport read FReport write FReport;
    property GetValueFunction:TCGGetValue
      read FGetValue write FGetValue;
    property ValueAndDerivativeFunction:TCGGetValueAndDerivative
      read FGetValueAndDerivative write FGetValueAndDerivative;
    function Minimize(MaxIterations:Integer;MinVariation:TFloat;
      var Point,Deriv:TFloats):TFloat;
    function LineMinimization(var Point:TFloats;const Direction:TFloats):TFloat;
    procedure Free;
    constructor Create(Mins:TFloats=nil;Maxs:TFloats=nil;Precision:TFloat=1e-8);
  end;

implementation

uses miscutils;

{ TCGMinimiser }

constructor TCGMinimiser.Create(Mins:TFloats=nil; Maxs: TFloats=nil;Precision:TFloat=1e-8);
begin
  inherited Create;
  FMins:=Copy(Mins,0,Length(Mins));
  FMaxs:=Copy(Maxs,0,Length(Mins));
  FPrecision:=Precision;
end;

procedure TCGMinimiser.Free;
begin
  inherited Free;
end;

function TCGMinimiser.LineMinimization(var Point: TFloats;
  const Direction: TFloats): TFloat;

const
     GOLD=1.61834;
     GLIMIT=100.0;
     TINY=1.0e-20;
     R=0.6180399;
     CR=1-R;
     TOL=1.0e-8;
     FIRSTSTEP=0.00001;

var ax,bx,cx:TFloat;
    xmin:TFloat;
    fa,fb,fc:TFloat;
    b:boolean;
    CurrentPoint:TFloats;
    Lmin,Lmax:TFloat;

procedure GetLineLimits;

var f:Integer;

begin
  if FMaxs<>nil then
    begin
    LMin:=1e30;
    LMax:=-1e30;
    for f:=0 to High(FMaxs) do
      if Abs(Direction[f])>1e-30 then
        begin
        LMin:=Min(LMin,(FMins[f]-Point[f])/Direction[f]);
        LMax:=Max(LMax,(FMAxs[f]-Point[f])/Direction[f]);
        end;
    end;
end;

function Evaluate(var x:TFloat):TFloat;

var f:Integer;
begin
  Inc(FLineCounts);
  if FMaxs<>nil then
    begin
    if x>LMax then
      begin
      x:=LMax;
      Result:=1e30;
      end
    else if x<LMin then
      begin
      x:=Lmin;
      Result:=1e30;
      end
    else
      begin
      for f:=0 to High(Point) do
        CurrentPoint[f]:=Point[f]+x*Direction[f];
      Result:=FGetValue(CurrentPoint);
      end
    end
    else
      begin
      for f:=0 to High(Point) do
        CurrentPoint[f]:=Point[f]+x*Direction[f];
      Result:=FGetValue(CurrentPoint);
      end

end;

function Bracket(var ax,bx:TFloat; out cx:TFloat;
                 var fa,fb:TFLoat; out fc:TFloat):Boolean;

var ulim,u,r,q,fu,dum:TFloat;

begin
     if (fb>fa) then
        begin
        dum:=ax;
        ax:=bx;
        bx:=dum;
        dum:=fa;
        fa:=fb;
        fb:=dum;
        end;
     cx:=bx+GOLD*(bx-ax);
     fc:=Evaluate(cx);
     while fb>fc do
           begin
           r:=(bx-ax)*(fb-fc);
           q:=(bx-cx)*(fb-fa);
           dum:=2*abs(q-r);
           if dum<TINY then dum:=TINY;
           u:=bx-((bx-cx)*q-(bx-ax)*r)/dum;
           ulim:=bx+GLIMIT*(cx-bx);
           if (bx-u)*(u-cx)>0 then
              begin
              fu:=Evaluate(u);
              if fu<fc then
                 begin
                 ax:=bx;
                 bx:=u;
                 fa:=fb;
                 fb:=fu;
                 break;
                 end
              else if fu>fb then
                   begin
                   cx:=u;
                   fc:=fu;
                   break;
                   end;
              u:=cx+GOLD*(cx-bx);
              fu:=Evaluate(u);
              end
           else if (cx-u)*(u-ulim)>0 then
                begin
                fu:=Evaluate(u);
                if fu<fc then
                   begin
                   bx:=cx;
                   cx:=u;
                   u:=cx+GOLD*(cx-bx);
                   fb:=fc;
                   fc:=fu;
                   fu:=Evaluate(u);
                   end;
                end
           else if (u-ulim)*(ulim-cx)>=0 then
                begin
                u:=ulim;
                fu:=Evaluate(u);
                end
           else
               begin
               u:=cx+GOLD*(cx-bx);
               fu:=Evaluate(u);
               end;
           ax:=bx;
           bx:=cx;
           cx:=u;
           fa:=fb;
           fb:=fc;
           fc:=fu;
           end;
    Result:=true;
end;

procedure Brent;

const ZEPS=10e-10;

var a,b,d,etemp,fu,fv,fw,fx,p,q,r,u,v,w,x,xm,e,tol1,tol2:TFloat;
    iter:integer;

begin
     e:=0;
     d:=0;
     if ax<cx then a:=ax
        else a:=cx;
     if ax>cx then b:=ax
        else b:=cx;
     x:=bx;
     w:=bx;
     v:=bx;
     fw:=fb;
     fv:=fb;
     fx:=fb;
     for iter:=1 to 100 do
         begin
         xm:=0.5*(a+b);
         tol1:=FPrecision*Abs(x)+ZEPS;
         tol2:=2.0*tol1;
         if abs(x-xm)<tol2-0.5*(b-a) then
            begin
            xmin:=x;
            break;
            end;
         if abs(e)>tol1 then
            begin
            r:=(x-w)*(fx-fv);
            q:=(x-v)*(fx-fw);
            p:=(x-v)*q-(x-w)*r;
            q:=2*(q-r);
            if q>0 then p:=-p;
            q:=abs(q);
            etemp:=e;
            e:=d;
            if (abs(p)>=abs(0.5*q*etemp)) or (p<=q*(a-x)) or
               (p>=q*(b-x)) then
               begin
               if x>=xm then e:=a-x
               else e:=b-x;
               d:=CR*e;
               end
            else begin
                 d:=p/q;
                 u:=x+d;
                 if (u-a<tol2) or (b-u<tol2) then
                    begin
                    if xm-x>0 then d:=tol1
                    else d:=-tol1;
                    end;
                 end;
            end
            else
                begin
                if x>=xm then e:=a-x else e:=b-x;
                d:=CR*e;
                end;
            if abs(d)>tol1 then
               u:=x+d
            else
                begin
                if d>0 then u:=x+tol1
                else u:=x-tol1;
                end;
         fu:=Evaluate(u);
         if fu<=fx then
            begin
            if u>=x then a:=x else b:=x;
            v:=w;
            w:=x;
            x:=u;
            fv:=fw;
            fw:=fx;
            fx:=fu;
            end
         else
             begin
             if u<x then a:=u else b:=u;
             if (fu<=fw) or (w=x) then
                begin
                v:=w;
                w:=u;
                fv:=fw;
                fw:=fu;
                end
             else if (fu<=fv) or (v=x) or (v=w) then
                  begin
                  v:=u;
                  fv:=fu;
                  end;
             end;
         end;
     xmin:=x;
     fa:=fx;
end;

procedure SetPoints(x:TFloat);

var f:Integer;

begin
  for f:=0 to High(Point) do
    Point[f]:=Point[f]+x*Direction[f];

end;

begin
  GetLineLimits;
  FLineCounts:=0;
  SetLength(CurrentPoint,Length(Point));
  ax:=0;
  fa:=Evaluate(ax);
  xmin:=ax;
  bx:=FPrecision;
  repeat
    fb:=Evaluate(bx);
    if (fb=fa) and (bx<FPrecision*1e4) then bx:=bx*10;
  until (fb<>fa) or (bx>=FPrecision*1e4);
  b:=Bracket(ax,bx,cx,fa,fb,fc);
  if b and (bx<>cx) then Brent;
  Result:=fa;
  SetPoints(xmin);
  CurrentPoint:=nil;
end;

function TCGMinimiser.Minimize(MaxIterations: Integer;
  MinVariation: TFloat; var Point, Deriv: TFloats): TFloat;

var f,j:integer;
    g,h:TFloats;
    dgg,dg,gam,oldr:TFloat;
begin
  SetLength(g,Length(Point));
  SetLength(h,Length(Point));
  Result:=FGetValue(Point);
  for f:=0 to High(Point) do
    begin
    g[f]:=-Deriv[f];
    Deriv[f]:=g[f];
    h[f]:=g[f];
    end;
  for f:=1 to MaxIterations do
    begin
    oldr:=Result;
    LineMinimization(Point,Deriv);
    Result:=FGetValueAndDerivative(Point,Deriv);
    if Assigned(FReport) then FReport(f,Result);
    dgg:=0;
    dg:=0;
    for j:=0 to High(Point) do
      begin
      dg:=dg+sqr(g[j]);
      dgg:=dgg+(Deriv[j]+g[j])*Deriv[j];
      end;
    if (dg=0) or (oldr-Result<=MinVariation) then Break
    else
      begin
      gam:=dgg/dg;
      for j:=0 to High(Point) do
        begin
        g[j]:=-Deriv[j];
        Deriv[j]:=g[j]+gam*h[j];
        h[j]:=Deriv[j];
        end;
      end;
    end;
  g:=nil;
  h:=nil;
end;


end.

