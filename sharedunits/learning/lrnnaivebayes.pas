{*******************************************************************************
This source file is public domain. It is provided as is, with no explicit
or implied warranties regarding anything whatsoever.
********************************************************************************
Author: Ludwig Krippahl
Date: 30.12.2013
Purpose:
  Naive Bayes classifiers
Requirements:
Revisions:
To do:
*******************************************************************************}

unit lrnnaivebayes;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, lrnhistogram, basetypes;



function NBClassify(ClassHistograms:THistMatrix;
                      //a matrix of Num Classes * Num Features of TFloatHistogram
                      //this is the trained NB classifier,
                      //assumed to have normalized log probabilities
                    Data:TMatrix;
                      //a matrix of Num Features * Num Points of TFLoat
                      //These are the features for classification, assumed all to be TFloat
                    out Scores:TMatrix
                      //Num Classes * Num Points matrix of the scores obtained for each class
                      //(scores are the non-normalized log probabilities)
                    ):TIntegers;
                      //Num Points array with the class maximizing the score

implementation


function NBClassify(ClassHistograms:THistMatrix; Data:TMatrix; out Scores:TMatrix):TIntegers;

var
  tot,f,classix,cl,feat:Integer;
  score:TFloat;
  tmpvals:TFloats;
  priorlogp:TFloats;

begin
  if Data<>nil then
    begin
    SetLength(Result,Length(Data[0]));
    SetLength(Scores,Length(ClassHistograms));
    SetLength(priorlogp,Length(ClassHistograms));

    //compute priors from number of data points in histograms
    //assumes all histograms of same class have same number of points
    tot:=0;
    for cl:=0 to High(priorlogp) do
      begin
      tot:=tot+ClassHistograms[cl,0].Count;
      priorlogp[cl]:=ClassHistograms[cl,0].Count;
      end;
    for cl:=0 to High(priorlogp) do
      priorlogp[cl]:=Ln(priorlogp[cl]/tot);
    for cl:=0 to High(priorlogp) do
      begin
      Scores[cl]:=FilledFloats(Length(Data[0]),priorlogp[cl]);
      for feat:=0 to High(Data) do
        begin
        tmpvals:=TFloatHistogram(ClassHistograms[cl,feat]).BinVals(Data[feat]);
        Scores[cl]:=Sum(Scores[cl],tmpvals);
        end;
      end;
    for f:=0 to High(Scores[0]) do
      begin
      classix:=0;
      score:=Scores[0,f];
      for cl:=1 to High(Scores) do
        if score<Scores[cl,f] then
          begin
          classix:=cl;
          score:=Scores[cl,f];
          end;
      Result[f]:=classix;
      end;
    end
  else
    begin
    Result:=nil;
    Scores:=nil;
    end;
end;

end.

