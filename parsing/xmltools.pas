unit xmltools;

{
  Help parse XML files

}

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, DOM;


 // These functions return empty strings if the child or attribute is not found
function GetNodeAttribute(node: TDOMNode; atname: string): string;
function GetChildText(node:TDOMNOde; childname:string):string;

implementation

function GetNodeAttribute(node: TDOMNode; atname: string): string;

var pn:TDOMNode;

begin
  Result:='';
  pn:=node.Attributes.GetNamedItem(atname);
  if pn<>nil then
    Result:=pn.FirstChild.NodeValue;
end;

function GetChildText(node: TDOMNOde; childname: string): string;

var pn:TDOMNode;

begin
  Result:='';
  pn:=node.FindNode(childname);
  if pn<>nil then
    Result:=pn.FirstChild.NodeValue;
end;

end.

