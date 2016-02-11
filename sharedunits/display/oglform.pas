{*******************************************************************************
This source file is public domain. It is provided as is, with no explicit
or implied warranties regarding anything whatsoever.
********************************************************************************
Author: Ludwig Krippahl
Date: 9.1.2011
Purpose:
  OpenGL form based on the IBase3DDisplay interface, to instantiate a display
  using OpenGL. This form should be accessed via the interface methods.
Requirements:
Revisions:
To do: Comments
*******************************************************************************}
unit oglform;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, FileUtil, Forms, Controls, Graphics, Dialogs,
  OpenGLContext,gl,glu,glext,base3ddisplay,basetypes,geomutils,LCLProc, types;

type
  //color format compatible with glu
  TGluRGBA=array[0..3] of GLFloat;
  TGluMatrix=array[0..15] of GLFLoat;

  { TOpenGLForm }

  TOpenGLForm = class(TForm,IBase3DDisplay)
    OGLPb: TOpenGLControl;
    procedure FormDropFiles(Sender: TObject; const FileNames: array of String);
    procedure FormResize(Sender: TObject);
    procedure OGLPbMouseDown(Sender: TObject; Button: TMouseButton;
      Shift: TShiftState; X, Y: Integer);
    procedure OGLPbMouseMove(Sender: TObject; Shift: TShiftState; X, Y: Integer
      );
    procedure OGLPbMouseWheel(Sender: TObject; Shift: TShiftState;
      WheelDelta: Integer; MousePos: TPoint; var Handled: Boolean);
    procedure OGLPbMouseWheelDown(Sender: TObject; Shift: TShiftState;
      MousePos: TPoint; var Handled: Boolean);
    procedure OGLPbPaint(Sender: TObject);
    procedure OGLPbResize(Sender: TObject);

  private
    { private declarations }
    procedure ResetMouseRotation;

  protected
      FInitialized:Boolean;
      FDistanceToCenter:TFLoat;
      FBaseSphere,FBaseCube:TQuadMesh;
      FDisplayLists:array of GLUInt;
      FObjectLists: array of T3DObjectList;
      FDisplayColors: array of TGluRGBA; //colors for displaylists
      FDisplayShine: array of TGLFloat; //colors for displaylists
      FBackground:TGluRGBA;

      FCameraX,FCameraY,FCameraZ:TFloat;

      //mouse tracking
      FCX,FCY:Integer;                   //current mouse position
      FOX,FOY:Integer;                   //old mouse position
      FModelMat:TGluMatrix;              //rotation matrix
      FModelViewMat:TGluMatrix;         //Full model view matrix
      //Light attributes
      FLAmbient,FLDiffuse,FLPosition,FLSpecular: TGluRGBA;

      //for testing only
      FAngle:Single;
      procedure ClearDisplayLists;
      procedure CompileList(ObjectList: T3DObjectList; ListName: GLUInt);
      procedure SetMaterial(Material:T3DMaterial);
      procedure InitOGL;
    public
      property ModelViewMat:TGluMatrix read FModelViewMat;
      property ZDist:TFloat read FDistanceToCenter;
      constructor Create(TheOwner: TComponent); override;
      //IBase3DInterface functions
      procedure AddObjectList(ObjectList:T3DObjectList);
      procedure ClearObjectLists;
      procedure Compile;
      procedure ResizeForm(NewWidth,NewHeight:Longint);
      procedure Refresh;
      procedure SetQuality(Quality:Integer);
      function GetImage:TBitMap;
      procedure RotateMatrix(Horizontal,Vertical:TFLoat);
      procedure CameraPos(out X,Y,Z:TFLoat);
  end;


implementation

{$R *.lfm}

function RGBAtoGlu(C:TRGBAColor):TGluRGBA;

//to convert from base3ddisplay color format
//No typecasting to avoid conflict with different float formats

begin
 Result[0]:=C[0];
 Result[1]:=C[1];
 Result[2]:=C[2];
 Result[3]:=C[3];
end;

function GluColor(R,G,B,A:TGLFloat):TGluRGBA;

begin
 Result[0]:=R;
 Result[1]:=G;
 Result[2]:=B;
 Result[3]:=A;
end;


{ TOpenGLForm }

constructor TOpenGLForm.Create(TheOwner: TComponent);
begin
  inherited Create(TheOwner);
  FDistanceToCenter:=100;
  OGLPb.MakeCurrent;
  SetQuality(3);
  FBackground:=ColorWhite;    // sets background color

end;

procedure TOpenGLForm.OGLPbPaint(Sender: TObject);
begin
  if not FInitialized then
    InitOGL;
  Refresh;
end;

procedure TOpenGLForm.OGLPbResize(Sender: TObject);
begin
  if (FInitialized) and OGLPb.MakeCurrent then
    begin
    glViewport (0, 0, OGLPb.Width, OGLPb.Height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(30.0, OGLPb.Width/OGLPb.Height, 1, 600.0);
    glMatrixMode(GL_MODELVIEW);
    end;
end;

procedure TOpenGLForm.RotateMatrix(Horizontal, Vertical: TFLoat);
begin
  //Set rotation and reset X Y
  glLoadIdentity;
  glRotatef(Horizontal,0,1,0);
  glRotatef(Vertical,1,0,0);
  glMultMatrixf(FModelMat);
  glGetFloatv(GL_MODELVIEW_MATRIX,FModelMat);
end;

procedure TOpenGLForm.CameraPos(out X, Y, Z: TFLoat);

begin
  X:=FCameraX;
  Y:=FCameraY;
  Z:=FCameraZ;
end;

procedure TOpenGLForm.ResetMouseRotation;
begin
  //Set rotation and reset X Y
  RotateMatrix(0.5*(FCX-FOX),0.5*(FCY-FOY));
  FOX:=FCX;
  FOY:=FCY;
end;

procedure TOpenGLForm.FormResize(Sender: TObject);
begin
  OGlPb.Top:=3;
  OGlPb.Left:=3;
  OGlPb.Width:=Width-6;
  OGlPb.Height:=Height-6;

end;

procedure TOpenGLForm.FormDropFiles(Sender: TObject;
  const FileNames: array of String);
begin

end;

procedure TOpenGLForm.OGLPbMouseDown(Sender: TObject; Button: TMouseButton;
  Shift: TShiftState; X, Y: Integer);
begin
  //set old and current mouse position to this position
  FOX:=X;
  FOY:=Y;
  FCX:=X;
  FCY:=Y;
end;

procedure TOpenGLForm.OGLPbMouseMove(Sender: TObject; Shift: TShiftState; X,
  Y: Integer);
begin
  if Shift=[ssLeft] then
    begin
    FCX:=X;
    FCY:=Y;
    OGLPb.Invalidate;
    end;
end;

procedure TOpenGLForm.OGLPbMouseWheel(Sender: TObject; Shift: TShiftState;
  WheelDelta: Integer; MousePos: TPoint; var Handled: Boolean);
begin
  FDistanceToCenter:=FDistanceToCenter-WheelDelta/20;
  OGLPb.Invalidate;
end;

procedure TOpenGLForm.OGLPbMouseWheelDown(Sender: TObject; Shift: TShiftState;
  MousePos: TPoint; var Handled: Boolean);
begin

end;

procedure TOpenGLForm.ClearDisplayLists;

var f:Integer;

begin
  for f:=0 to High(FDisplayLists) do
    glDeleteLists(FDisplayLists[f],1);
  FDisplayLists:=nil;
end;

procedure TOpenGLForm.CompileList(ObjectList: T3DObjectList; ListName: GLUInt);

procedure DoSphere(Rad:TFloat;sphC:TCoord);

var
  f,g:Integer;
  c0,c:TCoord;

begin
  with FBaseSphere do
    for f:=0 to High(Faces) do
      for g:=0 to 3 do
        begin
        c0:=Points[Faces[f,g]];
        glNormal3f(c0[0],c0[1],c0[2]);
        c:=Multiply(c0,Rad);
        c:=Add(c,sphC);
        //TO DO:if UseTextures then glTexCoord2f(,);
        glVertex3f(c[0],c[1],c[2]);
        end;
end;

procedure DoCuboid(TL,BR:TCoord);

var
  f,g:Integer;
  c:TCoord;

begin
  with FBaseCube do
    begin
    Points[0]:=BR;
    Points[1]:=Coord(BR[0],BR[1],TL[2]);
    Points[2]:=Coord(BR[0],TL[1],BR[2]);
    Points[3]:=Coord(BR[0],TL[1],TL[2]);
    Points[4]:=Coord(TL[0],BR[1],BR[2]);
    Points[5]:=Coord(TL[0],BR[1],TL[2]);
    Points[6]:=Coord(TL[0],TL[1],BR[2]);
    Points[7]:=TL;
    for f:=0 to High(Faces) do
      for g:=0 to 3 do
        begin
        c:=Points[Faces[f,g]];
        glNormal3f(Normals[f,0],Normals[f,1],Normals[f,2]);
        //TO DO:if UseTextures then glTexCoord2f(,);
        glVertex3f(c[0],c[1],c[2]);
        end;
    end;
end;

var
  f:Integer;

begin
  glNewList(ListName, GL_COMPILE);
  SetMaterial(ObjectList.Material);
  glBegin(GL_QUADS);

  glEnable (GL_BLEND);
  //glBlendFunc (GL_SRC_ALPHA, GL_DST_ALPHA);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  //glBlendFunc(GL_ONE_MINUS_DST_ALPHA,GL_DST_ALPHA);

  for f:=0 to High(ObjectList.Objects) do
    with ObjectList.Objects[f] do
    case ObjectType of
      otSphere:DoSphere(Rad,sphC);
      {otCilinder:(cylC1,cylC2:TCoord;Rad1,Rad2:TFLoat);}
      otCuboid:DoCuboid(cubTopLeft,cubBotRight);
    end;
  glEnd;
  // do {otLine:(linC1,linC2:TCoord); here, and others not quads
  glEndList;
end;

procedure TOpenGLForm.SetMaterial(Material: T3DMaterial);

var
  colix:Integer;

begin
  with Material do
    begin

    //TO DO: Textures, need handles
    //if Texture<>'' then glEnable(GL_TEXTURE_2D)
    //else glDisable(GL_TEXTURE_2D);

    colix:=Length(FDisplayColors);
    SetLength(FDisplayColors,colix+3);
    //store and set each color, for sending pointers
    //TO DO: check if better way to do this...
    FDisplayColors[colix]:=RGBAToGlu(Material.Ambient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, FDisplayColors[colix]);
    FDisplayColors[colix+1]:=RGBAToGlu(Material.Diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, FDisplayColors[colix+1]);
    FDisplayColors[colix+2]:=RGBAToGlu(Material.Specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, FDisplayColors[colix+2]);
    SetLength(FDisplayShine,Length(FDisplayShine)+1);
    FDisplayShine[High(FDisplayShine)]:=Material.Shine;
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, @FDisplayShine[High(FDisplayShine)]);
    end;
end;

procedure TOpenGLForm.InitOGL;

procedure InitLights;

begin
  FLAmbient:=GluColor(1, 1, 1, 1 );
  FLDiffuse:=GluColor( 1, 1, 1, 1 );
  FLSpecular:=GluColor( 1, 1, 1, 1 );
  FLPosition:=GluColor( 1.0, 1.0, 1, 0 );
end;

begin
  if FInitialized then exit;
  FInitialized:=true;
  InitLights;
  {setting lighting conditions}
  glEnable(GL_LIGHT0);
  glLightfv(GL_LIGHT0, GL_AMBIENT, FLAmbient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, FLDiffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, FLSpecular);
  glLightfv(GL_LIGHT0, GL_POSITION, FLPosition);

  glMaterialfv(GL_FRONT, GL_DIFFUSE, ColorLtGray);
  glMaterialfv(GL_FRONT, GL_AMBIENT, ColorLtGray);

  glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);

  glCullFace(GL_BACK);
  glEnable(GL_CULL_FACE);

  //Fog
  //glEnable (GL_FOG);
  glFogi (GL_FOG_MODE,GL_LINEAR);
  glFogfv (GL_FOG_COLOR, FBackground);
  glFogf (GL_FOG_DENSITY, 0.01);
  glFogf(GL_FOG_START, 20);
  glFogf(GL_FOG_END, 100);
  glHint (GL_FOG_HINT, GL_NICEST);

  //textures
  //TO DO:  FTextures.BindTextures;

  glClearColor(FBackground[0],
               FBackground[1],
               FBackground[2],
               FBackground[3]);

  glClearDepth(1.0);
  glDepthFunc(GL_LEQUAL);           // the type of depth test to do
  glEnable(GL_DEPTH_TEST);          // enables depth testing
  glShadeModel(GL_SMOOTH);          // enables smooth color shading
  glEnable(GL_LIGHTING);


  glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
  glHint(GL_POLYGON_SMOOTH_HINT,GL_NICEST);
  glHint(GL_PERSPECTIVE_CORRECTION_HINT,GL_NICEST);

  glMatrixMode (GL_PROJECTION);    { prepare for and then }
  glLoadIdentity ();               { define the projection }
  glFrustum (-2.0, 2.0, -2.0, 2.0, 5, 200); { transformation }
  glMatrixMode (GL_MODELVIEW);  { back to modelview matrix }
  glViewport (0, 0, OGLPb.Width, OGLPb.Height);
  OGLPbResize(nil);

  //set the matrix for rotation
  glLoadIdentity;
  glGetFloatv(GL_MODELVIEW_MATRIX,FModelMat);



end;

procedure TOpenGLForm.AddObjectList(ObjectList: T3DObjectList);
begin
  SetLength(FObjectLists,Length(FObjectLists)+1);
  FObjectLists[High(FObjectLists)]:=ObjectList;
  //copy the objects array so that caller can reuse the original
  FObjectLists[High(FObjectLists)].Objects:=
    Copy(ObjectList.Objects,0,Length(Objectlist.Objects));
end;

procedure TOpenGLForm.ClearObjectLists;
begin
  FObjectLists:=nil;
  FDisplayColors:=nil;
  FDisplayShine:=nil;
end;

procedure TOpenGLForm.Compile;
// create display lists for all objects
//TO DO: should check which object lists changed, delete those display lists
// and reuse them (that's what FDisplayLists is for...)

var
  f:Integer;

begin
  ClearDisplayLists;
  SetLength(FDisplayLists,Length(FObjectLists));
  for f:=0 to High(FObjectLists) do
    begin
    FDisplayLists[f]:=f+1;
    CompileList(FObjectLists[f],f+1);
    end;
end;


procedure TOpenGLForm.Refresh;
  //This is the drawing function, called by DisplayWindow, which is called by
  //the Open GL

procedure RecomputeCamera;
var
  viewport:array[0..3] of Integer;
  model,proj:array [0..15] of GLDouble;
  glX,glY,glZ:PGLDouble;

begin
  New(glX);
  New(glY);
  New(glZ);
  glGetIntegerv(GL_VIEWPORT, viewport);
  glGetDoublev(GL_MODELVIEW_MATRIX,model);
  glGetFloatv(GL_MODELVIEW_MATRIX,FModelViewMat);
  glGetDoublev(GL_PROJECTION_MATRIX,proj);
  gluUnProject((viewport[2]-viewport[0])/2 , (viewport[3]-viewport[1])/2,0.0,
        model, proj, viewport, glX,glY,glZ);
  FCameraX:=glX^;
  FCameraY:=glY^;
  FCameraZ:=glZ^;

  Dispose(glX);
  Dispose(glY);
  Dispose(glZ);
end;

var
  f:Integer;

begin
  glClearColor(FBackground[0],
               FBackground[1],
               FBackground[2],
               FBackground[3]);

  glClear(GL_COLOR_BUFFER_BIT or GL_DEPTH_BUFFER_BIT);
  glEnable(GL_DEPTH_TEST);


  glPushMatrix;



  //set objects in place and rotation
  ResetMouseRotation;
  glLoadIdentity;
  glTranslatef(0,0,-FDistanceToCenter);
  glMultMatrixf(FModelMat);

  for f:=0 to High(FDisplayLists) do
    glCallList(FDisplayLists[f]);

  RecomputeCamera;

  glPopMatrix;
  OGLPb.SwapBuffers;
  glFinish()

end;

procedure TOpenGLForm.SetQuality(Quality: Integer);

begin
  FBaseSphere:=QuadSphere(Quality);
  FBaseCube:=QuadCube;
end;

function TOpenGLForm.GetImage: TBitMap;


var
  buffer,p:PByte;
  r,g,b,a:Byte;
  y,x:Integer;

begin
  Result:=TBitMap.Create;
  Result.PixelFormat := pf32bit;
  Result.Width:=OGLPb.Width;
  Result.Height:=OGLPb.Height;
  GetMem( buffer, OGLPb.Width*OGLPb.Height*SizeOf(DWord));
  Refresh;
  glReadBuffer(GL_Front);
  glPixelStorei(GL_PACK_ALIGNMENT, 4);
  glReadPixels( 0, 0, OGLPb.Width,OGLPb.Height, GL_RGBA, GL_UNSIGNED_BYTE, buffer);
  p:=buffer;
  Result.BeginUpdate(True);
  for y:=Result.Height-1 downto 0 do
    for x:=0 to Result.Width-1 do
      begin
      r:=p^;Inc(p);
      g:=p^;Inc(p);
      b:=p^;Inc(p);
      Result.Canvas.Pixels[x,y]:=r+g shl 8 +b shl 16+$02000000;
      Inc(p);
      end;
  Result.EndUpdate(False);
  FreeMem(buffer);
  writeln('ok');
end;


procedure TOpenGLForm.ResizeForm(NewWidth, NewHeight: Longint);
begin
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(     15.0,   // Field-of-view angle
                   Width/Height,   // Aspect ratio of view volume
                         0.1,   // Distance to near clipping plane
                    100.0 ); // Distance to far clipping plane
  glViewport(0, 0, Width, Height);
end;

end.
