unit oglform;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, FileUtil, Forms, Controls, Graphics, Dialogs,
  OpenGLContext,gl,glu,base3ddisplay,basetypes,LCLProc;

type
  //color format compatible with glu
  TGluRGBA=array[0..3] of GLFloat;

  { TOpenGLForm }

  TOpenGLForm = class(TForm,IBase3DDisplay)
    OGLPb: TOpenGLControl;
  private
    { private declarations }
  protected
      FTop,FLeft,FHeight,FWidth:Integer;
      FForm:TOpenGLForm;
      FBaseSphere,FBaseCube:TQuadMesh;
      FDisplayLists:array of GLUInt;
      FObjectLists: array of T3DObjectList;
      FDisplayColors: array of TGluRGBA; //colors for displaylists
      //Light attributes
      FLAmbient,FLDiffuse,FLPosition,FLSpecular: TGluRGBA;

      //for testing only
      FAngle:Single;
      procedure ClearDisplayLists;
      function CompileList(ObjectList:T3DObjectList):GLUInt;
      procedure SetMaterial(Material:T3DMaterial);
    public
      constructor Create(TheOwner: TComponent); override;
      //IBase3DInterface functions
      procedure AddObjectList(ObjectList:T3DObjectList);
      procedure ClearObjectLists;
      procedure Compile;
      procedure Resize(NewWidth,NewHeight:Integer);
      procedure Refresh;
      procedure SetQuality(Quality:Integer);

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
  OGLPb.MakeCurrent;
    SetQuality(3);

  FLAmbient:=GluColor(0.2, 0.2, 0.2, 1 );
  FLDiffuse:=GluColor( 0.7, 0.7, 0.7, 1 );
  FLSpecular:=GluColor( 1, 1, 1, 1 );
  FLPosition:=GluColor( 0.0, 100.0, 0, 20 );

  FTop:=Top;
  FLeft:=Left;
  FHeight:=Height;
  FWidth:=Width;

  glClearColor(0.0, 1.0, 0.0, 1);
  glClearDepth(1.0);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
  glHint(GL_PERSPECTIVE_CORRECTION_HINT,GL_NICEST);


  //Default lights and Rendering options

  glShadeModel(GL_SMOOTH);
  glLightfv(GL_LIGHT0, GL_SPECULAR,FLSpecular);
  glLightfv(GL_LIGHT0, GL_POSITION,FLPosition);
  glLightfv(GL_LIGHT0, GL_AMBIENT,FLAmbient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE,FLDiffuse);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHTING);

  glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);

  glCullFace(GL_BACK);
  glEnable(GL_CULL_FACE);
end;


procedure TOpenGLForm.ClearDisplayLists;
begin
  FDisplayLists:=nil;
end;

function TOpenGLForm.CompileList(ObjectList: T3DObjectList): GLUInt;

procedure DoSphere(Rad:TOCFloat;sphC:TOCCoord);

var
  f,g:Integer;
  c0,c:TOCCoord;

begin
  with FBaseSphere do
    for f:=0 to High(Faces) do
      for g:=0 to 3 do
        begin
        c0:=Points[Faces[f,g]];
        glNormal3f(c0[0],c0[1],c0[2]);
        c:=ScaleVector(c0,Rad);
        c:=AddVectors(c,sphC);
        //TODO:if UseTextures then glTexCoord2f(,);
        glVertex3f(c[0],c[1],c[2]);
        end;
end;

var f:Integer;

begin
  //setup
  glNewList(Result, GL_COMPILE);
  SetMaterial(ObjectList.Material);
  glBegin(GL_QUADS);
  for f:=0 to High(ObjectList.Objects) do
    with ObjectList.Objects[f] do
    case ObjecType of
      otSphere:DoSphere(Rad,sphC);
      {otCilinder:(cylC1,cylC2:TOCCoord;Rad1,Rad2:TOCFLoat);
      otCuboid:(cubTopLeft,cubBotRight:TOCCoord);}
    end;
  glEnd;
  // do {otLine:(linC1,linC2:TOCCoord); here, and others not quads
  glEndList;
end;

procedure TOpenGLForm.SetMaterial(Material: T3DMaterial);

var colix:Integer;

begin
  with Material do
    begin

    //TODO: Textures, need handles
    //if Texture<>'' then glEnable(GL_TEXTURE_2D)
    //else glDisable(GL_TEXTURE_2D);

    colix:=Length(FDisplayColors);
    SetLength(FDisplayColors,colix+3);
    //store and set each color, for sending pointers
    //TODO: check if better way to do this...
    FDisplayColors[colix]:=RGBAToGlu(Material.Ambient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, @FDisplayColors[colix]);
    FDisplayColors[colix+1]:=RGBAToGlu(Material.Diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, @FDisplayColors[colix+1]);
    FDisplayColors[colix+2]:=RGBAToGlu(Material.Specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, @FDisplayColors[colix+2]);
    end;
end;

procedure TOpenGLForm.AddObjectList(ObjectList: T3DObjectList);
begin
  SetLength(FObjectLists,Length(FObjectLists)+1);
  FObjectLists[High(FObjectLists)]:=ObjectList;
end;

procedure TOpenGLForm.ClearObjectLists;
begin
  FObjectLists:=nil;
end;

procedure TOpenGLForm.Compile;
begin

end;


procedure TOpenGLForm.Refresh;
  //This is the drawing function, called by DisplayWindow, which is called by
  //the Open GL

var
  f:Integer;
  col:TGluRGBA;

begin
  glClearColor(0.0, 0, 0.0, 1);
  glClear(GL_COLOR_BUFFER_BIT or GL_DEPTH_BUFFER_BIT);
  glEnable(GL_DEPTH_TEST);
  glPushMatrix;
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity;

  //Set rotation here
  //TODO:translation needs variable parameter
  glTranslatef(0,0,-20);

  for f:=0 to High(FDisplayLists) do
    glCallList(FDisplayLists[f]);

  glPopMatrix;
  OGLPb.SwapBuffers;

end;

procedure TOpenGLForm.SetQuality(Quality: Integer);

begin
  FBaseSphere:=QuadSphere(Quality);
  FBaseCube:=QuadCube;
end;

procedure TOpenGLForm.Resize(NewWidth, NewHeight: Integer);
begin
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(     30.0,   // Field-of-view angle
                   Width/Height,   // Aspect ratio of view volume
                         0.1,   // Distance to near clipping plane
                    100.0 ); // Distance to far clipping plane
  glViewport(0, 0, Width, Height);
  DebugLn('Resize');

end;

end.
