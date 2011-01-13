{*******************************************************************************
This file is part of the Open Chemera Library.
This work is public domain. Enjoy.
********************************************************************************
Author: Ludwig Krippahl
Date: 9.1.2011
Purpose:
  Open GL window manager, using GLUT.

  NOTE: this unit initializes the GLUT library in the initialization section.

Requirements:
  glut32.dll
  Open gl library. Either installed, or dll obtained from www.opengl.org
  or Nate Robins, http://www.xmission.com/~nate/glut.html

  Currently only supports using one TGLUTWindow object in the application

Revisions:
To do:
  Fix the requirement form only one TGLUTWindow object, by implementing
  functions to update the CurrentGlutWindow variable.
  Delete quadric
*******************************************************************************}
unit glutwindow;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, base3ddisplay, basetypes, glut, gl, glu, LCLProc;

type

    { TGlutWindow }

    TGlutWindow=class(TInterfacedObject,IBase3DDisplay)
    protected
      FTop,FLeft,FHeight,FWidth:Integer;
      FMainWindow:Integer; //OpenGL Id of main window
      FQuadric:PGLUquadric;

      //Light attributes
      FLAmbient,FLDiffuse,FLPosition,FLSpecular: TRGBAColor;

      //for testing only
      FAngle:Single;

    public
      constructor Create(Top,Left,Height,Width:Integer);

      //IBase3DInterface functions
      function AddSphere(Coords:TOCCoords;Rad:TOCFloat;MatId:Integer;Quality:Integer=-1):Integer;
      function AddCylinder(C1,C2:TOCCoords;R1,R2:TOCFloat;MatId:Integer;Quality:Integer=-1):Integer;
      function AddLine(C1,C2:TOCCoords;R1,R2:TOCFloat;MatId:Integer):Integer;
      function AddMaterial(Material:T3DMaterial):Integer;
      procedure ClearObjects;
      procedure ClearMaterials;
      procedure Refresh;
      procedure Resize(Width,Height:Integer);
    end;

implementation

var
   //Keeps track of current manager object to use for OGL callback functions
   //(which cannot be of object)
   //Should never need to be changed or used UNLESS using more than one
   //TGLUTWindowObject
   CurrentGlutWindow:TGlutWindow;

procedure InitializeGLUT;
//initializes GLUT library
//adapted from Lazarus OpenGL tutorials

var
   CmdCount, I: Integer;
   Cmd:array of pchar;
begin
   CmdCount := ParamCount + 1;
   SetLength(Cmd, CmdCount);
   for I := 0 to CmdCount - 1 do
     Cmd[I] := PChar(ParamStr(I));
   glutInit(@CmdCount, @Cmd);
 end;

procedure DisplayWindow; cdecl;
begin
  CurrentGlutWindow.Refresh;
end;

procedure ReshapeWindow(Width, Height: Integer); cdecl;
begin
  CurrentGLUTWindow.Resize(Width,Height);
end;

{ TGlutWindow }

constructor TGlutWindow.Create(Top, Left, Height, Width: Integer);
begin
  inherited Create;
  CurrentGLUTWindow:=Self;
  FLAmbient:=RGBAColor(0.2, 0.2, 0.2, 1 );
  FLDiffuse:=RGBAColor( 0.7, 0.7, 0.7, 1 );
  FLSpecular:=RGBAColor( 1, 1, 1, 1 );
  FLPosition:=RGBAColor( 0.0, 100.0, 0, 20 );

  FTop:=Top;
  FLeft:=Left;
  FHeight:=Height;
  FWidth:=Width;
  glutInitDisplayMode(GLUT_RGBA or GLUT_DOUBLE or GLUT_DEPTH);

  glClearColor(0.0, 1.0, 0.0, 1);
  glClearDepth(1.0);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
  glHint(GL_PERSPECTIVE_CORRECTION_HINT,GL_NICEST);

  //Lighting

  //Create quadric for cylinders and spheres
  FQuadric:=gluNewQuadric();
  gluQuadricOrientation(FQuadric,GLU_OUTSIDE);
  gluQuadricNormals(FQuadric, GLU_SMOOTH);


  glutInitWindowPosition(Left, Top);
  glutInitWindowSize(Width, Height);
  FMainWindow:=glutCreateWindow('Display');
  glutDisplayFunc(@DisplayWindow);
  glutReshapeFunc(@ReshapeWindow);

  //Testing Onlyu
  FAngle:=0;
  glutIdleFunc(@DisplayWindow);

  glutMainLoop;
end;

function TGlutWindow.AddSphere(Coords: TOCCoords; Rad: TOCFloat;
  MatId: Integer; Quality: Integer): Integer;
begin

end;

function TGlutWindow.AddCylinder(C1, C2: TOCCoords; R1, R2: TOCFloat;
  MatId: Integer; Quality: Integer): Integer;
begin

end;

function TGlutWindow.AddLine(C1, C2: TOCCoords; R1, R2: TOCFloat; MatId: Integer
  ): Integer;
begin

end;

function TGlutWindow.AddMaterial(Material: T3DMaterial): Integer;
begin

end;

procedure TGlutWindow.ClearObjects;
begin

end;

procedure TGlutWindow.ClearMaterials;
begin

end;

procedure TGlutWindow.Refresh;
  //This is the drawing function, called by DisplayWindow, which is called by
  //the Open GL

var
  f:Integer;
  col:TRGBAColor;

begin
  glClearColor(0.0, 0, 0.0, 1);
  glClear(GL_COLOR_BUFFER_BIT or GL_DEPTH_BUFFER_BIT);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
  glPushMatrix;
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity;

  glShadeModel(GL_SMOOTH);
  glLightfv(GL_LIGHT0, GL_SPECULAR,FLSpecular);
  glLightfv(GL_LIGHT0, GL_POSITION,FLPosition);
  glLightfv(GL_LIGHT0, GL_AMBIENT,FLAmbient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE,FLDiffuse);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHTING);
  //glEnable(GL_COLOR_MATERIAL);


  //glEnable(GL_NORMALIZE);
  col:=RGBAColor(0.8,1,0.2,0);
  //glColor4fv(col);
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, col);
  glMaterialfv(GL_FRONT, GL_SPECULAR, col);
  glMateriali(GL_FRONT, GL_SHININESS, 60);

  glTranslatef(0,0,-20);
  glRotatef(FAngle,1,1,0);

  //gluCylinder(FQuadric,3,3,6,31,31);
  glutSolidTeapot(2);
  {glColor3f(0.0, 0.0, 1.0); // blue reflective properties
  for f:=1 to 100 do
      begin
      glLoadIdentity;
      glRotatef(Random(360),Random(2),Random(2),0);
      glTranslatef(Random(20)-10,Random(20)-10,Random(20)-50);
      gluSphere(FQuadric,1.3,32,32);
  //    gluCylinder(CurrentGLUTWindow.FQuadric,1.0,1.0,10.0,32,32);

      end;}
  DebugLn('Refresh');
  glPopMatrix;
  glutSwapBuffers;
  FAngle:=FAngle+0.01;
end;

procedure TGlutWindow.Resize(Width, Height: Integer);
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

initialization
  InitializeGLUT;
end.



