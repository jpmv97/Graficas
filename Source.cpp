#include <windows.h>
#include <stdio.h>
#include <math.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/gl.h>
#include "glut.h"

#include "CG_Tec.h"

using namespace std;

// Global variables
// ------------------------------------
// Light Components
static float globAmb = 0.3;
static float globDif = 0.4;
static float globSpe = 0.7;

// Material Components
static float globMatAmb = 0.5;
static float globMatDif = 0.5;
static float globMatSpe = 0.5;
static float globExp = 4.0;

static bool extrudeORsor = false;

static float xSpeed = 0.0, ySpeed = 0.0, xAngle = 0.0, yAngle = 0.0;

static int m = 20; // number of divisions in the solid of revolution
static int n = 24; // number of points in the profile NOTE!!! ADDED 4 MORE POINTS!!

#define PI 3.1416

Point3 profile[24];
Point3 normals[24];
Point3 totalPoints[20][24];
Point3 totalNormals[20][24];
Affine4 transform;

RGBpixmap myTex;
GLuint theTexture[1];


//Proyectos
Mesh firstobj;
void createSORPoints()
{
	float angle = -2.0*PI / (float)m;

	transform.setIdentityMatrix();

	// calculate transform
	if (extrudeORsor)
	{
		// Need Translation
		Point3 translate(0.0, 0.0, 0.3);
		transform.ApplyTranslate(translate);
	}
	else
	{
		// Need Rotation Y
		transform.ApplyRotateY(angle);
	}

	// calculate all points
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i == 0)
			{
				totalPoints[i][j].set(profile[j]);
				totalNormals[i][j].set(normals[j]);
			}
			else
			{
				totalPoints[i][j].set(transform.MultiplyPoint(totalPoints[i - 1][j]));
				totalNormals[i][j].set(transform.MultiplyPoint(totalNormals[i - 1][j]));
			}
		}
	}
}

void updateLights(void)
{
	// Update LIGHT properties
	GLfloat myAmb[] = { globAmb, globAmb, globAmb, 1.0 };
	GLfloat myDif[] = { globDif, globDif, globDif, 1.0 };
	GLfloat mySpe[] = { globSpe, globSpe, globSpe, 1.0 };
	glLightfv(GL_LIGHT0, GL_AMBIENT, myAmb);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, myDif);
	glLightfv(GL_LIGHT0, GL_SPECULAR, mySpe);
}

void init(void)
{
	glClearColor(0.0, 0.0, 0.0, 0.0);

	//glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	GLfloat mat_spec[] = { 0.393548, 0.271906, 0.166721, 1.0 };
	GLfloat mat_shiny[] = { 0.2 };
	GLfloat mat_amb[] = { 0.2125, 0.1275, 0.054, 1.0 };
	GLfloat mat_dif[] = { 0.714, 0.4284, 0.18144, 1.0 };

	//SET material properties
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_spec);
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shiny);
	glMaterialfv(GL_FRONT, GL_AMBIENT, mat_amb);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_dif);

	//DEFINE light properties (for lights in scene)
	GLfloat while_light[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_post[] = { 1.0, 1.0, 1.0, 0.0 };
	//SET light properties
	glLightfv(GL_LIGHT0, GL_POSITION, light_post);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, while_light);
	glLightfv(GL_LIGHT0, GL_SPECULAR, while_light);

	//Activate (enable) lights
	glEnable(GL_LIGHT0);

	//Enable lighting
	glEnable(GL_LIGHTING);

	//Enable depth testing (for hidden surface removal)
	glEnable(GL_DEPTH_TEST);

	// Update LIGHT properties
	//updateLights();

	// Initialize the profile
	// NOTE!!! NEW!!! Make a Square Base
	
	char objname[100] = "arbol.obj";
	firstobj.readOBJ(objname);
	//char texturaName[10] = "itesm.bmp";
	//InitializeTexture(&myTex, theTexture, texturaName);

}

void setMaterials(void)
{
	GLfloat matAmbi[] = { globMatAmb, globMatAmb, globMatAmb, 1.0 };
	GLfloat matDiff[] = { globMatDif, globMatDif, globMatDif, 1.0 };
	GLfloat matSpec[] = { globMatSpe, globMatSpe, globMatSpe, 1.0 };

	glMaterialfv(GL_FRONT, GL_AMBIENT, matAmbi);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, matDiff);
	glMaterialfv(GL_FRONT, GL_SPECULAR, matSpec);
	glMaterialfv(GL_FRONT, GL_SHININESS, &globExp);
}

void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);

	// Set material properties
	//setMaterials();
	glPushMatrix();
	//arbol1
	firstobj.Draw();
	firstobj.DrawBoundingBox();

	glTranslated(40, 0, 0);
	
	firstobj.Draw();
	firstobj.DrawBoundingBox();

	glPopMatrix();
	glutSwapBuffers();
}

void spinner(void)
{
	// alter angles by small amount
	xAngle += xSpeed;
	yAngle += ySpeed;
	display();
}

void reshape(int w, int h)
{
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);
	float aspect = (((float)w) / (float)h);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60.0f, aspect, 0.1f, 80.0f);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(7.0, 15.0, 70.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0);
}

// Use the arrow keys (in OpenGL) for controlling cylinder
void processSpecialKeys(int key, int x, int y) {

	switch (key) {
	case GLUT_KEY_UP:
		m += 1;
		if (m > 20)
			m = 20;
		break;
	case GLUT_KEY_DOWN:
		m -= 1;
		if (m < 2)
			m = 2;
		break;
	case GLUT_KEY_LEFT:
		n += 1;
		if (n > 20)
			n = 20;
		break;
	case GLUT_KEY_RIGHT:
		n -= 1;
		if (n < 1)
			n = 1;
		break;
	}
	glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y)
{
	switch (key) {

	case 't':
		if (extrudeORsor)
			extrudeORsor = false;
		else
			extrudeORsor = true;

		createSORPoints();
		glutPostRedisplay();
		break;

	case '*':
		xSpeed += 0.1;
		break;
	case '+':
		xSpeed -= 0.1;
		break;
	case '-':
		ySpeed += 0.1;
		break;
	case '_':
		ySpeed -= 0.1;
		break;

		// Set Smooth Shading (Gouraud)
	case 's':
		glShadeModel(GL_SMOOTH);
		break;
		// Set Flat Shading
	case 'f':
		glShadeModel(GL_FLAT);
		break;

		// Increase Ambient Light component (light)
	case 'A':
		globAmb = min(1.0, globAmb + 0.1);
		break;
		// Decrease Ambient Light component (light)
	case 'a':
		globAmb = max(0.0, globAmb - 0.1);
		break;

		// Increase Diffuse Light component (light)
	case 'D':
		globDif = min(1.0, globDif + 0.1);
		break;
		// Decrease Diffuse Light component (light)
	case 'd':
		globDif = max(0.0, globDif - 0.1);
		break;

		// Increase Specular Light component (light)
	case 'P':
		globSpe = min(1.0, globSpe + 0.1);
		break;
		// Decrease Specular Light component (light)
	case 'p':
		globSpe = max(0.0, globSpe - 0.1);
		break;

		// Increase Ambient refelctivity component (material)
	case 'M':
		globMatAmb = min(1.0, globMatAmb + 0.1);
		break;
		// Decrease Ambient refelctivity component (material)
	case 'm':
		globMatAmb = max(0.0, globMatAmb - 0.1);
		break;

		// Increase Diffuse refelctivity component (material)
	case 'N':
		globMatDif = min(1.0, globMatDif + 0.1);
		break;
		// Decrease Diffuse refelctivity component (material)
	case 'n':
		globMatDif = max(0.0, globMatDif - 0.1);
		break;

		// Increase Specular refelctivity component (material)
	case 'B':
		globMatSpe = min(1.0, globMatSpe + 0.1);
		break;
		// Decrease Specular refelctivity component (material)
	case 'b':
		globMatSpe = max(0.0, globMatSpe - 0.1);
		break;

		// Increase Specular Exponent (material)
	case 'E':
		globExp = min(128.0, globExp + 2.0);
		break;
		// Decrease Specular Exponent (material)
	case 'e':
		globExp = max(0.0, globExp - 2.0);
		break;

		// Reset All Light Values
	case 'R':
		globAmb = 0.3;
		globDif = 0.4;
		globSpe = 0.7;
		globExp = 4.0;

		globMatAmb = 0.5;
		globMatDif = 0.5;
		globMatSpe = 0.5;

		break;
	case 'r':
		xAngle = 0.0;
		yAngle = 0.0;
		break;

		// Quit
	case 0x1B:
	case 'q':
	case 'Q':
		exit(0);
		break;
	default:
		break;
	}
}

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(640, 480);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Final Project");

	init();

	glutKeyboardFunc(keyboard);
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutIdleFunc(spinner);
	glutSpecialFunc(processSpecialKeys);

	glutMainLoop();
	return 0;
}