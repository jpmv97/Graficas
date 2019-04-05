#pragma once
#pragma warning(disable : 4996)
#include <windows.h>
#include <assert.h>
#include <math.h>

#include <string>
#include <iostream>
#include <fstream>
#include <strstream>

#include <stdlib.h>
#include <vector>

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/gl.h>
#include "glut.h"

using namespace std;

typedef unsigned short ushort;
typedef unsigned long ulong;
typedef unsigned char uchar;

// ###########################################################################################//
/// CG_Tec.h ////
// J Guerrero J
// Classes and meshes based on code examples from Computer Graphics Usign OpenGL, Hill & Kelley
// OBJ Reader based on opengl-tutorials.com
// ###########################################################################################//

//############################### class Point2 #######################
class Point2 {
public:
	float x, y;

	Point2(float xx, float yy)
	{
		x = xx; y = yy;
	}

	Point2()
	{
		x = y = 0;
	}

	void set(float xx, float yy)
	{
		x = xx; y = yy;
	}
};

//@@@@@@@@@@@@@@@@@@ Point3 class @@@@@@@@@@@@@@@@
class Point3 {
public:
	float x, y, z;

	Point3(float xx, float yy, float zz) { x = xx; y = yy; z = zz; }
	Point3() { x = y = z = 0; }

	void build4tuple(float v[])
	{
		// load 4-tuple with this color: v[3] = 1 for homogeneous
		v[0] = x; v[1] = y; v[2] = z; v[3] = 1.0f;
	}
	void set(float dx, float dy, float dz) { x = dx; y = dy; z = dz; }
	void set(Point3 p) { x = p.x; y = p.y; z = p.z; }
};

//@@@@@@@@@@@@@@@@@@ Vector3 class @@@@@@@@@@@@@@@@
class Vector3 {
public:
	float x, y, z;

	Vector3(float xx, float yy, float zz) { x = xx; y = yy; z = zz; }
	Vector3(Vector3& v) { x = v.x; y = v.y; z = v.z; }
	Vector3() { x = y = z = 0; } //default constructor
	Vector3 cross(Vector3 b) //return this cross b
	{
		Vector3 c(y*b.z - z * b.y, z*b.x - x * b.z, x*b.y - y * b.x);
		return c;
	}

	float dot(Vector3 b) // return this dotted with b
	{
		return x * b.x + y * b.y + z * b.z;
	}
	void flip() { x = -x; y = -y; z = -z; } // reverse this vector
	void normalize() //adjust this vector to unit length
	{
		double sizeSq = x * x + y * y + z * z;
		if (sizeSq < 0.0000001)
		{
			cerr << "\nnormalize() sees vector (0,0,0)!";
			return; // does nothing to zero vectors;
		}
		float scaleFactor = 1.0 / (float)sqrt(sizeSq);
		x *= scaleFactor; y *= scaleFactor; z *= scaleFactor;
	}
	void set(float dx, float dy, float dz) { x = dx; y = dy; z = dz; }
	void set(const Vector3& v) { x = v.x; y = v.y; z = v.z; }
	void setDiff(Point3& a, Point3& b)//set to difference a - b
	{
		x = a.x - b.x; y = a.y - b.y; z = a.z - b.z;
	}
};

// @@@@@@@@@@@@@@@@@@@@@ Color3 class @@@@@@@@@@@@@@@@
class Color3 { // holds an red,green,blue 3-tuple
public:
	float red, green, blue;

	Color3()
	{
		red = green = blue = 0;
	}

	Color3(float r, float g, float b)
	{
		red = r; green = g; blue = b;
	}

	Color3(Color3& c)
	{
		red = c.red; green = c.green; blue = c.blue;
	}

	void add(float r, float g, float b)
	{
		red += r; green += g; blue += b;
	}

	void add(Color3& src, Color3& refl)
	{ // add the product of source color and reflection coefficient
		red += src.red   * refl.red;
		green += src.green * refl.green;
		blue += src.blue  * refl.blue;
	}

	void add(Color3& colr)
	{ // add colr to this color
		red += colr.red; green += colr.green; blue += colr.blue;
	}

	void build4tuple(float v[])
	{// load 4-tuple with this color: v[3] = 1 for homogeneous
		v[0] = red; v[1] = green; v[2] = blue; v[3] = 1.0f;
	}

	void set(float r, float g, float b)
	{
		red = r; green = g; blue = b;
	}

	void set(Color3& c)
	{
		red = c.red; green = c.green; blue = c.blue;
	}
};

//@@@@@@@@@@@@@@@@@ Material class @@@@@@@@@@@@@@
class Material {
public:
	Color3 ambient, diffuse, specular, emissive;
	int illum;					// Illumination parameters, as per OBJ

	float specularExponent; //, reflectivity, transparency, speedOfLight;
	//float specularFraction, surfaceRoughness;	

	GLuint textureName[1];
	GLfloat textureParam;

	char mtlName[32]; // NOTE! Assume name will not be longer than 32 characters

	Material() {
		setDefault();
	}

	Material(const Material& m) {
		set((Material &)m);
	}

	~Material() {}

	void set(Material& m)
	{

		//textureType = m.textureType;
		//numParams = m.numParams;
		//for(int i = 0; i < numParams; i++)
		//	params[i] = m.params[i];
		illum = m.illum;
		//transparency = m.transparency;
		//speedOfLight = m.speedOfLight;
		//reflectivity = m.reflectivity;
		specularExponent = m.specularExponent;
		//specularFraction = m.specularFraction;
		//surfaceRoughness = m.surfaceRoughness;
		ambient.set(m.ambient);
		diffuse.set(m.diffuse);
		specular.set(m.specular);
		emissive.set(m.emissive);

		textureParam = m.textureParam;
		textureName[0] = m.textureName[0];

		//for(int i = 0; i < 32; i++)
		//	mtlName[i] = m.mtlName[i];
		memcpy(&mtlName[0], &m.mtlName[0], 32 * sizeof(char));
	}

	void setDefault()
	{
		//textureType = 0; // for none
		//numParams = 0;
		//reflectivity = transparency = 0.0;
		//speedOfLight = 
		specularExponent = 1.0;
		//specularFraction = 0.0;
		//surfaceRoughness = 1.0;
		ambient.set(0.1f, 0.1f, 0.1f);
		diffuse.set(0.8f, 0.8f, 0.8f);
		specular.set(0, 0, 0);
		emissive.set(0, 0, 0);

		textureParam = GL_MODULATE;
		textureName[0] = -1;
	}

	//void ResetAndKeepMtl
};

// @@@@@@@@@@@@@@@@@@@@@ Affine4 class @@@@@@@@@@@@@@@@
class Affine4 {
	// manages homogeneous affine transformations
	// including inverse transformations
	// and a stack to put them on
public:
	float m[16]; // hold a 4 by 4 matrix (row-major)

	Affine4()
	{
		// make identity transform
		// NOTE: Row-major
		m[0] = m[5] = m[10] = m[15] = 1.0;
		m[1] = m[2] = m[3] = m[4] = 0.0;
		m[6] = m[7] = m[8] = m[9] = 0.0;
		m[11] = m[12] = m[13] = m[14] = 0.0;
	}

	void preMult(Affine4 n)
	{	// premultiplies this with n
		float sum;
		Affine4 tmp;
		tmp.set(*this); // tmp copy
		// following mult's : this = tmp * n
		for (int c = 0; c < 4; c++)
			for (int r = 0; r < 4; r++)
			{
				sum = 0;
				for (int k = 0; k < 4; k++)
					sum += n.m[4 * k + r] * tmp.m[4 * c + k];
				m[4 * c + r] = sum;
			}// end of for loops
	}

	void postMult(Affine4 n)
	{// postmultiplies this with n
		float sum;
		Affine4 tmp;
		tmp.set(*this); // tmp copy
		for (int c = 0; c < 4; c++)// form this = tmp * n
			for (int r = 0; r < 4; r++)
			{
				sum = 0;
				for (int k = 0; k < 4; k++)
					sum += tmp.m[4 * k + r] * n.m[4 * c + k];
				m[4 * c + r] = sum;
			}// end of for loops
	}

	void set(Affine4 a)
	{	// set this matrix to a
		for (int i = 0; i < 16; i++)
			m[i] = a.m[i];
	}

	void setIdentityMatrix()
	{	// make identity transform
		m[0] = m[5] = m[10] = m[15] = 1.0;
		m[1] = m[2] = m[3] = m[4] = 0.0;
		m[6] = m[7] = m[8] = m[9] = 0.0;
		m[11] = m[12] = m[13] = m[14] = 0.0;
	}

	// Apply a Rotation Matrix and Post/Pre multiply with current m
	void ApplyRotateX(float angle, bool post = true)
	{
		// { 1,  0,  0,  0,
		//   0,  c, -s,  0,
		//   0,  s,  c,  0,
		//   0,  0,  0,  1 };
		Affine4 temp;				// Default is I
		float ca = cos(angle);
		float sa = sin(angle);
		temp.m[5] = ca;
		temp.m[6] = sa;
		temp.m[9] = -sa;
		temp.m[10] = ca;

		if (post)
			postMult(temp);
		else
			preMult(temp);
	}

	// Apply a Rotation Matrix and Post/Pre multiply with current m
	void ApplyRotateY(float angle, bool post = true)
	{
		// { c,  0,  s,  0,
		//   0,  1,  0,  0,
		//  -s,  0,  c,  0,
		//   0,  0,  0,  1 };
		Affine4 temp;				// Default is I
		float ca = cos(angle);
		float sa = sin(angle);
		temp.m[0] = ca;
		temp.m[2] = -sa;
		temp.m[8] = sa;
		temp.m[10] = ca;

		if (post)
			postMult(temp);
		else
			preMult(temp);
	}

	// Apply a Rotation Matrix and Post/Pre multiply with current m
	void ApplyRotateZ(float angle, bool post = true)
	{
		// { c, -s,  0,  0,
		//   s,  c,  0,  0,
		//   0,  0,  1,  0,
		//   0,  0,  0,  1 };
		Affine4 temp;				// Default is I
		float ca = cos(angle);
		float sa = sin(angle);
		temp.m[0] = ca;
		temp.m[1] = sa;
		temp.m[4] = -sa;
		temp.m[5] = ca;

		if (post)
			postMult(temp);
		else
			preMult(temp);
	}

	// Apply a Translation Matrix and Post/Pre multiply with current m
	void ApplyTranslate(Point3 translate, bool post = true)
	{
		// { 1,  0,  0, tX,
		//   0,  1,  0, tY,
		//   0,  0,  1, tZ,
		//   0,  0,  0,  1 };
		Affine4 temp;				// Default is I
		temp.m[12] = translate.x;
		temp.m[13] = translate.y;
		temp.m[14] = translate.z;

		if (post)
			postMult(temp);
		else
			preMult(temp);
	}

	// Apply a Scaling Matrix and Post/Pre multiply with current m
	void ApplyScale(Point3 scale, bool post = true)
	{
		// {sX,  0,  0,  0,
		//   0, sY,  0,  0,
		//   0,  0, sZ,  0,
		//   0,  0,  0,  1 };
		Affine4 temp;				// Default is I
		temp.m[0] = scale.x;
		temp.m[5] = scale.y;
		temp.m[10] = scale.z;

		if (post)
			postMult(temp);
		else
			preMult(temp);
	}

	Point3 MultiplyPoint(Point3 point)
	{
		Point3 retval;
		/// do matrix by vector mulitiplication
		retval.x = m[0] * point.x + m[4] * point.y + m[8] * point.z + m[12];
		retval.y = m[1] * point.x + m[5] * point.y + m[9] * point.z + m[13];
		retval.z = m[2] * point.x + m[6] * point.y + m[10] * point.z + m[14];
		// Technically, we should also have a fourth line here - but nowhere to store it in a Point3!

		return retval;
	}

}; // end of Affine4 class

//  --- ------- ---- Classes for Texture Reading --- ---- --- ----
// @@@@@@@@@@@@@@@@@@@@@ mRGB class @@@@@@@@@@@@@@@@
class mRGB {
public:
	unsigned char r, g, b;

	mRGB() { r = g = b = 0; }
	mRGB(mRGB& p) { r = p.r; g = p.g; b = p.b; }
	mRGB(uchar rr, uchar gg, uchar bb) { r = rr; g = gg; b = bb; }
	void set(uchar rr, uchar gg, uchar bb) { r = rr; g = gg; b = bb; }
};

class RGBpixmap {
private:
	fstream inf;
	ushort getShort() //helper function
	{ //BMP format uses little-endian integer types
		// get a 2-byte integer stored in little-endian form
		char ic;
		ushort ip;
		inf.get(ic); ip = ic;  //first byte is little one
		inf.get(ic);  ip |= ((ushort)ic << 8); // or in high order byte
		return ip;
	}
	ulong getLong() //helper function
	{  //BMP format uses little-endian integer types
		// get a 4-byte integer stored in little-endian form
		ulong ip = 0;
		char ic = 0;
		unsigned char uc = ic;
		inf.get(ic); uc = ic; ip = uc;
		inf.get(ic); uc = ic; ip |= ((ulong)uc << 8);
		inf.get(ic); uc = ic; ip |= ((ulong)uc << 16);
		inf.get(ic); uc = ic; ip |= ((ulong)uc << 24);
		return ip;
	}

public:
	int nRows, nCols;
	mRGB *pixel;

	// Read into memory an mRGB image from an uncompressed BMP file.
	// return false on failure, true on success
	bool readBMPFile(char *fname)
	{
		//read binary char's
		inf.open(fname, ios::in | ios::binary);
		if (!inf)
		{
			cout << " can't open file: " << fname << endl;
			return false;
		}

		int k, row, col, numPadBytes, nBytesInRow;

		// read the file header information
		char ch1, ch2;
		inf.get(ch1); inf.get(ch2); //type: always 'BM'
		ulong fileSize = getLong();
		ushort reserved1 = getShort();    // always 0
		ushort reserved2 = getShort();     // always 0
		ulong offBits = getLong(); // offset to image - unreliable
		ulong headerSize = getLong();     // always 40
		ulong numCols = getLong(); // number of columns in image
		ulong numRows = getLong(); // number of rows in image
		ushort planes = getShort();      // always 1
		ushort bitsPerPixel = getShort();    //8 or 24; allow 24 here
		ulong compression = getLong();      // must be 0 for uncompressed
		ulong imageSize = getLong();       // total bytes in image
		ulong xPels = getLong();    // always 0
		ulong yPels = getLong();    // always 0
		ulong numLUTentries = getLong();    // 256 for 8 bit, otherwise 0
		ulong impColors = getLong();       // always 0

		if (bitsPerPixel != 24)
		{ // error - must be a 24 bit uncompressed image
			cout << "not a 24 bit/pixelimage (RGB), or is compressed!\n";
			inf.close();
			return false;
		}

		//add bytes at end of each row so total # is a multiple of 4
		// round up 3*numCols to next mult. of 4
		nBytesInRow = ((3 * numCols + 3) / 4) * 4;
		numPadBytes = nBytesInRow - 3 * numCols; // need this many
		nRows = numRows; // set class's data members
		nCols = numCols;
		pixel = new mRGB[nRows * nCols]; //make space for array
		if (!pixel) {
			cout << "Out of memory!\n";
			return false; // out of memory!
		}

		long count = 0;
		char dum;
		for (row = 0; row < nRows; row++) // read pixel values
		{
			for (col = 0; col < nCols; col++)
			{
				char r, g, b;
				inf.get(b); inf.get(g); inf.get(r); //read bytes
				pixel[count].r = r; //place them in colors
				pixel[count].g = g;
				pixel[count++].b = b;
			}
			for (k = 0; k < numPadBytes; k++) //skip pad bytes at row's end
				inf >> dum;
		}
		inf.close();

		return true; // success
	}
};

void InitializeTexture(RGBpixmap *textureObject, GLuint *texturename, char *filename)
{
	// Get OpenGL to automatically generate the texture "names"
	glGenTextures(1, texturename);

	// Create the texture data
	textureObject->readBMPFile(filename);

	// Assign the texture - Link the information in pixel to OpenGL
	// Use "Bind" to let OpenGL know what texture we are using
	glBindTexture(GL_TEXTURE_2D, texturename[0]);

	// Set wrapping and interpolation
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	// Set how the texture will be mapped to our polygon(s)
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	// Specify how OpenGL will read our pixel values (by byte)
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	// "Create" the OpenGL texture - link the actual texture (from name) to the image data (pixel)
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, textureObject->nCols, textureObject->nRows, 0, GL_RGB, GL_UNSIGNED_BYTE, textureObject->pixel);

	// Actually enable the texture
	glEnable(GL_TEXTURE_2D);
};
// --- ------- ---- End Classes for Texture Reading --- -- --- ----

//################## class FaceID ################
class FaceID {
public:
	std::vector< int > vertIndex;
	std::vector< int > uvIndex;
	std::vector< int > normIndex;
};

//################## class FaceList ################
class FaceList {
public:
	int numVerts;					// How many vertices (data) per face
	Material mtrl;
	std::vector< FaceID > faces;	// The data

	FaceList() { numVerts = 0; }
	~FaceList() {}
	FaceList(const FaceList &obj) {
		for (unsigned int i = 0; i < obj.faces.size(); i++)
			faces.push_back(obj.faces[i]);
		numVerts = obj.numVerts;
		mtrl = obj.mtrl;
	}
};

//@$@$@$@$@$@$@$@$@$@ Mesh class @$@$@$@$@$@$@$@$@$
// TO DO (?): Mesh Subdivision
class Mesh
{
public:
	enum RenderMode {
		MODE_WIRE = 0,
		MODE_SOLID,
		MODE_WIRE_SOLID
	};

private:
	RenderMode mode;
	Material lastmtrl;
	Affine4 transf, invTransf;
	std::vector< Material > materials;
	bool debug;
	Point3 min;
	Point3 max;

public:
	std::vector< Point3 > verts;
	std::vector< Point2 > uvs;
	std::vector< Point3 > norms;
	std::vector< FaceList > meshFaces;

private:
	void ProcessBounds(Point3 point)
	{
		// Check if this point is greater than max
		if (point.x > max.x)
			max.x = point.x;

		if (point.y > max.y)
			max.y = point.y;

		if (point.z > max.z)
			max.z = point.z;

		// Check if this point is less than min
		if (point.x < min.x)
			min.x = point.x;

		if (point.y < min.y)
			min.y = point.y;

		if (point.z < min.z)
			min.z = point.z;
	}

	bool readTexture(char * path, Material * currentMaterial)
	{
		RGBpixmap theTexture;
		if (!theTexture.readBMPFile(path))
		{
			printf("Impossible to open the Texture file (via Mesh) !\n");
			return false;
		}

		// Generate OpenGL Texture Name
		glGenTextures(1, &currentMaterial->textureName[0]);

		// Assign the texture - Link the information in pixel to OpenGL
		// Use "Bind" to let OpenGL know what texture we are using
		glBindTexture(GL_TEXTURE_2D, currentMaterial->textureName[0]);

		// Set (my default) wrapping and interpolation
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

		// Set how the texture will be mapped to our polygon(s)
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, currentMaterial->textureParam);

		// Specify how OpenGL will read our pixel values (by byte)
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

		// "Create" the OpenGL texture - link the actual texture (from name) to the image data (pixel)
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, theTexture.nCols, theTexture.nRows, 0, GL_RGB, GL_UNSIGNED_BYTE, theTexture.pixel);

		// Enable Texture done when drawing

		return true;
	}

	bool readMTL(const char * path)
	{
		int countMtls = 0;
		Material currentMaterial;

		FILE * mtlfile = fopen(path, "r");
		if (mtlfile == NULL) {
			printf("Impossible to open the mtl file !\n");
			return false;
		}

		// Parse...
		while (1)
		{
			char mtlLineHeader[128];

			// read the first word of the line
			int res = fscanf(mtlfile, "%s", mtlLineHeader);

			if (res == EOF)
				break;
			// else
			if (strcmp(mtlLineHeader, "newmtl") == 0)
			{
				// We have a new material definition
				// Check to see if we already have a material defined
				// Assume our array/vector is empty if this is the first one; otherwise, push back previous mtl
				if (countMtls != 0)
					materials.push_back(currentMaterial);

				// Reset material values
				currentMaterial.setDefault();
				fscanf(mtlfile, "%s", &currentMaterial.mtlName[0]);
				countMtls++;
			}
			else if (strcmp(mtlLineHeader, "Kd") == 0)
			{
				fscanf(mtlfile, "%f %f %f\n", &currentMaterial.diffuse.red, &currentMaterial.diffuse.green, &currentMaterial.diffuse.blue);
			}
			else if (strcmp(mtlLineHeader, "Ka") == 0)
			{
				fscanf(mtlfile, "%f %f %f\n", &currentMaterial.ambient.red, &currentMaterial.ambient.green, &currentMaterial.ambient.blue);
			}
			else if (strcmp(mtlLineHeader, "Ks") == 0)
			{
				fscanf(mtlfile, "%f %f %f\n", &currentMaterial.specular.red, &currentMaterial.specular.green, &currentMaterial.specular.blue);
			}
			else if (strcmp(mtlLineHeader, "Ns") == 0)
			{
				fscanf(mtlfile, "%f\n", &currentMaterial.specularExponent);
			}
			else if (strcmp(mtlLineHeader, "Ke") == 0)
			{
				fscanf(mtlfile, "%f %f %f\n", &currentMaterial.emissive.red, &currentMaterial.emissive.green, &currentMaterial.emissive.blue);
			}
			else if (strcmp(mtlLineHeader, "d") == 0)
			{
				float d[1];
				fscanf(mtlfile, "%f\n", &d[0]);
				//currentMaterial.transparency = 1.0 - d[0];
			}
			else if (strcmp(mtlLineHeader, "Tr") == 0)
			{
				float d[1];
				fscanf(mtlfile, "%f\n", &d[0]);
				//currentMaterial.transparency = d[0];
				////fscanf(mtlfile, "%f\n", &currentMaterial.transparency);
			}
			else if (strcmp(mtlLineHeader, "illum") == 0)
			{
				//Illumination parameters - various options
				//0. Color on and Ambient off
				//1. Color on and Ambient on
				//2. Highlight on
				//3. Reflection on and Ray trace on
				//4. Transparency: Glass on, Reflection: Ray trace on
				//5. Reflection: Fresnel on and Ray trace on
				//6. Transparency: Refraction on, Reflection: Fresnel off and Ray trace on
				//7. Transparency: Refraction on, Reflection: Fresnel on and Ray trace on
				//8. Reflection on and Ray trace off
				//9. Transparency: Glass on, Reflection: Ray trace off
				//10. Casts shadows onto invisible surfaces
				fscanf(mtlfile, "%i\n", &currentMaterial.illum);
			}
			// Texture map - AMBIENT - NOTE: Basic OpenGL uses only one at a time
			else if (strcmp(mtlLineHeader, "map_Ka") == 0)
			{
				// Will ignore - use only diffuse
				char name[32];
				fscanf(mtlfile, "%s\n", &name[0]);
			}
			// Texture map - DIFFUSE - NOTE: Basic OpenGL uses only one at a time
			else if (strcmp(mtlLineHeader, "map_Kd") == 0)
			{
				// Load image and prepare texture
				char name[32];
				fscanf(mtlfile, "%s\n", &name[0]);

				readTexture(name, &currentMaterial);
			}
			// Texture map - SPECULAR - NOTE: Basic OpenGL uses only one at a time
			else if (strcmp(mtlLineHeader, "map_Ks") == 0)
			{
				// Will ignore - use only diffuse
				char name[32];
				fscanf(mtlfile, "%s\n", &name[0]);
			}
		}

		// Add last (or only) material to array
		materials.push_back(currentMaterial);

		return true;
	}

	void tellMaterialsGL(int i)
	{
		glEnable(GL_LIGHTING);

		// Set material properties, as defined previously
		glMaterialfv(GL_FRONT, GL_SPECULAR, (GLfloat*)&meshFaces[i].mtrl.specular);
		glMaterialfv(GL_FRONT, GL_SHININESS, (GLfloat*)&meshFaces[i].mtrl.specularExponent);
		glMaterialfv(GL_FRONT, GL_AMBIENT, (GLfloat*)&meshFaces[i].mtrl.ambient);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, (GLfloat*)&meshFaces[i].mtrl.diffuse);
		glMaterialfv(GL_FRONT, GL_EMISSION, (GLfloat*)&meshFaces[i].mtrl.emissive);

		// Set Texture properties, if applicable
		if (meshFaces[i].mtrl.textureName > 0) {
			glEnable(GL_TEXTURE_2D);
			glBindTexture(GL_TEXTURE_2D, meshFaces[i].mtrl.textureName[0]);
			glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, meshFaces[i].mtrl.textureParam);
		}
		else {
			glDisable(GL_TEXTURE_2D);
		}
	}

	void DrawFace(int mesh, int face, int verticesPerFace)
	{
		bool uvData = true;
		bool normalData = true;
		int normIndex, uvIndex, vertIndex;

		if (meshFaces[mesh].faces[face].uvIndex.size() == 0)
			uvData = false;

		if (meshFaces[mesh].faces[face].normIndex.size() == 0)
			normalData = false;

		for (int k = 0; k < verticesPerFace; k++) {

			if (normalData)
				normIndex = meshFaces[mesh].faces[face].normIndex[k];
			if (uvData)
				uvIndex = meshFaces[mesh].faces[face].uvIndex[k];
			vertIndex = meshFaces[mesh].faces[face].vertIndex[k];

			if (debug)
				printf("%d %d %d\n", normIndex, uvIndex, vertIndex);

			if (uvData)
				glTexCoord2fv((GLfloat*)&uvs[uvIndex - 1]);
			if (normalData)
				glNormal3fv((GLfloat*)&norms[normIndex - 1]);
			glVertex3fv((GLfloat*)&verts[vertIndex - 1]);
		}
	}

	void DrawFaces(int i)
	{
		// Find out how many vertices(data) per face
		int numVerts = meshFaces[i].numVerts;

		if (numVerts == 3)
			glBegin(GL_TRIANGLES);
		else if (numVerts == 4)
			glBegin(GL_QUADS);

		// pass index and data to render
		for (unsigned int j = 0; j < meshFaces[i].faces.size(); j++) {
			DrawFace(i, j, numVerts);
			if (debug)
				printf("%d %d %d\n", i, j, numVerts);
		}
		glEnd();
	}

	void DrawEdges(int i)
	{
		// Find out how many vertices(data) per face
		int numVerts = meshFaces[i].numVerts;

		// pass index and data to render
		for (unsigned int j = 0; j < meshFaces[i].faces.size(); j++) {
			glBegin(GL_LINE_LOOP);
			DrawFace(i, j, numVerts);
			glEnd();
		}
	}

public:
	Mesh() {
		mode = MODE_SOLID;
		debug = false;
	}

	Mesh(const char * path) {
		mode = MODE_SOLID;
		readOBJ(path);
		debug = false;
	}

	// Read data from OBJ file 
	bool readOBJ(const char * path)
	{
		bool negIndex = false;	// Used to determine if we have negative indices
		bool noUVs = false;		// Used if UV data not included
		bool noNormals = false;	// Used if normal data not included
		std::vector< unsigned int > vertexIndices, uvIndices, normalIndices;
		FaceList currentFaceList;

		FILE * file = fopen(path, "r");
		if (file == NULL) {
			printf("Impossible to open the OBJ file !\n");
			return false;
		}

		// OBJ headers (some of them...)
		// v	vertex (x, y, z [,w])
		// vt	texture coords (u, [v ,w])
		// vn	normales (x,y,z)
		// vp	parameter space ( u [,v] [,w] )
		// f	faces (vertex/texture/normal) - Base 1 (i.e. start at 1)
		// l	line element
		//
		// mtllib [external .mtl file name]
		// usemtl [material name]
		//
		// o [object name]
		// g [group name]
		//
		// s [option]	option = "1" or "off"
		while (1)
		{
			char lineHeader[128];

			// read the first word of the line
			int res = fscanf(file, "%s", lineHeader);

			if (res == EOF)
				break; // EOF = End Of File. Quit the loop.
			// else : parse lineHeader
			//
			// TO DO : 
			// Parse "o", "g", "s", "l", "vp"
			if (strcmp(lineHeader, "#") == 0) {
				char lineFromFile[128];
				fscanf(file, "%[^\n]%*c", lineFromFile);
				if (debug)
					printf("Comment in OBJ file: %s\n", lineFromFile);
			}
			else if (strcmp(lineHeader, "v") == 0) {
				Point3 vertex;
				fscanf(file, "%f %f %f\n", &vertex.x, &vertex.y, &vertex.z);
				verts.push_back(vertex);
				ProcessBounds(vertex);
			}
			else if (strcmp(lineHeader, "vt") == 0) {
				Point2 uv;
				fscanf(file, "%f %f\n", &uv.x, &uv.y);
				uvs.push_back(uv);
			}
			else if (strcmp(lineHeader, "vn") == 0) {
				Point3 normal;
				fscanf(file, "%f %f %f\n", &normal.x, &normal.y, &normal.z);
				norms.push_back(normal);
			}
			else if (strcmp(lineHeader, "f") == 0) {
				// Read full line - we do not know how many vertices per face (nor if data complete v/vt/vn)
				char lineFromFile[128];

				fscanf(file, "%[^\n]%*c", lineFromFile);
				if (debug)
					printf("%s\n", lineFromFile);

				// We have to read characters and determine what the format is, and how many data points we have.
				// "v/vt/vn" or "v//vn" or "v/vt"
				// all of the above either 3 or 4 (or 5 or more?) data points
				// first character can also be "-" for negative face indices
				//bool negIndex = false;
				int i = 0;
				char c[2] = "\n";		// size = 2 to include a terminating null-character 
				char cPrev[2] = "\n";

				int from = 0;
				int face = 0;
				int slashes = 0;
				bool prevCharIsSlash = false;
				char lineToParse[32];	// assume our biggest data point v/t/n will not be larger than this
				char nullChars[32];		// pointer with null values to reset 
				for (int i = 0; i < 32; i++)
					nullChars[i] = 0;

				// Possiblity of up to 6 vertices per face!
				unsigned int vertexIndex[6], uvIndex[6], normalIndex[6];
				int matches = 0;

				while (lineFromFile[i])
				{
					c[0] = lineFromFile[i];

					// We reached the end of the string
					if (strcmp(c, "\n") == 0) {
						break;
					}
					else if (strcmp(c, "-") == 0) {
						negIndex = true;
					}
					else if (strcmp(c, "/") == 0) {
						slashes++;
						if (strcmp(cPrev, "/") == 0)
							prevCharIsSlash = true;
					}
					//else if ( strcmp( c, " " ) == 0 )
					else if (isspace(c[0])) {
						// spaces are our separators
						// can be: between the "f" and first data point
						//         between each data point
						//         between the last data point and "\n"
						// Slashes separate data types
						// make sure it is not the first sapce
						if (slashes != 0) {
							// Fill array with zeros...
							memcpy(&lineToParse[0], &nullChars[0], (32) * sizeof(char));
							// ...then memcpy
							memcpy(&lineToParse[0], &lineFromFile[from], (i - from) * sizeof(char));

							// "v/vt/vn" or "v//vn" or "v/vt"
							if (slashes == 1) {
								matches += sscanf(&lineToParse[0], "%d/%d", &vertexIndex[face], &uvIndex[face]);
								noNormals = true;
							}
							else if (slashes == 2) {
								if (prevCharIsSlash) {
									matches += sscanf(&lineToParse[0], "%d//%d", &vertexIndex[face], &normalIndex[face]);
									noUVs = true;
									prevCharIsSlash = false;
								}
								else {
									matches += sscanf(&lineToParse[0], "%d/%d/%d", &vertexIndex[face], &uvIndex[face], &normalIndex[face]);
								}
							}
							slashes = 0;
							face++;
							from = i;
						}
					}

					if (debug)
						printf("%c", c[0]);

					cPrev[0] = c[0];
					i++;
				}

				// Handle the last tuple of data
				if (slashes != 0) {
					// Fill array with zeros...
					memcpy(&lineToParse[0], &nullChars[0], (32) * sizeof(char));
					// ...then memcpy
					memcpy(&lineToParse[0], &lineFromFile[from], (i - from) * sizeof(char));

					// "v/vt/vn" or "v//vn" or "v/vt"
					if (slashes == 1) {
						matches += sscanf(&lineToParse[0], "%d/%d", &vertexIndex[face], &uvIndex[face]);
						noNormals = true;
					}
					else if (slashes == 2) {
						if (prevCharIsSlash) {
							matches += sscanf(&lineToParse[0], "%d//%d", &vertexIndex[face], &normalIndex[face]);
							noUVs = true;
							prevCharIsSlash = false;
						}
						else {
							matches += sscanf(&lineToParse[0], "%d/%d/%d", &vertexIndex[face], &uvIndex[face], &normalIndex[face]);
						}
					}
					slashes = 0;
					face++;
				}

				if (debug)
					printf("\n");

				// Check to see if there was a change in number of vertices per face
				if (currentFaceList.numVerts != face) {
					// Push previous MeshList Object on to vector
					if (currentFaceList.faces.size() != 0)
						meshFaces.push_back(currentFaceList);

					// Reset (Create) new MeshList Object
					currentFaceList.mtrl.set(lastmtrl);
					currentFaceList.numVerts = face;
					currentFaceList.faces.clear();
				}

				if ((matches != 3 * currentFaceList.numVerts && !noUVs && !noNormals) ||
					(matches != 2 * currentFaceList.numVerts && (noUVs || noNormals))) {
					printf("File can't be read by our simple parser : ( Try exporting with other options\n");
					break;
				}

				FaceID thisFace;
				for (int ii = 0; ii < currentFaceList.numVerts; ii++) {
					thisFace.vertIndex.push_back(vertexIndex[ii]);
					if (!noUVs)
						thisFace.uvIndex.push_back(uvIndex[ii]);
					if (!noNormals)
						thisFace.normIndex.push_back(normalIndex[ii]);
				}
				currentFaceList.faces.push_back(thisFace);

			}
			else if (strcmp(lineHeader, "mtllib") == 0) {
				char mtlFileName[128];
				fscanf(file, "%s", mtlFileName);

				readMTL(mtlFileName);
			}
			else if (strcmp(lineHeader, "usemtl") == 0) {
				char mtlName[128];
				fscanf(file, "%s", mtlName);

				if (materials.size() == 0)
					printf("No materials loaded! Cannot assign material.\n");
				else
				{
					for (unsigned int i = 0; i < materials.size(); i++)
					{
						if (strcmp(&mtlName[0], &materials[i].mtlName[0]) == 0)
						{
							lastmtrl.set(materials[i]);
							currentFaceList.mtrl.set(lastmtrl);
							printf("Material %s assigned to current Facelist.\n", &mtlName[0]);
							break;
						}
					}
				}


			}
		}

		// Push the last faceList on to vector
		if (currentFaceList.faces.size() != 0)
			meshFaces.push_back(currentFaceList);

		// Handle negative indices
		if (negIndex)
		{
			int normIndex, uvIndex, vertIndex;

			int numVerts = verts.size();
			int numUVs = uvs.size();
			int numNorms = norms.size();
			int numMeshes = meshFaces.size();

			for (int i = 0; i < numMeshes; i++)
			{
				int numVertsPerFace = meshFaces[i].numVerts;
				int numFaces = meshFaces[i].faces.size();

				for (int j = 0; j < numFaces; j++) {
					for (int k = 0; k < numVertsPerFace; k++) {
						if (!noNormals)
							normIndex = meshFaces[i].faces[j].normIndex[k];
						if (!noUVs)
							uvIndex = meshFaces[i].faces[j].uvIndex[k];
						vertIndex = meshFaces[i].faces[j].vertIndex[k];

						// Invert the indices, so that they are positive (and ONE-based)
						if (!noNormals)
							meshFaces[i].faces[j].normIndex[k] = numNorms + normIndex + 1;
						if (!noUVs)
							meshFaces[i].faces[j].uvIndex[k] = numUVs + uvIndex + 1;
						meshFaces[i].faces[j].vertIndex[k] = numVerts + vertIndex + 1;
					}
				}
			}
		}

		return true;
	}

	// Free up memory used by this mesh.
	void freeMesh() {
		makeEmpty();
		verts.shrink_to_fit();
		uvs.shrink_to_fit();
		norms.shrink_to_fit();
		meshFaces.shrink_to_fit();
	}

	// Check to see if we have data in this mesh
	int isEmpty()
	{
		int whatWeHave = 0;

		for (int i = 0; i < meshFaces.size(); i++)
		{
			whatWeHave += meshFaces[i].faces.size();
		}

		whatWeHave += verts.size();
		whatWeHave += uvs.size();
		whatWeHave += norms.size();

		return whatWeHave;
	}

	// Make mesh empty (without freeing memory)
	void makeEmpty()
	{
		verts.clear();
		uvs.clear();
		norms.clear();
		meshFaces.clear();
	}

	//<<<<<<<<<<<<<<<<<<<<<<<<<<<< setRenderMode >>>>>>>>>>>>>>>>>>>>>>>>
	void setRenderMode(RenderMode m)
	{
		mode = m;
	}

	void ApplyTransform(Affine4 m)
	{
		transf.postMult(m);
	}

	void ComputeBounds()
	{
		int numVerts = verts.size();
		for (int i = 0; i < numVerts; i++)
			ProcessBounds(verts[i]);
	}

	Point3 GetMinimumBounds()
	{
		return min;
	}

	Point3 GetMaximumBounds()
	{
		return max;
	}

	// Scale/Translate mesh so that it is centered at origin and fits within (but not necessarily is) a 1x1x1 cube
	// NOTE: only affine transformation is used, data is not re-written
	// NOTE: overwrites any previous transformation
	void NormalizeMesh()
	{
		// Translate mesh based on average min+max, all in negative direction
		Point3 translate;
		translate.x = -(min.x + max.x) / 2.0;
		translate.y = -(min.y + max.y) / 2.0;
		translate.z = -(min.z + max.z) / 2.0;

		// Scale mesh based on the size; assume always positive
		// Scale everything by the biggest dimension, in order to keep proportions
		Point3 scaleMesh;
		scaleMesh.x = max.x - min.x;
		scaleMesh.y = max.y - min.y;
		scaleMesh.z = max.z - min.z;

		float scale = scaleMesh.x;
		if (scale < scaleMesh.y)
			scale = scaleMesh.y;
		if (scale < scaleMesh.z)
			scale = scaleMesh.z;

		scale = 1.0 / scale;

		Affine4 m;

		m.m[12] = translate.x / scaleMesh.x;
		m.m[13] = translate.y / scaleMesh.y;
		m.m[14] = translate.z / scaleMesh.z;
		m.m[0] = scale;
		m.m[5] = scale;
		m.m[10] = scale;

		transf.set(m);
	}

	void Draw()
	{
		glPushMatrix();
		glMultMatrixf(transf.m);

		// Start render loop
		for (unsigned int i = 0; i < meshFaces.size(); i++)
		{
			tellMaterialsGL(i);

			switch (mode) {
			case MODE_WIRE:
				DrawEdges(i);
				break;
			case MODE_SOLID:
				DrawFaces(i);
				break;
			default:
				DrawEdges(i);
				DrawFaces(i);
				break;
			}
		}
		glPopMatrix();
	}

	void DrawBoundingBox()
	{
		glPushMatrix();
		glMultMatrixf(transf.m);

		glDisable(GL_LIGHTING);
		glColor3f(1.0f, 1.0f, 1.0f);

		glLineWidth(1.0);
		glBegin(GL_LINE_LOOP);
		glVertex3f(min.x, min.y, min.z);
		glVertex3f(max.x, min.y, min.z);
		glVertex3f(max.x, max.y, min.z);
		glVertex3f(min.x, max.y, min.z);
		glVertex3f(min.x, max.y, max.z);
		glVertex3f(max.x, max.y, max.z);
		glVertex3f(max.x, min.y, max.z);
		glVertex3f(min.x, min.y, max.z);
		glEnd();

		glBegin(GL_LINES);
		glVertex3f(min.x, min.y, max.z);
		glVertex3f(min.x, max.y, max.z);
		glEnd();

		glBegin(GL_LINES);
		glVertex3f(max.x, min.y, min.z);
		glVertex3f(max.x, min.y, max.z);
		glEnd();

		glBegin(GL_LINES);
		glVertex3f(max.x, max.y, min.z);
		glVertex3f(max.x, max.y, max.z);
		glEnd();

		glBegin(GL_LINES);
		glVertex3f(min.x, min.y, min.z);
		glVertex3f(min.x, max.y, min.z);
		glEnd();

		glPopMatrix();
	}
};
