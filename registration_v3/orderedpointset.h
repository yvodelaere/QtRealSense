////////////////////////////////////////////////////////////////////
//
// This class can read an abs file (ordered point set of 3D coordinates)
// into a structure and manipulate the points
//
// This class provides methods for reading and writing pointsets
// (reading ABS, writing as VRML, ABS and 3D points file)
// furthermore, methods are provided for accessing points, testing validness
// of points, obtaining basic statistics (boundingbox, average stepsize).
// 
// Also projection planes vplane,hplane are provided to project points on
// horizontal or vertical planes
//
// Author: Luuk Spreeuwers, University of Twente
// Date: 2005-2009
// Last update: 29-07-2009
// version: 2.0
//
// 29-5-2013: added (limited) support for reading VRML files
////////////////////////////////////////////////////////////////////////

#ifndef ORDEREDPOINTSET_DEFINED
#define ORDEREDPOINTSET_DEFINED

#include "point3d.h"
#include "plane3d.h"

namespace utw3dface {

class Polygon
{
public:
	int pointindex[4];

	int &operator[](int index) { return pointindex[index]; }
	Polygon(int p1,int p2, int p3, int p4)
	{
		pointindex[0]=p1;
		pointindex[1]=p2;
		pointindex[2]=p3;
		pointindex[3]=p4;
	}
	Polygon() {}

	Polygon& operator=(const Polygon p)
	{
		pointindex[0]=p.pointindex[0];
		pointindex[1]=p.pointindex[1];
		pointindex[2]=p.pointindex[2];
		pointindex[3]=p.pointindex[3];
	}

	Point3D Centre(Point3D *p)
	{
		double x=0,y=0,z=0;
		int i;
		int n=(pointindex[3]==-1) ? 3:4;
		for (i=0; i<n; i++)
		{
			x+=p[pointindex[i]].x;
			y+=p[pointindex[i]].y;
			z+=p[pointindex[i]].z;
		}
		x/=n;
		y/=n;
		z/=n;

		return Point3D(x,y,z);
	}
};

class OrderedPointSet
{
public:
	int width,height;
	int *flags;
	double *X;
	double *Y;
	double *Z;

	Polygon *polygon;
	int npolygons;

	int nvalidpoints;
	double minX,maxX;
	double minY,maxY;
	double minZ,maxZ;
	double dx,dy;

	Plane3D vplane;
	Plane3D hplane;

	OrderedPointSet();
	~OrderedPointSet();

	int Allocate(int width,int height);
	int Free();

	int GetStats();

	int GetPoint(int x,int y,double *point);
	Point3D GetPoint(int x,int y);
	int ValidPoint(int x,int y);

	int Read(const char *file);
	int ReadAbs(const char *file);
	int Read3DPoints(const char *file);
	int ReadVRML(const char *file);
	int WriteAbs(const char *file);
	int WriteVRML(const char *file);
	int WriteVRMLsurface(const char *filename);
	int Write3DPoints(const char *file);

	int Resample(double res);
};

}
#endif // ORDEREDPOINTSETDEFINED
