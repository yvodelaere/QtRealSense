////////////////////////////////////////////////////////////////////
//
// A class for 3D points
// Methods are provided for all common operations on 3D points
// like addition, subtraction, distance, normalisation, norm etc.
// 
//
// Author: Luuk Spreeuwers, University of Twente
// Date: 2005-2009
// Last update: 29-07-2009
// version: 2.0
//
////////////////////////////////////////////////////////////////////////

#ifndef POINT3D_DEFINED
#define POINT3D_DEFINED

#include <math.h>

namespace utw3dface {

class Point3D
{
public:
	double x,y,z;

	Point3D(double x,double y,double z): x(x),y(y),z(z) {}
	Point3D(const Point3D &p): x(p.x),y(p.y),z(p.z) {}
	Point3D(const double *p): x(p[0]),y(p[1]),z(p[2]) {}
	Point3D(): x(0),y(0),z(0) {}
	~Point3D() {}
	Point3D &operator=(const Point3D &p) { x=p.x; y=p.y; z=p.z; return *this; }
	Point3D &operator+=(const Point3D &p) { x+=p.x; y+=p.y; z+=p.z; return *this; }
	Point3D &operator-=(const Point3D &p) { x-=p.x; y-=p.y; z-=p.z; return *this; }
	Point3D &operator/=(const double k) { x/=k; y/=k; z/=k; return *this; }
	Point3D &operator*=(const double k) { x*=k; y*=k; z*=k; return *this; }
	double Distance(const Point3D &p) { return sqrt((x-p.x)*(x-p.x)+(y-p.y)*(y-p.y)+(z-p.z)*(z-p.z)); }
	double Norm() { return sqrt(x*x+y*y+z*z); }
	void Normalise() { *this/=(Norm()); } 
	double Inp(const Point3D &p) { return p.x*x+p.y*y+p.z*z; }
};

Point3D operator-(const Point3D &p,const Point3D &q);
Point3D operator+(const Point3D &p,const Point3D &q);
double operator*(const Point3D &p,const Point3D &q);
Point3D operator*(double d,const Point3D &p);
Point3D operator*(const Point3D &p,double d);
Point3D operator/(double d,const Point3D &p);
Point3D operator/(const Point3D &p,double d);
Point3D operator^(const Point3D &p,const Point3D &q);

}
#endif // POINT3D_DEFINED
