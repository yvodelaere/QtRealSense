////////////////////////////////////////////////////////////////////
//
// A class that defines a plane in 3D using a normal n and a constant k
// Points (x,y,z) in the plane are defined by: x*n.x + y*n.y + z*n.z = k
// 
// The class provides a methods to determine the distance of a 3D point 
// to the plane and to project a 3D point onto the plane
//
// Author: Luuk Spreeuwers, University of Twente
// Date: 2005-2009
// Last update: 29-07-2009
// version: 2.0
//
////////////////////////////////////////////////////////////////////////

#ifndef PLANE3D_DEFINED
#define PLANE3D_DEFINED

namespace utw3dface {

class Plane3D
{
public:
	Point3D n;
	double k;

	Plane3D() {}
	Plane3D(const Point3D &n,double k): n(n),k(k) {}
	~Plane3D() {}
	
	double x(double y,double z) { return (n.x==0) ? 0 : (k-n.y*y-n.z*z)/n.x; }
	double y(double x,double z) { return (n.y==0) ? 0 : (k-n.x*x-n.z*z)/n.y; }
	double z(double x,double y) { return (n.z==0) ? 0 : (k-n.x*x-n.y*y)/n.z; }

	double SignedDistance(const Point3D &p) { return k-p.x*n.x-p.y*n.y-p.z*n.z; }
	Point3D ProjectPoint(const Point3D &p) 
	{ 
		Point3D q=n;
        	q*=SignedDistance(p);
        	q+=p;
		return q;
	}
};

}
#endif // PLANE3D_DEFINED
