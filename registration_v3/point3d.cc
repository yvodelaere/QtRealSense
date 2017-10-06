#include <stdio.h>
#include <stdlib.h>
//#include <unistd.h>   
#ifdef linux
#include <malloc.h>
#endif
#include <string.h>
//#include <values.h>
#include <math.h>

#include "point3d.h"

namespace utw3dface {

Point3D operator-(const Point3D &p,const Point3D &q) { return Point3D(p.x-q.x,p.y-q.y,p.z-q.z); }
Point3D operator+(const Point3D &p,const Point3D &q) { return Point3D(p.x+q.x,p.y+q.y,p.z+q.z); }
double operator*(const Point3D &p,const Point3D &q) { return p.x*q.x+p.y*q.y+p.z*q.z; }
Point3D operator*(double d,const Point3D &p) { return Point3D(d*p.x,d*p.y,d*p.z); }
Point3D operator*(const Point3D &p,double d) { return Point3D(d*p.x,d*p.y,d*p.z); }
Point3D operator/(double d,const Point3D &p) { return Point3D(p.x/d,p.y/d,p.z/d); }
Point3D operator/(const Point3D &p,double d) { return Point3D(p.x/d,p.y/d,p.z/d); }
Point3D operator^(const Point3D &p,const Point3D &q) { return Point3D(p.y*q.z-p.z*q.y,p.z*q.x-p.x*q.z,p.x*q.y-p.y*q.x); }
}
