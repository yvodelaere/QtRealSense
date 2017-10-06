////////////////////////////////////////////////////////////////////
//
// This class represents an unordered set of 3D points
//
// Methods are provided to add a point, to access a point and to copy
// from another unordered point set. Next methods are available for
// sorting of points, writing to VRML, getting basic statistics (bounding box)
// and the centre of gravity.
// 
// Author: Luuk Spreeuwers, University of Twente
// Date: 2005-2009
// Last update: 29-07-2009
// version: 2.0
//
////////////////////////////////////////////////////////////////////////


#ifndef UNORDEREDPOINTSET_INCLUDED
#define UNORDEREDPOINTSET_INCLUDED

#include "point3d.h"

namespace utw3dface {

class UnorderedPointSet
{
public:
	int allocationstride;
	int npointsallocated;
	int npoints;
	Point3D *point;
	Point3D cog;

	double minX,maxX;
	double minY,maxY;
	double minZ,maxZ;
	double dx,dy;

	UnorderedPointSet(int n=-1);
	~UnorderedPointSet();

	int Reallocate(int n);
	int Free();

	int AddPoint(const Point3D &p);
	Point3D &operator[](int i);
	void CopyFrom(UnorderedPointSet &ups);

	Point3D &CentreOfGravity();

	void SortXYZ();

	int WriteVRML(const char *file);

	int GetStats();
};

}
#endif // UNORDEREDPOINTSET_INCLUDED
