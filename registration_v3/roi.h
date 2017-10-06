////////////////////////////////////////////////////////////////////
//
// Class to determine 3D ROI for the face surface (FaceRoi)
// Five helper classes are defined: UpsCell and ListOfUpsCells and Cylinder, CylCell and ListOfCylCells.
// The FaceRoi class contains a ListOfUpsCells and a ListOfCylCells.
// An UpsCell is a cell containing a number of 3D points.
// The approach is as follows:
//	1. The unordered pointset describing the facial surface is divided into a course grid in xy with size cellsize [20 mm]
//	2. For each grid point an unordered point set is created with the points falling within this grid box
//	3. These unordered point lists are added to the ListOfUpsCells (lups), if at least minpointspercell are present [20]
//	4. The UpsCells in the ListOfOpsCells are split into cells with a depth range less than 20 mm
//	5. Again cells with less than minpointspercell are dropped
//	6. For each cell the mean position and normal (using eigenvalue analysis) are determined
//	7. For a good surface patch, the largest eigenvalue represents the normal and should have a large magnitude (quality)
//	8. Cells with poor quality (<minnormalquality [0.3]) are dropped
//	9. Next cylinders are fitted between pairs of UpsCells using their positions and normals to define the cylinder
//	10. Only cylinders are considered defined by Upscells with:
//		* a max distance in y less than cylpoints_maxydist [50 mm]
//		* a max distance in x less than cylpoints_maxxdist [150 mm]
//		* a radius between facemaxr [75 mm] and faceminr 100 [mm]
//		* a max tilt (deviation from vertical) of less than maxtilt [45 deg]
//	11. For each cylinder the support is calculated, i.e. the number of UpsCells that are close enough and have normals 
//		that point in more or less the same direction as the normals of the cylinder
//	12. All cylinders and their support are put into a list of CylCells
//	13. Cylinders that are close to each other are merged
//	14. The cylinder with the best support is selected
//	15. Finally all points in the original unordered points set with a distance larger than maxdisttocylinder [75mm] are dropped
//
// Author: Luuk Spreeuwers, University of Twente
// Date: 2005-2009
// Last update: 31-07-2009
// version: 2.0
//
////////////////////////////////////////////////////////////////////////

#ifndef ROI_INCLUDED
#define ROI_INCLUDED

#include <cv.h>
#include "unorderedpointset.h"
#include "rangeimage.h"

namespace utw3dface {

//////////////////////////////////////////////////////////////////
// helper classes UpsCell and ListOfUpsCells
//////////////////////////////////////////////////////////////////

class UpsCell:public UnorderedPointSet
{
public:
	Point3D mean;
	Point3D normal;
	double quality;

	UpsCell *pnext;
	UpsCell *pprev;

	UpsCell();
	UpsCell(UnorderedPointSet &ups);
	~UpsCell();

	Point3D &GetNormal();
};

class ListOfUpsCells
{
public:
	UpsCell *begin;
	UpsCell *current;

	int ncells;

	ListOfUpsCells();
	UpsCell *Remove();
	~ListOfUpsCells();

	UpsCell *Append(UnorderedPointSet &ups);
	UpsCell *Insert(UnorderedPointSet &ups);

	UpsCell *Current();
	UpsCell *Next();
	UpsCell *Prev();
	UpsCell *Begin();
};

class Cylinder
{
public:
	Point3D c;	// a point on the axis
	Point3D a;	// direction of the axis
	double r;	// radius of the cylinder

	double support;

	double maxdist;
	double maxangle;
	double faceheight;
	double faceminr,facemaxr;

	void Init();
	Cylinder();
	Cylinder(const Cylinder &cyl);
	Cylinder(const Cylinder *pcyl);
	Cylinder(const Point3D &p1,const Point3D &p2,const Point3D &n1,const Point3D &n2);
	Cylinder &operator=(const Cylinder &cyl)
	{
		c=cyl.c;
		a=cyl.a;
		r=cyl.r;
		support=cyl.support;
		maxdist=cyl.maxdist;
		maxangle=cyl.maxangle;
		faceheight=cyl.faceheight;
		faceminr=cyl.faceminr;
		facemaxr=cyl.facemaxr;
	}

	double Distance(const Point3D &p);
	double AngleDiff(const Point3D &p,const Point3D &n);
	int Support(const Point3D &p,const Point3D &n);
	double Support(UnorderedPointSet &ups,UnorderedPointSet &upsnormals);
	double Support(UnorderedPointSet &ups,double &dist);
};

class CylCell:public Cylinder
{
public:
	CylCell *pnext;
	CylCell *pprev;

	CylCell(const Point3D &p1,const Point3D &p2,const Point3D &n1,const Point3D &n2);
	CylCell(const Cylinder &cyl);
	CylCell(const Cylinder *pcyl);
	~CylCell();
};

class ListOfCylCells
{
public:
	CylCell *begin;
	CylCell *current;

	int ncells;

	ListOfCylCells();
	CylCell *Remove();

	~ListOfCylCells();

	CylCell *Append(Cylinder &cyl);
	CylCell *Insert(Cylinder &cyl);

	CylCell *Current();
	CylCell *Next();
	CylCell *Prev();
	CylCell *Begin();
};

class FaceRoi
{
public:
	ListOfUpsCells lups;
	ListOfCylCells lcyl;
	UnorderedPointSet rups,rupsn;	// reduced ups and normals
	int w,h;
	double cellsize;
	int minpointspercell;
	double maxtilt;
	double minnormalquality;
	double cylpoints_maxydist;
	double cylpoints_maxxdist;
	double maxdisttocylinder;
	UnorderedPointSet upsroi;
	Cylinder cylinder;

	void Init();

	void GetRupsAndNormals(UnorderedPointSet &ups);	// get reduced ups and normals

	void DumpRups(Cylinder &cyl,const char *filename,Point3D &p,int onlysupported=0);

	FaceRoi();

	FaceRoi(UnorderedPointSet &ups);

	~FaceRoi();
};

}

#endif // ROI_INCLUDED
