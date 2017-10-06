////////////////////////////////////////////////////////////////////
//
// Advanced class for representation of range images and construction of 
//range images from 3D point clouds
//
// A range image is defined as a plane in 3D (defined by directions u,v and origin o) 
// onto which points are projected. Attached to the plane is a grid defined by 
// gridspacing du,dv and origin ou,ov and width,height. A range image is created 
// from a point cloud by projecting the points onto the grid and accumulating the 
// points in each grid point like for a 2D histogram. Each grid point is a bin 
// of the histogram. The average depth (distance of point to the plane of the range image) 
// is set to the pixel value, the flag contains the number of points accumulated and 
// the stddev is the standard deviation of depth values.
//
// In addition there is some code to determine if all points projected to a grid 
// point lay on the same surface or if there might be points from the front and 
// the back of a 3D object. The approach followed using back and backflag probably 
// needs improvement.
//
// There is also a method for studying the effect of compression/decompression
// on the range image: CompressDecompress().
//
// Author: Luuk Spreeuwers, University of Twente
// Date: 2005-2009
// Last update: 31-07-2009
// version: 2.0
//
////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifdef linux
#include <malloc.h>
#endif
#include "compat.h"

#ifndef RANGEIMAGE_DEFINED
#define RANGEIMAGE_DEFINED

#include "unorderedpointset.h"

namespace utw3dface {

class RangeImage
{
public:
	int width;
	int height;
	double *pixel;
	double *stddev;
	double *flag;
	double *front;
	double *back;
	double *frontflag;
	double *backflag;

	unsigned char *jpeg_buf;
	int jpeg_buf_size;

	Point3D o;
	Point3D u;
	Point3D v;

	double du;	// size of pixel grid in u and v directions
	double dv;
	double ou;	// offset relative to origin of grid in pixels
	double ov;

	// transformation matrix to transform 3D points to points on the plane defined by u,v,o and on the grid defined by du,dv,ou,ov
	double T[3][4];	

	UnorderedPointSet *ups;	// optional poiter to corresponding point set

	//////////////////////////////////////////////////////////////////////
	// helper methods for construction and destruction
	/////////////////////////////////////////////////////////////////////

	void calculateT();	// calculate transformatio matrix T
	void Clear(); 		// set all allocated memory to zeros
	void Free();		// free all dynamically allocated space

	// Initialise range image: 
	// 	3D plane defined by u,v,o and grid du,dv,ou,ov and size width,height. 
	// 	Allocate and initialise all required memory 
	// 	Calculate transformation matrix for for projection of 3D points to range image
	void Init(const Point3D &o,const Point3D &u,const Point3D &v,double du,double dv,double ou,double ov,int width,int height);

	/////////////////////////////////////////////////////////////////////
	// Constructors
	/////////////////////////////////////////////////////////////////////

	// Create a range image and call Init
	RangeImage(const Point3D &o,const Point3D &u,const Point3D &v,double du,double dv,double ou,double ov,int width,int height);

	RangeImage();	// create an empty range image to be initialised later

	// Create and initialise a range image with:
	//	size defined by width, height
	// 	origin (0,0,0)
	//	u coinsiding with the x-axis and v coinsiding with the inverted y-axis (v=-y)
	//	grid spacing du=dv=1
	// 	grid origin ou=ov=0
	RangeImage(int width,int height);

	///////////////////////////////////////////////////////////////////
	// Destructor
	///////////////////////////////////////////////////////////////////

	~RangeImage();	// release all resources

	///////////////////////////////////////////////////////////////////
	// Other methods
	///////////////////////////////////////////////////////////////////

	void CopyDataFrom(RangeImage &ri);	// copy data from other range image, if necessary, the size of the grid is adjusted

	// projects point p on the range image plane and returns u and v coordinates 
	// and distance to range image (x,y,z resp. of returned point)
	Point3D TransformPoint(const Point3D &p);

	// set depth value in range image corresponding to point p; the flag field is set to 255 (high value)
	void SetDepth(const Point3D &p);

	// projects point p on the range image, sets/adjusts the distance of p to the range image plane 
	// and increments the count of the corresponding pixel in the range image (flag field)
	// i.e. the range image is actually a histogram
	// there is a quick fix for determining front and back of the face (only front should be counted) using the backflag field
	// the spread of the entries is recorded in the stddev field
	void AccumulateDepth(const Point3D &p);

	// Determines the range image for a complete unordered point set as present in the field ups
	// it calls AccumulateDepth() for each individual point and handles swapping back and front if necessary
	// the minflag parameter is the minimum number of points in a 'bin' of the grid of the range image to
	// accept it as reliable
	void AccumulateDepth(double minflag=0);

	// Determines the range image for a complete unordered point set ups
	// it calls AccumulateDepth(double minflag=0)
	void AccumulateDepth(UnorderedPointSet &ups,double minflag=0);

	// simple way to set depth in range image by calling SetDepth() for each point of the unordered point set in the ups field
	void SetDepth();

	// simple way to set depth in range image by calling SetDepth() for each point of the unordered point set ups 
	void SetDepth(UnorderedPointSet &ups);

	double GetDepth(double dx,double dy,double &flag);	// get depth and flag (# entries in bin) at position in range image
	double GetStddev(double dx,double dy);			// get stddev of depth at position in range image
	double GetDepth(double dx,double dy);			// get depth at position in range image

	//////////////////////////////////////////////////////////////////////
	// file I/O
	//////////////////////////////////////////////////////////////////////
	int WritePGM(const char *filename,int stretch=0);	// write range image as pgm file (only depth)
	int ReadPNM(const char *filename);			// read range image from pnm file (only depth)
	int WriteSFI(const char *filename);			// write range image as sfi file (only depth)
	int ReadSFI(const char *fname);				// read range image from sfi file (only depth)

	int WriteFlagPGM(const char *filename);			// write flags as pgm file (mainly for debugging)
	int WriteStddevPGM(const char *filename,int stretch=0);	// write std dev as pgm file (mainly for debugging)

	////////////////////////////////////////////////////////////////////
	// quick access to elements without range checking etc.
	////////////////////////////////////////////////////////////////////
	double *operator[](int y)
	{
		return pixel+y*width;
	}

	double &Pixel(int x,int y)
	{
		return pixel[y*width+x];
	}

	double &Stddev(int x,int y)
	{
		return stddev[y*width+x];
	}

	double &Flag(int x,int y)
	{
		return flag[y*width+x];
	}

	void CompressDecompress(int quality=90);	// special function to study effect of jpeg compression/decompression
};

}
#endif // RANGEIMAGE_DEFINED
