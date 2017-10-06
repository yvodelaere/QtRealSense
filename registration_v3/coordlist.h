////////////////////////////////////////////////////////////////////
//
// This class represents a list of coordinates of any dimension
//
// Memory is allocated in chunks in order to make handling efficient
// CoordLists can be read and written to files, copied and dynamically extended/shrunken
//
// Author: Luuk Spreeuwers, University of Twente
// Date: 2011
// Last update: 12 July 2011
// version: 1.0
//
////////////////////////////////////////////////////////////////////////

#ifndef COORDLIST_INCLUDED
#define COORDLIST_INCLUDED

namespace utw3dface {

class CoordList
{
public:
	// dimension of the points in the list
	int dim;
	// number of points in the list
	int npoints;
	// coordinates of the points stored in sequence, first all coordinates of point 0, then all of point 1 etc.
	float *coords;
private:
	int allocated_points;
	int allocation_incr;

public:
	// constructor and destructor
	CoordList(int dim=0); // create an empty CoordList with dimension dim
	CoordList(CoordList &cl); // create a CoordList from another CoordList
	~CoordList(); // destroy a CoordList

	int Empty(); // set npoints to 0, allocated space remains
	int Resize(int npoints); // resize allocated space and number of points to npoints
	int Realloc(int allocatedpoints); // resizes allocated space; if allocatedpoints<npoints, then resize to npoints
	float* Get(int index); // get coords of point at index
	int Set(int index,float *coords); // set coords of point at index
	float* operator[](int index) { return Get(index); }
	int Add(float* p); // at a point at the end of the list
	int Insert(int i,float* p); // insert a point at index i
	int Remove(int i); // remove a point at index i
	int Add(int n,float* p); // add n points at the end of the list
	int Insert(int i,int n,float* p); // insert n points at the end of the list
	int Remove(int i,int n); // remove n points at index i
	operator float*() { return coords; } 
	CoordList& operator=(CoordList &cl); // copy a coordlist
	int Read(char *file); // read a coordlist from file
	int Write(char *file); // write a coordlist to file
};

}

#endif // COORDLIST_INCLUDED
