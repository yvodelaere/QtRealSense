////////////////////////////////////////////////////////////////////
//
// This class represents a spline in an arbitrary dimensional space
//
// It uses the CoordList class to represent the list of points
// It uses the spline source of the gnu plotutils library
//
// Author: Luuk Spreeuwers, University of Twente
// Date: 2011
// Last update: 13 July 2011
// version: 1.0
//
////////////////////////////////////////////////////////////////////////

#ifndef SPLINE_INCLUDED
#define SPLINE_INCLUDED

#include <coordlist.h>

namespace utw3dface {

class Spline
{
public:
	CoordList *points;
	CoordList *abscissa;

	// signals if spline is closed, default 0 (=open)
	int closed;
	// tension 0->piecewise cubic, +Inf->polygonal, default 0.0
	float tension;
	// determine boundaries of spline 0->natural spline, default 1.0; not used in closed spline.
	float boundary_condition;

private:
	// these are used by the Fit() method
	// list of points separated after dimension
	double *y; 
	// list of derivatives at points separated after dimension
	double *z;
	// list of abscissa
	double *t;

public:
	// constructor & destructor
	Spline(int dim);
	Spline(CoordList &points,int closed=0); // also calls Fit()
	~Spline();

	int AddPoint(float *coords,float t); // add point at absciss t, must call Fit() afterwards
	int SetAbscissaGeoDist(int npoints); // set abscissa using npoints between each point on spline
	int SetAbscissa(float *t=NULL); // set abscissa, must be monotonously increasing, default from 0..1
	int Fit(); // Fit spline; if necessary, first genereate valid abscissa; required for Interpolate()
	int Interpolate(float t,float *coords); // get point on spline for absciss t
	CoordList *InterpolateRange(float t0,float t1,int n); // interpolate spline between t0 and t1 with n points
	int Normal2D(float t,float dt,float *n); // calculate normal on 2D spline at t: (-dy/dt,dx/dt); length of n is set to 1
	float Curvature2D(float t,float dt); // calculate curvature at t using step dt
};

}

#endif // SPLINE_DEFINED
