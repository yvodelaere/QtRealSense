#include <stdio.h>

#include "opt.h"
#include "point3d.h"
#include "unorderedpointset.h"
#include "rangeimage.h"

#ifndef NOSE_INCLUDED
#define NOSE_INCLUDED

namespace utw3dface
{

class Nose
{
public:
	double length;
	double height;
	double tiplength;
	double mouthlength;
	double foreheadlength;
	
	double x;
	double y;
	double phi;
	double sinphi,cosphi;

	void SetSize(double length,double height,double tiplength,double mouthlength,double foreheadlength)
	{
		this->length=length;
		this->height=height;
		this->tiplength=tiplength;
		this->mouthlength=mouthlength;
		this->foreheadlength=foreheadlength;
	}

	Nose(double x,double y,double phi): x(x), y(y), phi(phi) 
	{ 
		SetSize(30,30,5,10,20); 
		sinphi=sin(phi);
		cosphi=cos(phi);
	}

	~Nose() {}

	double Start() { return -mouthlength; }
	double End() { return length+foreheadlength; }

	int PartOf(double x)
	{
		return (x>-mouthlength && x<length+foreheadlength);
	}

	double NY(double nx)
	{
		if (nx<0 || nx>length)
			return 0;
		if (nx<tiplength)
			return height;
		return (1-(nx-tiplength)/(length-tiplength))*height;
	}

	double NX(const Point3D &p)
	{
		 return (p.x-x)*cosphi+(p.y-y)*sinphi;
	}

	double NY(const Point3D &p)
	{
		 return (p.y-y)*cosphi-(p.x-x)*sinphi;
	}

	int PartOf(const Point3D &p)
	{
		double nx=NX(p);
		return (nx>-mouthlength && nx<length+foreheadlength);
	}

	double Dist(const Point3D &p)
	{
		double nx=NX(p);
		double ny=NY(p);

		double dy=fabs(ny-NY(nx));
		if (ny>0 && ny<height && fabs(nx)<10)
		{
			double dx=fabs(nx);
			return (dx<dy) ? dx : dy;
		}

		return dy;
	}

	double Dist(UnorderedPointSet &ups)
	{
		double sum=0;
		int i,n=0,nn=0;
		double minnx=length,maxnx=0;
		for (i=0; i<ups.npoints; i++)
		{
			if (!PartOf(ups[i]))
				continue;
			sum+=Dist(ups[i]);
			n++;
			double nx=NX(ups[i]);
			if (nx<minnx)
				minnx=nx;
			if (nx>maxnx)
				maxnx=nx;
			if (nx>=0 && nn<length)
				nn++;
		}

		if (minnx>5 || maxnx<length-5 || nn<5)
			sum=100000;

		return (n>0) ? sum/n : 100000;
	}

	double Corr(UnorderedPointSet &ups)
	{
		double nosemean=height*(tiplength+0.5*(length-tiplength))/(mouthlength+length+foreheadlength);
		double profmean=0;
		int n=0,i;
		for (i=0; i<ups.npoints; i++)
		{
			if (!PartOf(ups[i]))
				continue;
			double ny=NY(ups[i]);
			profmean+=ny;
			n++;
		}
		if (n<1)
			return 0;
		profmean/=n;

		double corr=0,nosevar=0,profvar=0;
		for (i=0; i<ups.npoints; i++)
		{
			if (!PartOf(ups[i]))
				continue;
			double profy=NY(ups[i])-profmean;
			double ny=NY(NX(ups[i]))-nosemean;
			corr+=ny*profy;
			nosevar+=ny*ny;
			profvar+=profy*profy;
		}

		 return corr/sqrt(nosevar*profvar);
	}

	Point3D XY(double nx)
	{
		double ny=NY(nx);
		return Point3D(x+nx*cos(phi)-ny*sin(phi),y+nx*sin(phi)+ny*cos(phi),0);
	}
};

double nosemodel(Opt *opt);

double fitnose_rough(RangeImage &rir,UnorderedPointSet &ups,int minflag,double theta_opt,double phi_opt,double dx_opt,double &gamma_opt,int &nose_x,int &nose_y,double &nose_z,int LR=0,int nearfrontal=0);
double fitnose2profile(UnorderedPointSet &profile,double noselength,double &nx,double &ny,double &nphi,double &nl,double &nh,double &nt,int nosefitmethod);

}

#endif // NOSE_INCLUDED

