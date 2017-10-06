#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef linux
#include <values.h>
#include <malloc.h>
#endif
#include <math.h>
#include "compat.h"

#include "unorderedpointset.h"

namespace utw3dface {

int UnorderedPointSet::Reallocate(int n)
{
	Point3D *newpoint=new Point3D[n];
	if (newpoint==NULL)
	{
		fprintf(stderr,"UnorderedPointSet::Reallocate - failed to allocate memory\n");
		exit(1);
	}
	
	if (npoints>0)
	{
		if (npoints<n)
			memcpy(newpoint,point,npoints*sizeof(Point3D));
		else 
		{
			memcpy(newpoint,point,n*sizeof(Point3D));
			npoints=n;
		}
	}

	delete [] point;
	point=newpoint;
	npointsallocated=n;

	return 1;
}

int UnorderedPointSet::Free()
{
	if (point!=NULL)
		delete [] point;
	point=NULL;
	
	return 1;
}

UnorderedPointSet::UnorderedPointSet(int n)
{
	allocationstride=10000;
	point=NULL;
	npoints=0;
	Reallocate((n<0) ? allocationstride : n);
}

UnorderedPointSet::~UnorderedPointSet()
{
	Free();
}

int UnorderedPointSet::AddPoint(const Point3D &p)
{
	if (npoints>=npointsallocated)
		Reallocate(npoints+allocationstride);
	point[npoints++]=p;

	return 1;
}

Point3D &UnorderedPointSet::operator[](int index)
{
	//fprintf(stderr, "index = %d, npoints = %d \n", index, npoints);
	if (index<0 || index>=npoints)
	{
		fprintf(stderr,"UnorderedPointSet::operator[] - index out of range\n");
		fprintf(stderr, "index = %d, npoints = %d \n", index, npoints);
		exit(1);
	}

	return point[index];
}

static int cmp_double(const void *pa,const void *pb)
{
	double a=*((double*)pa);
	double b=*((double*)pb);

	return (a<b) ? 1 : (a>b) ? -1 : 0;
}

Point3D &UnorderedPointSet::CentreOfGravity()
{
	// iterative outlayer removal
	int niter=3;
	int j,n;

	double *dist=new double[npoints];
	double *dists=new double[npoints];

	memset(dist,0,npoints*sizeof(double));

	// determine the c.o.g.
	cog=Point3D(0,0,0);
	for (j=0; j<npoints; j++)
		cog+=point[j];
	cog/=npoints;
//fprintf(stderr,"cog=(%g %g %g)\n",cog.x,cog.y,cog.z);
		
	for (j=0; j<npoints; j++)
		dists[j]=dist[j]=cog.Distance(point[j]);

	qsort(dists,npoints,sizeof(double),cmp_double);

	// discard 10% points furthest away and recalculate
	double maxdist=dists[npoints/10];
	cog=Point3D(0,0,0);
	for (n=0,j=0; j<npoints; j++)
	{
		if (dist[j]<maxdist)
		{
			cog+=point[j];
			n++;
		}
	}
	cog/=n;
//fprintf(stderr,"cog=(%g %g %g)\n",cog.x,cog.y,cog.z);

	delete [] dist;
	delete [] dists;

	return cog;
}

static int compareXYZ(const void *vp,const void *vq)
{
	Point3D *p=(Point3D*)vp;
	Point3D *q=(Point3D*)vq;

	return (p->x>q->x)? 1 : (p->x<q->x) ? -1 : (p->y>q->y) ? 1 : (p->y<q->y) ? -1 : (p->z>q->z) ? 1 : (p->z<q->z) ? -1 : 0;
}

void UnorderedPointSet::SortXYZ()
{
	qsort(point,npoints,sizeof(Point3D),compareXYZ);
}

int UnorderedPointSet::WriteVRML(const char *file)
{
	FILE *fvrml;
	fvrml=fopen(file,"w");
	if (fvrml==NULL)
	{
		fprintf(stderr,"UnorderedPointSet::WriteVRML - cannot open %s for writing\n",file);
		exit(1);
	}

	fprintf(fvrml,"#VRML V2.0 utf8\n");
	fprintf(fvrml,"Transform {\n");
	fprintf(fvrml,"  rotation 1 0 0 3.14159265\n");
	fprintf(fvrml,"    children [\n");
	fprintf(fvrml,"      Shape {\n");
	fprintf(fvrml,"        appearance Appearance {\n");
	fprintf(fvrml,"        }\n");
	fprintf(fvrml,"        geometry PointSet {\n");
	fprintf(fvrml,"          coord Coordinate {\n");
	fprintf(fvrml,"            point [\n");
	fprintf(fvrml,"# SCAN COORIDNATES BEGIN\n");
	
	int i;
	for (i=0; i<npoints; i++)
		fprintf(fvrml,"%g %g %g\n",point[i].x,point[i].y,point[i].z);

	fprintf(fvrml,"# SCAN COORIDNATES END\n");
	fprintf(fvrml,"            ]\n");
	fprintf(fvrml,"          }\n");
	fprintf(fvrml,"          color Color {\n");
	fprintf(fvrml,"            color [\n");
	fprintf(fvrml,"# COLOR VALUES BEGIN\n");

	for (i=0; i<npoints; i++)
        	fprintf(fvrml,"%g %g %g\n",0.5,0.5,0.5);

	fprintf(fvrml,"# COLOR VALUES END\n");
	fprintf(fvrml,"            ]\n");
	fprintf(fvrml,"          }\n");
	fprintf(fvrml,"        }\n");
	fprintf(fvrml,"      }\n");
	fprintf(fvrml,"    ]\n");
	fprintf(fvrml,"}\n");

	fclose(fvrml);

	return 1;
}

int UnorderedPointSet::GetStats()
{
	if (npoints==0)
		return 0;

	minX=maxX=point[0].x;
	minY=maxY=point[0].y;
	minZ=maxZ=point[0].z;

	int i;
	for (i=0; i<npoints; i++)
	{
		Point3D p=point[i];
		if (p.x<minX) minX=p.x; else if (p.x>maxX) maxX=p.x;
		if (p.y<minY) minY=p.y; else if (p.y>maxY) maxY=p.y;
		if (p.z<minZ) minZ=p.z; else if (p.z>maxZ) maxZ=p.z;
	}

	// for calculation of dx and dy we assume a homegeneous distribution of
	// the points in the area minX..MaxX,minY..maxY
	dx=dy=(maxX-minX)*(maxY-minY)/npoints;

	if (dx>50 || dy>50)
	{
		fprintf(stderr,"UnorderedPointSet::GetStats(): too few points or too big face\n");
		exit(1);
	}

	return 1;
}

void UnorderedPointSet::CopyFrom(UnorderedPointSet &ups)
{
	npoints=0;
	int i;
	for (i=0; i<ups.npoints; i++)
		AddPoint(ups[i]);
}
}
