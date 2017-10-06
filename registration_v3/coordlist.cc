#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef linux
#include <values.h>
#include <malloc.h>
#endif
#include <math.h>
#include "compat.h"

#include "coordlist.h"

namespace utw3dface {

CoordList::CoordList(int dim)
{
	if (dim<0)
	{
		fprintf(stderr,"CoordList: attempt to create CoordList with dim<0\n");
		exit(1);
	}

	coords=NULL;
	allocated_points=0;
	allocation_incr=10;
	npoints=0;
	this->dim=dim;
}

CoordList::CoordList(CoordList &cl)
{
	coords=NULL;
	allocated_points=0;
	allocation_incr=10;
	npoints=0;
	this->dim=cl.dim;
	Add(cl.npoints,cl);
}

CoordList::~CoordList()
{
	if (coords != NULL)
		free(coords);
}

int CoordList::Realloc(int npoints)
{
	float *coords;
	
	if (npoints<this->npoints)
		npoints=this->npoints;

	if (npoints==0)
		coords=NULL;
	else
	{
		coords=(float*)malloc(sizeof(float)*dim*npoints);
		if (coords==NULL)
			return 0;
	}

	if (this->coords != NULL)
	{
		memcpy(coords,this->coords,sizeof(float)*dim*this->npoints);
		free(this->coords);
	}

	this->coords=coords;
	this->allocated_points=npoints;

	return 1;
}

int CoordList::Resize(int npoints)
{
	if (npoints<0)
		return 0;
	if (npoints>this->allocated_points)
		Realloc(npoints);

	this->npoints=npoints;

	return 1;
}

float* CoordList::Get(int index)
{
	if (index<0 || index>=this->npoints)
		return NULL;

	return this->coords+index*this->dim;
}

int CoordList::Set(int index,float* coords)
{
	if (index<0 || index>=this->npoints)
		return 0;

	memcpy(this->coords+index*this->dim,coords,this->dim*sizeof(float));

	return 1;
}

int CoordList::Add(float* coords)
{
	// note: we need a copy if we copy from a point in this list!
	float *copy;

	copy=(float*)malloc(sizeof(float)*this->dim);
	if (copy==NULL)
		return 0;
	memcpy(copy,coords,this->dim*sizeof(float));

	if (this->allocated_points <= this->npoints)
		if (!Realloc(this->npoints+this->allocation_incr))
		{
			free(copy);
			return 0;
		}
	
	memcpy(this->coords+this->npoints*this->dim,copy,this->dim*sizeof(float));

	this->npoints++;
	free(copy);
		
	return 1;
}

int CoordList::Insert(int index,float* coords)
{
	// note: we need the copy if we copy from a point in this list!
	float *copy;

	if (index<0 || index>this->npoints)
		return 0;

	copy=(float*)malloc(sizeof(float)*this->dim);
	if (copy==NULL)
		return 0;
	memcpy(copy,coords,this->dim*sizeof(float));

	if (this->allocated_points <= this->npoints)
		if (!Realloc(this->npoints+this->allocation_incr))
		{
			free(copy);
			return 0;
		}
	
	if (index<this->npoints)
		memmove(this->coords+(index+1)*this->dim,this->coords+index*this->dim,(this->npoints-index)*this->dim*sizeof(float));

	memcpy(this->coords+index*this->dim,copy,this->dim*sizeof(float));

	this->npoints++;
	free(copy);
		
	return 1;
}
	
int CoordList::Remove(int index)
{
	if (index<0 || index>=this->npoints)
		return 0;

	memmove(this->coords+index*this->dim,this->coords+(index+1)*this->dim,(this->npoints-index)*this->dim*sizeof(float));

	this->npoints--;
		
	return 1;
}

int CoordList::Add(int n,float* coords)
{
	if (this->allocated_points < this->npoints+n || n<1)
		if (!Realloc(this->npoints+n+this->allocation_incr))
			return 0;
	
	memcpy(this->coords+this->npoints*this->dim,coords,this->dim*sizeof(float)*n);

	this->npoints+=n;
		
	return 1;
}


int CoordList::Insert(int index,int n,float* coords)
{
	if (index<0 || index>this->npoints || n<1)
		return 0;

	if (this->allocated_points < this->npoints+n)
		if (!Realloc(this->npoints+n+this->allocation_incr))
			return 0;
	
	if (index<this->npoints)
		memmove(this->coords+(index+n)*this->dim,this->coords+index*this->dim,(this->npoints-index)*this->dim*sizeof(float));

	memcpy(this->coords+index*this->dim,coords,n*this->dim*sizeof(float));

	this->npoints+=n;
		
	return 1;
}

int CoordList::Remove(int index,int n)
{
	if (index<0 || index+n>this->npoints || n<1)
		return 0;

	if (index+n<this->npoints)
		memmove(this->coords+index*this->dim,this->coords+(index+n)*this->dim,(this->npoints-index-n)*this->dim*sizeof(float));

	this->npoints-=n;
		
	return 1;
}

int CoordList::Empty()
{
	this->npoints=0;

	return 1;
}

#define ISICOORDLISTID	"# isiCoordList"

int CoordList::Read(char *file)
{
	FILE *f;

	f=fopen(file,"r");
	if (f==NULL)
		return 0;

	int dim=0,npoints=0;
	int d,i;
	char s[1024];
	float *val;

	while (fgets(s,1023,f) != NULL && s[0] == '#' || s[0]=='\n');
	if (s == NULL)
	{
		fclose(f);
		return 0;
	}

	if (sscanf(s,"%d%d",&npoints,&dim) != 2) // in case of empty line
		fscanf(f,"%d%d",&npoints,&dim);

	if (npoints>0 && dim>0)
	{
		Empty();
		this->dim=dim;
		Resize(npoints);
	}

	val=(float*)malloc(sizeof(float)*dim);
	if (val==NULL)
	{
		fclose(f);
		return 0;
	}
	for (i=0; i<npoints; i++)
	{
		for (d=0; d<dim; d++)
			fscanf(f,"%f",&val[d]);
		Set(i,val);
	}

	free(val);

	fclose(f);
}

int CoordList::Write(char *file)
{
	FILE *f;

	f=fopen(file,"w");
	if (f==NULL)
		return 0;

	int i,d;

	fprintf(f,"%s\n",ISICOORDLISTID);
	fprintf(f,"# first line contains # points and dimension\n");
	fprintf(f,"%d %d\n",this->npoints,this->dim);
	// insert 2 empty lines for gnuplot, the data can be plotted with index 1
	fprintf(f,"\n\n");
	for (i=0; i<this->npoints; i++)
	{
		for (d=0; d<this->dim; d++)
			fprintf(f,"%g ",Get(i)[d]);
		fprintf(f,"\n");
	}

	fclose(f);
}

CoordList& CoordList::operator=(CoordList &cl)
{
	Empty();
	this->dim=cl.dim;
	Realloc(cl.npoints);
	Add(cl.npoints,cl);
}

} 

