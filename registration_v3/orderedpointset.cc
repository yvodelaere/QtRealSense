#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#ifdef linux
#include <malloc.h>
#include <values.h>
#endif
#include <math.h>
#include <float.h>
#include "compat.h"

#include "orderedpointset.h"

namespace utw3dface {

extern int verbose;

OrderedPointSet::OrderedPointSet()
{
	npolygons=nvalidpoints=width=height=0;
	X=NULL;
	Y=NULL;
	Z=NULL;
	flags=NULL;
	polygon=NULL;
}

int OrderedPointSet::Free()
{
	width=height=nvalidpoints=0;
	if (X!=NULL)
		delete [] X;
	if (Y!=NULL)
		delete [] Y;
	if (Z!=NULL)
		delete [] Z;
	if (flags!=NULL)
		delete [] flags;
	if (polygon!=NULL)
		delete [] polygon;

	return 1;
}

OrderedPointSet::~OrderedPointSet()
{
	Free();
}

int OrderedPointSet::Allocate(int width,int height)
{
	this->width=width;
	this->height=height;
	X=new double [width*height];
	Y=new double [width*height];
	Z=new double [width*height];
	flags=new int [width*height];

	memset(flags,0,width*height*sizeof(int));

	if (X==NULL || Y==NULL || Z==NULL || flags==NULL)
	{
		fprintf(stderr,"OrderedPointSet::Allocate - failed to allocate memory\n");
		exit(1);
	}

	return 1;
}

int comparepoints(const void *p1,const void *p2)
{
	Point3D *point1=(Point3D*)p1;
	Point3D *point2=(Point3D*)p2;

	if (point1->y>point2->y)
		return 1;
	if (point1->y<point2->y)
		return -1;
	if (point1->x>point2->x)
		return 1;
	if (point1->x<point2->x)
		return -1;
	
	return 0;
}

int OrderedPointSet::Read(const char *file)
{
	if (this->ReadAbs(file))
		return 1;
	if (this->ReadVRML(file))
		return 1;
	if (this->Read3DPoints(file))
		return 1;

	return 0;
}

// return double in d and/or next non white space
int readdoublefromfile(FILE *f,double &d)
{
	char s[256];
	int c,i;
	for (;;)
	{
		c=getc(f);
		if (isspace(c) || c==',')
			continue;
		if (c==EOF)
			return -1;
		if (c=='#') // skip to eol
		{
			for (;;)
			{
				c=fgetc(f);
				if (c==EOF)
					return -1;
				if (c=='\r' || c=='\n')
					break;
			}
		}
		if (isdigit(c) || c=='.' || c=='-')
			break;
		else
			return c;
	}

	int dots=0;
	if (c=='.')
		dots++;
	int exps=0;
	s[0]=c;
	for (i=1; i<256; i++)
	{
		c=getc(f);
		if (isdigit(c) || c=='.' || c=='-' || c=='+' || c=='e' || c=='E')
			s[i]=c;
		else
			break;
		if (c=='.')
			dots++;
		if (c=='e' || c=='E')
			exps++;
	}
	s[i]='\0';

	if (dots>1 || exps>1)
	{
		fprintf(stderr,"Error reading floating point number: multiple dots or exp\n");
		exit(1);
	}
	d=atof(s);

	return c;
}

// very limited point extractor from VRML
int OrderedPointSet::ReadVRML(const char *file)
{
	Point3D *p=new Point3D[1000020]; // max 1 million points!
	FILE *f;
	f=fopen(file,"r");
	if (f==NULL)
	{
		fprintf(stderr,"OrderedPointSet::ReadVRML - failed to open %s for input\n",file);
		exit(1);
	}

	// read header
	char s[256];
	fscanf(f,"%s",s);
	if (strcmp(s,"#VRML")!=0)
	{
		// not a VRML file
		fclose(f);
		return 0;
	}

	// scan for "point" and [
	for (;;)
	{
		int n=fscanf(f,"%s",s);
		if (n==EOF)
		{
			fprintf(stderr,"OrderedPointSet::ReadVRML - failed to read %s\n",file);
			exit(1);
		}
		if (strcmp(s,"point[")==0)
			break;
		if (strcmp(s,"point")==0)
		{
			fscanf(f,"%s",s);
			if (strcmp(s,"[")==0)
				break;
		}
	}
		

	int i;
	double zmin=0, zmax=0;

	for (i=0;;i++)
	{
		double x,y,z;

		if (i>=1000000)
		{
			fprintf(stderr,"too many points in OrderedPointSet::ReadVRML!\n");
			exit(1);
		}

		int c;
		c=readdoublefromfile(f,x);
		if (c==']')
			break;
		c=readdoublefromfile(f,y);
		c=readdoublefromfile(f,z);

		p[i]=Point3D(x,y,z);
		if (i==0)
			zmax=zmin=z;
		else 
		{
			if (z<zmin)
				zmin=z;
			else if (z>zmax)
				zmax=z;
		}       
	}

	int npoints=i;

/* code is too shaky at the moment
	// check if there are polygons given as well
	// scan for "coordIndex" and [
	int haspolygons=0;
	for (;;)
	{
		int n=fscanf(f,"%s",s);
		if (n==EOF)
			break;
		if (strcmp(s,"coordIndex[")==0)
		{
			haspolygons=1;
			break;
		}
		if (strcmp(s,"coordIndex")==0)
		{
			fscanf(f,"%s",s);
			if (strcmp(s,"[")==0)
			{
				haspolygons=1;
				break;
			}
		}
	}
	// skip rest of line
	fgets(s,255,f);

	if (haspolygons)
	{
		npolygons=0;
		Polygon *pol=new Polygon[100000];

		int stop=0;
		for (i=0;;i++)
		{
			double x,y,z;
	
			if (i>=100000)
			{
				fprintf(stderr,"too many polygons in OrderedPointSet::ReadVRML!\n");
				exit(1);
			}
			if (fgets(s,255,f)==NULL)
			{
				fprintf(stderr,"failed reading polygons in OrderedPointSet::ReadVRML!\n");
			}
			int j;
			for (j=0; j<strlen(s); j++)
			{
				if (s[j]==',')
					s[j]=' ';
				if (s[j]==']')
					stop=1;
			}
			if (stop)
				break;
			int p1,p2,p3,p4,p5;
			int k=sscanf(s,"%d%d%d%d",&p1,&p2,&p3,&p4);
			if (k<4)
			{
				fprintf(stderr,"failed to read polygons!\n");
				exit(1);
			}
			if (p4!=-1)
			{
				k=sscanf(s,"%d",&p5);
			}

			if (p4!=-1 && p5!=-1)
			{
				fprintf(stderr,"only supports 3/4 point polygons!\n");
				exit(1);
			}
			pol[i]=Polygon(p1,p2,p3,p4);
		}

		npolygons=i;
		polygon=new Polygon[npolygons];

		for (i=0; i<npolygons; i++)
			polygon[i]=pol[i];
	
		if (verbose)
			fprintf(stderr,"OrderedPointSet::ReadVRML: read %d polygons\n",npolygons);

		delete pol;
					
	}
*/

	fclose(f);

/* this is a method to generate more, interpolated points
	for (i=0; i<npolygons; i++)
	{
		Point3D c=polygon[i].Centre(p);
		p[npoints++]=c;
	}
*/
	qsort(p,npoints,sizeof(Point3D),comparepoints);

	width=20;
	height=npoints/width+1;
	if (height<20)
	{
		fprintf(stderr,"too few points in OrderedPointSet::Read3DPoints()\n");
		exit(1);
	}
	Allocate(width,height);

	// for padding last line of data
	for (i=npoints; i<npoints+20; i++)
		p[i]=p[0];

	// read data
	int x,y;
	for (y=0; y<height; y++)
	for (x=0; x<width; x++)
	{
		int i=y*width+x;
		flags[i]=(i<npoints) ? 1 : 0;
		X[i]=p[i].x;
		Y[i]=p[i].y;
		Z[i]=p[i].z;
	}

	delete [] p;

	return 1;
}

int OrderedPointSet::Read3DPoints(const char *file)
{
	Point3D *p=new Point3D[1000020]; // max 1 million points!
	FILE *f;
	f=fopen(file,"r");
	if (f==NULL)
	{
		fprintf(stderr,"OrderedPointSet::Read3DPoints - failed to open %s for input\n",file);
		exit(1);
	}

	int i;
	char s[1024];
	double zmin=0, zmax=0;

	for (i=0;;i++)
	{
		double x,y,z;

		if (i>=1000000)
		{
			fprintf(stderr,"too many points in OrderedPointSet::Read3DPoints!\n");
			exit(1);
		}
		if (fgets(s,1023,f)==NULL)
			break;
		if (sscanf(s,"%lf%lf%lf",&x,&y,&z) != 3)
		{
			fprintf(stderr,"failed to read coordinates!\n");
			exit(1);
		}
		p[i]=Point3D(x,y,z);
		if (i==0)
			zmax=zmin=z;
		else 
		{
			if (z<zmin)
				zmin=z;
			else if (z>zmax)
				zmax=z;
		}       
	}

	fclose(f);

	int npoints=i;

	qsort(p,npoints,sizeof(Point3D),comparepoints);

	width=20;
	height=npoints/width+1;
	if (height<20)
	{
		fprintf(stderr,"too few points in OrderedPointSet::Read3DPoints()\n");
		exit(1);
	}
	Allocate(width,height);

	// for padding last line of data
	for (i=npoints; i<npoints+20; i++)
		p[i]=p[0];

	// read data
	int x,y;
	for (y=0; y<height; y++)
	for (x=0; x<width; x++)
	{
		int i=y*width+x;
		flags[i]=(i<npoints) ? 1 : 0;
		X[i]=p[i].x;
		Y[i]=p[i].y;
		Z[i]=p[i].z;
	}

	delete [] p;

	return 1;
}

int OrderedPointSet::ReadAbs(const char *file)
{
	FILE *f;
	f=fopen(file,"r");
	if (f==NULL)
	{
		fprintf(stderr,"OrderedPointSet::ReadAbs - failed to open %s for input\n",file);
		exit(1);
	}

	// read header
	char s[256],magic1[256],magic2[256],magic3[256];
	fgets(s,255,f);
	sscanf(s,"%d%s",&height,magic1);
	fgets(s,255,f);
	sscanf(s,"%d%s",&width,magic2);
	fgets(s,255,f);
	sscanf(s,"%s",magic3);

	if (strcasecmp(magic1,"rows")!=0 || strcasecmp(magic2,"columns")!=0 || strcasecmp(magic3,"pixels")!=0)
	{
		// this is not an abs file; try reading a pts file
		fclose(f);
		return 0;
	}

	if (width<1 || width>100000 || height<1 || height>100000)
	{
		fprintf(stderr,"OrderedPointSet::ReadAbs - dimensions out of bounds!\n");
		exit(1);
	}

	Allocate(width,height);

	// read data
	int x,y;
	for (y=0; y<height; y++)
	for (x=0; x<width; x++)
	{
		int i=y*width+x;
		if (fscanf(f,"%d",&flags[i])!= 1)
		{
			fprintf(stderr,"OrderedPointSet::ReadAbs - file truncated!\n");
			exit(1);
		}
	}
	for (y=0; y<height; y++)
	for (x=0; x<width; x++)
	{
		int i=y*width+x;
		if (fscanf(f,"%lf",&X[i]) != 1)
		{
			fprintf(stderr,"OrderedPointSet::ReadAbs - file truncated!\n");
			exit(1);
		}
	}
	for (y=0; y<height; y++)
	for (x=0; x<width; x++)
	{
		int i=y*width+x;
		if (fscanf(f,"%lf",&Y[i]) != 1)
		{
			fprintf(stderr,"OrderedPointSet::ReadAbs - file truncated!\n");
			exit(1);
		}
	}
	for (y=0; y<height; y++)
	for (x=0; x<width; x++)
	{
		int i=y*width+x;
		if (fscanf(f,"%lf",&Z[i]) != 1)
		{
			fprintf(stderr,"OrderedPointSet::ReadAbs - file truncated!\n");
			exit(1);
		}
	}

	fclose(f);

	return 1;
}

int OrderedPointSet::GetStats()
{
	if (flags==NULL || X==NULL || Y==NULL || Z==NULL)
		return 0;

	nvalidpoints=0;
	minX=minY=minZ=MAXFLOAT;
	maxX=maxY=maxZ=-MAXFLOAT;

	dx=dy=0;
	int ndx=0,ndy=0;

	int x,y;
	for (y=0; y<height; y++)
	for (x=0; x<width; x++)
	{
		int i=y*width+x;
		if (flags[i]==1)
		{
			nvalidpoints++;
			if (X[i]<minX) minX=X[i]; else if (X[i]>maxX) maxX=X[i];
			if (Y[i]<minY) minY=Y[i]; else if (Y[i]>maxY) maxY=Y[i];
			if (Z[i]<minZ) minZ=Z[i]; else if (Z[i]>maxZ) maxZ=Z[i];
			if (i>0 && flags[i-1]==1)
			{
				dx+=fabs(X[i]-X[i-1]);
				ndx++;
			}
			if (i-width>0 && flags[i-width]==1)
			{
				dy+=fabs(Y[i]-Y[i-width]);
				ndy++;
			}
		}	
	}
	if (ndx>0)
		dx/=ndx;
	if (ndy>0)
		dy/=ndy;

	return 1;
}

int OrderedPointSet::Write3DPoints(const char *file)
{
	FILE *f;
	f=fopen(file,"w");
	if (f==NULL)
	{
		fprintf(stderr,"OrderedPointSet::Write3DPoints - cannot open %s for writing\n",file);
		exit(1);
	}

	int x,y;
	for (y=0; y<height; y++)
	for (x=0; x<width; x++)
	{
		int i=y*width+x;
		if (flags[i]==1)
			fprintf(f,"%g %g %g\n",X[i],Y[i],Z[i]);
	}

	fclose(f);

	return 1;
}

int OrderedPointSet::WriteVRML(const char *file)
{
	FILE *fvrml;
	fvrml=fopen(file,"w");
	if (fvrml==NULL)
	{
		fprintf(stderr,"OrderedPointSet::WriteVRML - cannot open %s for writing\n",file);
		exit(1);
	}

	fprintf(fvrml,"#VRML V2.0 utf8\n");
	fprintf(fvrml,"Transform {\n");
	fprintf(fvrml,"  rotation 1 0 0 3.14159265\n");
	fprintf(fvrml,"    children [\n");
	fprintf(fvrml,"      Shape {\n");
	fprintf(fvrml,"	appearance Appearance {\n");
	fprintf(fvrml,"	}\n");
	fprintf(fvrml,"	geometry PointSet {\n");
	fprintf(fvrml,"	  coord Coordinate {\n");
	fprintf(fvrml,"	    point [\n");
	fprintf(fvrml,"# SCAN COORIDNATES BEGIN\n");
	
	int x,y;
	for (y=0; y<height; y++)
	for (x=0; x<width; x++)
	{
		int i=y*width+x;
		if (flags[i]==1)
			fprintf(fvrml,"%g %g %g\n",X[i],Y[i],Z[i]);
	}

	fprintf(fvrml,"# SCAN COORIDNATES END\n");
	fprintf(fvrml,"	    ]\n");
	fprintf(fvrml,"	  }\n");
	fprintf(fvrml,"	  color Color {\n");
	fprintf(fvrml,"	    color [\n");
	fprintf(fvrml,"# COLOR VALUES BEGIN\n");

	for (y=0; y<height; y++)
	for (x=0; x<width; x++)
	{
		int i=y*width+x;
		if (flags[i]==1)
			fprintf(fvrml,"%g %g %g\n",0.5,0.5,0.5);
	}

	fprintf(fvrml,"# COLOR VALUES END\n");
	fprintf(fvrml,"	    ]\n");
	fprintf(fvrml,"	  }\n");
	fprintf(fvrml,"	}\n");
	fprintf(fvrml,"      }\n");
/*
	fprintf(fvrml,"      Shape {\n");
	fprintf(fvrml,"	appearance Appearance {\n");
	fprintf(fvrml,"	}\n");
	fprintf(fvrml,"	geometry PointSet {\n");
	fprintf(fvrml,"	  coord Coordinate {\n");
	fprintf(fvrml,"	    point [\n");
	fprintf(fvrml,"# SCAN COORIDNATES BEGIN\n");
	
	double fx,fy,fz;
	for (fy=minY; fy<maxY; fy+=2)
	for (fz=minZ; fz<maxZ; fz+=2)
		fprintf(fvrml,"%g %g %g\n",vplane.x(fy,fz),fy,fz);

	for (fx=minX; fx<maxX; fx+=2)
	for (fz=minZ; fz<maxZ; fz+=2)
		fprintf(fvrml,"%g %g %g\n",fx,hplane.y(fy,fz),fz);

	fprintf(fvrml,"# SCAN COORIDNATES END\n");
	fprintf(fvrml,"	    ]\n");
	fprintf(fvrml,"	  }\n");
	fprintf(fvrml,"	  color Color {\n");
	fprintf(fvrml,"	    color [\n");
	fprintf(fvrml,"# COLOR VALUES BEGIN\n");

	for (fy=minY; fy<maxY; fy+=2)
	for (fz=minZ; fz<maxZ; fz+=2)
		fprintf(fvrml,"%g %g %g\n",0.5,0.5,0.5);

	for (fx=minX; fx<maxX; fx+=2)
	for (fz=minZ; fz<maxZ; fz+=2)
		fprintf(fvrml,"%g %g %g\n",0.5,0.5,0.5);

	fprintf(fvrml,"# COLOR VALUES END\n");
	fprintf(fvrml,"	    ]\n");
	fprintf(fvrml,"	  }\n");
	fprintf(fvrml,"	}\n");
	fprintf(fvrml,"      }\n");
*/
	fprintf(fvrml,"    ]\n");
	fprintf(fvrml,"}\n");

	fclose(fvrml);

	return 1;
}

int OrderedPointSet::WriteVRMLsurface(const char *filename)
{
	FILE *fvrml=fopen(filename,"w");
	if (fvrml==NULL)
	{
		fprintf(stderr,"failed to open %s for writing\n",filename);
		exit(1);
	}

	fprintf(fvrml,"#VRML V2.0 utf8\n");
	fprintf(fvrml,"Transform {\n");
	fprintf(fvrml,"  rotation 1 0 0 3.14159265\n");
	fprintf(fvrml,"    children [\n");
	fprintf(fvrml,"Shape {\n\tappearance Appearance {\n\tmaterial Material {\n\t\tdiffuseColor 0.8 0.8 0.8\n\t\tambientIntensity 0.2\n\t\temissiveColor 0.0 0.0 0.0\n\t\tspecularColor 0.0 0.0 0.0\n\t\tshininess 0.2\n\t\ttransparency 0.0\n\t\t}\n\t}geometry IndexedFaceSet {\n\tcoord Coordinate {\n\t\t point [\n"); 
	int x,y;
	for (y=0; y<height; y++)
	for (x=0; x<width; x++)
	{
		int i=y*width+x;
		fprintf(fvrml,"%g %g %g\n",X[i],Y[i],Z[i]);
	}

	fprintf(fvrml,"\t\t\t]\n\t\t}\n\t\tccw TRUE\n\t\tconvex TRUE\n\t\tcreaseAngle 0\n\t\tsolid TRUE\n\t\tcoordIndex [\n");

	if (polygon!=NULL)
	{
		if (verbose)
			fprintf(stderr,"writing original polygons\n");
		int i,j;
		for (i=0; i<npolygons; i++)
		{
			for (j=0; j<4; j++)
			{
				if (polygon[i][j]>0)
					fprintf(fvrml,"%d ",polygon[i][j]);	
			}
			fprintf(fvrml,"-1\n");
		}
	}
	else
	{
		for (y=0; y<height; y++)
		for (x=0; x<width; x++)
		{
			if ( flags[y*width+x]==1 && flags[y*width+x+1]==1 && flags[(y+1)*width+x]==1 && flags[(y+1)*width+x+1]==1)
			{
				fprintf(fvrml,"%d %d %d -1\n",y*width+x,y*width+x+1,(y+1)*width+x);
				fprintf(fvrml,"%d %d %d -1\n",(y+1)*width+x,y*width+x+1,(y+1)*width+x+1);
			}
		}
	}

	fprintf(fvrml,"\t\t\t]\n\t\tcolorIndex [ ]\n\t}\n}");
	fprintf(fvrml,"    ]\n");
	fprintf(fvrml,"}\n");

	fclose(fvrml);

	return 1;
}

int OrderedPointSet::GetPoint(int x,int y,double *point)
{
	if (x<0 || x>=width || y<0 || y>=height)
		return 0;
	int offset=y*width+x;
	if (flags[offset]==0)
		return 0;

	point[0]=X[offset];
	point[1]=Y[offset];
	point[2]=Z[offset];

	return 1;
}

Point3D OrderedPointSet::GetPoint(int x,int y)
{
	if (x<0 || x>=width || y<0 || y>=height)
		return Point3D(0,0,0);
	int offset=y*width+x;
	if (flags[offset]==0)
		return Point3D(0,0,0);

	return Point3D(X[offset],Y[offset],Z[offset]);
}

int OrderedPointSet::ValidPoint(int x,int y)
{
	if (x<0 || x>=width || y<0 || y>=height)
		return 0;

	return (flags[y*width+x]!=0);
}

int OrderedPointSet::WriteAbs(const char *file)
{
	FILE *f;
	f=fopen(file,"w");
	if (f==NULL)
	{
		fprintf(stderr,"OrderedPointSet::WriteAbs - failed to open %s for output\n",file);
		exit(1);
	}

	// write header
	fprintf(f,"%d rows\n%d columns\npixels (flag X Y Z):\n",height,width);

	// write data
	int x,y;
	for (y=0; y<height; y++)
	for (x=0; x<width; x++)
		fprintf(f,"%d ",flags[y*width+x]);
	fprintf(f,"\n");
	for (y=0; y<height; y++)
	for (x=0; x<width; x++)
		fprintf(f,"%g ",X[y*width+x]);
	fprintf(f,"\n");
	for (y=0; y<height; y++)
	for (x=0; x<width; x++)
		fprintf(f,"%g ",Y[y*width+x]);
	fprintf(f,"\n");
	for (y=0; y<height; y++)
	for (x=0; x<width; x++)
		fprintf(f,"%g ",Z[y*width+x]);
	fprintf(f,"\n");

	fclose(f);

	return 1;
}

int OrderedPointSet::Resample(double res)
{
	GetStats();
	int gridwidth=int((maxX-minX)/res);
	int gridheight=int((maxY-minY)/res);

	if (gridwidth>1024 || gridheight>1024)
	{
		if (verbose)
			fprintf(stderr,"OrderedPointSet::Resample - Cannot resample to grid>1024*1024, skipping resampling\n");
		return 0;
	}

	if (verbose)
		fprintf(stderr,"OrderedPointSet::Resample - attempting to resample to %g mm grid of %dx%d\n",res,gridwidth,gridheight);

	int i,j;
	int grid[gridheight*gridwidth];
	for (j=0; j<gridheight*gridwidth; j++)
		grid[j]=-1;

	int dropped=0;
	for (i=0; i<width*height; i++)
	{
		if (flags[i]==0)
			continue;
		j=int((X[i]-minX)/res)+gridwidth*int((Y[i]-minY)/res);
		if (grid[j]==-1)
			grid[j]=i;
 		else if (Z[i]<Z[grid[j]])
		{
			dropped++;
			grid[j]=i;
		}
	}

	double *Xnew=new double[gridheight*gridwidth];
	double *Ynew=new double[gridheight*gridwidth];
	double *Znew=new double[gridheight*gridwidth];
	int *Fnew=new int[gridheight*gridwidth];

	int empty0=0;
	for (i=0; i<gridwidth*gridheight; i++)
	{
		if (grid[i]==-1)
		{
			// try to interpolate
			empty0++;
		}
		else 
		{
			Xnew[i]=X[grid[i]];
			Ynew[i]=Y[grid[i]];
			Znew[i]=Z[grid[i]];
			Fnew[i]=1;
		}
	}

	int iter,empty=0;
	for (iter=1; iter<5; iter++)
	{
		empty=0;
		int x,y;
		for (y=1; y<gridheight-1; y++)
		for (x=1; x<gridwidth-1; x++)
		{
			i=x+y*gridwidth;
			if (grid[i]==-1)
			{
				// try to interpolate
				empty++;
			}
			else 
			{
				int dir;
				for (dir=0; dir<4; dir++)
				{
					switch (dir)
					{
						case 0: j=i-1; break;
						case 1: j=i+1; break;
						case 2: j=i-gridwidth; break;
						case 3: j=i+gridwidth; break;
					}
					if (grid[j]==-1)
					{
						grid[j]-=iter;
						Fnew[j]=1;
						Xnew[j]=Xnew[i];
						Ynew[j]=Ynew[i];
						Znew[j]=Znew[i];
					}
					else if (grid[j]<-1)
					{
						if (grid[i]>=0)
						{
							// this should never happen
							grid[j]-=iter;
							Fnew[j]=1;
							Xnew[j]=Xnew[i];
							Ynew[j]=Ynew[i];
							Znew[j]=Znew[i];
						}
						else
						{
							Xnew[j]=(Xnew[j]/grid[j]+Xnew[i]/grid[i])*(grid[i]+grid[j]);
							Ynew[j]=(Ynew[j]/grid[j]+Ynew[i]/grid[i])*(grid[i]+grid[j]);
							Znew[j]=(Znew[j]/grid[j]+Znew[i]/grid[i])*(grid[i]+grid[j]);
						}
					}
				}
			}
		}
		fprintf(stderr,"%d\n",empty);
	}
	
	if (verbose)
		fprintf(stderr,"dropped %d; empty0=%d; empty %d\n",dropped,empty0,empty);

	delete [] Xnew;
	delete [] Ynew;
	delete [] Znew;
	delete [] Fnew;

	return 1;
}
}
