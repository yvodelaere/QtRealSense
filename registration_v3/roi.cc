#include <cv.h>
#include "unorderedpointset.h"
#include "rangeimage.h"
#include "roi.h"

namespace utw3dface {

extern int verbose;
extern int debug;

UpsCell::UpsCell()
{
	pnext=NULL;
	pprev=NULL;
}

UpsCell::UpsCell(UnorderedPointSet &ups)
{
	CopyFrom(ups);
	pnext=NULL;
	pprev=NULL;
}

UpsCell::~UpsCell()
{
}

Point3D &UpsCell::GetNormal()
{
	if (npoints<3)
	{
		quality=0.0;
		normal=Point3D(1,0,0);
		return normal;
	}

	int i;
	mean=Point3D(0,0,0);
	for (i=0; i<npoints; i++)
		mean+=point[i];
	mean/=npoints;

	CvMat *zeromean=cvCreateMat(3,npoints,CV_64F);
	for (i=0; i<npoints; i++)
	{
		cvSetReal2D(zeromean,0,i,point[i].x-mean.x);
		cvSetReal2D(zeromean,1,i,point[i].y-mean.y);
		cvSetReal2D(zeromean,2,i,point[i].z-mean.z);
	}

	CvMat *eigenvalues=cvCreateMat(3,1,CV_64F);
	CvMat *eigenvectors=cvCreateMat(3,3,CV_64F);

	cvSVD(zeromean,eigenvalues,eigenvectors,NULL,CV_SVD_MODIFY_A | CV_SVD_U_T);
	
	// make sure all eigenvectors are always represented in the same way, i.e. set x-coordinate to positive
	if (cvGetReal2D(eigenvectors,0,0)<0)
	{
		for (i=0; i<3; i++)
			cvSetReal2D(eigenvectors,0,i,-cvGetReal2D(eigenvectors,0,i));
	}
	if (cvGetReal2D(eigenvectors,1,0)<0)
	{
		for (i=0; i<3; i++)
			cvSetReal2D(eigenvectors,1,i,-cvGetReal2D(eigenvectors,1,i));
	}
	if (cvGetReal2D(eigenvectors,2,0)<0)
	{
		for (i=0; i<3; i++)
			cvSetReal2D(eigenvectors,2,i,-cvGetReal2D(eigenvectors,2,i));
	}

	quality=cvGetReal1D(eigenvalues,2)/(cvGetReal1D(eigenvalues,0)+cvGetReal1D(eigenvalues,1));

	normal=Point3D(cvGetReal2D(eigenvectors,2,0),cvGetReal2D(eigenvectors,2,1),cvGetReal2D(eigenvectors,2,2));

	if (verbose>2)
		fprintf(stderr,"npoints: %d eigenvalues: %g %g %g quality=%g normal: %g %g %g\n",npoints,cvGetReal1D(eigenvalues,0),cvGetReal1D(eigenvalues,1),cvGetReal1D(eigenvalues,2),quality,normal.x,normal.y,normal.z);

	cvReleaseMat(&eigenvalues);
	cvReleaseMat(&eigenvectors);
	cvReleaseMat(&zeromean);

	return normal;
}

ListOfUpsCells::ListOfUpsCells()
{
	current=begin=NULL;
	ncells=0;
}

UpsCell *ListOfUpsCells::Remove()
{
	if (ncells==0)
		return NULL;

	UpsCell *cell=current;

	current=(cell->pnext!=NULL) ? cell->pnext : cell->pprev;

	if (cell->pprev!=NULL)
		cell->pprev->pnext=cell->pnext;
	if (begin==cell)
		begin=cell->pnext;
	if (cell->pnext!=NULL)
		cell->pnext->pprev=cell->pprev;

	delete cell;
	ncells--;

	return current;
}

ListOfUpsCells::~ListOfUpsCells()
{
	while (Remove()!=NULL);
}

UpsCell *ListOfUpsCells::Append(UnorderedPointSet &ups)
{
	UpsCell *cell=new UpsCell(ups);

	if (ncells==0)
	{
		current=begin=cell;
		cell->pprev=NULL;
		cell->pnext=NULL;
	}
	else
	{
		cell->pprev=current;
		cell->pnext=current->pnext;
		if (current->pnext!=NULL)
			current->pnext->pprev=cell;
		current->pnext=cell;
	}

	current=cell;
	ncells++;

	return current;
}

UpsCell *ListOfUpsCells::Insert(UnorderedPointSet &ups)
{
	UpsCell *cell=new UpsCell(ups);

	if (ncells==0)
	{
		current=begin=cell;
		cell->pprev=NULL;
		cell->pnext=NULL;
	}
	else
	{
		cell->pnext=current;
		cell->pprev=current->pprev;
		if (current->pprev!=NULL)
			current->pprev->pnext=cell;
		if (begin==current)
			begin=cell;
		current->pprev=cell;
	}

	current=cell;
	ncells++;

	return current;
}

UpsCell *ListOfUpsCells::Current()
{
	return current;
}

UpsCell *ListOfUpsCells::Next()
{
	if (current!=NULL && current->pnext!=NULL)
	{
		current=current->pnext;
		return current;
	}
	else
		return NULL;
}

UpsCell *ListOfUpsCells::Prev()
{
	if (current!=NULL && current->pprev!=NULL)
	{
		current=current->pprev;
		return current;
	}
	else
		return NULL;
}

UpsCell *ListOfUpsCells::Begin()
{
	current=begin;
	return current;
}

void Cylinder::Init()
{
	maxdist=20; // [mm]
	maxangle=M_PI/4;
	faceheight=200;
	faceminr=75;
	facemaxr=100;
	support=0;
}

Cylinder::Cylinder()
{
	Init();
}

Cylinder::Cylinder(const Cylinder &cyl)
{
	Init();
	a=cyl.a;
	c=cyl.c;
	r=cyl.r;
	support=cyl.support;
}

Cylinder::Cylinder(const Cylinder *pcyl)
{
	Init();
	a=pcyl->a;
	c=pcyl->c;
	r=pcyl->r;
	support=pcyl->support;
}

Cylinder::Cylinder(const Point3D &p1,const Point3D &p2,const Point3D &n1,const Point3D &n2)
{
	Init();

	// axis of cylinder is perpendicular to both n1 and n2: a=n1Xn2
	a=n1^n2;
	double norm=a.Norm();
	a/=(norm==0) ? 1 : norm;

	// make sure axis always points upwards
	if (a.y<0)
		a*=-1;

	// get a plane through p1 and perpendicular to a and n1: (aXn1)^T*x = (aXn1)^T*p1
	Point3D an1=a^n1;

	// the intersection of the plane with the line through p2 with direction n2 is on the axis of the cylinder
	// (aXn1)^T*(p2+lambda*n2) = (aXn1)^T*p1
	double lambda=an1*(p1-p2)/(an1*n2);
	c=p2+lambda*n2;

	r=fabs(lambda);
}

double Cylinder::Distance(const Point3D &p)
{
	// s = c + a^T*(p-c)*a is the intersection of a plane through p perpendicular to a and the line c+lambda*a
	Point3D s=c+(a*(p-c))*a;
	Point3D ps=p-s;
	Point3D cs=c-s;

	double dr=fabs(ps.Norm()-r);
	if (ps.z<-r/2)
		dr=100000;
	double da=cs.Norm();

	return (da<0.5*faceheight) ? dr : dr+da-0.5*faceheight;
}

double Cylinder::AngleDiff(const Point3D &p,const Point3D &n)
{
	double angle;

	// s = c + a^T*(p-c)*a is the intersection of a plane through p perpendicular to a and the line c+lambda*a
	Point3D ps=(p-c)-(a*(p-c))*a;
	double norm=ps.Norm();	// p is on the line c+lambda*a => we cannot calculate an angle!
	if (norm==0)
		angle=M_PI/2;
	else
	{
		ps/=norm;
		double angle=acos(ps*n);
		if (angle>M_PI/2)
			angle=M_PI-angle;
	}

	return angle;
}

int Cylinder::Support(const Point3D &p,const Point3D &n)
{
	double angle; 

	// s = c + a^T*(p-c)*a is the intersection of a plane through p perpendicular to a and the line c+lambda*a
	Point3D s=c+(a*(p-c))*a;
	double sc=(s.y-c.y)/a.y;
	if (fabs(sc)>0.5*faceheight)
		return 0;
	Point3D ps=p-s;
	if (ps.z<0)
		return 0;
	double norm=ps.Norm();
	if (norm==0)
		angle=M_PI/2;
	else
	{
		ps/=norm;
		angle=acos(ps*n);
		if (angle>M_PI/2)
			angle=M_PI-angle;
	}
	
	if (norm>r+maxdist && norm-r<100)
		return -1;
	return (fabs(norm-r)<maxdist && angle<maxangle);
}

double Cylinder::Support(UnorderedPointSet &ups,UnorderedPointSet &upsnormals)
{
	support=0;
	int k;
	for (k=0; k<ups.npoints; k++)
		support+=Support(ups[k],upsnormals[k]);

	//support/=r;

	return support;
}

double Cylinder::Support(UnorderedPointSet &ups,double &dist)
{
	int i;
	double n=0;
	dist=0;
	for (i=0; i<ups.npoints; i++)
	{
		double d=Distance(ups[i]);
		if (d<maxdist)
		{
			n++;
			dist+=d;
		}
	}

	dist/=(n>0) ? n : 1;

	return n;
}

CylCell::CylCell(const Point3D &p1,const Point3D &p2,const Point3D &n1,const Point3D &n2):Cylinder(p1,p2,n1,n2)
{
	pnext=NULL;
	pprev=NULL;
}

CylCell::CylCell(const Cylinder &cyl): Cylinder(cyl)
{
	pnext=NULL;
	pprev=NULL;
}

CylCell::CylCell(const Cylinder *pcyl): Cylinder(pcyl)
{
	pnext=NULL;
	pprev=NULL;
}

CylCell::~CylCell()
{
}

ListOfCylCells::ListOfCylCells()
{
	current=begin=NULL;
	ncells=0;
}

CylCell *ListOfCylCells::Remove()
{
	if (ncells==0)
		return NULL;

	CylCell *cell=current;

	current=(cell->pnext!=NULL) ? cell->pnext : cell->pprev;

	if (cell->pprev!=NULL)
		cell->pprev->pnext=cell->pnext;
	if (begin==cell)
		begin=cell->pnext;
	if (cell->pnext!=NULL)
		cell->pnext->pprev=cell->pprev;

	delete cell;
	ncells--;

	return current;
}

ListOfCylCells::~ListOfCylCells()
{
	while (Remove()!=NULL);
}

CylCell *ListOfCylCells::Append(Cylinder &cyl)
{
	CylCell *cell=new CylCell(cyl);

	if (ncells==0)
	{
		current=begin=cell;
		cell->pprev=NULL;
		cell->pnext=NULL;
	}
	else
	{
		cell->pprev=current;
		cell->pnext=current->pnext;
		if (current->pnext!=NULL)
			current->pnext->pprev=cell;
		current->pnext=cell;
	}

	current=cell;
	ncells++;

	return current;
}

CylCell *ListOfCylCells::Insert(Cylinder &cyl)
{
	CylCell *cell=new CylCell(cyl);

	if (ncells==0)
	{
		current=begin=cell;
		cell->pprev=NULL;
		cell->pnext=NULL;
	}
	else
	{
		cell->pnext=current;
		cell->pprev=current->pprev;
		if (current->pprev!=NULL)
			current->pprev->pnext=cell;
		if (begin==current)
			begin=cell;
		current->pprev=cell;
	}

	current=cell;
	ncells++;

	return current;
}

CylCell *ListOfCylCells::Current()
{
	return current;
}

CylCell *ListOfCylCells::Next()
{
	if (current!=NULL && current->pnext!=NULL)
	{
		current=current->pnext;
		return current;
	}
	else
		return NULL;
}

CylCell *ListOfCylCells::Prev()
{
	if (current!=NULL && current->pprev!=NULL)
	{
		current=current->pprev;
		return current;
	}
	else
		return NULL;
}

CylCell *ListOfCylCells::Begin()
{
	current=begin;
	return current;
}

///////////////////////////////////////////////////////////////////
// methods for FaceRoi class
///////////////////////////////////////////////////////////////////

void FaceRoi::Init()
{
	cellsize=20;			// cell size in mm
	minpointspercell=20;
	maxtilt=cos(45*M_PI/180);	// allow deviation of 45 degrees from vertical axis
	minnormalquality=0.3;		// the largest eigenvalue must have magnitude >0.3
	cylpoints_maxxdist=150;		// max distance in x-dir between cells defining a cylinder
	cylpoints_maxydist=50;		// max distance in y-dir between cells defining a cylinder
	maxdisttocylinder=75;		// max dist of 3D points to cylinder
}

void FaceRoi::GetRupsAndNormals(UnorderedPointSet &ups)	// get reduced ups and normals
{
	ups.GetStats();
	w=int((ups.maxX-ups.minX)/cellsize)+1;
	h=int((ups.maxY-ups.minY)/cellsize)+1;

	if (verbose)
		fprintf(stderr,"FaceRoi: created grid with %d x %d cells\n",w,h);

	UnorderedPointSet *grid=new UnorderedPointSet[w*h];

	int i,x,y;
	for (i=0; i<ups.npoints; i++)
	{
		x=int((ups[i].x-ups.minX)/cellsize);
		y=int((ups[i].y-ups.minY)/cellsize);
		grid[y*w+x].AddPoint(ups[i]);
	}
	
	for (i=0; i<w*h; i++)
	{
		if (grid[i].npoints>=minpointspercell)
			lups.Append(grid[i]);
	}

	delete [] grid;

	if (verbose)
		fprintf(stderr,"FaceRoi: %d cells with >= %d points\n",lups.ncells,minpointspercell);

	UpsCell *pups;

	for (pups=lups.Begin(); pups!=NULL; pups=lups.Next())
		pups->GetStats();

	for (pups=lups.Begin(); pups!=NULL; pups=lups.Next())
	{
		if (pups->npoints<minpointspercell)
		{
			lups.Remove();
			lups.Prev();
			continue;
		}
		while (pups->maxZ-pups->minZ>cellsize)
		{
			UnorderedPointSet ups;
			UnorderedPointSet leftover;
			double zmin=pups->minZ;
			int i;
			for (i=0; i<pups->npoints; i++)
			{
				if ((*pups)[i].z-zmin<cellsize)
					ups.AddPoint((*pups)[i]);
				else
					leftover.AddPoint((*pups)[i]);
			}
			if (ups.npoints>minpointspercell)
				lups.Insert(ups);
			lups.current->GetStats();
			pups->CopyFrom(leftover);
			pups->GetStats();
			if (pups->npoints<minpointspercell)
				break;
		}
	}

	if (verbose)
		fprintf(stderr,"%d cells with > %d points and size %gx%gx%g [mm]\n",lups.ncells,minpointspercell,cellsize,cellsize,cellsize);

	for (pups=lups.Begin(); pups!=NULL; pups=lups.Next())
		pups->GetNormal();

	for (pups=lups.Begin(); pups!=NULL; pups=lups.Next())
	{
		if (pups->quality<minnormalquality)
		{
			rups.AddPoint(pups->mean);
			rupsn.AddPoint(pups->normal);
		}
	}

	if (verbose)
		fprintf(stderr,"FaceRoi: %d cells with good normals\n",rups.npoints);
}

void FaceRoi::DumpRups(Cylinder &cyl,const char *filename,Point3D &p,int onlysupported)
{
	UnorderedPointSet upsc;

	upsc.CopyFrom(rups);

	// add the axis
	double lambda;
	for (lambda=-300; lambda<300; lambda++)
		upsc.AddPoint(cyl.c+lambda*cyl.a);

	Point3D e1=cyl.a^(cyl.a+Point3D(-1,0,0));
	e1.Normalise();
	Point3D e2=cyl.a^e1;
	double phi;
	double h2=cyl.faceheight*0.5;
	for (phi=0; phi<M_PI; phi+=0.01)
	{
		upsc.AddPoint(cyl.c+h2*cyl.a+cyl.r*(sin(phi)*e1+cos(phi)*e2));
		upsc.AddPoint(cyl.c+cyl.r*(sin(phi)*e1+cos(phi)*e2));
		upsc.AddPoint(cyl.c-h2*cyl.a+cyl.r*(sin(phi)*e1+cos(phi)*e2));

		upsc.AddPoint(cyl.c+h2*cyl.a+(100+cyl.r)*(sin(phi)*e1+cos(phi)*e2));
		upsc.AddPoint(cyl.c-h2*cyl.a+(100+cyl.r)*(sin(phi)*e1+cos(phi)*e2));
	}
	for (lambda=-h2; lambda<h2; lambda++)
	{
		upsc.AddPoint(cyl.c+lambda*cyl.a+cyl.r*(sin(0.0)*e1+cos(0.0)*e2));
		upsc.AddPoint(cyl.c+lambda*cyl.a+cyl.r*(sin(M_PI/2)*e1+cos(M_PI/2)*e2));
		upsc.AddPoint(cyl.c+lambda*cyl.a+cyl.r*(sin(M_PI)*e1+cos(M_PI)*e2));
	}

	int k;
	for (k=0; k<rups.npoints; k++)
	{
		if (onlysupported || cyl.Support(rups[k],rupsn[k])==0)
			continue;

		int j;
		for (j=0; j<10; j++)
			upsc.AddPoint(rups[k]+j*rupsn[k]);
	}

	Point3D s=cyl.c+(cyl.a*(p-cyl.c))*cyl.a;
	Point3D n=(p-s);
	if (n.Norm()>0.001)
	{
		n.Normalise();
		Point3D t=n^cyl.a;
		for (lambda=-10; lambda<10; lambda++)
		{
			upsc.AddPoint(p+lambda*cyl.a);
			upsc.AddPoint(p+lambda*t);
		}
	}

	upsc.WriteVRML(filename);
}

FaceRoi::FaceRoi()
{
	Init();
}

FaceRoi::FaceRoi(UnorderedPointSet &ups)
{
	Init();

	GetRupsAndNormals(ups);

	int i,j;
	for (i=0; i<rups.npoints; i++)
	for (j=0; j<rups.npoints; j++)
	{
		if (i==j)
			continue;
		if (fabs(rups[i].y-rups[j].y)>cylpoints_maxydist)	// maximum y-distance between points is 50 [mm]
			continue;
		if (rups[j].x-rups[i].x<50 || rups[j].x-rups[i].x>cylpoints_maxxdist)	// x-distance is 50-150 [mm]
			continue;
		Cylinder cyl(rups[i],rups[j],rupsn[i],rupsn[j]);
		if (cyl.r>cyl.facemaxr || cyl.r<cyl.faceminr)
			continue;
		if (fabs(cyl.a.y)<maxtilt)	// face must be more or less vertical
			continue;
		double support=cyl.Support(rups,rupsn);
		
		CylCell *pcyl;
		for (pcyl=lcyl.Begin(); pcyl!=NULL; pcyl=lcyl.Next())
		{
			if (cyl.support>pcyl->support)
				break;
		}
		if (pcyl!=NULL)
			lcyl.Insert(cyl);
		else
			lcyl.Append(cyl);

	}

	if (verbose)
	{
		fprintf(stderr,"roi: Checked %d cylinders\n",lcyl.ncells);
		if (lcyl.ncells>0)
			fprintf(stderr,"roi: max support: %g\n",lcyl.begin->support);
	}

	if (lcyl.ncells==0)
	{
		for (i=0; i<ups.npoints; i++)
		{
			upsroi.AddPoint(ups[i]);
		}
		return;
	}

	// merge cylinders close to eachother
	ListOfCylCells lcylr;
	double maxdphi=30*M_PI/180;
	double maxdc=40;
	while (lcyl.ncells>0)
	{
		CylCell *pcylr=lcylr.Append(*lcyl.Begin());
		lcyl.Remove();
		CylCell *pcyl;
		for (pcyl=lcyl.Begin(); pcyl!=NULL; pcyl=lcyl.Next())
		{
			double dcosphi=acos(pcyl->a*pcylr->a);
			if (dcosphi>M_PI/2)
				dcosphi=M_PI-dcosphi;
			Point3D cc=(pcyl->c-pcylr->c)-(pcylr->a*(pcyl->c-pcylr->c))*pcylr->a;
			double dc=cc.Norm();	
			if (dcosphi<maxdphi && dc<maxdc)
			{
				lcyl.Remove();
				lcyl.Prev();
			}
		}
	} 

	if (debug)
	{
		printf("# cylinder fitting stuff\n");
		CylCell *pcyl;
		for (pcyl=lcylr.Begin(); pcyl!=NULL; pcyl=lcylr.Next())
			printf("%g %g %g %g %g %g %g %g\n",pcyl->support,pcyl->c.x,pcyl->c.y,pcyl->c.z,pcyl->a.x,pcyl->a.y,pcyl->a.z,pcyl->r);
		printf("\n\n");
	}

	lcylr.Begin();
	Cylinder cyl(lcylr.current);
/*	
	double polar[10][36];
	int flag[10][36],iphi,iy;
	for (iy=0; iy<10; iy++)
	for (iphi=0; iphi<36; iphi++)
		flag[iy][iphi]=0;
	for (i=0; i<rups.npoints; i++)
	{
		Point3D p=rups[i];
		Point3D s=cyl.c+(cyl.a*(p-cyl.c))*cyl.a;
		double r=s.Distance(p);
		if (fabs(r-cyl.r)>30)
			continue;
		Point3D ps=(p-s)/r;
		double y=(p.y-cyl.c.y)/cyl.a.y;
		Point3D nz=Point3D(1,0,0)^cyl.a;
		Point3D nx=cyl.a^nz;
		double phi=atan2(nz*ps,nx*ps);
		iphi=(phi<0) ? int(18*(2+phi/M_PI)) : int(18*phi/M_PI);
		iy=int(y/cellsize + 0.5*cyl.faceheight/cellsize);
		if (iy>=0 && iy<10 && iphi>=0 && iphi<36)
		{
			if (flag[iy][iphi])
			{
				if (polar[iy][iphi]>r)
					polar[iy][iphi]=r;
			}
			else
			{
				polar[iy][iphi]=r;
				flag[iy][iphi]=1;
			}
		}
	}

	int mphi;
	for (mphi=0; mphi<36; mphi++)
	{
		double d=0;
		int nd=0;

		for (iphi=0; iphi<36; iphi++)
		for (iy=0; iy<10; iy++)
		{
			if (!flag[iy][iphi])
				continue;
			int jphi=2*mphi-iphi;
			while (jphi<0)
				jphi+=36;
			while (jphi>=36)
				jphi-=36;
			if (!flag[iy][jphi])
				continue;
			d+=fabs(polar[iy][iphi]-polar[iy][jphi]);
			nd++;
		}

		printf("%d %g %d\n",mphi,d,nd);
	}


for (iy=0; iy<10; iy++)
for (iphi=0; iphi<36; iphi++)
fprintf(stderr,"%d %d %g %d\n",iphi,iy,polar[iy][iphi],flag[iy][iphi]);
*/
	if (debug)	// create some nice graphs of the cylinder and normals
	{
		DumpRups(cyl,"rups.wrl",cyl.c,0);
		DumpRups(cyl,"upsc.wrl",cyl.c,1);
/*
DumpRups(*lcylr.Next(),"rups1.wrl",cyl.c,0);
DumpRups(*lcylr.Next(),"rups2.wrl",cyl.c,0);
DumpRups(*lcylr.Next(),"rups3.wrl",cyl.c,0);
DumpRups(*lcylr.Next(),"rups4.wrl",cyl.c,0);
DumpRups(*lcylr.Next(),"rups5.wrl",cyl.c,0);
*/
	}

	for (i=0; i<ups.npoints; i++)
	{
		if (cyl.Distance(ups[i])>maxdisttocylinder)
			continue;
		upsroi.AddPoint(ups[i]);
	}

	cylinder=cyl;
}

FaceRoi::~FaceRoi()
{
}

}
