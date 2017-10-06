#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#ifdef linux
#include <unistd.h>
#include <malloc.h>
#endif
#include <string.h>
#include <math.h>
#include <time.h>
#include <libgen.h>
#include <X11/Xarch.h>

#include <cv.h>
#include <highgui.h>

#include "matrixio.h"
#include "lda.h"

#include "point3d.h"
#include "plane3d.h"
#include "orderedpointset.h"
#include "unorderedpointset.h"
#include "rangeimage.h"
#include "nose.h"
#include "filter.h"
#include "holes.h"
#include "symmetry.h"
#include "register.h"
#include "roi.h"

#include "faceut3d.h"

void readups(char *file, UnorderedPointSet &ups)
{
	OrderedPointSet ops;
	ops.Read(file);

	int x,y;
	for (y=0; y<ops.height; y++)
	for (x=0; x<ops.width; x++)
	{
		if (ops.ValidPoint(x,y))
			ups.AddPoint(ops.GetPoint(x,y));
	}
	
}

int make1D(RangeImage &ri,IplImage *mask,double *v)
{
	int x,y,k=0;
	for (y=0; y<mask->height; y++)
	{
		unsigned char *pmask=cvPtr2D(mask,y,0);
		for (x=0; x<mask->width; x++)
		{
			if (pmask[x])
				v[k++]=ri.Pixel(x,y);
		}
	}

	return 1;
}

int make1D(Register &reg,IplImage *mask,double *v)
{
	return make1D(reg.ri,mask,v);
}

void transform_im(RangeImage &ri,RangeImage &ro,double tilt_angle,double vertical_shift,double phi)
{
	ro.Init(ri.o,ri.u,ri.v,ri.du,ri.dv,ri.ou,ri.ov,ri.width,ri.height);

	double ny=85.0/130.0*(ri.height-1);
	double nx=0.5*(ri.width-1);
	double cosa=cos(tilt_angle);
	double sina=sin(tilt_angle);

	int x,y;
	int y1f[ri.height];
	for (x=0; x<ri.width; x++)
	{
		for (y=0; y<ri.height; y++)
			y1f[y]=0;
		for (y=0; y<ri.height; y++)
		{
			double z=ri.Pixel(x,y);
			int y1=int((y-ny-vertical_shift)*cosa-z*sina+ny);
			double z1=z*cosa+(y-ny-vertical_shift)*sina;
	
			if (y1>=0 && y1<ri.height)
			{
				ro.Pixel(x,y1)=z1;
				y1f[y1]=1;
			}
		}
		for (y=0; y<ri.height; y++)
			if (y1f[y]==0)
			{
				double val=0;
				int d=0,n=0;
				for (d=1; d<ri.height; d++)
				{
					if ((y-d)>=0 && y1f[y-d]==1)
					{
						val+=ro.Pixel(x,y-d);
						n++;
					}
					if ((y+d)<ri.height && y1f[y+d]==1)
					{
						val+=ro.Pixel(x,y+d);
						n++;
					}
					if (n>0)
						break;
				}
				ro.Pixel(x,y)=val/n;
			}
	}
	if (phi!=0)
	{
		RangeImage rtmp;
		rtmp.CopyDataFrom(ro);
		double cosphi=cos(phi);
		double sinphi=sin(phi);
		for (y=0; y<ro.height; y++)
		for (x=0; x<ro.width; x++)
		{
			double x1=cosphi*(x-nx)-sinphi*(y-ny);
			double y1=sinphi*(x-nx)-cosphi*(y-ny);
			ro.Pixel(x,y)=rtmp.GetDepth(x1,y1);
		}
	}
}

void *getscorethread(void *p)
{
	ThreadedVotes *tv=(ThreadedVotes*)p;

	int i,j;
	for (i=tv->start; i<tv->end; i++)
		tv->votes[i]=tv->mfr->Votes(*(tv->probe),(*(tv->gallery))[i]);
}

int checksfi(char *filename)
{
	char s[8];
	FILE *f=fopen(filename,"r");
	if (f==NULL)
	{
		fprintf(stderr,"registerprobe: failed to open %s for reading\n",filename);
		exit(1);
	}
	fgets(s,8,f);
	fclose(f);

	return (strcmp(s,"CSU_SFI")==0) ? 1 : 0;
}

void registerprobesfi(Classifier &mfr,char *probe_sfi1,char *probe_sfi2,Features &probe)
{
	double templatevector[mfr.nfeatures];

	RangeImage regp[2];
	regp[0].ReadSFI(probe_sfi1);
	regp[1].ReadSFI(probe_sfi2);

	RangeImage regpf[35][2];
	double tilt,shift;
	int k=0;
	for (tilt=-0.04; tilt<=0.04; tilt+=0.02)
	for (shift=-3; shift<=3; shift+=1)
	{
		transform_im(regp[0],regpf[k][0],tilt,shift);
		transform_im(regp[1],regpf[k++][1],tilt,shift);
	}
	
	int i,j;
	for (i=0; i<35; i++)
	{
		mfr.GetFeatures(regpf[i],templatevector);
		probe.Add(templatevector);
	}
}

void registerprobe(Classifier &mfr,char *probe_abs,Features &probe)
{
	double templatevector[mfr.nfeatures];

	UnorderedPointSet probeups;
	readups(probe_abs,probeups);
	Register regprobe(probeups,resolution,holefilling,spikeremoval,ellipticalmask,reflectionremoval,backgroundremoval,nosefitmethod1,LR,symmetrize,motion_threshold,maxshift,nearfrontal);
	RangeImage regp[2];
	regp[0].CopyDataFrom(regprobe.ri);
	
/* this code is much slower than using image transformation, but probably more accurate
	// for fine registration
	RangeImage regpf[35][2];
	Point3D o=regprobe.ri.o;
	Point3D v=regprobe.ri.v;
	Point3D u=regprobe.ri.u;
	double ou=regprobe.ri.ou;
	double ov=regprobe.ri.ov;
	double tilt,shift;
	int k=0;
	for (tilt=-0.04; tilt<=0.04; tilt+=0.02)
	for (shift=-3; shift<=3; shift+=1)
	{
		regprobe.ri.o=o-shift*v;
		regprobe.ri.v=cos(tilt)*v+sin(tilt)*(u^v);
		regprobe.ri.calculateT();
		regprobe.ri.AccumulateDepth(regprobe.ups);
		regprobe.PostProc();
		regpf[k++][0].CopyDataFrom(regprobe.ri);
	}
*/
	regprobe.FitNose(2);
	regprobe.PostProc();
	regp[1].CopyDataFrom(regprobe.ri);
/*
	o=regprobe.ri.o;
	v=regprobe.ri.v;
	u=regprobe.ri.u;
	ou=regprobe.ri.ou;
	ov=regprobe.ri.ov;
	k=0;
	for (tilt=-0.04; tilt<=0.04; tilt+=0.02)
	for (shift=-3; shift<=3; shift+=1)
	{
		regprobe.ri.o=o-shift*v;
		regprobe.ri.v=cos(tilt)*v+sin(tilt)*(u^v);
		regprobe.ri.calculateT();
		regprobe.ri.AccumulateDepth(regprobe.ups);
		regprobe.PostProc();
		regpf[k++][1].CopyDataFrom(regprobe.ri);
	}
*/

	RangeImage regpf[35][2];
	double tilt,shift;
	int k=0;
	for (tilt=-0.04; tilt<=0.04; tilt+=0.02)
	for (shift=-3; shift<=3; shift+=1)
	{
		transform_im(regp[0],regpf[k][0],tilt,shift);
		transform_im(regp[1],regpf[k++][1],tilt,shift);
	}
	
	int i,j;
	for (i=0; i<35; i++)
	{
//char s[256];
//sprintf(s,"finereg%d.sfi",i);
//regpf[i][0].WriteSFI(s);
		mfr.GetFeatures(regpf[i],templatevector);
		probe.Add(templatevector);
	}
}


double *Features::operator[](int i)
{
	return features[i];
}

Features::Features(int nfeatures)
{
	nitems=0;
	this->nfeatures=nfeatures;
	features=NULL;
}


int Features::Add(double *v)
{
	int i;
	double **p=new double*[nitems+1];
	p[nitems]=new double[nfeatures];
	for (i=0; i<nfeatures; i++)
		p[nitems][i]=v[i];

	if (nitems>0)
	{
		for (i=0; i<nitems; i++)
			p[i]=features[i];
		delete features;
	}
	features=p;
	nitems++;

	return nitems;
}

int Features::Write(FILE *f)
{
	int i,j;
	for (i=0; i<nitems; i++)
	{
		for (j=0; j<nfeatures; j++)
			fprintf(f,"%g ",features[i][j]);
		fprintf(f,"\n");
	}

	return 1;
}

int Features::Write(char *filename)
{
	FILE *f=fopen(filename,"w");
	if (f==NULL)
	{
		fprintf(stderr,"Features::Write(): failed to open %s for output\n",filename);
		exit(1);
	}
	Write(f);
	fclose(f);
}

int Features::ReadMax1000(FILE *f)
{
	char s[nfeatures*20],*status;
	double *p[1000],l; 
	int i,j;
	for (i=0;i<1000;i++)
	{
		do { status=fgets(s,nfeatures*20-1,f); } while (s[0]=='#');
		if (status==NULL)
			break;
		p[i]=new double[nfeatures];
		char *pthis=s, *pnext=NULL;
		for (j=0; j<nfeatures; j++)
		{
			p[i][j]=strtod(pthis,&pnext);
			if (pthis[0]=='\0' || pthis==pnext)
			{
				fprintf(stderr,"Features::Read(): failure reading featurs\n");
				exit(1);
			}
			pthis=pnext;
		}
	}

	if (i>0)
	{
		double **pold=features;
		features=new double*[i+nitems];
		for (j=0; j<nitems; j++)
			features[j]=pold[j];
		for (j=0; j<i; j++)
			features[j+nitems]=p[j];
		
		nitems+=i;
	}

	return i;
}

int Features::Read(FILE *f)
{
	int i,j=0;
	do
	{
		i=ReadMax1000(f);
		j+=i;
	} while (i>0);

	return j;
}

int Features::Read(char *filename)
{
	FILE *f=fopen(filename,"r");
	if (f==NULL)
	{
		fprintf(stderr,"Features::Read(): failed to open %s for input\n",filename);
		exit(1);
	}
	Read(f);
	fclose(f);
}

void Features::Empty()
{
	if (nitems>0)
	{
		int i;
		for (i=0; i<nitems; i++)
			delete [] features[i];
		delete [] features;
	}
	nitems=0;
}

Features::~Features()
{
	Empty();
}

// make a templatevector of a set of registered images for gallery
double *Classifier::Enroll(RangeImage *reg,double *templatevector)
{
	int i;
	double v[maxdim];
	double *pg=templatevector;
	for (i=0; i<nclassifiers; i++)
	{
		make1D(reg[regmode[i]],maskim[i],v);
		lda[i].enroll(v,pg);
		pg+=nlda[i];
	}

	return templatevector;
}

// get templatevector of set of registered images for probe
double *Classifier::GetFeatures(RangeImage *reg,double *templatevector)
{
	int i;
	double v[maxdim];
	double *pp=templatevector;
	for (i=0; i<nclassifiers; i++)
	{
		make1D(reg[regmode[i]],maskim[i],v);
		lda[i].features(v,pp);
		pp+=nlda[i];
	}

	return templatevector;
}

double *Classifier::Scores(double *probe,double *gallery)
{
	double *pp=probe;
	double *pg=gallery;

	int i;
	for (i=0; i<nclassifiers; i++)
	{
		score[i]=lda[i].likelihoodratio(pg,pp);
		pp+=nlda[i];
		pg+=nlda[i];
	}

	return score;
}

int Classifier::Votes(double *scores)
{
	int i,nvotes=0;
	for (i=0; i<nclassifiers; i++)
	{
		if (scores[i]>threshold[i])
			nvotes++;
	}

	return nvotes;
}

int Classifier::Votes(double *probe,double *gallery)
{
	double *pp=probe;
	double *pg=gallery;
	double score;

	int i,nvotes=0;
	for (i=0; i<nclassifiers; i++)
	{
		score=lda[i].likelihoodratio(pg,pp);
		pp+=nlda[i];
		pg+=nlda[i];
		if (score>threshold[i])
			nvotes++;
	}

	return nvotes;
}

// for fine registration
int Classifier::Votes(Features &probe,double *gallery)
{
	double *pp;
	double *pg=gallery;
	double bestscore[nclassifiers];

	int i,j;
	for (i=0; i<nclassifiers; i++)
		bestscore[i]=-100000;

	for (j=0; j<probe.nitems; j++)
	{
		pp=probe[j];
		pg=gallery;
		for (i=0; i<nclassifiers; i++)
		{
			double score=lda[i].likelihoodratio(pg,pp);
			if (score>bestscore[i])
				bestscore[i]=score;
			pp+=nlda[i];
			pg+=nlda[i];
		}
	}

	int nvotes=0;
	for (i=0; i<nclassifiers; i++)
	{
		if (bestscore[i]>threshold[i])
			nvotes++;
	}

	return nvotes;
}

// write classifier as a c-file
int Classifier::WriteC(char *file)
{
	FILE *f=fopen(file,"w");
	if (f==NULL)
	{
		fprintf(stderr,"Classifier::WriteC(): Could not open %s for output\n",file);
		exit(1);
	}

	fprintf(f,"#include \"faceut3d.h\"\n\n");

	int n,i,j;
	for (n=0; n<nclassifiers; n++)
	{
		fprintf(f,"char *mfrid%d=(char *) \"%s\";\n",n,id[n]);
		fprintf(f,"int mfrregmode%d=%d;\n",n,regmode[n]);

		// start reading classifier n
		fprintf(f,"double mfrlda%dmeanZ[%d]={%g",n,dim[n],cvGetReal1D(lda[n].meanZ,0));
		for (i=1; i<lda[n].meanZ->rows; i++)
			fprintf(f,",%g",cvGetReal1D(lda[n].meanZ,i));
		fprintf(f,"};\n\n");

		fprintf(f,"double mfrlda%dRho1[%d]={%g",n,lda[n].Rho1->rows,cvGetReal2D(lda[n].Rho1,0,0));
		for (i=1; i<lda[n].Rho1->rows; i++)
			fprintf(f,",%g",cvGetReal2D(lda[n].Rho1,i,i));
		fprintf(f,"};\n\n");
		
		fprintf(f,"double mfrlda%dDdiag[%d]={%g",n,lda[n].Ddiag->rows,cvGetReal2D(lda[n].Ddiag,0,0));
		for (i=1; i<lda[n].Ddiag->rows; i++)
			fprintf(f,",%g",cvGetReal1D(lda[n].Ddiag,i));
		fprintf(f,"};\n\n");
		
		fprintf(f,"\tdouble mfrlda%dT[%d]={%g",n,lda[n].T->rows*lda[n].T->cols,cvGetReal2D(lda[n].T,0,0));
		for (i=0; i<lda[n].T->rows; i++)
		for (j=0; j<lda[n].T->cols; j++)
		{
			if (j==0 && i==0)
				continue;
			fprintf(f,",%g",cvGetReal2D(lda[n].T,i,j));
		}
		fprintf(f,"};\n\n");

		if (maskim[n]!=NULL)
		{
			fprintf(f,"char *mfrmaskim%d=(char *)\"",n);
			int i,j;
			for (i=0; i<maskim[n]->height; i++)
			for (j=0; j<maskim[n]->width; j++)
				fprintf(f,(cvGetReal2D(maskim[n],i,j)==0) ? "0":"1");
			fprintf(f,"\";\n\n");
		}
	}

	fprintf(f,"void InitBuiltInClassifier(Classifier &mfr)\n");
	fprintf(f,"{\n");
	fprintf(f,"\tint i,j;\n");
	fprintf(f,"\tmfr.nclassifiers=%d;\n",nclassifiers);
	fprintf(f,"\tmfr.lda=new LDA[%d];\n",nclassifiers);
	fprintf(f,"\tmfr.dim=new int[%d];\n",nclassifiers);
	fprintf(f,"\tmfr.nlda=new int[%d];\n",nclassifiers);
	fprintf(f,"\tmfr.maskim=new IplImage*[%d];\n",nclassifiers);
	fprintf(f,"\tmfr.threshold=new double[%d];\n",nclassifiers);
	fprintf(f,"\tmfr.score=new double[%d];\n",nclassifiers);
	fprintf(f,"\tmfr.id=new char*[%d];\n",nclassifiers);
	fprintf(f,"\tmfr.regmode=new int[%d];\n",nclassifiers);

	fprintf(f,"\tmfr.maxdim=%d;\n",maxdim);
	fprintf(f,"\tmfr.nfeatures=%d;\n",nfeatures);

	for (n=0; n<nclassifiers; n++)
	{
		fprintf(f,"\n// Setting classifier %d\n",n);
		fprintf(f,"\tmfr.id[%d]=mfrid%d;\n",n,n);
		fprintf(f,"\tmfr.regmode[%d]=mfrregmode%d;\n",n,n);

		// start reading classifier n
		fprintf(f,"\tmfr.lda[%d].cleanup();\n",n);
		fprintf(f,"\tmfr.lda[%d].meanZ=cvCreateMatHeader(%d,1,CV_64F);\n",n,dim[n]);
		fprintf(f,"\tcvSetData(mfr.lda[%d].meanZ,mfrlda%dmeanZ,%d);\n",n,n,dim[n]*sizeof(double));

		fprintf(f,"\tmfr.lda[%d].Rho1=cvCreateMat(%d,%d,CV_64F);\n",n,lda[n].Rho1->rows,lda[n].Rho1->cols);
		fprintf(f,"\tcvSetZero(mfr.lda[%d].Rho1);\n",n);
		fprintf(f,"\tfor (i=0; i<%d; i++)\n",lda[n].Rho1->rows);
			fprintf(f,"\t\tcvSetReal2D(mfr.lda[%d].Rho1,i,i,mfrlda%dRho1[i]);\n",n,n);
		
		fprintf(f,"\tmfr.lda[%d].Ddiag=cvCreateMat(%d,1,CV_64F);\n",n,lda[n].Ddiag->rows);
		fprintf(f,"\tfor (i=0; i<%d; i++)\n",lda[n].Ddiag->rows);
			fprintf(f,"\t\tcvSetReal1D(mfr.lda[%d].Ddiag,i,mfrlda%dDdiag[i]);\n",n,n);
		
		fprintf(f,"\tmfr.lda[%d].T=cvCreateMatHeader(%d,%d,CV_64F);\n",n,lda[n].T->rows,lda[n].T->cols);
		fprintf(f,"\tcvSetData(mfr.lda[%d].T,mfrlda%dT,%d);\n",n,n,lda[n].T->cols*sizeof(double));

		fprintf(f,"\tmfr.lda[%d].T1=cvCreateMatHeader(%d,%d,CV_64F);\n",n,lda[n].T1->rows,lda[n].T1->cols);
       		fprintf(f,"\tmfr.lda[%d].S1diag=cvCreateMatHeader(%d,1,CV_64F);\n",n,lda[n].S1diag->rows);

       		fprintf(f,"\tmfr.lda[%d].Rho1_T=cvCreateMat(%d,%d,CV_64F);\n",n,lda[n].Rho1_T->rows,lda[n].Rho1_T->cols);
       		fprintf(f,"\tcvMatMul(mfr.lda[%d].Rho1,mfr.lda[%d].T,mfr.lda[%d].Rho1_T);\n",n,n,n);
		// end reading classifier n

		fprintf(f,"\tmfr.dim[%d]=%d;\n",n,dim[n]);
		fprintf(f,"\tmfr.nlda[%d]=%d;\n",n,nlda[n]);
		fprintf(f,"\tmfr.threshold[%d]=%g;\n",n,threshold[n]);

		if (maskim[n]==NULL)
			fprintf(f,"\tmfr.maskim[%d]=NULL);\n");
		else
		{
			fprintf(f,"\tmfr.maskim[%d]=cvCreateImage(cvSize(%d,%d),IPL_DEPTH_8U,1);\n",n,maskim[n]->width,maskim[n]->height);
			fprintf(f,"\tfor (i=0; i<mfr.maskim[%d]->height; i++)\n",n);
			fprintf(f,"\tfor (j=0; j<mfr.maskim[%d]->width; j++)\n",n);
				fprintf(f,"\t\tcvSetReal2D(mfr.maskim[%d],i,j,(mfrmaskim%d[i*mfr.maskim[%d]->width+j]=='0')?0:1);\n",n,n,n);
		}
	}

	fprintf(f,"}\n");

	fclose(f);
}

int Classifier::Read(char *file)
{
	FILE *f=fopen(file,"r");
	if (f==NULL)
	{
		fprintf(stderr,"Classifier::Read(): Could not open %s for input\n",file);
		exit(1);
	}

	char s[4096],sid[256];
	do { fgets(s,4095,f); } while (s[0]=='#');
	sscanf(s,"%d",&nclassifiers);

    if (localverbose2)
		fprintf(stderr,"Classifier::Read(): Reading %d classifiers\n",nclassifiers);

	lda=new LDA[nclassifiers];
	dim=new int[nclassifiers];
	nlda=new int[nclassifiers];
	maskim=new IplImage*[nclassifiers];
	threshold=new double[nclassifiers];
	score=new double[nclassifiers];
	id=new char*[nclassifiers];
	regmode=new int[nclassifiers];

	int n;
	nfeatures=0;
	for (n=0; n<nclassifiers; n++)
	{
		do { fgets(s,4095,f); } while (s[0]=='#');
		sscanf(s,"%s",sid);
		id[n]=strdup(sid);
		do { fgets(s,4095,f); } while (s[0]=='#');
		sscanf(s,"%d",&regmode[n]);

        if (localverbose2>1)
		{
			fprintf(stderr,"classifier id=%s\n",id[n]);
			fprintf(stderr,"regmode =%d\n",regmode[n]);
		}

		lda[n].load(f);		
		dim[n]=lda[n].T->cols;
		if (dim[n]>maxdim)
			maxdim=dim[n];
		nlda[n]=lda[n].T->rows;
		nfeatures+=nlda[n];

		int width,height;
		do { fgets(s,4095,f); } while (s[0]=='#');
		sscanf(s,"%d%d\n",&width,&height);

        if (localverbose2>1)
			fprintf(stderr,"Classifier::Read(): reading mask image %dx%d\n",width,height);
		if (width==0 && height==0)
			maskim[n]=NULL;
		else
		{
			maskim[n]=cvCreateImage(cvSize(width,height),IPL_DEPTH_8U,1);
			int i,j;
			for (i=0; i<height; i++)
			{
				fgets(s,4095,f);
				for (j=0; j<width; j++)
					cvSetReal2D(maskim[n],i,j,(s[j]=='0')?0:1);
			}
		}
	}
    if (localverbose2)
		fprintf(stderr,"Classifier::Read(): nfeatures=%d\n",nfeatures);

	// read thresholds
	do { fgets(s,4095,f); } while (s[0]=='#');
	char *pthis=s, *pnext=NULL;
	for (n=0; n<nclassifiers; n++)
	{
		threshold[n]=strtod(pthis,&pnext);
		if (pthis[0]=='\0' || pthis==pnext)
		{
			fprintf(stderr,"Classifier::Read(): failure reading thresholds\n");
			exit(1);
		}
		pthis=pnext;
	}
    if (localverbose2>1)
		fprintf(stderr,"Classifier::Read(): Read %d thresholds\n",nclassifiers);

	fclose(f);
}
