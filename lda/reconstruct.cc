#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

#include <cv.h>

#include "matrixio.h"
#include "lda.h"

using namespace utw3dface;

namespace utw3dface {
int verbose=0;
}

typedef struct {char a,b,c,d;} EL;
typedef union {float f; EL elem; } NUM;

class Roi
{
public:
	double left,right,bottom,top;

	Roi(double left,double right,double bottom,double top): left(left),right(right),bottom(bottom),top(top)
	{
		if (left>=right || bottom>=top)
		{
			fprintf(stderr,"illegal ROI\n");
			exit(1);
		}
	} 
};

CvMat *readsfi(char *fname)
{
	FILE *f;
	int channels,width,height;
	NUM junk;
	char firstline[1000];
	char ftype[1000];

	f = fopen( fname, "rb" );
	if ( !f ) { fprintf(stderr,"Can't open %s\n", fname); exit(1); }

	fgets(firstline,1000,f);
	sscanf(firstline,"%s %d %d %d", ftype,&width, &height, &channels);

	if( !(strcmp(ftype,"CSU_SFI") == 0) )
	{
		fclose(f);
		fprintf(stderr,"not a SFI image\n");
		exit(1);
	}

	if (channels!=1)
	{
		fprintf(stderr,"readsfi: can only read 1 channel images\n");
		exit(1);
	}

	if (verbose)
		fprintf(stderr,"reading %d x %d sfi image from %s\n",width,height,fname);

	CvMat *mat=cvCreateMat(height,width,CV_32F);

	int i,j;
	for(j=0;j<height;j++)
	for(i=0;i<width;i++)
	{
		/* read in the correct byte order for floating point values */
		if (__BYTE_ORDER == __LITTLE_ENDIAN)
		{
		    fread(&(junk.elem.d),1,1,f);
		    fread(&(junk.elem.c),1,1,f);
		    fread(&(junk.elem.b),1,1,f);
		    fread(&(junk.elem.a),1,1,f);
		}
		else{
		    fread(&(junk.elem.a),1,1,f);
		    fread(&(junk.elem.b),1,1,f);
		    fread(&(junk.elem.c),1,1,f);
		    fread(&(junk.elem.d),1,1,f);
		}

		cvSetReal2D(mat,j,i,junk.f);
	}

	fclose(f);

	return mat;
}

void writesfi(char *fname,CvMat *mat)
{
	FILE *f;
	int channels=1,width=mat->cols,height=mat->rows;
	NUM junk;

	f = fopen( fname, "wb" );
	if ( !f ) { fprintf(stderr,"Can't open %s\n", fname); exit(1); }

	fprintf(f,"%s %d %d %d\n","CSU_SFI",width,height,1);

	if (verbose)
		fprintf(stderr,"writing %d x %d sfi image to %s\n",width,height,fname);

	int i,j;
	for(j=0;j<height;j++)
	for(i=0;i<width;i++)
	{
		junk.f=cvGetReal2D(mat,j,i);

		/* write in the correct byte order for floating point values */
		if (__BYTE_ORDER == __LITTLE_ENDIAN)
		{
		    fwrite(&(junk.elem.d),1,1,f);
		    fwrite(&(junk.elem.c),1,1,f);
		    fwrite(&(junk.elem.b),1,1,f);
		    fwrite(&(junk.elem.a),1,1,f);
		}
		else{
		    fwrite(&(junk.elem.a),1,1,f);
		    fwrite(&(junk.elem.b),1,1,f);
		    fwrite(&(junk.elem.c),1,1,f);
		    fwrite(&(junk.elem.d),1,1,f);
		}
	}

	fclose(f);
}

CvMat *make1D(CvMat *mat,Roi &roi,CvMat *v=NULL)
{
	int width=mat->cols;
	int h_min=int(mat->rows*roi.bottom);
	int h_max=int(mat->rows*roi.top);
	int w_min=int(mat->cols*roi.left);
	int w_max=int(mat->cols*roi.right);

	if (v==NULL)
		v=cvCreateMat((h_max-h_min)*(w_max-w_min),1,CV_64F);

	if (v->cols*v->rows!=(h_max-h_min)*(w_max-w_min))
	{
		fprintf(stderr,"feature vector size mismatch\n");
		exit(1);
	}

	int i,j;
	for (i=h_min; i<h_max; i++)
	for (j=w_min; j<w_max; j++)
		cvSetReal1D(v,(i-h_min)*(w_max-w_min)+j-w_min,cvGetReal2D(mat,i,j));

	return v;	
}

void make2D(CvMat *mat,Roi &roi,CvMat *v)
{
	int width=mat->cols;
	int h_min=int(mat->rows*roi.bottom);
	int h_max=int(mat->rows*roi.top);
	int w_min=int(mat->cols*roi.left);
	int w_max=int(mat->cols*roi.right);

	if (v->cols*v->rows!=(h_max-h_min)*(w_max-w_min))
	{
		fprintf(stderr,"feature vector size mismatch\n");
		exit(1);
	}

	int i,j;
	for (i=h_min; i<h_max; i++)
	for (j=w_min; j<w_max; j++)
		cvSetReal2D(mat,i,j,cvGetReal1D(v,(i-h_min)*(w_max-w_min)+j-w_min));
}

CvMat *readsfi1D(char *fname,Roi &roi,CvMat *v=NULL)
{
	CvMat *mat=readsfi(fname);
	v=make1D(mat,roi,v);

	cvReleaseMat(&mat);

	return v;	
}

void train(LDA &lda,int nPCA,int nLDA,char *flist,Roi &roi,char *fidlist=NULL)
{
	// first determine how long the list is
	FILE *f;
	f=fopen(flist,"r");
	if (f==NULL)
	{
		fprintf(stderr,"failed to read list of filenames\n");
		exit(1);
	}
	
	FILE *fid=NULL;
	if (fidlist!=NULL)
	{
		fid=fopen(fidlist,"r");
		if (fid==NULL)
		{
			fprintf(stderr,"failed to read list of id's\n");
			exit(1);
		}
	}

	char fname[1000];
	int i,j;
	for (i=0; ; i++)
	{
		int status=fscanf(f,"%s",fname);
		if (status==EOF || status==0)
			break;
	}
	int ll=i;

	if (verbose)
		fprintf(stderr,"processing %d images in list\n",ll);

	rewind(f);

	// read first image to find out size
	fscanf(f,"%s",fname);
	CvMat *im=readsfi1D(fname,roi);
	CvMat *Z=cvCreateMat(im->rows,ll,CV_64F);
	CvMat *iZ=cvCreateMat(1,ll,CV_64F);
	cvReleaseMat(&im);
	rewind(f);

	int *idlist=new int[100000];
        for (i=0; i<ll; i++)
                idlist[i]=0;
	int id,nid=0;

	CvMat col=cvMat(Z->rows,1,CV_64F);
	for (i=0; i<ll; i++)
	{
		if (verbose)
			fprintf(stderr,"%d ",i);

		fscanf(f,"%s",fname);
		cvGetCol(Z,&col,i);
		readsfi1D(fname,roi,&col);

		if (fid!=NULL)
		{
			if (fscanf(fid,"%d",&id)!=1)
			{
				fprintf(stderr, "could not read id %d from id file\n",i);
				exit(1);
			}
		}
		else	// determine id from filename, not this only works properly for FRGC scans!
		{
			char *bname=basename(fname);
			sscanf(bname,"%d",&id);
		}
		if (idlist[id]==0)
			idlist[id]=++nid;
		cvSetReal1D(iZ,i,idlist[id]);
	}
	fclose(f);
	if (fid!=NULL)
		fclose(fid);

	delete [] idlist;

	if (verbose)
		fprintf(stderr,"%d identities, %d vectors, dim=%d\n",nid,ll,Z->rows);

	lda.train(Z,iZ,nPCA,nLDA);

	cvReleaseMat(&Z);
	cvReleaseMat(&iZ);
}

CvMat *getscorematrix(LDA &lda,char *flist,Roi &roi)
{
	// first determine how long the list is
	FILE *f;
	f=fopen(flist,"r");
	if (f==NULL)
	{
		fprintf(stderr,"failed to read list of filenames\n");
		exit(1);
	}

	char fname[1000];
	int i,j;
	for (i=0; ; i++)
	{
		int status=fscanf(f,"%s",fname);
		if (status==EOF || status==0)
			break;
	}
	int ll=i;

	if (verbose)
		fprintf(stderr,"processing %d images in list\n",ll);

	rewind(f);

	CvMat **gallery=new CvMat*[ll];

	for (i=0; i<ll; i++)
	{
		if (verbose)
			fprintf(stderr,"%d ",i);
		fscanf(f,"%s",fname);
		CvMat *v=readsfi1D(fname,roi);
		gallery[i]=lda.enroll(v);
		cvReleaseMat(&v);
	}

	CvMat *scorematrix=cvCreateMat(ll,ll,CV_64F);

	if (verbose)
		fprintf(stderr,"generating %d x %d score matrix\n",ll,ll);

	rewind(f);

	CvMat *probe;
	for (i=0; i<ll; i++)
	{
		if (verbose)
			fprintf(stderr,"%d ",i);
		fscanf(f,"%s",fname);
		CvMat *v=readsfi1D(fname,roi);
		probe=lda.features(v);
		if (verbose>1)
		{
			for (j=0; j<probe->rows; j++)
				printf("%g ",cvGetReal1D(probe,j));
			printf("\n");
		}
		cvReleaseMat(&v);
		for (j=0; j<ll; j++)
			cvSetReal2D(scorematrix,i,j,lda.likelihoodratio(gallery[j],probe));
		cvReleaseMat(&probe);
	}

	for (i=0; i<ll; i++)
		cvReleaseMat(&gallery[i]);

	delete [] gallery;

	fclose(f);

	return scorematrix;
}

const char *usage="usage: standardLDA [OPTIONS]\n\
OPTIONS:\n\
\t-I file     input image\n\
\t-O file     output image\n\
\t-i file     list of id's of training images as numbers from 1..99999 on one long line\n\
\t-l int      set number of LDA components (25)\n\
\t-p int      set number of PCA components (100)\n\
\t-r basename read meanZ, T,Rho1 and D matrices from files <basename>_meanZ/T/Rho/D.ascii\n\
\t-t file     list of training images (first 5 digits=id or use -i for id-file)\n\
\t-v	      increase verbosity\n\
\t-w basename write meanZ, T,Rho1 and D matrices in <basename>_meanZ/T/Rho/D.ascii\n\
\t-B float    bottom of roi of image [0..top> [0.0]\n\
\t-L float    left of roi of image [0..right> [0.0]\n\
\t-R float    right of roi of image <left..1] [1.0]\n\
\t-T float    top of roi of image <bottom..1] [1.0]\n\
";

main(int argc,char **argv)
{
	char c;
	char *traininglist=NULL;
	char *idlist=NULL;
	char *ofile=NULL;
	char *ifile=NULL;
	char *savebasename=NULL;
	char *loadbasename=NULL;
	int nPCA=100;
	int nLDA=25;
	double top=1.0;
	double bottom=0.0;
	double left=0.0;
	double right=1.0;

	while ((c=getopt(argc,argv,"I:f:O:r:w:t:vp:l:B:T:L:R:i:"))!=-1)
	{
		switch (c)
		{
			case 'i': idlist=optarg; break;
			case 't': traininglist=optarg; break;
			case 'v': verbose++; break;
			case 'L': sscanf(optarg,"%lf",&left); break;
			case 'R': sscanf(optarg,"%lf",&right); break;
			case 'T': sscanf(optarg,"%lf",&top); break;
			case 'B': sscanf(optarg,"%lf",&bottom); break;
			case 'p': sscanf(optarg,"%d",&nPCA); break;
			case 'l': sscanf(optarg,"%d",&nLDA); break;
			case 'O': ofile=optarg; break;
			case 'I': ifile=optarg; break;
			case 'r': loadbasename=optarg; break;
			case 'w': savebasename=optarg; break;
			default:
				fprintf(stderr,usage);
				exit(1);
		}
	}

	Roi roi(left,right,bottom,top);

	if (loadbasename==NULL && traininglist==NULL)
	{
		fprintf(stderr,"must provide either -L or -t option!\n");
		exit(1);
	}

	LDA lda;
	if (loadbasename!=NULL)
		lda.load(loadbasename);
	else
	{
		if (idlist==NULL)
			fprintf(stderr,"WARNING: extracting id's from filenames, only correct for FRGC data!\n");
		train(lda,nPCA,nLDA,traininglist,roi,idlist);
	}

	if (savebasename!=NULL)
		lda.save(savebasename);

	if (ifile==NULL)
		exit(0);

	CvMat *mat=readsfi(ifile);
	CvMat *v=make1D(mat,roi,NULL);

{
	int i;
	for (i=0; i<v->rows/2; i++)
		cvSetReal1D(v,i,cvGetReal1D(lda.meanZ,i));
}
{
	int i;
	for (i=0; i<v->rows; i++)
		if (cvGetReal1D(v,i)==0)
			cvSetReal1D(v,i,cvGetReal1D(lda.meanZ,i));
}

        CvMat *v_zeromean=cvCreateMat(lda.T1->cols,1,CV_64F);
        cvSub(v,lda.meanZ,v_zeromean);

        CvMat *v_pca=cvCreateMat(lda.T1->rows,1,CV_64F);
        cvMatMul(lda.T1,v_zeromean,v_pca);

{
	int i;
	for (i=0; i<v_pca->rows; i++)
		printf("%d %g\n",i,cvGetReal1D(v_pca,i)/cvGetReal1D(lda.S1diag,i));
}
	double norm=cvNorm(v_pca);
	fprintf(stderr,"norm=%g\n",norm);

	CvMat *T1T=cvCreateMat(lda.T1->cols,lda.T1->rows,CV_64F);
	cvTranspose(lda.T1,T1T);
        CvMat *vrec_zeromean=cvCreateMat(lda.T1->cols,1,CV_64F);
	cvMatMul(T1T,v_pca,vrec_zeromean);
        CvMat *vrec=cvCreateMat(lda.T1->cols,1,CV_64F);
        cvAdd(lda.meanZ,vrec_zeromean,vrec);

        CvMat *vdiff=cvCreateMat(lda.T1->cols,1,CV_64F);
        cvSub(v,vrec,vdiff);
/*
	CvScalar avg_s=cvAvg(vdiff);
	double avg=avg_s.val[0];
	int i;
	for (i=vrec->rows/2; i<vrec->rows; i++)
	{
		cvSetReal1D(vrec,i,cvGetReal1D(v,i));
	}
*/

	cvSet(mat,cvScalar(0));
	make2D(mat,roi,vrec);

	if (ofile!=NULL)
		writesfi(ofile,mat);
}

