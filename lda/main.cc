/*
 * 04-01-2008: fixed a rather nasty bug, now reads id file for training with other than FRGC data
 * 15-04-2009: added difs threhold
 * 11-10-2012: added robust training using median (doesn't work very well)
 * 19-10-2012: added subtract average as option
 * 26-3-2013: added suppport for reading multiple registrations in 1 using -J option
 */

#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <libgen.h>
#include <X11/Xarch.h>

//#include <cv.h>
#include <highgui.h>

#include "matrixio.h"
#include "lda.h"

using namespace utw3dface;

namespace utw3dface {
int verbose=0;
}

typedef struct {char a,b,c,d;} EL;
typedef union {float f; EL elem; } NUM;

// flags that determine if derivative of input image is used

int dX=0;
int dY=0;
int subtractAverage=0;

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

int writesfi(const char *filename,CvMat *im)
{
        FILE *f;
        f=fopen(filename,"wb");
        if (f==NULL)
        {
                fprintf(stderr,"failed to open %s\n",filename);
                exit(1);
        }

        fprintf(f,"%s %d %d %d\n","CSU_SFI",im->cols,im->rows,1);

        NUM junk;
        int x,y;
        for(y=0;y<im->rows;y++)
        for(x=0;x<im->cols;x++)
        {
                junk.f=(float)cvGetReal2D(im,y,x);

                if (BYTE_ORDER == LITTLE_ENDIAN)
                {
                        fwrite(&(junk.elem.d),1,1,f);
                        fwrite(&(junk.elem.c),1,1,f);
                        fwrite(&(junk.elem.b),1,1,f);
                        fwrite(&(junk.elem.a),1,1,f);
                }
                else
                {
                        fwrite(&(junk.elem.a),1,1,f);
                        fwrite(&(junk.elem.b),1,1,f);
                        fwrite(&(junk.elem.c),1,1,f);
                        fwrite(&(junk.elem.d),1,1,f);
                }
        }

        fclose(f);

        return 1;
}

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
		if (BYTE_ORDER == LITTLE_ENDIAN)
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

	if (subtractAverage==1)
	{
		double average=0;
		for(j=0;j<height;j++)
		for(i=0;i<width;i++)
			average+=cvGetReal2D(mat,j,i);
		average/=height*width;
		for(j=0;j<height;j++)
		for(i=0;i<width;i++)
			cvSetReal2D(mat,j,i,cvGetReal2D(mat,j,i)-average);
	}

	if (dX)
	{
		CvMat *dX=cvCreateMat(height,width,CV_32F);
		for(j=1;j<height-1;j++)
		for(i=1;i<width-1;i++)
			cvSetReal2D(dX,j,i,cvGetReal2D(mat,j,i+1)-cvGetReal2D(mat,j,i-1));
		for(j=1;j<height-1;j++)
		for(i=1;i<width-1;i++)
			cvSetReal2D(mat,j,i,cvGetReal2D(dX,j,i));
		for(j=0;j<height;j++)
		{
			cvSetReal2D(mat,j,0,0.0);
			cvSetReal2D(mat,j,width-1,0.0);
		}
		for(i=0;i<width;i++)
		{
			cvSetReal2D(mat,0,i,0.0);
			cvSetReal2D(mat,height-1,i,0.0);
		}
		cvReleaseMat(&dX);
	}
	if (dY)
	{
		CvMat *dY=cvCreateMat(height,width,CV_32F);
		for(j=1;j<height-1;j++)
		for(i=1;i<width-1;i++)
			cvSetReal2D(dY,j,i,cvGetReal2D(mat,j+1,i)-cvGetReal2D(mat,j-1,i));
		for(j=1;j<height-1;j++)
		for(i=1;i<width-1;i++)
			cvSetReal2D(mat,j,i,cvGetReal2D(dY,j,i));
		for(j=0;j<height;j++)
		{
			cvSetReal2D(mat,j,0,0.0);
			cvSetReal2D(mat,j,width-1,0.0);
		}
		for(i=0;i<width;i++)
		{
			cvSetReal2D(mat,0,i,0.0);
			cvSetReal2D(mat,height-1,i,0.0);
		}
	}

	return mat;
}

CvMat *transform_im(CvMat *im,double tilt_angle,double vertical_shift)
{
	double ny=85.0/130.0*(im->rows-1);
	double cosa=cos(tilt_angle);
	double sina=sin(tilt_angle);

	CvMat *result=cvCreateMat(im->rows,im->cols,CV_32F);
	
	int x,y;
	int y1f[im->rows];
	for (x=0; x<im->cols; x++)
	{
		for (y=0; y<im->rows; y++)
			y1f[y]=0;
		for (y=0; y<im->rows; y++)
		{
			double z=cvGetReal2D(im,y,x);
			int y1=int((y-ny-vertical_shift)*cosa-z*sina+ny);
			double z1=z*cosa+(y-ny-vertical_shift)*sina;
	
			if (y1>=0 && y1<im->rows)
			{
				cvSetReal2D(result,y1,x,z1);
				y1f[y1]=1;
			}
		}
		for (y=0; y<im->rows; y++)
			if (y1f[y]==0)
			{
				double val=0;
				int d=0,n=0;
				for (d=1; d<im->rows; d++)
				{
					if ((y-d)>=0 && y1f[y-d]==1)
					{
						val+=cvGetReal2D(result,y-d,x);
						n++;
					}
					if ((y+d)<im->rows && y1f[y+d]==1)
					{
						val+=cvGetReal2D(result,y+d,x);
						n++;
					}
					if (n>0)
						break;
				}
				cvSetReal2D(result,y,x,val/n);
			}
	}

	return result;
}

CvMat *rotate_z(CvMat *im,double angle)
{
	CvMat *result=cvCreateMat(im->rows,im->cols,CV_32F);
	CvMat* rot=cvCreateMat(2,3,CV_32F);
	cv2DRotationMatrix(cvPoint2D32f(0.5*(im->cols-1),0.5*(im->rows-1)),angle,1.0,rot);
	cvWarpAffine(im,result,rot);

	cvReleaseMat(&rot);

	return result;
}

CvMat *rotate_x(CvMat *im,double max_dz)
{
	CvMat *result=cvCreateMat(im->rows,im->cols,CV_32F);
	int i,j;
	for (i=0; i<im->rows; i++)
	{
		double dz=(i-0.5*im->rows)/(0.5*im->rows)*max_dz;
		for (j=0; j<im->cols; j++)
			cvSetReal2D(result,i,j,cvGetReal2D(im,i,j)+dz);
	}

	return result;
}

CvMat *rotate_y(CvMat *im,double max_dz)
{
	CvMat *result=cvCreateMat(im->rows,im->cols,CV_32F);
	int i,j;
	for (j=0; j<im->cols; j++)
	{
		double dz=(i-0.5*im->cols)/(0.5*im->cols)*max_dz;
		for (i=0; i<im->rows; i++)
			cvSetReal2D(result,i,j,cvGetReal2D(im,i,j)+dz);
	}

	return result;
}

CvMat *translate_z(CvMat *im,double dz)
{
	CvMat *result=cvCreateMat(im->rows,im->cols,CV_32F);
	int i,j;
	for (i=0; i<im->rows; i++)
	for (j=0; j<im->cols; j++)
		cvSetReal2D(result,i,j,cvGetReal2D(im,i,j)+dz);

	return result;
}

CvMat *make1D(CvMat *mat,Roi &roi,IplImage *maskim,CvMat *v=NULL)
{
	int width=mat->cols;
	int h_min=int(mat->rows*roi.bottom);
	int h_max=int(mat->rows*roi.top);
	int w_min=int(mat->cols*roi.left);
	int w_max=int(mat->cols*roi.right);
	double *d=new double[mat->cols*mat->rows];

	if (maskim!=NULL && (mat->cols!=maskim->width || mat->rows!=maskim->height))
	{
		fprintf(stderr,"mask does not match image size\n");
		exit(1);
	}

	double mean=0,l;
	int i,j,k=0;

	if (maskim!=NULL)
	{
		for (i=h_min; i<h_max; i++)
		{
			unsigned char *pmask=cvPtr2D(maskim,i,w_min);
			float *pmat=(float *) cvPtr2D(mat,i,w_min);
			for (j=w_min; j<w_max; j++)
			{
	/*
				if (maskim!=NULL && cvGetReal2D(maskim,i,j)==0)
					continue;
				d[k++]=cvGetReal2D(mat,i,j);
	*/
	
				if (pmask[j])
				{
					d[k++]=pmat[j];
					mean+=pmat[j];
				}
			}
		}
	}
	else
	{
		for (i=h_min; i<h_max; i++)
		{
			float *pmat=(float *) cvPtr2D(mat,i,w_min);
			for (j=w_min; j<w_max; j++)
			{
				mean+=pmat[j];
				d[k++]=pmat[j];
			}
		}
	}
	mean/=(k>0) ? k : 1;

	if (subtractAverage==2)
	{
		for (i=0; i<k; i++)
			d[k]-=mean;
	}

	if (v==NULL)
	{
		if (maskim)
			v=cvCreateMat(k,1,CV_64F);
		else	
			v=cvCreateMat((h_max-h_min)*(w_max-w_min),1,CV_64F);
	}

	if (maskim!=0)
	{
		if (k!=v->cols*v->rows)
		{
			fprintf(stderr,"feature vector size mismatch\n");
			exit(1);
		}
	}
	else if (v->cols*v->rows!=(h_max-h_min)*(w_max-w_min))
	{
		fprintf(stderr,"feature vector size mismatch\n");
		exit(1);
	}

/*
	// for some reason this doesn't work if v is not NULL as passed to function!
	double *pv=v->data.db;
	for (i=0; i<k; i++)
		pv[i]=d[i];
*/
	//double *pv=(double *) cvPtr1D(v,0);
	//for (i=0; i<k; i++)
		//pv[i]=d[i];
	for (i=0; i<k; i++)
		cvSetReal1D(v,i,d[i]);
	
	delete [] d;

	return v;	
}

CvMat *readsfi1D(char *fname,Roi &roi,IplImage *maskim,CvMat *v=NULL)
{
	CvMat *mat=readsfi(fname);
	v=make1D(mat,roi,maskim,v);

	cvReleaseMat(&mat);

	return v;	
}

void train(LDA &lda,int nPCA,int nLDA,char *flist,Roi &roi,IplImage *maskim,char *fidlist=NULL,int robust=0,int retrainlda=0)
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
	CvMat *im=readsfi1D(fname,roi,maskim);
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
		readsfi1D(fname,roi,maskim,&col);

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

	if (retrainlda)
		lda.retrainLDA(Z,iZ,nPCA,nLDA);
	else if (robust)
		lda.trainR(Z,iZ,nPCA,nLDA);
	else
		lda.train(Z,iZ,nPCA,nLDA);

	cvReleaseMat(&Z);
	cvReleaseMat(&iZ);
}

CvMat *getscorematrix(LDA &lda,char *flist,Roi &roi,IplImage *maskim,int use_difs_threshold=0,double difs_threshold=0.0,int joinmultiple=0)
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
	int *gallery_difs_accept=new int[ll];

	for (i=0; i<ll; i++)
	{
		if (verbose)
			fprintf(stderr,"%d ",i);
		fscanf(f,"%s",fname);
		CvMat *v;
		if (joinmultiple)
		{
			CvMat *mat=readsfi(fname);
			int w=mat->cols/5;
			int h=mat->rows/5;
			CvMat *submat=cvCreateMat(h,w,CV_32F);
			CvRect rect=cvRect(2*w,2*h,w,h);
			cvGetSubRect(mat,submat,rect);
			v=make1D(submat,roi,maskim);
			cvReleaseMat(&mat);
			cvReleaseMat(&submat);
		}
		else
			v=readsfi1D(fname,roi,maskim);

		gallery[i]=lda.enroll(v);
		if (use_difs_threshold)
			gallery_difs_accept[i]=(lda.DIFS(v)<difs_threshold);
		cvReleaseMat(&v);
	}

	CvMat *scorematrix=cvCreateMat(ll,ll,CV_64F);
	if (joinmultiple)
	{
		for (i=0; i<ll; i++)
		for (j=0; j<ll; j++)
			cvSetReal2D(scorematrix,i,j,-100000.0);
	}

	if (verbose)
		fprintf(stderr,"generating %d x %d score matrix\n",ll,ll);

	rewind(f);

	CvMat *probe;
	int probe_difs_accept;
	for (i=0; i<ll; i++)
	{
		if (verbose)
			fprintf(stderr,"%d ",i);
		fscanf(f,"%s",fname);

		if (joinmultiple)
		{
			CvMat *mat=readsfi(fname);
			int w=mat->cols/5;
			int h=mat->rows/5;
			CvMat *submat=cvCreateMat(h,w,CV_32F);

			int ih,iw;
			for (ih=0; ih<5; ih++)
			for (iw=0; iw<5; iw++)
			{
				CvRect rect=cvRect(iw*w,ih*h,w,h);
				cvGetSubRect(mat,submat,rect);
				CvMat *v=make1D(submat,roi,maskim);
	
				CvMat *probe=lda.features(v);
				cvReleaseMat(&v);

				for (j=0; j<ll; j++)
				{
					double likelihoodratio=lda.likelihoodratio(gallery[j]->data.db,probe->data.db);
					double maxlikelihoodratio=cvGetReal2D(scorematrix,i,j);
					if (likelihoodratio>maxlikelihoodratio)
						cvSetReal2D(scorematrix,i,j,likelihoodratio);
				}
				cvReleaseMat(&probe);
			}

			cvReleaseMat(&mat);
			cvReleaseMat(&submat);
		}
		else
		{
			CvMat *v=readsfi1D(fname,roi,maskim);

			probe=lda.features(v);
			if (use_difs_threshold)
			probe_difs_accept=(lda.DIFS(v)<difs_threshold);
			if (verbose>1)
			{
				for (j=0; j<probe->rows; j++)
					printf("%g ",cvGetReal1D(probe,j));
				printf("\n");
			}
			cvReleaseMat(&v);
			if (use_difs_threshold)
			{
				double likelihoodratio;
				for (j=0; j<ll; j++)
				{
					if (probe_difs_accept && gallery_difs_accept[j])
	//					likelihoodratio=lda.likelihoodratio(gallery[j],probe);
						likelihoodratio=lda.likelihoodratio(gallery[j]->data.db,probe->data.db);
					else
						likelihoodratio=-10000;
					cvSetReal2D(scorematrix,i,j,likelihoodratio);
				}
			}
			else
			{
				for (j=0; j<ll; j++)
				{
	//				cvSetReal2D(scorematrix,i,j,lda.likelihoodratio(gallery[j],probe));
					cvSetReal2D(scorematrix,i,j,lda.likelihoodratio(gallery[j]->data.db,probe->data.db));
	
/*
// test code using cos
			int k;
			double score=0,na=0,nb=0,ak,bk,ed=0;
			for (k=0; k<lda.T->rows; k++)
			{
				ak=gallery[j]->data.db[k];
				bk=probe->data.db[k];
				score+=ak*bk;
				na+=ak*ak;
				nb+=bk*bk;
			}
			score/=sqrt(na*nb);
			cvSetReal2D(scorematrix,i,j,score);
*/
				}
			}
			cvReleaseMat(&probe);
		}
	}

	for (i=0; i<ll; i++)
		cvReleaseMat(&gallery[i]);

	delete [] gallery;
	delete [] gallery_difs_accept;

	fclose(f);

	return scorematrix;
}

CvMat *getscorematrix_finereg(LDA &lda,char *flist,Roi &roi,IplImage *maskim,int use_difs_threshold=0,double difs_threshold=0.0)
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
	int *gallery_difs_accept=new int[ll];

	for (i=0; i<ll; i++)
	{
		if (verbose)
			fprintf(stderr,"%d ",i);
		fscanf(f,"%s",fname);
		CvMat *v=readsfi1D(fname,roi,maskim);
		gallery[i]=lda.enroll(v);
		if (use_difs_threshold)
			gallery_difs_accept[i]=(lda.DIFS(v)<difs_threshold);
		cvReleaseMat(&v);
	}

	CvMat *scorematrix=cvCreateMat(ll,ll,CV_64F);
	for (i=0; i<ll; i++)
	for (j=0; j<ll; j++)
			cvSetReal2D(scorematrix,i,j,-100000.0);

	if (verbose)
		fprintf(stderr,"generating %d x %d score matrix\n",ll,ll);

	rewind(f);

	CvMat *probe;
	int probe_difs_accept;
	for (i=0; i<ll; i++)
	{
		if (verbose)
			fprintf(stderr,"%d ",i);
		fscanf(f,"%s",fname);
		CvMat* imin=readsfi(fname);
		double tilt=0,shift=0;
		for (tilt=-0.04; tilt<=0.04; tilt+=0.02)
		for (shift=-3; shift<=3; shift+=1)
		{
			CvMat *imout=transform_im(imin,tilt,shift);
			CvMat *v=make1D(imout,roi,maskim);

			probe=lda.features(v);
			cvReleaseMat(&imout);

			if (use_difs_threshold)
				probe_difs_accept=(lda.DIFS(v)<difs_threshold);
			cvReleaseMat(&v);

			double likelihoodratio;
			if (use_difs_threshold)
			{
				for (j=0; j<ll; j++)
				{
					if (probe_difs_accept && gallery_difs_accept[j])
						likelihoodratio=lda.likelihoodratio(gallery[j]->data.db,probe->data.db);
					else
						likelihoodratio=-10000;
					if (likelihoodratio>cvGetReal2D(scorematrix,i,j))
						cvSetReal2D(scorematrix,i,j,likelihoodratio);
				}
			}
			else
			{
				for (j=0; j<ll; j++)
				{
					likelihoodratio=lda.likelihoodratio(gallery[j]->data.db,probe->data.db);
					if (likelihoodratio>cvGetReal2D(scorematrix,i,j))
						cvSetReal2D(scorematrix,i,j,likelihoodratio);
				}
			}
			cvReleaseMat(&probe);
		}
		cvReleaseMat(&imin);
	}

	for (i=0; i<ll; i++)
		cvReleaseMat(&gallery[i]);

	delete [] gallery;
	delete [] gallery_difs_accept;

	fclose(f);

	return scorematrix;
}


// just compare 2 images and get the score; first in enroll
double getscore(LDA &lda,char *im1,char *im2,Roi &roi,IplImage *maskim,int joinmultiple=0)
{
	double likelihoodratio;

	if (joinmultiple)
	{
		double maxlikelihoodratio=-1000000;

		// get gallery
	        CvMat *mat=readsfi(im1);
		int w=mat->cols/5;
		int h=mat->rows/5;
		CvMat *submat=cvCreateMat(h,w,CV_32F);
		CvRect rect=cvRect(2*w,2*h,w,h);
		cvGetSubRect(mat,submat,rect);
        	CvMat *v=make1D(submat,roi,maskim);
		CvMat *gallery=lda.enroll(v);
		cvReleaseMat(&v);
        	cvReleaseMat(&mat);

		// get multiple probes and pick one with highest likelihoodratio
	        mat=readsfi(im2);

		int ih,iw;
		for (ih=0; ih<5; ih++)
		for (iw=0; iw<5; iw++)
		{
			CvRect rect=cvRect(iw*w,ih*h,w,h);
			cvGetSubRect(mat,submat,rect);
        		v=make1D(submat,roi,maskim);
	
			CvMat *probe=lda.features(v);
			cvReleaseMat(&v);
	
			likelihoodratio=lda.likelihoodratio(gallery,probe);
			if (likelihoodratio>maxlikelihoodratio)
				maxlikelihoodratio=likelihoodratio;

fprintf(stderr,"%d %d %g\n",ih,iw,likelihoodratio);
			cvReleaseMat(&probe);
		}
        	cvReleaseMat(&submat);
        	cvReleaseMat(&mat);
		cvReleaseMat(&gallery);

		likelihoodratio=maxlikelihoodratio;
	}
	else
	{
		CvMat *v=readsfi1D(im1,roi,maskim);
		CvMat *gallery=lda.enroll(v);
		cvReleaseMat(&v);
	
	
		v=readsfi1D(im2,roi,maskim);
		CvMat *probe=lda.features(v);
		cvReleaseMat(&v);
	
		likelihoodratio=lda.likelihoodratio(gallery,probe);

		cvReleaseMat(&probe);
		cvReleaseMat(&gallery);
	}

	return likelihoodratio;
}

// just compare 2 images and get the score; first in enroll and do fineregistration on 2nd image
double getscore_finereg(LDA &lda,char *im1,char *im2,Roi &roi,IplImage *maskim)
{
	CvMat *v=readsfi1D(im1,roi,maskim);
	CvMat *gallery=lda.enroll(v);
	cvReleaseMat(&v);

	double maxlr=-1000000;
	CvMat* imin=readsfi(im2);
	double tilt=0,shift=0,d=0;
	for (tilt=-0.04; tilt<=0.04; tilt+=0.02)
	for (shift=-3; shift<=3; shift+=1)
	//for (d=-3; d<=3; d+=1)
	{
		CvMat *imout=transform_im(imin,tilt,shift);
		v=make1D(imout,roi,maskim,v);
	//	int i;
	//	for (i=0; i<v->rows; i++)
	//		cvSetReal1D(v,i,cvGetReal1D(v,i)+d);

		CvMat *probe=lda.features(v);
		cvReleaseMat(&v);
		cvReleaseMat(&imout);

		double likelihoodratio=lda.likelihoodratio(gallery,probe);
		if (maxlr<likelihoodratio)
			maxlr=likelihoodratio;
		fprintf(stdout,"%g %g %g\n",tilt,shift,likelihoodratio);
	}

	cvReleaseMat(&imin);

	return maxlr;
}

void dumpfeatures(char *ffeatures,LDA &lda,char *flist,Roi &roi,IplImage *maskim)
{
	FILE *ff;
	ff=fopen(ffeatures,"w");
	if (ff==NULL)
	{
		fprintf(stderr,"cannot open %s for writing\n",ffeatures);
		exit(1);
	}

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

	CvMat *probe;
	for (i=0; i<ll; i++)
	{
		if (verbose)
			fprintf(stderr,"%d ",i);
		fscanf(f,"%s",fname);
		CvMat *v=readsfi1D(fname,roi,maskim);
		probe=lda.features(v);
		fprintf(ff,"%s ",fname);
		for (j=0; j<probe->rows; j++)
			fprintf(ff,"%g ",cvGetReal1D(probe,j));
		fprintf(ff,"\n");
		cvReleaseMat(&v);
		cvReleaseMat(&probe);
	}

	fclose(ff);
	fclose(f);
}

void dumpdifsdffs(char *fdifsdffs,LDA &lda,char *flist,Roi &roi,IplImage *maskim)
{
	FILE *ff;
	ff=fopen(fdifsdffs,"w");
	if (ff==NULL)
	{
		fprintf(stderr,"cannot open %s for writing\n",fdifsdffs);
		exit(1);
	}

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

	CvMat *probe;
	for (i=0; i<ll; i++)
	{
		if (verbose)
			fprintf(stderr,"%d ",i);
		fscanf(f,"%s",fname);
		CvMat *v=readsfi1D(fname,roi,maskim);
		double difs=lda.DIFS(v);
		double dffs=lda.DFFS(v);
		fprintf(ff,"%s %g %g\n",fname,difs,dffs);
		cvReleaseMat(&v);
	}

	fclose(ff);
	fclose(f);
}

const char *usage="usage: lda [OPTIONS] [probe-im target-im] \n\
OPTIONS:\n\
\t-A int      subtract average from input image (0=not; 1=average over whole image; 2=over masked part) [0]\n\
\t-a file     list of target images (first five digets=id)\n\
\t-b          use robust training\n\
\t-c file     load binary classifier file\n\
\t-d file     dump difs and dffs to file\n\
\t-D float    set a threshold to difs (not during training) [default none]\n\
\t-f file     dump features to file\n\
\t-F          perform fine registration by shifting and tilting probe image\n\
\t-i file     list of id's of training images as numbers from 1..99999 on one long line\n\
\t-J          assume input image is a 5x5 mosaic of input images with fine variations of registration\n\
\t-l int      set number of LDA components (25)\n\
\t-m file     read a mask file to use as ROI in PGM format \n\
\t-o file     file to store scorematrix\n\
\t-p int      set number of PCA components (100)\n\
\t-r basename read meanZ, T,Rho1 and D matrices from files <basename>_meanZ/T/Rho/D.ascii\n\
\t-s file     save classifier to binary file\n\
\t-t file     list of training images (first 5 digits=id or use -i for id-file)\n\
\t-v	      increase verbosity\n\
\t-w basename write meanZ, T,Rho1 and D matrices in <basename>_meanZ/T/Rho/D.ascii\n\
\t-z          retrain lda (use together with -r and -t)\n\
\t-B float    bottom of roi of image [0..top> [0.0]\n\
\t-L float    left of roi of image [0..right> [0.0]\n\
\t-R float    right of roi of image <left..1] [1.0]\n\
\t-T float    top of roi of image <bottom..1] [1.0]\n\
\t-X          take the x-derivative of the input image\n\
\t-Y          take the y-derivative of the input image\n\
";

main(int argc,char **argv)
{
	char c;
	char *targetlist=NULL;
	char *traininglist=NULL;
	char *idlist=NULL;
	char *ofile=NULL;
	char *ffeatures=NULL;
	char *fdifsdffs=NULL;
	char *savebasename=NULL;
	char *savebinfile=NULL;
	char *loadbinfile=NULL;
	char *loadbasename=NULL;
	int nPCA=100;
	int nLDA=25;
	double top=1.0;
	double bottom=0.0;
	double left=0.0;
	double right=1.0;
	char *maskfile=NULL;
	IplImage *maskim=NULL;
	int use_difs_threshold=0;
	double difs_threshold=1.0;
	int robust=0;
	int retrainlda=0;
	int fineregistration=0;
	int joinmultiple=0;

	while ((c=getopt(argc,argv,"JFA:zbXYD:d:m:f:o:r:w:t:a:vp:l:B:T:L:R:i:s:c:"))!=-1)
	{
		switch (c)
		{
			case 'A': sscanf(optarg,"%d",&subtractAverage); break;
			case 'B': sscanf(optarg,"%lf",&bottom); break;
			case 'D': use_difs_threshold=1; sscanf(optarg,"%lf",&difs_threshold); break;
			case 'F': fineregistration=1; break;
			case 'J': joinmultiple=1; break;
			case 'L': sscanf(optarg,"%lf",&left); break;
			case 'R': sscanf(optarg,"%lf",&right); break;
			case 'T': sscanf(optarg,"%lf",&top); break;
			case 'X': dX=1; break;
			case 'Y': dY=1; break;
			case 'a': targetlist=optarg; break;
			case 'b': robust=1; break;
			case 'c': loadbinfile=optarg; break;
			case 'd': fdifsdffs=optarg; break;
			case 'f': ffeatures=optarg; break;
			case 'i': idlist=optarg; break;
			case 'l': sscanf(optarg,"%d",&nLDA); break;
			case 'm': maskfile=optarg; break;
			case 'o': ofile=optarg; break;
			case 'p': sscanf(optarg,"%d",&nPCA); break;
			case 'r': loadbasename=optarg; break;
			case 's': savebinfile=optarg; break;
			case 't': traininglist=optarg; break;
			case 'v': verbose++; break;
			case 'w': savebasename=optarg; break;
			case 'z': retrainlda=1; break;
			default:
				fprintf(stderr,usage);
				exit(1);
		}
	}

	Roi roi(left,right,bottom,top);

	if (loadbinfile==NULL && loadbasename==NULL && traininglist==NULL)
	{
		fprintf(stderr,"must provide either -r or -t option!\n");
		exit(1);
	}

	if (maskfile!=NULL)
	{
		maskim=cvLoadImage(maskfile,CV_LOAD_IMAGE_GRAYSCALE);
		if (maskim==NULL)
		{
			fprintf(stderr,"failed to read mask image %s\n",maskfile);
			exit(1);
		}
	}

	LDA lda;
	if (loadbasename!=NULL)
		lda.load(loadbasename);

	if (loadbinfile!=NULL)
	{
		FILE *f=fopen(loadbinfile,"r");
		if (f==NULL)
		{
			fprintf(stderr,"Could not open %s for input\n",loadbinfile);
			exit(1);
		}
		lda.load(f);
		nPCA=lda.T1->rows;
		nLDA=lda.T->rows;
		char s[4096];
		int width,height;
		int i,j;
		do 
		{ 
			fgets(s,4095,f);
		} while (s[0]=='#');
		sscanf(s,"%d%d\n",&width,&height);

		if (verbose)
			fprintf(stderr,"reading mask image %dx%d\n",width,height);
		if (width==0 && height==0)
			maskim=NULL;
		else
		{
			maskim=cvCreateImage(cvSize(width,height),IPL_DEPTH_8U,1);
			for (i=0; i<height; i++)
			{
				fgets(s,4095,f);
				for (j=0; j<width; j++)
					cvSetReal2D(maskim,i,j,(s[j]=='0')?0:1);
			}
		}
		fclose(f);
	}

	if ((loadbasename==NULL && loadbinfile==NULL) || retrainlda)
	{
		if (idlist==NULL)
			fprintf(stderr,"WARNING: extracting id's from filenames, only correct for FRGC data!\n");
		train(lda,nPCA,nLDA,traininglist,roi,maskim,idlist,robust,retrainlda);
	}

	if (savebasename!=NULL)
		lda.save(savebasename);

	if (savebinfile!=NULL)
	{
		FILE *f=fopen(savebinfile,"w");
		if (f==NULL)
		{
			fprintf(stderr,"Could not open %s for output\n",savebinfile);
			exit(1);
		}
		lda.save(f);

		fprintf(f,"# Mask image\n");
		if (maskim==NULL)
			fprintf(f,"%d %d\n",0,0);
		else
		{
			fprintf(f,"%d %d\n",maskim->width,maskim->height);
			int i,j;
			for (i=0; i<maskim->height; i++)
			{
				for (j=0; j<maskim->width; j++)
					fprintf(f,"%d",(cvGetReal2D(maskim,i,j)!=0));
				fprintf(f,"\n");
			}
			fclose(f);
		}
	}

	if (ffeatures!=NULL)
		dumpfeatures(ffeatures,lda,targetlist,roi,maskim);

	if (fdifsdffs!=NULL)
		dumpdifsdffs(fdifsdffs,lda,targetlist,roi,maskim);

	if (ofile!=NULL && targetlist!=NULL)
	{
		CvMat *ScoreMatrix;
		clock_t t0=clock();
		if (fineregistration)
			ScoreMatrix=getscorematrix_finereg(lda,targetlist,roi,maskim,use_difs_threshold,difs_threshold);
		else
			ScoreMatrix=getscorematrix(lda,targetlist,roi,maskim,use_difs_threshold,difs_threshold,joinmultiple);
		clock_t t1=clock();
		if (verbose)
			fprintf(stderr,"calculation of score matrix (%d scores) took %.6f seconds\n",ScoreMatrix->rows*ScoreMatrix->cols,double(t1-t0)/CLOCKS_PER_SEC);
		write_beematrix(ofile,ScoreMatrix,(char *) "output/FRGC_Exp_2.0.3_Target.bxml.tmp",(char *) "output/FRGC_Exp_2.0.3_Target.bxml.tmp");
	}

	if (optind+1<argc)	// two more arg => just compare and report score
	{
		double score;

		if (fineregistration)
			score=getscore_finereg(lda,argv[optind],argv[optind+1],roi,maskim);
		else
			score=getscore(lda,argv[optind],argv[optind+1],roi,maskim,joinmultiple);
		if (verbose)
			fprintf(stderr,"%s %s %g\n",argv[optind],argv[optind+1],score);
		else
			fprintf(stderr,"%g\n",score);
	}
	
}

