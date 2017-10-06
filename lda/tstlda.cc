#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <libgen.h>
#include <X11/Xarch.h>

#include <cv.h>
#include <highgui.h>

#include "matrixio.h"
#include "lda.h"

using namespace utw3dface;

namespace utw3dface {
int verbose=0;
}

main(int argc,char **argv)
{
	int l=4;
	int p=7;

	int dim=10;
	int nsamples=1000;

	int i,j;

	CvMat *Z=cvCreateMat(dim,nsamples,CV_64F);
	CvMat col=cvMat(Z->rows,1,CV_64F);
	FILE *f=fopen("trainingset.txt","r");
	double x;
	for (i=0; i<nsamples; i++)
	{
		cvGetCol(Z,&col,i);
		for (j=0; j<dim; j++)
		{
			fscanf(f,"%lf",&x);
			cvSetReal1D(&col,j,x);
		}
	}
	fclose(f);

	CvMat *iZ=cvCreateMat(1,nsamples,CV_64F);
	f=fopen("id.txt","r");
	int id;
	for (i=0; i<nsamples; i++)
	{
		fscanf(f,"%d",&id);
		cvSetReal1D(iZ,i,id);
	}
	fclose(f);
	
	LDA lda;
	lda.train(Z,iZ,p,l);

	CvMat *Zt=cvCreateMat(dim,nsamples,CV_64F);
	f=fopen("testset.txt","r");
	for (i=0; i<nsamples; i++)
	{
		cvGetCol(Zt,&col,i);
		for (j=0; j<dim; j++)
		{
			fscanf(f,"%lf",&x);
			cvSetReal1D(&col,j,x);
		}
	}
	fclose(f);

	for (i=0; i<nsamples; i++)
	{
		cvGetCol(Zt,&col,i);
		CvMat *gallery=lda.enroll(&col);
		for (j=0; j<nsamples; j++)
		{
			cvGetCol(Zt,&col,j);
			CvMat *probe=lda.features(&col);
	
			double likelihoodratio=lda.likelihoodratio(gallery,probe);
			printf("%g ",likelihoodratio);
			cvReleaseMat(&probe);
		}
		printf("\n");
		cvReleaseMat(&gallery);
	}
}

