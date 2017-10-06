#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#ifdef linux
#include <unistd.h>
#endif
#include <time.h>

#include <cv.h>

#include "matrixio.h"
#include "lda.h"

namespace utw3dface {

extern int verbose;

#define limit_dim_by_classes 0

LDA::LDA()
{
	meanZ=NULL;
	T=NULL;
	Tinv=NULL;
	Rho1=NULL;
	Ddiag=NULL;
	Rho1_T=NULL;
	T1=NULL;
	S1diag=NULL;
}

void LDA::cleanup()
{
	if (meanZ!=NULL)
		cvReleaseMat(&meanZ);
	if (T!=NULL)
		cvReleaseMat(&T);
	if (Tinv!=NULL)
		cvReleaseMat(&Tinv);
	if (Rho1!=NULL)
		cvReleaseMat(&Rho1);
	if (Ddiag!=NULL)
		cvReleaseMat(&Ddiag);
	if (Rho1_T!=NULL)
		cvReleaseMat(&Rho1_T);
	if (T1!=NULL)
		cvReleaseMat(&T1);
	if (S1diag!=NULL)
		cvReleaseMat(&S1diag);
}

LDA::~LDA()
{
	cleanup();
}

int cmp(const void *pd1,const void *pd2)
{
	double d1=*((double *)pd1);
	double d2=*((double *)pd2);

	return (d1>d2) ? 1 : (d1<d2) ? -1 : 0;
}

CvScalar Median(CvMat *A,CvMat *mask=NULL)
{
	CvScalar median=cvScalar(0.0);

	int i,j;
	int n=A->cols;
	double *a=new double[n];

	j=0;
	for (i=0; i<n; i++)
	{
		if (mask==NULL || cvGetReal1D(mask,i)!=0)
			a[j++]=cvGetReal1D(A,i);
	}

	if (j<3)
		median=cvAvg(A,mask);
	else
	{
		qsort(a,j,sizeof(double),cmp);
		int l=int((j-2)*0.5);
		int r=int((j+1)*0.5);

		median=cvScalar(0.5*(a[l]+a[r]));
	}

	delete [] a;

	return median;
}

// Robust covariance matrix determination based on MAD (Median Average Deviation
int MadCov(CvMat *Z,CvMat *C)
{
	int d=Z->rows;
	int n=Z->cols;
	int i,j,k;
	double copy[n];
	int l=int((n-2)*0.5);
	int r=int((n+1)*0.5);

	if (n<3)
	{
		fprintf(stderr,"MadCov needs >3 samples!\n");
		exit(1);
	}

FILE *f=fopen("3050.txt","w");
for (j=0; j<n; j++)
	fprintf(f,"%d %g %g\n",j,cvGetReal2D(Z,30,j),cvGetReal2D(Z,50,j));
fclose(f);
for (j=0; j<n; j++)
	copy[j]=cvGetReal2D(Z,d-1,j);
double min=copy[0];
double max=copy[0];
for (j=0; j<n; j++)
	if (copy[j]>max) max=copy[j]; else if (copy[j]<min) min=copy[j];
int hist[100];
for (j=0; j<100; j++)
	hist[j]=0;
for (j=0; j<n; j++)
{
	int index=int((copy[j]-min)/(max-min)*99);
	hist[index]++;
}
for (j=0; j<100; j++)
	printf("%d %g %g\n",j,j*(max-min)/99+min,hist[j]*100/((max-min)*n));

	for (i=0; i<d; i++)
	{
		for (j=0; j<n; j++)
			copy[j]=fabs(cvGetReal2D(Z,i,j));
		qsort(copy,n,sizeof(double),cmp);
		double m=0.5*1.4826*(copy[l]+copy[r]);
		cvSetReal2D(C,i,i,m*m);
	}

	for (i=0; i<d; i++)
	for (k=i+1; k<d; k++)
	{
		for (j=0; j<n; j++)
			copy[j]=fabs(cvGetReal2D(Z,i,j)+cvGetReal2D(Z,k,j));
		qsort(copy,n,sizeof(double),cmp);
		double c1=0.5*1.4826*(copy[l]+copy[r]);

		for (j=0; j<n; j++)
			copy[j]=fabs(cvGetReal2D(Z,i,j)-cvGetReal2D(Z,k,j));
		qsort(copy,n,sizeof(double),cmp);
		double c2=0.5*1.4826*(copy[l]+copy[r]);

		double c=0.25*(c1*c1-c2*c2);
		cvSetReal2D(C,i,k,c);
		cvSetReal2D(C,k,i,c);
	}

	return 1;

}

int Cov(CvMat *Z,CvMat *C)
{
	int d=Z->rows;
	int n=Z->cols;
	int i,j,k;

	for (i=0; i<d; i++)
	for (k=i; k<d; k++)
	{
		double c=0;
		for (j=0; j<n; j++)
			c+=cvGetReal2D(Z,i,j)*cvGetReal2D(Z,k,j);
		c/=n-1.0;
		cvSetReal2D(C,i,k,c);
		cvSetReal2D(C,k,i,c);
	}

	return 1;

}

// robust version using median instead of mean
int LDA::get_means_and_variationsR(CvMat *Z,CvMat *iZ,
	CvMat *&meanZ,CvMat *&meanClasses,CvMat *&Zzeromean,CvMat *&Zvariations)
{
	int dZ=Z->rows;
	int nZ=Z->cols;

	if (verbose)
		fprintf(stderr,"lda: getting means and variations\n");

	if (iZ->rows!=1 || iZ->cols!=nZ)
	{
		fprintf(stderr,"wrong dimensions for training data idvector and datamatrix, try -r\n");
		exit(1);
	}

	int i,c;

	int nC=0;
	double min,max;
	cvMinMaxLoc(iZ,&min,&max);
	nC=int(max);

	meanZ=cvCreateMat(dZ,1,CV_64F);
	meanClasses=cvCreateMat(dZ,nC,CV_64F);
	Zzeromean=cvCreateMat(dZ,nZ,CV_64F);
	Zvariations=cvCreateMat(dZ,nZ,CV_64F);
	
	CvMat Zrowi=cvMat(1,nZ,CV_64F);
	for (i=0; i<dZ; i++)
	{
		cvGetRow(Z,&Zrowi,i);
		cvSet1D(meanZ,i,Median(&Zrowi));
	}	

	CvMat *meanZrep=cvCreateMat(dZ,nZ,CV_64F);
	cvRepeat(meanZ,meanZrep);
	cvSub(Z,meanZrep,Zzeromean);
	cvReleaseMat(&meanZrep);

	CvMat *ix=cvCreateMat(1,nZ,CV_8U);
	CvMat *meanCrep=cvCreateMat(dZ,nZ,CV_64F);
	CvMat meanCreprowi=cvMat(1,nZ,CV_64F);
	for (c=0; c<nC; c++)
	{
		cvCmpS(iZ,c+1,ix,CV_CMP_EQ);
		for (i=0; i<dZ; i++)
		{
			cvGetRow(Z,&Zrowi,i);
			CvScalar mean=Median(&Zrowi,ix);
			//CvScalar mean=cvAvg(&Zrowi,ix);
			cvSet2D(meanClasses,i,c,mean);
			cvGetRow(meanCrep,&meanCreprowi,i);
			cvSet(&meanCreprowi,mean,ix);
		}
	}

	cvSub(Z,meanCrep,Zvariations);
	
	cvReleaseMat(&ix);
	cvReleaseMat(&meanCrep);

	return 1;
}

// classical version using non-robust mean
int LDA::get_means_and_variations(CvMat *Z,CvMat *iZ,
	CvMat *&meanZ,CvMat *&meanClasses,CvMat *&Zzeromean,CvMat *&Zvariations)
{
	int dZ=Z->rows;
	int nZ=Z->cols;

	if (verbose)
		fprintf(stderr,"lda: getting means and variations\n");

	if (iZ->rows!=1 || iZ->cols!=nZ)
	{
		fprintf(stderr,"wrong dimensions for training data idvector and datamatrix, try -r\n");
		exit(1);
	}

	int i,c;

	int nC=0;
	double min,max;
	cvMinMaxLoc(iZ,&min,&max);
	nC=int(max);

	meanZ=cvCreateMat(dZ,1,CV_64F);
	meanClasses=cvCreateMat(dZ,nC,CV_64F);
	Zzeromean=cvCreateMat(dZ,nZ,CV_64F);
	Zvariations=cvCreateMat(dZ,nZ,CV_64F);
	
	CvMat Zrowi=cvMat(1,nZ,CV_64F);
	for (i=0; i<dZ; i++)
	{
		cvGetRow(Z,&Zrowi,i);
		cvSet1D(meanZ,i,cvAvg(&Zrowi));
	}	

	CvMat *meanZrep=cvCreateMat(dZ,nZ,CV_64F);
	cvRepeat(meanZ,meanZrep);
	cvSub(Z,meanZrep,Zzeromean);
	cvReleaseMat(&meanZrep);

	CvMat *ix=cvCreateMat(1,nZ,CV_8U);
	CvMat *meanCrep=cvCreateMat(dZ,nZ,CV_64F);
	CvMat meanCreprowi=cvMat(1,nZ,CV_64F);
	for (c=0; c<nC; c++)
	{
		cvCmpS(iZ,c+1,ix,CV_CMP_EQ);
		for (i=0; i<dZ; i++)
		{
			cvGetRow(Z,&Zrowi,i);
			CvScalar mean=cvAvg(&Zrowi,ix);
			cvSet2D(meanClasses,i,c,mean);
			cvGetRow(meanCrep,&meanCreprowi,i);
			cvSet(&meanCreprowi,mean,ix);
		}
	}

	cvSub(Z,meanCrep,Zvariations);
	
	cvReleaseMat(&ix);
	cvReleaseMat(&meanCrep);

	return 1;
}

int LDA::load(const char *loadbasename)
{
	cleanup();
	if (verbose)
		fprintf(stderr,"lda: reading meanZ, T, Rho1 and D matrices from files\n");
	char s[256];
	sprintf(s,"%s_meanZ.ascii",loadbasename);
	meanZ=read_matrix(s);
	sprintf(s,"%s_T.ascii",loadbasename);
	T=read_matrix(s);
	sprintf(s,"%s_Rho1.ascii",loadbasename);
	Rho1=read_matrix(s);
	sprintf(s,"%s_D.ascii",loadbasename);
	Ddiag=read_matrix(s);
	sprintf(s,"%s_T1.ascii",loadbasename);
	T1=read_matrix(s);
	sprintf(s,"%s_S1diag.ascii",loadbasename);
	S1diag=read_matrix(s);

	if (meanZ==NULL || T==NULL || Rho1==NULL || Ddiag==NULL)
	{
		fprintf(stderr,"failed to read meanZ, T, Rho1 and D matrices\n");
		exit(1);
	}
	if (T1==NULL || S1diag==NULL)
	{
		fprintf(stderr,"warning: failed to read T1 and/or S1diag matrix\n");
		exit(1);
	}

	Rho1_T=cvCreateMat(Rho1->rows,T->cols,CV_64F);
	cvMatMul(Rho1,T,Rho1_T);

	return 1;
}


// Robust version using median instead of mean and median of absolute deviation estimate for covariance
int LDA::trainR(CvMat *Z,CvMat *iZ,int nPCA,int nLDA)
{
	cleanup();

	int i;

	CvMat *Zzeromean,*meanClasses,*Zvariations;

	get_means_and_variationsR(Z,iZ,meanZ,meanClasses,Zzeromean,Zvariations);

	int nZ=Z->cols;
	int dZ=Z->rows;

	// determine the covariance matrix from Z using the MAD (median average deviation)
	if (verbose)
		fprintf(stderr,"lda: determine robust Covariance\n");
	CvMat *C=cvCreateMat(dZ,dZ,CV_64F);	
	CvMat *C1=cvCreateMat(dZ,dZ,CV_64F);	
	MadCov(Zzeromean,C);
	write_matrix("MC.ascii",C);
	//Cov(Zzeromean,C);
	//write_matrix("C.ascii",C);

	// PCA
	if (verbose)
		fprintf(stderr,"lda: starting PCA step\n");

	int nS1=dZ;
	CvMat *U1T;
	U1T=cvCreateMat(dZ,dZ,CV_64F);	
	S1diag=cvCreateMat(dZ,1,CV_64F);
	cvSVD(C,S1diag,U1T,NULL,CV_SVD_MODIFY_A | CV_SVD_U_T);
	// to make singular values from the eigenvalues again:
	for (i=0; i<nS1; i++)
		cvSetReal1D(S1diag,i,sqrt(cvGetReal1D(S1diag,i)*(nZ-1.0)));
	
	if (nPCA>S1diag->rows)
	{
		fprintf(stderr,"# PCA components should be <= # training vectors(%d) and # dimensions (%d)!\n",nZ,dZ);
		exit(1);
	}

	T1=cvCreateMat(nPCA,dZ,CV_64F);
	cvGetSubRect(U1T,T1,cvRect(0,0,dZ,nPCA)); // note: cvRect takes first width then height!

	// LDA
	if (verbose)
		fprintf(stderr,"lda: starting LDA step\n");

	double k=sqrt(nZ-1.0);
	CvMat *tmp1=cvCreateMat(nPCA,nPCA,CV_64F);
	cvSetZero(tmp1);
	for (i=0; i<nPCA; i++)
		cvSetReal2D(tmp1,i,i,k/cvGetReal1D(S1diag,i));
	CvMat *tmp2=cvCreateMat(nPCA,dZ,CV_64F);
	cvMatMul(tmp1,T1,tmp2);
	cvReleaseMat(&tmp1);
	CvMat *Yvariations=cvCreateMat(nPCA,nZ,CV_64F);
	cvMatMul(tmp2,Zvariations,Yvariations);
	cvReleaseMat(&Zvariations);
//!!! hier ook Mad!	
	CvMat *LC=cvCreateMat(nPCA,nPCA,CV_64F);	
	MadCov(Yvariations,LC);
	write_matrix("MLC.ascii",C);
	//Cov(Yvariations,LC);
	//write_matrix("LC.ascii",C);
	
	CvMat *Sigma2Ydiag;
	CvMat *U2T;
	CvMat *S2diag;
	if (nZ>=nPCA)
	{
		U2T=cvCreateMat(nPCA,nPCA,CV_64F);
		S2diag=cvCreateMat(nPCA,1,CV_64F);
		cvSVD(LC,S2diag,U2T,NULL,CV_SVD_MODIFY_A | CV_SVD_U_T);
		for (i=0; i<nPCA; i++)
			cvSetReal1D(S2diag,i,sqrt(cvGetReal1D(S2diag,i)*(nZ-1.0)));

		Sigma2Ydiag=cvCreateMat(nPCA,1,CV_64F);
		cvMul(S2diag,S2diag,Sigma2Ydiag,1.0/(nZ-1.0));
	}
	else
	{
		fprintf(stderr,"nPCA must be >= nZ\n");
		exit(1);
	}
	cvReleaseMat(&Yvariations);

	CvMat T2=cvMat(nLDA,nPCA,CV_64F);
	cvGetSubRect(U2T,&T2,cvRect(0,nPCA-nLDA,nPCA,nLDA));

	T=cvCreateMat(nLDA,dZ,CV_64F);
	cvMatMul(&T2,tmp2,T);
	cvReleaseMat(&U2T);
	cvReleaseMat(&tmp2);

	Rho1=cvCreateMat(nLDA,nLDA,CV_64F);
	cvSetZero(Rho1);
	for (i=0; i<nLDA; i++)
		cvSetReal2D(Rho1,i,i,1-cvGetReal1D(Sigma2Ydiag,i+nPCA-nLDA));
	
	Ddiag=cvCreateMat(nLDA,1,CV_64F);
	for (i=0; i<nLDA; i++)
	{
		double el=cvGetReal1D(Sigma2Ydiag,i+nPCA-nLDA);
		cvSetReal1D(Ddiag,i,1/((2-el)*el));
	}

	cvReleaseMat(&S2diag);
	cvReleaseMat(&Sigma2Ydiag);

	Rho1_T=cvCreateMat(nLDA,dZ,CV_64F);
	cvMatMul(Rho1,T,Rho1_T);

	return 1;
}		

// Classical version using non robust estimators for mean etc.
int LDA::train(CvMat *Z,CvMat *iZ,int nPCA,int nLDA)
{
	cleanup();

	int i;

	CvMat *Zzeromean,*meanClasses,*Zvariations;

	get_means_and_variations(Z,iZ,meanZ,meanClasses,Zzeromean,Zvariations);

	int nZ=Z->cols;
	int dZ=Z->rows;
	int nC=meanClasses->cols;

	// PCA
	if (verbose)
		fprintf(stderr,"lda: starting PCA step %d dimensions, %d samples, %d classes\n",dZ,nZ,nC);

	int nS1=(nZ>dZ) ? dZ : nZ;
	CvMat *U1T;
	if (nZ>dZ)	// more samples than dimensions
	{
		U1T=cvCreateMat(dZ,dZ,CV_64F);	// note nZ x dZ gives weird results!
		S1diag=cvCreateMat(dZ,1,CV_64F);
		cvSVD(Zzeromean,S1diag,U1T,NULL,CV_SVD_MODIFY_A | CV_SVD_U_T);
/*
		CvMat *ZzeromeanT=cvCreateMat(nZ,dZ,CV_64F);
		cvTranspose(Zzeromean,ZzeromeanT);
		U1T=cvCreateMat(dZ,dZ,CV_64F);
		S1diag=cvCreateMat(dZ,1,CV_64F);
		cvSVD(ZzeromeanT,S1diag,NULL,U1T,CV_SVD_MODIFY_A | CV_SVD_V_T);
		cvReleaseMat(&ZzeromeanT);
*/
	}
	else	// more dimensions than samples
	{
		U1T=cvCreateMat(nZ,dZ,CV_64F); // note: dZ x dZ gives weird results!
		S1diag=cvCreateMat(nZ,1,CV_64F);
		cvSVD(Zzeromean,S1diag,U1T,NULL,CV_SVD_MODIFY_A|CV_SVD_U_T);
	}
	//write_matrix("foo.ascii",U1T);
	
	if (nPCA>S1diag->rows)
	{
		fprintf(stderr,"# PCA components should be <= # training vectors(%d) and # dimensions (%d)!\n",nZ,dZ);
		exit(1);
	}

/*
	//regularisation
	double smallest_ne0=cvGetReal1D(S1diag,0);
	for (i=1; i<S1diag->rows; i++)
	{
		if (cvGetReal1D(S1diag,i)!=0)
			smallest_ne0=cvGetReal1D(S1diag,i);
		else
			cvSetReal1D(S1diag,i,smallest_ne0/1000);
	}
*/

	T1=cvCreateMat(nPCA,dZ,CV_64F);
	cvGetSubRect(U1T,T1,cvRect(0,0,dZ,nPCA)); // note: cvRect takes first width then height!

	// LDA
	if (verbose)
		fprintf(stderr,"lda: starting LDA step\n");

	double k=sqrt(nZ-1.0);
	CvMat *tmp1=cvCreateMat(nPCA,nPCA,CV_64F);
	cvSetZero(tmp1);
	for (i=0; i<nPCA; i++)
		cvSetReal2D(tmp1,i,i,k/cvGetReal1D(S1diag,i));
	CvMat *tmp2=cvCreateMat(nPCA,dZ,CV_64F);
	cvMatMul(tmp1,T1,tmp2);
	cvReleaseMat(&tmp1);
	CvMat *Yvariations=cvCreateMat(nPCA,nZ,CV_64F);
	cvMatMul(tmp2,Zvariations,Yvariations);
	cvReleaseMat(&Zvariations);
	
	CvMat *Sigma2Ydiag;
	CvMat *U2T;
	CvMat *S2diag;
	if (nZ>=nPCA)
	{
		U2T=cvCreateMat(nPCA,nPCA,CV_64F);
		S2diag=cvCreateMat(nPCA,1,CV_64F);
		cvSVD(Yvariations,S2diag,U2T,NULL,CV_SVD_MODIFY_A | CV_SVD_U_T);

		Sigma2Ydiag=cvCreateMat(nPCA,1,CV_64F);
		cvMul(S2diag,S2diag,Sigma2Ydiag,1.0/(nZ-1.0));
	}
	else
	{
		fprintf(stderr,"nPCA must be >= nZ\n");
		exit(1);
	}
	cvReleaseMat(&Yvariations);

	CvMat T2=cvMat(nLDA,nPCA,CV_64F);
	if (limit_dim_by_classes)
		//cvGetSubRect(U2T,&T2,cvRect(0,nC-nLDA-1,nPCA,nLDA));
		cvGetSubRect(U2T,&T2,cvRect(0,0,nPCA,nLDA));
	else
		cvGetSubRect(U2T,&T2,cvRect(0,nPCA-nLDA,nPCA,nLDA));

	T=cvCreateMat(nLDA,dZ,CV_64F);
	cvMatMul(&T2,tmp2,T);
	cvReleaseMat(&U2T);
	cvReleaseMat(&tmp2);

	Rho1=cvCreateMat(nLDA,nLDA,CV_64F);
	cvSetZero(Rho1);
	for (i=0; i<nLDA; i++)
		if (limit_dim_by_classes)
			//cvSetReal2D(Rho1,i,i,1-cvGetReal1D(Sigma2Ydiag,i+nC-1-nLDA));
			cvSetReal2D(Rho1,i,i,1-cvGetReal1D(Sigma2Ydiag,i));
		else
			cvSetReal2D(Rho1,i,i,1-cvGetReal1D(Sigma2Ydiag,i+nPCA-nLDA));
	
	Ddiag=cvCreateMat(nLDA,1,CV_64F);

	double el;
	for (i=0; i<nLDA; i++)
	{
		if (limit_dim_by_classes)
			//el=cvGetReal1D(Sigma2Ydiag,i+nC-1-nLDA);
			el=cvGetReal1D(Sigma2Ydiag,i);
		else
			el=cvGetReal1D(Sigma2Ydiag,i+nPCA-nLDA);
		cvSetReal1D(Ddiag,i,1/((2-el)*el));
	}

	cvReleaseMat(&S2diag);
	cvReleaseMat(&Sigma2Ydiag);

	Rho1_T=cvCreateMat(nLDA,dZ,CV_64F);
	cvMatMul(Rho1,T,Rho1_T);

	return 1;
}		

int LDA::retrainLDA(CvMat *Z,CvMat *iZ,int nPCA,int nLDA)
{
	int i;

	if (T1==NULL)
	{
		fprintf(stderr,"Cannot retrain if not trained before; use with -r option to read pretrained classifier\n");
		exit(1);
	}
	if (nPCA!=T1->rows)
	{
		fprintf(stderr,"For retraining use same #PCA and #LDA as for original training!\n");
		exit(1);
	}

	CvMat *Zzeromean,*meanClasses,*Zvariations;

	get_means_and_variations(Z,iZ,meanZ,meanClasses,Zzeromean,Zvariations);

	int nZ=Z->cols;
	int dZ=Z->rows;
	int nC=meanClasses->cols;

	// LDA
	if (verbose)
		fprintf(stderr,"lda: starting LDA step\n");

	double k=sqrt(nZ-1.0);
	CvMat *tmp1=cvCreateMat(nPCA,nPCA,CV_64F);
	cvSetZero(tmp1);
	for (i=0; i<nPCA; i++)
		cvSetReal2D(tmp1,i,i,k/cvGetReal1D(S1diag,i));
	CvMat *tmp2=cvCreateMat(nPCA,dZ,CV_64F);
	cvMatMul(tmp1,T1,tmp2);
	cvReleaseMat(&tmp1);
	CvMat *Yvariations=cvCreateMat(nPCA,nZ,CV_64F);
	cvMatMul(tmp2,Zvariations,Yvariations);
	cvReleaseMat(&Zvariations);
	
	CvMat *Sigma2Ydiag;
	CvMat *U2T;
	CvMat *S2diag;
	if (nZ>=nPCA)
	{
		U2T=cvCreateMat(nPCA,nPCA,CV_64F);
		S2diag=cvCreateMat(nPCA,1,CV_64F);
		cvSVD(Yvariations,S2diag,U2T,NULL,CV_SVD_MODIFY_A | CV_SVD_U_T);

		Sigma2Ydiag=cvCreateMat(nPCA,1,CV_64F);
		cvMul(S2diag,S2diag,Sigma2Ydiag,1.0/(nZ-1.0));
	}
	else
	{
		fprintf(stderr,"nPCA must be >= nZ\n");
		exit(1);
	}
	cvReleaseMat(&Yvariations);

	CvMat T2=cvMat(nLDA,nPCA,CV_64F);
	if (limit_dim_by_classes)
		//cvGetSubRect(U2T,&T2,cvRect(0,nC-nLDA,nPCA,nLDA));
		cvGetSubRect(U2T,&T2,cvRect(0,0,nPCA,nLDA));
	else
		cvGetSubRect(U2T,&T2,cvRect(0,nPCA-nLDA,nPCA,nLDA));


	if (T==NULL)
		T=cvCreateMat(nLDA,dZ,CV_64F);
	if (nLDA!=T->rows)
	{
		fprintf(stderr,"For retraining use same #PCA and #LDA as for original training!\n");
		exit(1);
	}

	cvMatMul(&T2,tmp2,T);
	cvReleaseMat(&U2T);
	cvReleaseMat(&tmp2);

	if (Rho1==NULL)
		Rho1=cvCreateMat(nLDA,nLDA,CV_64F);
	cvSetZero(Rho1);
	for (i=0; i<nLDA; i++)
		if (limit_dim_by_classes)
			//cvSetReal2D(Rho1,i,i,1-cvGetReal1D(Sigma2Ydiag,i+nC-nLDA));
			cvSetReal2D(Rho1,i,i,1-cvGetReal1D(Sigma2Ydiag,i));
		else
			cvSetReal2D(Rho1,i,i,1-cvGetReal1D(Sigma2Ydiag,i+nPCA-nLDA));
	
	if (Ddiag==NULL)
		Ddiag=cvCreateMat(nLDA,1,CV_64F);

	double el;
	for (i=0; i<nLDA; i++)
	{
		if (limit_dim_by_classes)
			//el=cvGetReal1D(Sigma2Ydiag,i+nC-nLDA);
			el=cvGetReal1D(Sigma2Ydiag,i);
		else
			el=cvGetReal1D(Sigma2Ydiag,i+nPCA-nLDA);
		cvSetReal1D(Ddiag,i,1/((2-el)*el));
	}

	cvReleaseMat(&S2diag);
	cvReleaseMat(&Sigma2Ydiag);

	if (Rho1_T==NULL)
		Rho1_T=cvCreateMat(nLDA,dZ,CV_64F);
	cvMatMul(Rho1,T,Rho1_T);

	return 1;
}		

int LDA::save(const char *savebasename)
{
	if (verbose)
		fprintf(stderr,"lda: writing meanZ, T, Rho and D matrices\n");

	char s[256];
	sprintf(s,"%s_meanZ.ascii",savebasename);
	write_matrix(s,meanZ);
	sprintf(s,"%s_T.ascii",savebasename);
	write_matrix(s,T);
	sprintf(s,"%s_Rho1.ascii",savebasename);
	write_matrix(s,Rho1);
	sprintf(s,"%s_D.ascii",savebasename);
	write_matrix(s,Ddiag);
	sprintf(s,"%s_T1.ascii",savebasename);
	write_matrix(s,T1);
	sprintf(s,"%s_S1diag.ascii",savebasename);
	write_matrix(s,S1diag);

	return 1;
}

double readval(FILE *f)
{
	char s[256];
	double val;

	do
	{
		if (fscanf(f,"%s",s)==EOF)
		{
			fprintf(stderr,"lda: failed reading classifier from truncated file\n");
			exit(1);
		}

		if (s[0]=='#')
		{
			// skip to end of line
			int c;
			for (;;)
			{
				c=fgetc(f);
				if (c=='\n')
					break;
				if (c==EOF)
				{
					fprintf(stderr,"lda: failed reading classifier from truncated file\n");
					exit(1);
				}
			}
			
		}
	} while (s[0]=='#');

	if (sscanf(s,"%lf",&val) != 1)
	{
		fprintf(stderr,"lda: failed reading classifier from file: expected floating point number\n");
		exit(1);
	}

	return val;
}

int LDA::load(FILE *f)
{
	char s[256];
	double val;
	int dZ,nPCA,nLDA,i,j;

	do { fgets(s,255,f); } while (s[0]=='#'); 
	sscanf(s,"%d%d%d",&dZ,&nPCA,&nLDA);
	if (dZ<1 || dZ>1000000 || nPCA<1 || nPCA>dZ || nLDA<1 || nLDA>nPCA)
	{
		fprintf(stderr,"Failed reading classifier from file: illegal dZ, nPCA or nLDA\n");
		exit(1);
	}

	cleanup();
	if (verbose)
		fprintf(stderr,"lda: reading meanZ, T, Rho1 and D matrices from file\n");

	meanZ=cvCreateMat(dZ,1,CV_64F);
	for (i=0; i<dZ; i++)
		cvSetReal1D(meanZ,i,readval(f));

	Rho1=cvCreateMat(nLDA,nLDA,CV_64F);
	cvSetZero(Rho1);
	for (i=0; i<nLDA; i++)
		cvSetReal2D(Rho1,i,i,readval(f));

	Ddiag=cvCreateMat(nLDA,1,CV_64F);
	for (i=0; i<nLDA; i++)
		cvSetReal1D(Ddiag,i,readval(f));

	T=cvCreateMat(nLDA,dZ,CV_64F);
	for (i=0; i<T->rows; i++)
	for (j=0; j<T->cols; j++)
		cvSetReal2D(T,i,j,readval(f));

	// skip to endofline
	while (fgetc(f)!='\n');

	// create dummy matrices
	T1=cvCreateMatHeader(nPCA,dZ,CV_64F);
	S1diag=cvCreateMatHeader(dZ,1,CV_64F);

	Rho1_T=cvCreateMat(Rho1->rows,T->cols,CV_64F);
	cvMatMul(Rho1,T,Rho1_T);

	return 1;
}

// save classifier to 1 file, but only mean, rho1, D and T
int LDA::save(FILE *f)
{
	int i,j;

	fprintf(f,"# classifier definition\n"); 
	fprintf(f,"# dimension input vector nPCA and nLDA\n");
	fprintf(f,"%d %d %d\n",meanZ->rows,T1->rows,T->rows);

	fprintf(f,"# meanZ\n");
	for (i=0; i<meanZ->rows; i++)
		fprintf(f,"%g ",cvGetReal1D(meanZ,i));
	fprintf(f,"\n");

	fprintf(f,"# Rho1\n");
	for (i=0; i<Rho1->rows; i++)
		fprintf(f,"%g ",cvGetReal2D(Rho1,i,i));
	fprintf(f,"\n");

	fprintf(f,"# Ddiag\n");
	for (i=0; i<Ddiag->rows; i++)
		fprintf(f,"%g ",cvGetReal1D(Ddiag,i));
	fprintf(f,"\n");

	fprintf(f,"# T\n");
	for (i=0; i<T->rows; i++)
	{
		for (j=0; j<T->cols; j++)
			fprintf(f,"%g ",cvGetReal2D(T,i,j));
		fprintf(f,"\n");
	}

	return 1;
}

// fast enroll not using opencv stuff; v_enrolled must be allocated beforehand
void LDA::enroll(double *v,double *v_enrolled)
{
	int nLDA=T->rows;
	int dZ=meanZ->rows;
	double v_zeromean[dZ];
	double *pmean=(double *)cvPtr1D(meanZ,0);

	int i,j;
	for (i=0; i<dZ; i++)
		v_zeromean[i]=v[i]-pmean[i];

	for (i=0; i<nLDA; i++)
	{
		double sum=0;
		double *pRho1_T=(double *)cvPtr2D(Rho1_T,i,0);
		for (j=0; j<dZ; j++)
			sum+=pRho1_T[j]*v_zeromean[j];
		v_enrolled[i]=sum;
	}
}

CvMat *LDA::enroll(CvMat *v)
{
	int nLDA=T->rows;
	int dZ=meanZ->rows;
	if (v->rows!=dZ)
	{
		fprintf(stderr,"lda: trying to enroll vector with wrong dimension\n");
		exit(1);
	}

	CvMat *v_zeromean=cvCreateMat(dZ,1,CV_64F);
	cvSub(v,meanZ,v_zeromean);

	CvMat *v_enrolled=cvCreateMat(nLDA,1,CV_64F);
	cvMatMul(Rho1_T,v_zeromean,v_enrolled);

	cvReleaseMat(&v_zeromean);

	return v_enrolled;
}

// fast features not using opencv stuff; v_features must be allocated beforehand
void LDA::features(double *v,double *v_enrolled)
{
	int nLDA=T->rows;
	int dZ=meanZ->rows;
	double v_zeromean[dZ];
	double *pmean=(double *)cvPtr1D(meanZ,0);

	int i,j;
	for (i=0; i<dZ; i++)
		v_zeromean[i]=v[i]-pmean[i];

	for (i=0; i<nLDA; i++)
	{
		double sum=0;
		double *pT=(double *)cvPtr2D(T,i,0);
		for (j=0; j<dZ; j++)
			sum+=pT[j]*v_zeromean[j];
		v_enrolled[i]=sum;
	}
}

CvMat *LDA::features(CvMat *v)
{
	int nLDA=T->rows;
	int dZ=meanZ->rows;
	if (v->rows!=dZ)
	{
		fprintf(stderr,"lda: trying to calculate features of vector with wrong dimension\n");
		exit(1);
	}

	CvMat *v_zeromean=cvCreateMat(dZ,1,CV_64F);
	cvSub(v,meanZ,v_zeromean);

	CvMat *v_features=cvCreateMat(nLDA,1,CV_64F);
	cvMatMul(T,v_zeromean,v_features);

	cvReleaseMat(&v_zeromean);

	return v_features;
}

double LDA::likelihoodratio(double *reference,double *test)
{
	int i;
	int nLDA=T->rows;
	double *pD=Ddiag->data.db,e,score=0;

	for (i=0; i<nLDA; i++)
	{
		e=test[i]-reference[i];
		score+=test[i]*test[i]-e*pD[i]*e;
	}

	return score;
}

double LDA::likelihoodratio(CvMat *reference,CvMat *test)
{
	int nLDA=T->rows;
	int dZ=meanZ->rows;
	CvMat *enrolled,*test_features;

	if (reference->rows==nLDA)
		enrolled=reference;
	else
		enrolled=enroll(reference);

	if (test->rows==nLDA)
		test_features=test;
	else
		test_features=features(test);

	CvMat *e=cvCreateMat(nLDA,1,CV_64F);
	CvMat *De=cvCreateMat(nLDA,1,CV_64F);

	cvSub(test_features,enrolled,e);
	cvMul(Ddiag,e,De);
	double score=-cvDotProduct(e,De)+cvDotProduct(test_features,test_features);

	if (enrolled!=reference)
		cvReleaseMat(&enrolled);
	if (test_features!=test)
		cvReleaseMat(&test_features);
	cvReleaseMat(&e);
	cvReleaseMat(&De);

	return score;
}

double LDA::DIFS(CvMat *v)
{
	int nLDA=T->rows;

	CvMat *v_features;
	
	if (v->rows==nLDA)
		v_features=v;
	else
		v_features=features(v);

	double difs=cvNorm(v_features)/nLDA;

	return difs;
}

double LDA::DFFS(CvMat *v)
{
	int nLDA=T->rows;
	int dZ=meanZ->rows;
	if (v->rows!=dZ)
	{
		fprintf(stderr,"lda: trying to calculate features of vector with wrong dimension\n");
		exit(1);
	}

	CvMat *v_features;
	v_features=features(v);

	if (Tinv==NULL)
	{
		Tinv=cvCreateMat(dZ,nLDA,CV_64F);
		cvInvert(T,Tinv,CV_SVD);
	}

	CvMat *v_zeromean=cvCreateMat(dZ,1,CV_64F);
	cvSub(v,meanZ,v_zeromean);

	CvMat *v_reconstructed=cvCreateMat(dZ,1,CV_64F);
	cvMatMul(Tinv,v_features,v_reconstructed);

	double dffs=cvNorm(v,v_reconstructed)/dZ;

	return dffs;
}
}
