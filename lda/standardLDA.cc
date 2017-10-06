#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

#include <cv.h>

#include "matrixio.h"

int verbose=0;

int preprocess(CvMat *Z,double a)
{
	int dZ=Z->rows;
	int nZ=Z->cols;

	CvMat *Zmax=cvCreateMat(1,nZ,CV_64F);

	int i;
	double min,max;
	CvMat Zcoli=cvMat(dZ,1,CV_64F);
	for (i=0; i<nZ; i++)
	{
		cvGetCol(Z,&Zcoli,i);
		cvMinMaxLoc(&Zcoli,&min,&max);
		cvSetReal1D(Zmax,i,max);
	}

	CvMat *Zmaxrep=cvCreateMat(dZ,nZ,CV_64F);
	cvRepeat(Zmax,Zmaxrep);
	CvMat *Zmaxis0=cvCreateMat(dZ,nZ,CV_64F);
	cvSub(Z,Zmaxrep,Zmaxis0);

	cvReleaseMat(&Zmax);
	cvReleaseMat(&Zmaxrep);

	CvMat *aminZ=cvCreateMat(dZ,nZ,CV_64F);
	cvSubRS(Zmaxis0,cvScalar(a),aminZ);
	cvReleaseMat(&Zmaxis0);
	
	cvLog(aminZ,Z);

	return 1;
}

int get_means_and_variations(CvMat *Z,CvMat *iZ,
	CvMat *&meanZ,CvMat *&meanClasses,CvMat *&Zzeromean,CvMat *&Zvariations)
{
	int dZ=Z->rows;
	int nZ=Z->cols;

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

char *usage="usage: standardLDA [OPTIONS]\n\
OPTIONS:\n\
\t-a file     target data file (MxN matrix, with M dim, N # samples)\n\
\t-A file     target id file (1xN matrix with class numbers staring from 1)\n\
\t-L basename load meanZ, T,Rho1 and D matrices from files <basename>_meanZ/T/Rho/D.ascii\n\
\t-l int      set number of LDA components (25)\n\
\t-o file     file to store scorematrix\n\
\t-P float    preprocess data using log(a-(Z-maxZ)) where a is the float (0=no preprocessing)\n\
\t-p int      set number of PCA components (100)\n\
\t-r          transpose the target and training data files\n\
\t-S basename save meanZ, T,Rho1 and D matrices in <basename>_meanZ/T/Rho/D.ascii\n\
\t-t file     training data file (MxN matrix, with M dim, N # samples)\n\
\t-T file     training id file (1xN matrix with class numbers staring from 1)\n\
\t-v          increase verbosity\n\
";

main(int argc,char **argv)
{
	char c;
	char *targetdatafile=NULL;
	char *targetidfile=NULL;
	char *trainingdatafile=NULL;
	char *trainingidfile=NULL;
	char *ofile=NULL;
	char *savebasename=NULL;
	char *loadbasename=NULL;
	int nPCA=100;
	int nLDA=25;
	int transpose=0;
	int time0=time(NULL);
	double aprep=0;

	while ((c=getopt(argc,argv,"P:o:L:S:rt:T:a:A:vp:l:"))!=-1)
	{
		switch (c)
		{
			case 'P': sscanf(optarg,"%lf",&aprep); break;
			case 't': trainingdatafile=optarg; break;
			case 'T': trainingidfile=optarg; break;
			case 'a': targetdatafile=optarg; break;
			case 'A': targetidfile=optarg; break;
			case 'v': verbose++; break;
			case 'p': sscanf(optarg,"%d",&nPCA); break;
			case 'l': sscanf(optarg,"%d",&nLDA); break;
			case 'r': transpose=1; break;
			case 'o': ofile=optarg; break;
			case 'L': loadbasename=optarg; break;
			case 'S': savebasename=optarg; break;
			default:
				fprintf(stderr,usage);
				exit(1);
		}
	}

	if (loadbasename==NULL && (trainingdatafile==NULL || trainingidfile==NULL))
	{
		fprintf(stderr,"must provide either -L or -t and -T options!\n");
		exit(1);
	}

	CvMat *meanZ,*T,*Rho1,*Ddiag;	// these are the final results of the training

	if (loadbasename!=NULL)
	{
		if (verbose)
			fprintf(stderr,"%d: reading meanZ, T, Rho1 and D matrices from files\n",time(NULL)-time0);
		char s[256];
		sprintf(s,"%s_meanZ.ascii",loadbasename);
        meanZ=utw3dface::read_matrix(s);
		sprintf(s,"%s_T.ascii",loadbasename);
        T=utw3dface::read_matrix(s);
		sprintf(s,"%s_Rho1.ascii",loadbasename);
        Rho1=utw3dface::read_matrix(s);
		sprintf(s,"%s_D.ascii",loadbasename);
        Ddiag=utw3dface::read_matrix(s);

		if (meanZ==NULL || T==NULL || Rho1==NULL || Ddiag==NULL)
		{
			fprintf(stderr,"failed to read meanZ, T, Rho1 and D matrices\n");
			exit(1);
		}

		nLDA=Ddiag->rows;
		if (verbose)
			fprintf(stderr,"%d: setting nLDA to %d\n",time(NULL)-time0,nLDA);
	}
	else // calculate from training set 
	{
		if (verbose)
			fprintf(stderr,"%d: reading training data\n",time(NULL)-time0);
        CvMat *iZ=utw3dface::read_matrix(trainingidfile);
        CvMat *Z=utw3dface::read_matrix(trainingdatafile,transpose);

		if (aprep>0)
		{
			if (verbose)
				fprintf(stderr,"%d: preprocessing training data\n",time(NULL)-time0);
			preprocess(Z,aprep);
		}

		int i;

		CvMat *Zzeromean,*meanClasses,*Zvariations;

		if (verbose)
			fprintf(stderr,"%d: starting get_means_and_variations ... ",time(NULL)-time0);
		get_means_and_variations(Z,iZ,meanZ,meanClasses,Zzeromean,Zvariations);
		if (verbose)
			fprintf(stderr,"finshed\n");

		int nZ=Z->cols;
		int dZ=Z->rows;

		cvReleaseMat(&Z);

		// PCA
		if (verbose)
			fprintf(stderr,"%d: starting PCA step ...",time(NULL)-time0);
	
		int nS1=(nZ>dZ) ? dZ : nZ;

		CvMat *U1T=cvCreateMat(nZ,dZ,CV_64F);
		CvMat *S1diag=cvCreateMat(nZ,1,CV_64F);
		cvSVD(Zzeromean,S1diag,U1T,NULL,CV_SVD_MODIFY_A | CV_SVD_U_T);
	
		CvMat T1=cvMat(nPCA,dZ,CV_64F);
		cvGetSubRect(U1T,&T1,cvRect(0,0,dZ,nPCA));

		if (verbose)
			fprintf(stderr,"finshed\n");

		// LDA
		if (verbose)
			fprintf(stderr,"%d: starting LDA step ...",time(NULL)-time0);

		double k=sqrt(nZ-1.0);
		CvMat *tmp1=cvCreateMat(nPCA,nPCA,CV_64F);
		cvSetZero(tmp1);
		for (i=0; i<nPCA; i++)
			cvSetReal2D(tmp1,i,i,k/cvGetReal1D(S1diag,i));
		CvMat *tmp2=cvCreateMat(nPCA,dZ,CV_64F);
		cvMatMul(tmp1,&T1,tmp2);
		cvReleaseMat(&tmp1);
		CvMat *Yvariations=cvCreateMat(nPCA,nZ,CV_64F);
		cvMatMul(tmp2,Zvariations,Yvariations);
		cvReleaseMat(&Zvariations);
	
		CvMat *U2T=cvCreateMat(nPCA,nPCA,CV_64F);
		CvMat *S2diag=cvCreateMat(nZ,1,CV_64F);
		cvSVD(Yvariations,S2diag,U2T,NULL,CV_SVD_MODIFY_A | CV_SVD_U_T);
		cvReleaseMat(&Yvariations);

		CvMat *Sigma2Ydiag=cvCreateMat(nZ,1,CV_64F);
		cvMul(S2diag,S2diag,Sigma2Ydiag,1.0/(nZ-1.0));
	
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

		cvReleaseMat(&S1diag);
		cvReleaseMat(&U1T);
		
		if (verbose)
			fprintf(stderr,"finshed\n");
	}

	if (savebasename!=NULL)
	{
		if (verbose)
			fprintf(stderr,"%d: writing meanZ, T, Rho and D matrices\n",time(NULL)-time0);

		char s[256];
		sprintf(s,"%s_meanZ.ascii",savebasename);
        utw3dface::write_matrix(s,meanZ);
		sprintf(s,"%s_T.ascii",savebasename);
        utw3dface::write_matrix(s,T);
		sprintf(s,"%s_Rho1.ascii",savebasename);
        utw3dface::write_matrix(s,Rho1);
		sprintf(s,"%s_D.ascii",savebasename);
        utw3dface::write_matrix(s,Ddiag);
	}

	if (ofile!=NULL)
	{
		if (verbose)
			fprintf(stderr,"%d: reading testing data\n",time(NULL)-time0);
        CvMat *iZTest=utw3dface::read_matrix(targetidfile);
        CvMat *ZTest=utw3dface::read_matrix(targetdatafile,transpose);

		if (aprep>0)
		{
			if (verbose)
				fprintf(stderr,"%d: preprocessing testing data\n",time(NULL)-time0);
			preprocess(ZTest,aprep);
		}

		// generate score matrix
		if (verbose)
			fprintf(stderr,"%d: generating score matrix\n",time(NULL)-time0);
		
		int nZTest=ZTest->cols;
		int dZ=ZTest->rows;

		// save time by precalculating lots of things
		CvMat ZTestcoli=cvMat(dZ,1,CV_64F);
		CvMat *diff=cvCreateMat(dZ,1,CV_64F);
		int i;
		for (i=0; i<nZTest; i++)
		{
			cvGetCol(ZTest,&ZTestcoli,i);
			cvSub(&ZTestcoli,meanZ,diff);
			cvCopy(diff,&ZTestcoli);
		}
		cvReleaseMat(&diff);
/*
		CvMat *ZTestzeromean=cvCreateMat(dZ,nZTest,CV_64F);
		CvMat *meanZrep=cvCreateMat(dZ,nZTest,CV_64F);
		cvRepeat(meanZ,meanZrep);
		cvSub(ZTest,meanZrep,ZTestzeromean);
		cvReleaseMat(&meanZrep);
		cvReleaseMat(&ZTest);
*/
		CvMat *ZTestzeromean=ZTest;

		CvMat *T_ZTestzeromean=cvCreateMat(nLDA,nZTest,CV_64F);
		cvMatMul(T,ZTestzeromean,T_ZTestzeromean);
		CvMat *Rho1_T_ZTestzeromean=cvCreateMat(nLDA,nZTest,CV_64F);
		cvMatMul(Rho1,T_ZTestzeromean,Rho1_T_ZTestzeromean);
		cvReleaseMat(&ZTestzeromean);

		CvMat *ScoreMatrix=cvCreateMat(nZTest,nZTest,CV_64F);
		cvSetZero(ScoreMatrix);

		int rtest,itest;
		CvMat *x=cvCreateMat(nLDA,1,CV_64F);
		CvMat *e=cvCreateMat(nLDA,1,CV_64F);
		CvMat *z=cvCreateMat(nLDA,1,CV_64F);
		CvMat *De=cvCreateMat(nLDA,1,CV_64F);
		for (rtest=0; rtest<nZTest; rtest++)
		{
			cvGetCol(Rho1_T_ZTestzeromean,x,rtest);
	
			for (itest=0; itest<nZTest; itest++)
			{
				cvGetCol(T_ZTestzeromean,z,itest);

				cvSub(z,x,e);
				
				cvMul(Ddiag,e,De);
				cvSetReal2D(ScoreMatrix,rtest,itest,-cvDotProduct(e,De)+cvDotProduct(z,z));

			}

			if (verbose>1)
				fprintf(stderr,"%d: finished row %d of %d\n",time(NULL)-time0,rtest,nZTest);
		}

		if (verbose)
			fprintf(stderr," finshed\n");

        utw3dface::write_beematrix(ofile,ScoreMatrix,"output/FRGC_Exp_2.0.3_Target.bxml.tmp","output/FRGC_Exp_2.0.3_Target.bxml.tmp");
	}
}

