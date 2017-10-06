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

#include <opencv2/opencv.hpp>
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
#include "builtinclassifier.h"

int localverbose=0;
int builtinclassifier=1;

#ifndef CLOCK_MONOTONIC
#define CLOCK_MONOTONIC 1

#include <sys/time.h>
int clock_gettime(int mode, struct timespec *t)
{
	struct timeval tp;
	struct timezone tz;
	gettimeofday(&tp,&tz);
	t->tv_sec=tp.tv_sec;
	t->tv_nsec=1000*tp.tv_usec;
	
	return 0;
}


#endif

const char *usage="usage: faceut3d [OPTIONS]\n\
OPTIONS:\n\
-c file   read classifier from .mfr file\n\
-C file   write C definition of classifier to file \n\
-g file   read gallery .abs file\n\
-G file   read list of gallery .abs files\n\
-m int    multithreaded calculations [nthreads=1]\n\
-N        handle non frontal images as well (slower)\n\
-p file   read probe .abs file\n\
-P file   read list of probe .abs files\n\
-r file   read features of gallery from file\n\
-v        be verbose\n\
-w file   write features of gallery as file\n\
";

char *default_probe=(char *) "probe.abs";
char *default_gallery=(char *) "gallery.abs";
char *default_classifier=(char *) "classifier.mfr";

int nearfrontal=1;

int main(int argc,char **argv)
{
	char *probe_abs=default_probe;
	char *gallery_abs=default_gallery;
	char *gallery_abslist=NULL;
	char *probe_abslist=NULL;
	char *classifier=default_classifier;
	char *classifier_Cfile=NULL;
	char *gallery_ofile=NULL;
	char *gallery_ifile=NULL;

	verbose=0;
	debug=0;
	clock_t t0,t1;
	struct timespec start, finish;

	int c;
	int i,j;

	while ((c=getopt(argc,argv,"P:m:G:w:?vp:g:c:r:C:N"))!=-1)
	{
		switch (c)
		{
			case 'N': nearfrontal=0; break;
			case 'm': sscanf(optarg,"%d",&nthreads); break;
			case 'p': probe_abs=optarg; break;
			case 'g': gallery_abs=optarg; break;
			case 'G': gallery_abslist=optarg; break;
			case 'P': probe_abslist=optarg; break;
			case 'c': classifier=optarg; break;
			case 'C': classifier_Cfile=optarg; break;
			case 'v': localverbose++; break;
			case 'w': gallery_ofile=optarg; break;
			case 'r': gallery_ifile=optarg; break;
			case '?':
			default: fprintf(stderr,usage); exit(1);
		}
	}

	if (nthreads<1 || nthreads>128)
	{
		fprintf(stderr,"number of threads must be in range 1..128\n");
		exit(1);
	}

	t0=clock();
	Classifier mfr;
	//mfr.Read(classifier);
	InitBuiltInClassifier(mfr);
	if (builtinclassifier==0)
		mfr.Read(classifier);
	t1=clock();
	if (localverbose)
		fprintf(stderr,"loading classifier took %.6f seconds\n",double(t1-t0)/CLOCKS_PER_SEC);
	
	Features gallery(mfr.nfeatures);
	Features probe(mfr.nfeatures);
	double templatevector[mfr.nfeatures];

	clock_gettime(CLOCK_MONOTONIC, &start);
	if (gallery_ifile!=NULL)
		gallery.Read(gallery_ifile);
	else if (gallery_abslist!=NULL)
	{
		FILE *f=fopen(gallery_abslist,"r");
		if (f==NULL)
		{
			fprintf(stderr,"Failed to open %s for input\n",gallery_abslist);
			exit(1);
		}
		
		char s[256],filename[256];
		int n=gallery.nitems;
		for (;;)
		{
			struct timespec start, finish;
			clock_gettime(CLOCK_MONOTONIC, &start);
			if (fgets(s,254,f)==NULL)
				break;
		
			sscanf(s,"%s",filename);
			UnorderedPointSet galleryups;
			readups(filename,galleryups);
			Register reggallery(galleryups,resolution,holefilling,spikeremoval,ellipticalmask,reflectionremoval,backgroundremoval,nosefitmethod1,LR,symmetrize,motion_threshold,maxshift,nearfrontal);
			RangeImage regg[2];
			regg[0].CopyDataFrom(reggallery.ri);
			reggallery.FitNose(2);
			reggallery.PostProc();
			regg[1].CopyDataFrom(reggallery.ri);
			regg[0].WriteSFI("output_orig");
			regg[1].WriteSFI("output_alt");

			mfr.Enroll(regg,templatevector);
			gallery.Add(templatevector);
			clock_gettime(CLOCK_MONOTONIC, &finish);
			if (localverbose)
				fprintf(stderr,"enrolling %s took %.6f seconds\n",filename,finish.tv_sec-start.tv_sec+(finish.tv_nsec-start.tv_nsec)/1e9);
		}
		fclose(f);
		if (localverbose)
			fprintf(stderr,"added %d items to the gallery\n",gallery.nitems-n);
	}
	else
	{
		UnorderedPointSet galleryups;
		readups(gallery_abs,galleryups);
		Register reggallery(galleryups,resolution,holefilling,spikeremoval,ellipticalmask,reflectionremoval,backgroundremoval,nosefitmethod1,LR,symmetrize,motion_threshold,maxshift,nearfrontal);
		RangeImage regg[2];
		regg[0].CopyDataFrom(reggallery.ri);
		reggallery.FitNose(2);
		reggallery.PostProc();
		regg[1].CopyDataFrom(reggallery.ri);
		regg[0].WriteSFI("output_orig");
		regg[1].WriteSFI("output_alt");

		mfr.Enroll(regg,templatevector);
		gallery.Add(templatevector);
	}

	clock_gettime(CLOCK_MONOTONIC, &finish);
	if (localverbose)
		fprintf(stderr,"reading gallery took %.6f seconds\n",finish.tv_sec-start.tv_sec+(finish.tv_nsec-start.tv_nsec)/1e9);


	char s[256],filename[256];
	FILE *f;
	if (probe_abslist!=NULL)
	{
		f=fopen(probe_abslist,"r");
		if (f==NULL)
		{
			fprintf(stderr,"Failed to open %s for input\n",probe_abslist);
			exit(1);
		}
	}
		
	for (;;)
	{
		if (probe_abslist!=NULL)
		{
			if (fgets(s,254,f)==NULL)
				break;
			sscanf(s,"%s",filename);
			probe_abs=filename;
		}

		clock_gettime(CLOCK_MONOTONIC, &start);
		if (checksfi(probe_abs))
		{
			// this is likely a list of sfi file 2 per line (normal and alt registration)
			char sfi1[256],sfi2[256];
			sscanf(s,"%s%s",sfi1,sfi2);
			registerprobesfi(mfr,sfi1,sfi2,probe);
		}
		else
			registerprobe(mfr,probe_abs,probe);
		clock_gettime(CLOCK_MONOTONIC, &finish);
		if (localverbose)
			fprintf(stderr,"Registration of probe %s took %.6f seconds\n",probe_abs,finish.tv_sec-start.tv_sec+(finish.tv_nsec-start.tv_nsec)/1e9);

		if (localverbose)
			fprintf(stderr,"# gallery=%d #probe=%d\n",gallery.nitems,probe.nitems);

/*
	double bestscores[mfr.nclassifiers];
	for (i=0mei<gallery.nitems; i++)
	{
		for (k=0; k<mfr.nclassifiers; k++)
			bestscores[k]=-100000;
		for (j=0; j<probe.nitems; j++)
		{
			mfr.Scores(probe[j],gallery[i]);
			for (k=0; k<mfr.nclassifiers; k++)
			{
				if (mfr.score[k]>bestscores[k])
					bestscores[k]=mfr.score[k];
			}
		}
		printf("%d ",mfr.Votes(bestscores));
	}
	printf("\n");
*/
		clock_gettime(CLOCK_MONOTONIC, &start);
	
		int votes[gallery.nitems];
		
		if (nthreads>1 && gallery.nitems>nthreads*2)
		{
			pthread_t thread[nthreads];
			ThreadedVotes tv[nthreads];
		
			int nt=gallery.nitems/nthreads;
			for (i=0; i<nthreads-1; i++)
			{
				tv[i].Init(&mfr,&probe,&gallery,votes,i*nt,(i+1)*nt);
				pthread_create(&thread[i],NULL,getscorethread,(void*)&tv[i]);
			}
			tv[nthreads-1].Init(&mfr,&probe,&gallery,votes,(nthreads-1)*nt,gallery.nitems);
			pthread_create(&thread[i],NULL,getscorethread,(void*)&tv[nthreads-1]);
	
			for (i=0; i<nthreads; i++)
				pthread_join(thread[i],NULL);
		}
		else
		{
			for (i=0; i<gallery.nitems; i++)
				votes[i]=mfr.Votes(probe,gallery[i]);
		}
		clock_gettime(CLOCK_MONOTONIC, &finish);
		if (localverbose)
			fprintf(stderr,"Calculation of %d scores took %.6f seconds\n",gallery.nitems,finish.tv_sec-start.tv_sec+(finish.tv_nsec-start.tv_nsec)/1e9);

		for (i=0; i<gallery.nitems; i++)
			printf("%d ",votes[i]);
		printf("\n");

		if (probe_abslist==NULL)
			break;

		probe.Empty();
	}
	if (probe_abslist!=NULL)
		fclose(f);

	if (gallery_ofile)
	{
		gallery.Write(gallery_ofile);
	}

	if (classifier_Cfile!=NULL)
	{
		if (builtinclassifier==0)
			mfr.WriteC(classifier_Cfile);
		else
		{
			fprintf(stderr,"-C option is disabled in this version!\n");
		}
	}

	return 0;
}
