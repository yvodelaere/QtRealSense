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

#ifndef FACEUT3D_INCLUDED
#define FACEUT3D_INCLUDED

extern int localverbose;

// global constants for registration
// resolution
const double resolution=1.5;

// use mask
const int ellipticalmask=0;

// postprocessing
const int holefilling=1;
const int spikeremoval=1;
const int backgroundremoval=1;
const int reflectionremoval=0;

// nose fitting method1=tip+bridge; method2=point below nose+dent above nose
const int nosefitmethod1=0;
const int nosefitmethod2=2;
const int LR=0;

// motion compensation
const double symmetrize=4.5;
const double motion_threshold=0;
const double maxshift=7.5;

extern int nearfrontal;
// end of caonstants for registration

using namespace utw3dface;

namespace utw3dface {
extern int verbose;
extern int debug;
extern int nthreads;
}

void readups(char *file, UnorderedPointSet &ups);
int make1D(RangeImage &ri,IplImage *mask,double *v);
int make1D(Register &reg,IplImage *mask,double *v);
void transform_im(RangeImage &ri,RangeImage &ro,double tilt_angle,double vertical_shift,double phi=0);

class Features
{
public:
	int nitems;
	int nfeatures;
	double **features;

	double *operator[](int i);
	Features(int nfeatures);
	int Add(double *v);
	int Write(FILE *f);
	int Write(char *filename);
	int ReadMax1000(FILE *f);
	int Read(FILE *f);
	int Read(char *filename);
	void Empty();
	~Features();
};

class Classifier
{
public:
	int nclassifiers;
	LDA *lda;
	int *dim;
	int *nlda;
	IplImage **maskim;
	double *threshold;
	double *score;
	int *regmode;
	char **id;
	int nfeatures;
	int maxdim;

	// make a templatevector of a set of registered images for gallery
	double *Enroll(RangeImage *reg,double *templatevector);

	// get templatevector of set of registered images for probe
	double *GetFeatures(RangeImage *reg,double *templatevector);

	double *Scores(double *probe,double *gallery);
	int Votes(double *scores);
	int Votes(double *probe,double *gallery);

	// for fine registration
	int Votes(Features &probe,double *gallery);

	// write classifier as a c-file
	int WriteC(char *file);

	int Read(char *file);
};

class ThreadedVotes
{
public:
	Classifier *mfr;
	Features *probe;
	Features *gallery;
	int *votes;
	int start,end;

	void Init(Classifier *mfr,Features *probe,Features *gallery,int *votes,int start,int end)
	{
		this->mfr=mfr;
		this->probe=probe;
		this->gallery=gallery;
		this->votes=votes;
		this->start=start;
		this->end=end;
	}
};

void *getscorethread(void *p);
int checksfi(char *filename);
void registerprobesfi(Classifier &mfr,char *probe_sfi1,char *probe_sfi2,Features &probe);
void registerprobe(Classifier &mfr,char *probe_abs,Features &probe);

#endif // FACEUT3D_INCLUDED
