////////////////////////////////////////////////////////////////////////
//
// This class is the main registration class
// It takes a point cloud in the form of an UnorederdPointSet as input and delivers a registered, 
// resampled range image as output
// The regisration parameters: rotation matrix and translation can be obtained using the
// methods GetRotationMatrix and GetTranslationVector
// It also provides methods to obtain a profile.
// Furthermore it stores some quality indications of the registration process:
//	spikesQuality - fraction of surface points which are not spikes
//	holesQuality - fraction of surface points which are not small holes
//	registrationQuality - indication of quality of first estimate of symmetry axis and nose
//	noiseQuality - fraction of surface points which are not noise
//	completenessQuality - fraction of surface points which are not big holes
//	noseQuality - quality of nose detection
//	overallQuality - product of all above quality factors
// All of the quality numbers have a range of [0..1], although, because of poor normalisation, it is possible that
// sometimes values slightly outside this range are reported.
//
// Author: Luuk Spreeuwers, University of Twente
// Date: 2005-2011
// Last update: 28-08-2012
// version: 3.0
//
////////////////////////////////////////////////////////////////////////

#ifndef REGISTER_INCLUDED
#define REGISTER_INCLUDED

#include "unorderedpointset.h"
#include "rangeimage.h"

namespace utw3dface {

class Register
{
public:
	UnorderedPointSet ups;
	RangeImage ri;
	RangeImage rirroi;
	
	// parameters that control registration; these are initialised in Init()
	double course_registration_resolution;
	int course_registration_ri_width;
	double resolution;
	int holefilling;
	int spikeremoval;
	int ellipticalmask;
	int reflectionremoval;
	int backgroundremoval;
	int nosefitmethod;
	int LR; // fit nose left, right or both
	double symmetrize;
	double motion_threshold;
	double maxshift;
	int nearfrontal;
	
	// nose pars
	double nx,ny,nphi,nl,nh,nt;

	// methods to get profile store the results in here
	UnorderedPointSet profile;
	UnorderedPointSet maxprofile;

	// quality measures; set to 0 in Init()
	double spikesQuality;
	double holesQuality;
	double registrationQuality;
	double noiseQuality;
	double completenessQuality;
	double overallQuality;
	double noseQuality;

	// initialisation
	void Init();
	// actual registration
	void DoReg(UnorderedPointSet &pointcloud);
	// hole filling etc at the end (included in DoReg)
	void PostProc();
	// nosefitting (included in DoReg)
	void FitNose(int nosefitmethod=0);

	// constructors, calls Init()
	Register();

	// constructor calls Init() and DoReg()
	Register(UnorderedPointSet &ups,double resolution,int holefilling=1,int spikeremoval=1,int ellipticalmask=1,int reflectionremoval=1,int backgroundremoval=1,int nosefitmethod=0,int LR=0,double symmetrize=0.0,double motion_threshold=0.0,double maxshift=7.5,int nearfrontal=0);

	// destructor (empty)
	~Register() {};

	// get rotation matrix from registered range image ri
	int GetRotationMatrix(float *m)
	{
		RangeImage riO(ri.o,ri.u,-1*(ri.v),1,1,0,0,220,260);
		m[0]=(float)riO.T[0][0]; m[1]=(float)riO.T[0][1]; m[2]=(float)riO.T[0][2];
		m[3]=(float)riO.T[1][0]; m[4]=(float)riO.T[1][1]; m[5]=(float)riO.T[1][2];
		m[6]=(float)riO.T[2][0]; m[7]=(float)riO.T[2][1]; m[8]=(float)riO.T[2][2];

		return 1;
	}

	// get translation vector from registered range image ri
	int GetTranslationVector(float *t)
	{
		RangeImage riO(ri.o,ri.u,-1*(ri.v),1,1,0,0,220,260);
		t[0]=(float)riO.T[0][3]; t[1]=(float)riO.T[1][3]; t[2]=(float)riO.T[2][3];

		return 1;
	}

	// get a raw profile, i.e. all points of the unordered point set less than 5 mm from the xy plane
	int GetProfile(double maxd=5);

	// get a profile, and get the local max in a range of 1 mm that does not deviate 
	// more than 5 mm from the local average
	int GetMaxProfile();
};

}

#endif // REGISTER_INCLUDED
