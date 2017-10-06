/*
 * 25-4-2013: repaired small error in symmetrise that used negative x-coordinates for right shifts for first pixel
 */

#include <time.h>
#include <stdio.h>

#include "opt.h"

#include "rangeimage.h"
#include "register.h"
#include "plane3d.h"
#include "filter.h"
#include "nose.h"
#include "symmetry.h"
#include "holes.h"
#include "roi.h"

namespace utw3dface {

int nthreads=1;
int verbose=0;
int debug=0;
int Profile=0;
double profile_theta=0.5*M_PI;

void Register::Init()
{
	completenessQuality=0;
	holesQuality=0;
	spikesQuality=0;
	noiseQuality=0;
	registrationQuality=0;
	overallQuality=0;
	noseQuality=0;

	course_registration_resolution=5;	// mm
	course_registration_ri_width=200;	// mm
	resolution=0.67;			// mm
	holefilling=1;
	spikeremoval=1;
	ellipticalmask=1;
	reflectionremoval=0;
	backgroundremoval=1;
	nosefitmethod=0;
	LR=0;	// fit nose left&right
	symmetrize=0;
	motion_threshold=0.0;
	nearfrontal=0;
	maxshift=7.5;
}

Register::Register()
{
	Init();
}

Register::Register(UnorderedPointSet &pointcloud,double resolution,int holefilling,int spikeremoval,int ellipticalmask,int reflectionremoval,int backgroundremoval,int nosefitmethod,int LR,double symmetrize,double motion_threshold,double maxshift,int nearfrontal)
{
	Init();

	this->resolution=resolution;
	this->holefilling=holefilling;
	this->spikeremoval=spikeremoval;
	this->ellipticalmask=ellipticalmask;
	this->reflectionremoval=reflectionremoval;
	this->backgroundremoval=backgroundremoval;
	this->nosefitmethod=nosefitmethod;	
	this->LR=LR;
	this->symmetrize=symmetrize;
	this->motion_threshold=motion_threshold;
	this->maxshift=maxshift;
	this->nearfrontal=nearfrontal;

	DoReg(pointcloud);
}

void Register::DoReg(UnorderedPointSet &pointcloud)
{
	int i;
	int minflag=1;

	if (verbose)
		fprintf(stderr,"input point cloud #npoints=%d\n",pointcloud.npoints);

	// first find an roi
	FaceRoi faceroi(pointcloud);
	ups.CopyFrom(faceroi.upsroi);

        if (debug)
        {
fprintf(stderr,"faceroi.cylinder.c=%g %g %g\n",faceroi.cylinder.c.x,faceroi.cylinder.c.y,faceroi.cylinder.c.z);
                FaceRoi froi;
                int i;
                for (i=0; i<ups.npoints; i++)
                {
                        froi.rups.AddPoint(ups[i]);
                        froi.rupsn.AddPoint(Point3D(0,0,0));
                }

                froi.DumpRups(faceroi.cylinder,"roicyl2.wrl",faceroi.cylinder.c);
        }


	// get rid of shoulders that might upset the symmetry axis detection
	filtershoulders(ups);

	// get the extend of the remaining point cloud and average distance between points in x- and y-directions
	ups.GetStats();
	if (verbose)
	{
		fprintf(stderr,"statistics remaining point cloud:\n");
		fprintf(stderr,"npoints=%d\n",ups.npoints);
		fprintf(stderr,"minx=%g maxx=%g\n",ups.minX,ups.maxX);
		fprintf(stderr,"miny=%g maxy=%g\n",ups.minY,ups.maxY);
		fprintf(stderr,"minz=%g maxz=%g\n",ups.minZ,ups.maxZ);
		fprintf(stderr,"dx=%g dy=%g\n",ups.dx,ups.dy);
	}

	// get the centre of gravity of the point cloud
	Point3D cog=ups.CentreOfGravity();
	if (verbose)
		fprintf(stderr,"cog=(%g %g %g)\n",cog.x,cog.y,cog.z);

	// create a low resolution range image from point cloud by projection on xy plane taking
	// the centre of gravity as the origin, a width of course_registration_ri_width mm and a 
	// height of the extent of the point cloud
	double dx=course_registration_resolution;
	double dy=course_registration_resolution;
	int width=int(course_registration_ri_width/dx);
	int height=int((ups.maxY-ups.minY)/dy);
	double ox=0.5*width;
	double oy=(ups.maxY-cog.y)/dy;
	ri.Init(cog,Point3D(1,0,0),Point3D(0,-1,0),dx,dy,ox,oy,width,height);
	ri.AccumulateDepth(ups);

	if (debug)
	{
		ri.WritePGM("proj2xy.pgm",1);
		ri.WriteStddevPGM("proj2xy_stddev.pgm",1);
		ri.WriteFlagPGM("proj2xy_flag.pgm");
	}

	// drop serious outlayer points by low pass filtering the face and dropping points to far away from
	// the average
	spikesQuality=filter(ri,ups,100);

	if (backgroundremoval)
		backgroundfilter(ri,ups);

	clock_t t0=clock();

	// first find a rough estimate for the through plane rotation,in-plane rotation and hor. displacement of nose
	ArrayOfSymmetryScores symmetryscores;

	int nsymm;
	if (Profile)
	{
		//symmetryscores[0]=SymmetryScore(M_PI/2,0,-10,1);
		double theta;
		double theta_max=0.5*M_PI;
        	double theta_min=0.5*M_PI;
        	double theta_step=0.05*M_PI;
		theta_min=theta_max=profile_theta;
        	int ntheta=int((theta_max-theta_min)/theta_step)+1;
		double phi;
        	double phi_min=-0*M_PI;
        	double phi_max=0*M_PI;
        	double phi_step=0.05*M_PI;
        	int nphi=int((phi_max-phi_min)/phi_step)+1;
		int dx;
        	int dx_min=-3*ri.width/4;
        	int dx_max=3*ri.width/4;
        	int dx_step=1;
        	int ndx=int((dx_max-dx_min)/dx_step)+1;

		symmetryscores.Init(ntheta*nphi*ndx);
		nsymm=ntheta*nphi*ndx;
		
		int i=0;
		for (theta=theta_min; theta<=theta_max; theta+=theta_step)
		for (phi=phi_min; phi<=phi_max; phi+=phi_step)
		for (dx=dx_min; dx<=dx_max; dx+=dx_step)
			symmetryscores[i++]=SymmetryScore(theta,phi,dx,0);
		//symmetryscores[10]=SymmetryScore(M_PI/2,0,-13,1);
	}
	else
		nsymm=rough_symmetry(ri,ups,symmetryscores,nearfrontal);
	if (debug)
	{
		for (i=0; i<nsymm; i++)
			fprintf(stderr,"%g %g %g %g\n",symmetryscores[i].theta,symmetryscores[i].phi,symmetryscores[i].dx,symmetryscores[i].score);
	}

	double theta_opt=symmetryscores[0].theta;
	double dx_opt=symmetryscores[0].dx;
	double phi_opt=symmetryscores[0].phi;

	clock_t t1=clock();
	if (verbose)
		fprintf(stderr,"rough symmetry plane took %.6f seconds\n",double(t1-t0)/CLOCKS_PER_SEC);

	RangeImage rir(ri.o,Point3D(1,0,0),Point3D(0,-1,0),ri.du,ri.dv,ri.ou,ri.ov,ri.width,ri.height);

	if (debug)
	{
		ri.u=Point3D(cos(theta_opt),0,sin(theta_opt));
		ri.calculateT();
		ri.AccumulateDepth(minflag);
		ri.WritePGM("foomr.pgm",1);

		double cosphi=cos(phi_opt);
		double sinphi=sin(phi_opt);
		double tx=ox-ox*cosphi+oy*sinphi+dx_opt;
		double ty=oy-ox*sinphi-oy*cosphi;
		memset(rir.pixel,0,width*height*sizeof(double));
		int x,y;
		for (y=0; y<height; y++)
		for (x=0; x<width; x++)
		{
			int x1=int(x*cosphi-y*sinphi+tx);
			int y1=int(x*sinphi+y*cosphi+ty);
			if (x1<0 || x1>=width || y1<0 || y1>=height)
				continue;
			rir.Pixel(x,y)=ri.Pixel(x1,y1);
		}

		rir.WritePGM("foor.pgm",1);

		/* code to show the pixels with less than 10 mm distance */
		int n=0;
		for (y=0; y<rir.height; y++)
		for (x=0; x<rir.width/2; x++)
		{
			if (fabs(rir.Pixel(x,y)-rir.Pixel(rir.width-x-1,y))>10)
				rir.Pixel(rir.width-x-1,y)=rir.Pixel(x,y)=0;
			else
				n++;
		}
		rir.WritePGM("food.pgm",1);
		fprintf(stderr,"food-n=%d\n",n);
	}

	// roughly fit 3d nose model and get gamma_opt (tilt), nose_x, nose_y, nose_z and returns correlation
	int nose_x=rir.width/2,nose_y=rir.height/2;
	double gamma_opt=0,nose_z=0;

	double maxcorr=-1;
	int opt_i=0;

/*
	if (Profile)
	{
		
		// Profile experiments
		//symmetryscores[0]=SymmetryScore(M_PI/2,0,-10,1);
	}
	else
*/
	{
		for (i=0; i<nsymm; i++)
		{
			RangeImage rir2(rir.o,rir.u,rir.v,rir.du,rir.dv,rir.ou,rir.ov,rir.width,rir.height);
	if (debug)
	{
		ri.Init(rir.o,rir.u,rir.v,rir.du,rir.dv,rir.ou,rir.ov,rir.width,rir.height);
        	Point3D u=Point3D(cos(symmetryscores[i].theta),0,sin(symmetryscores[i].theta));
        	Point3D v=Point3D(0,-1,0);
        	ri.o.x+=symmetryscores[i].dx*ri.du;
        	Point3D u1=cos(symmetryscores[i].phi)*u+sin(symmetryscores[i].phi)*v;
        	Point3D v1=cos(symmetryscores[i].phi)*v-sin(symmetryscores[i].phi)*u;
        	ri.u=u1;
        	ri.v=v1;
		ri.calculateT();
		ri.AccumulateDepth(ups);
		GetProfile(1);
		int j;
		FILE *f;
		char s[64];
		sprintf(s,"prof_sym_%d.txt",i);
		f=fopen(s,"w");
		for (j=0; j<profile.npoints; j++)
			fprintf(f,"%g %g %g\n",profile[j].x,profile[j].y,profile[j].z);
		fclose(f);
		sprintf(s,"prof_sym_%d.pgm",i);
		ri.WritePGM(s,1);
	}
			int nx,ny;
			double gamma,nz;
			double corr=fitnose_rough(rir2,ups,minflag,symmetryscores[i].theta,symmetryscores[i].phi,symmetryscores[i].dx,gamma,nx,ny,nz,LR,nearfrontal);
			if (verbose)
				fprintf(stderr,"fitting nose for symmetry %d: %g\n",i,corr);
			if (corr>maxcorr)
			{
				if (maxcorr<0.6 || Profile)	// if the nose is good already, we choose for the best symmetry
				{
					maxcorr=corr;
					phi_opt=symmetryscores[i].phi;
					theta_opt=symmetryscores[i].theta;
					dx_opt=symmetryscores[i].dx;
					nose_x=nx;
					nose_y=ny;
					nose_z=nz;
					gamma_opt=gamma;
					opt_i=i;
				}
			}
	
		}
	
	}
	registrationQuality=maxcorr;

	if (verbose)
	{
		fprintf(stderr,"selected symmetry %d based on best nose and symmetry\n",opt_i);
	}

	if (maxcorr<=0)
	{
		overallQuality=0;
		ri.Init(Point3D(0,0,0),Point3D(1,0,0),Point3D(0,-1,0),0.67,0.67,82,127,165,195);
		ri.AccumulateDepth(ups);
		return;
	}

	Point3D u=Point3D(cos(theta_opt),0,sin(theta_opt));
	Point3D v=Point3D(0,-1,0);
	rir.ou-=dx_opt;
	Point3D u1=cos(phi_opt)*u+sin(phi_opt)*v;
	Point3D v1=cos(phi_opt)*v-sin(phi_opt)*u;
	rir.u=u1;
	rir.v=cos(gamma_opt)*v1+sin(gamma_opt)*(u1^v1);
	rir.calculateT();
	rir.AccumulateDepth(ups,minflag);

	if (debug==1)
	{
		rir.Pixel(nose_x,nose_y)+=100;
		rir.Pixel(nose_x+1,nose_y)+=100;
		rir.Pixel(nose_x-1,nose_y)+=100;
		rir.Pixel(nose_x,nose_y+1)+=100;
		rir.Pixel(nose_x,nose_y-1)+=100;
	}
	if (debug)
		rir.WritePGM("rough_reg.pgm",1);

	clock_t t2=clock();
	if (verbose)
	{
		fprintf(stderr,"rough registration took %.6f seconds\n",double(t2-t1)/CLOCKS_PER_SEC);
		fprintf(stderr,"gamma=%g theta=%g phi=%g dx=%g\n",gamma_opt,theta_opt,phi_opt,dx_opt);
		fprintf(stderr,"nose_x=%d nose_y=%d nose_z=%g\n",nose_x,nose_y,nose_z);
	}

	// define a 3-D ROI and create a new ups with only points inside the ROI
	UnorderedPointSet upsroi;
	double cx=nose_x;
	double cy=nose_y;//-20.0/rir.dv;
	double rx=55/rir.du;
	double ry=55/rir.dv;
	for (i=0; i<ups.npoints; i++)
	{
		Point3D p=rir.TransformPoint(ups[i]);
		if (fabs(p.z-nose_z)<50 && (p.x-cx)*(p.x-cx)/(rx*rx)+(p.y-cy)*(p.y-cy)/(ry*ry)<=1)
			upsroi.AddPoint(ups[i]);
	}

	Point3D o=rir.o+(((nose_x-rir.ou)*rir.du)*rir.u+((nose_y-rir.ov)*rir.dv)*rir.v+nose_z*(rir.u^rir.v));
	rirroi.Init(o,rir.u,rir.v,1,1,55,55,110,110);
	rirroi.AccumulateDepth(upsroi);

	if (verbose)
		fprintf(stderr,"%d points in roi\n",upsroi.npoints);

	if (upsroi.npoints<100)
	{
		overallQuality=0;
		ri.Init(rirroi.o,rirroi.u,rirroi.v,0.67,0.67,82,127,165,195);
		ri.AccumulateDepth(ups);
		return;
	}

	if (debug)
		rirroi.WritePGM("foorroi.pgm",1);

	if (debug)
		printf("\n\n");

/*
	UnorderedPointSet upsroit;
	for (i=0; i<upsroi.npoints; i++)
		upsroit.AddPoint(rirroi.TransformPoint(upsroi[i]));

	RangeImage rirroim(o,rir.u,rir.v,1,1,55,55,110,110);

	ox=rirroi.ou;
	phi=-0.1*M_PI;
	for (phi=-0.1*M_PI; phi<=0.1*M_PI; phi+=0.01*M_PI)
	for (ox=rirroi.ou-10; ox<=rirroi.ou+10; ox+=0.5)
	{
		Point3D n(cos(phi),sin(phi),0);
		Point3D o(ox,rirroi.ov,0);
		double k=n*o;
		double sum=0;
		int nhits=0;
		int i,flag;
		for (i=0; i<upsroit.npoints; i++)
		{
			Point3D q=upsroit[i];
			Point3D q1=q+2*(k-q*n)*n;

			q1.z=rirroi.GetDepth(q1.x,q1.y,flag);
//fprintf(stderr,"%g %g %g %g\n",q1.x,q1.y,q1.z,q.z);
			if (flag && fabs(q1.z-q.z)<10)
			{
				sum+=fabs(q1.z-q.z);
				int x=int(q1.x+0.5);
				int y=int(q1.y+0.5);
				rirroim.Pixel(x,y)=fabs(q1.z-q.z);

				nhits++;
			}
		}

		if (debug)
			printf("%g %g %d %g %g\n",phi,ox,nhits,sum,sum/nhits/nhits);
	}

	if (debug)
		rirroim.WritePGM("foorroim.pgm",0);
*/
	if (Profile)
	{
		fprintf(stderr,"Profile, skipping fine symmetry\n");
	}
	else
	{
		Opt opt(2,symmetryfunction);
		opt.SetPar(0,rirroi.ou,rirroi.ou-10,rirroi.ou+10,0.1);
		opt.SetPar(1,0,-0.1*M_PI,0.1*M_PI,0.001*M_PI);
		opt.SetData(&rirroi);
	
		Opt opt2(1,symmetryfunction2);
		opt2.SetPar(0,0,-0.1*M_PI,0.1*M_PI,0.001);
		opt2.SetData(&opt);

		if (debug)
		{
			printf("\n\n");
			double phi;
			for (phi=-0.1*M_PI; phi<=0.1*M_PI; phi+=0.001*M_PI)
				printf("%g %g\n",phi,opt2.GetFuncValue(&phi));
		}
	
		double phi,phi_min;
		for (phi=-0.1*M_PI; phi<=0.1*M_PI; phi+=0.01*M_PI)
		{
			double val=opt2.GetFuncValue(&phi);
			if (val<opt2.GetMinimum())
			{
				opt2.SetPar(0,phi,phi-0.05,phi+0.05,0.001);
				opt2.minimum=val;
				phi_min=phi;
			}
		}
	
		for (phi=phi_min-0.02*M_PI; phi<=phi_min+0.02*M_PI; phi+=0.001*M_PI)
		{
			double val=opt2.GetFuncValue(&phi);
			if (val<opt2.GetMinimum())
			{
				opt2.SetPar(0,phi,phi-0.05,phi+0.05,0.001);
				opt2.minimum=val;
				phi_min=phi;
			}
		}
		opt2.SetPar(0,phi_min,phi_min-0.05,phi_min+0.05,0.001);
	
		//opt2.ParabolicFit(5,5);
	
		if (verbose)
		{
			fprintf(stderr,"Minimum = %g\n",opt2.GetMinimum());
			fprintf(stderr,"par[0] = %g\n",opt2.GetParValue(0));
		}
	
		double theta=opt2.GetParValue(0);
		rirroi.u=cos(theta)*rirroi.u+sin(theta)*(rirroi.u^rirroi.v);
		rirroi.calculateT();
		rirroi.AccumulateDepth();
		if (debug)
			rirroi.WritePGM("final0.pgm");
	
		opt.SetPar(0,0.5*(opt.min[0]+opt.max[0]),opt.min[0],opt.max[0],opt.acc[0]);
		opt.SetPar(1,0.5*(opt.min[1]+opt.max[1]),opt.min[1],opt.max[1],opt.acc[1]);
		opt.ParabolicFit(20,3);

		if (verbose)
		{
			fprintf(stderr,"Minimum = %g\n",opt.GetMinimum());
			fprintf(stderr,"par[0] = %g\n",opt.GetParValue(0));
			fprintf(stderr,"par[1] = %g\n",opt.GetParValue(1));
		}
/*
double pp[2];
pp[0]=opt.GetParValue(0);
for (pp[1]=opt.min[1]; pp[1]<=opt.max[1]; pp[1]+=opt.acc[1])
{
	double val=opt.GetFuncValue(pp);
	printf("%g %g\n",pp[1],val);
}
*/

		rirroi.o += (opt.GetParValue(0)-rirroi.ou)*rirroi.u*rirroi.du;
		phi=opt.GetParValue(1);
		u=cos(phi)*rirroi.u+sin(phi)*rirroi.v;
		v=-sin(phi)*rirroi.u+cos(phi)*rirroi.v;
		rirroi.u=u;
		rirroi.v=v;
		rirroi.calculateT();
		rirroi.AccumulateDepth(ups);
	
		if (debug)
			rirroi.WritePGM("final1.pgm");
	
		clock_t t3=clock();
		if (verbose)
			fprintf(stderr,"fine symmetry detection took %.6f seconds\n",double(t3-t2)/CLOCKS_PER_SEC);
	
	}

	ri.Init(rirroi.o,rirroi.u,rirroi.v,rirroi.du,rirroi.dv,rirroi.ou,rirroi.ov,rirroi.width,rirroi.height);

	// fit a cylinder to get first estimate of up axis
	FILE *f;
	if (debug)
		f=fopen("cyl.txt","w");

	Cylinder cyl;
	cyl.maxdist=10;
	
	gamma_opt=0;
	double dz_opt=0;
	double min_cost=upsroi.npoints*10;
	double gamma;
	for (gamma=-M_PI/2; gamma<=M_PI/2; gamma+=0.01*M_PI)
	{
		cyl.r=100;
		cyl.a=cos(gamma)*ri.v+sin(gamma)*(ri.u^ri.v);
		cyl.a=-1*cyl.a;
		cyl.c=ri.o-cyl.r*(ri.u^cyl.a);
		int i;
		double dist=0;
		double n=cyl.Support(upsroi,dist);
		double cost=1.0/(1+n);
		if (cost<min_cost)
		{
			min_cost=cost;
			gamma_opt=gamma;
		}
		if (debug)
			fprintf(f,"%g %g %g %g\n",gamma,dist,n,min_cost);
	}

	if (debug)
		fclose(f);

	if (verbose)
		fprintf(stderr,"gamma_opt=%g\n",gamma_opt);
	cyl.a=(cos(gamma_opt)*ri.v+sin(gamma_opt)*(ri.u^ri.v));
	cyl.a=-1*cyl.a;
	cyl.c=ri.o-cyl.r*(ri.u^cyl.a);

	if (debug)
	{
		FaceRoi froi;
		int i;
		for (i=0; i<upsroi.npoints; i++)
		{
			froi.rups.AddPoint(upsroi[i]);
			froi.rupsn.AddPoint(Point3D(0,0,0));
		}
		froi.DumpRups(cyl,"cyl.wrl",ri.o);
	}

	ri.v=cos(gamma_opt)*ri.v+sin(gamma_opt)*(ri.u^ri.v);
	rirroi.v=ri.v;
	ri.o.z+=dz_opt;

	// finished fitting cylinder

	GetProfile();

	//double nx,ny,nphi,nl,nh,nt=5;
	nt=5;
	double mind=1000000;
	double noselength,best_noselength=50;
/*
	for (noselength=30; noselength<80; noselength+=10)
	{
		double d=fitnose2profile(profile,noselength,nx,ny,nphi,nl,nh,nt);
		if (debug)
			fprintf(stderr,"%g %g\n",noselength,d);
		if (d<mind)
		{
			mind=d;
			best_noselength=noselength;
		}
	}
*/
/*
	noseQuality=fitnose2profile(profile,best_noselength,nx,ny,nphi,nl,nh,nt,nosefitmethod);

	if (verbose)
		fprintf(stderr,"x=%g y=%g phi=%g\n",nx,ny,nphi);

	rirroi.o=rirroi.o-nx*rirroi.v-ny*(rirroi.u^rirroi.v);
	rirroi.v=cos(nphi)*rirroi.v+sin(nphi)*(rirroi.u^rirroi.v);
	rirroi.calculateT();
	rirroi.AccumulateDepth(ups);

	if (debug)
	{
		rirroi.WritePGM("final2.pgm");
		rirroi.WriteStddevPGM("final2stddev.pgm",0);
		rirroi.WriteFlagPGM("final2flag.pgm");
	}

	//ri.Init(rirroi.o,rirroi.u,rirroi.v,1,1,55,85,110,130);
	int rix0=int(55/resolution+0.5);
	int riy0=int(85/resolution+0.5);
	int riwidth=2*rix0+1;
	int riheight=int(130/resolution)+1;
	//ri.Init(rirroi.o,rirroi.u,rirroi.v,0.67,0.67,82,127,165,195);
	ri.Init(rirroi.o,rirroi.u,rirroi.v,resolution,resolution,rix0,riy0,riwidth,riheight);
	ri.AccumulateDepth(ups);
	if (debug)
		ri.WritePGM("final3.pgm");
*/
	FitNose(nosefitmethod);

	PostProc();
}

void Register::FitNose(int nosefitmethod)
{
	double best_noselength=50;
	noseQuality=fitnose2profile(profile,best_noselength,nx,ny,nphi,nl,nh,nt,nosefitmethod);

	if (verbose)
		fprintf(stderr,"x=%g y=%g phi=%g\n",nx,ny,nphi);

	ri.Init(rirroi.o,rirroi.u,rirroi.v,rirroi.du,rirroi.dv,rirroi.ou,rirroi.ov,rirroi.width,rirroi.height);

	ri.o=rirroi.o-nx*rirroi.v-ny*(rirroi.u^rirroi.v);
	ri.v=cos(nphi)*rirroi.v+sin(nphi)*(rirroi.u^rirroi.v);
	ri.calculateT();
	ri.AccumulateDepth(ups);

	if (debug)
	{
		ri.WritePGM("final2.pgm");
		ri.WriteStddevPGM("final2stddev.pgm",0);
		ri.WriteFlagPGM("final2flag.pgm");
	}

	//ri.Init(rirroi.o,rirroi.u,rirroi.v,1,1,55,85,110,130);
	int rix0=int(55/resolution+0.5);
	int riy0=int(85/resolution+0.5);
	int riwidth=2*rix0+1;
	int riheight=int(130/resolution)+1;
	//ri.Init(rirroi.o,rirroi.u,rirroi.v,0.67,0.67,82,127,165,195);
	ri.Init(ri.o,ri.u,ri.v,resolution,resolution,rix0,riy0,riwidth,riheight);
	ri.AccumulateDepth(ups);
	if (debug)
		ri.WritePGM("final3.pgm");
}

void Register::PostProc()
{
	UnorderedPointSet fups;
	fups.CopyFrom(ups);

	filterpickfrontal(ri,fups);
	if (debug)
		ri.WritePGM("final3frontal.pgm");

	if (spikeremoval==1)
		noiseQuality=filter(ri,fups,5);
	else if (spikeremoval>=2) // with other smoothing
		noiseQuality=filter(ri,fups,5,spikeremoval);

	if (debug)
		ri.WritePGM("spikesremoved.pgm");

	if (reflectionremoval)
		reflectionfilter(ri,fups);

	if (holefilling)
	{
		Nose nosef(0,0,0);
		//nosef.SetSize(nl,nh,5,0,0); 
		nosef.SetSize(0,0,0,0,0); 
		completenessQuality=1-fill_big_holes_using_symmetry(ri);
		holesQuality=1-fillholes(ri,nosef);

		if (debug)
			ri.WritePGM("holesfilled.pgm");
	}

/*
	filterhair(ri,fups);

	if (holefilling)
	{
		Nose nosef(0,0,0);
		nosef.SetSize(0,0,0,0,0); 
		//fill_big_holes_using_symmetry(ri);
		fillholes(ri,nosef);
	}
*/

/*
	motion compensation by shifting lines horizontally to get best symmetry
	is only applied if significant average shifts are detected
	the significant average shift is defined by the parameter motion_threshold
*/
	if (symmetrize>0)
	{
		double shifts[ri.height];

		int hs=int(symmetrize/(resolution*3));
		if (hs<1)
			hs=1;

		if (hs>ri.height/4)
		{
			fprintf(stderr,"symmetrize too big! should be less than image height/4.\n");
			exit(1);
		}

		int x,y,dy;
		int maxshift_pix=int(maxshift/resolution);

		for (y=hs; y<ri.height-hs; y++)
		{
			int shift,bestshift=0;
			double dmin=100000;
			int hw=ri.width/2;
	
			for (shift=-maxshift_pix; shift<=maxshift_pix; shift++)
			{
				double d=0;
				for (dy=-hs; dy<=hs; dy++)
				for (x=1; x<hw-maxshift_pix; x++)
					d+=fabs(ri.Pixel(hw+x+shift,y+dy)-ri.Pixel(hw-x+shift,y+dy));
				if (debug)
					fprintf(stdout,"%d %d %g\n",y,shift,d);
				if (d<dmin)
				{
					dmin=d;
					bestshift=shift;
				}
			}

			double dl=0;
			for (dy=-hs; dy<=hs; dy++)
			for (x=1; x<hw-maxshift_pix; x++)
				dl+=fabs(ri.Pixel(hw+x+bestshift-1,y+dy)-ri.Pixel(hw-x+bestshift-1,y+dy));
			double dr=0;
			for (dy=-hs; dy<=hs; dy++)
			for (x=1; x<hw-maxshift_pix; x++)
				dr+=fabs(ri.Pixel(hw+x+bestshift+1,y+dy)-ri.Pixel(hw-x+bestshift+1,y+dy));
	
			// fit parabola through dl,dmin,dr: d=ax^2+bx+c; set x=0 at bestshift
			double c=dmin;
			double a=(dr-c+dl-c)/2;
			double b=(dr-c)-a;
			// find minimum of parabola: 2ax+b=0 => x=-b/2a
			double fshift=(a==0) ? 0 : -b/(2*a);

			if (debug)
				fprintf(stdout,"#a=%g b=%g c=%g dl=%g dm=%g dr=%g shift=%g\n",a,b,c,dl,dmin,dr,bestshift+fshift);
			
			if (dl<dmin || dr<dmin)
				shifts[y]=-bestshift;
			else
				shifts[y]=-bestshift-fshift;

			if (debug)
				fprintf(stdout,"\n\n");
		}

		for (y=0; y<hs; y++)
			shifts[y]=shifts[hs];
		for (y=ri.height-hs; y<ri.height; y++)
			shifts[y]=shifts[ri.height-hs-1];

		// check if there is actually significant motion
		// use a big filter to get average shift
		int fw=ri.height/8;
		int motion=0;

		for (y=fw; y<ri.height-fw; y++)
		{
			int yy;
			double sum=0;
			for (yy=y-fw; yy<y+fw; yy++)
				sum+=shifts[yy];
			double ashift=sum/(2*fw+1);
			if (ashift>=motion_threshold)
				motion=1;
		}
		if (verbose)
		{
			if (motion)
				fprintf(stderr,"Detected motion: shift>%g, compensating ...\n",motion_threshold);
			else
				fprintf(stderr,"No significant motion detected. No compensation.\n");
		}

		if (motion)
		{
			RangeImage ris;
			ris.CopyDataFrom(ri);
	
			for (y=0; y<ri.height; y++)
			{
				int shift=int(floor(shifts[y]));
				double fshift=shifts[y]-shift;
	
				if (verbose)
					fprintf(stderr,"shift line %d=%g %d %g\n",y,shifts[y],shift,fshift);
	
				if (shift>=0)
				{
					for (x=0; x<shift+1; x++)
						ri.Pixel(x,y)=ri.Pixel(0,y);
					for (x=shift+1; x<ri.width; x++)
						ri.Pixel(x,y)=(1-fshift)*ris.Pixel(x-shift,y)+fshift*ris.Pixel(x-(shift+1),y);
				} else if (shift<0)
				{
					for (x=0; x<ri.width+shift; x++)
						ri.Pixel(x,y)=(1-fshift)*ris.Pixel(x-shift,y)+fshift*ris.Pixel(x-(shift+1),y);
					for (x=ri.width+shift; x<ri.width; x++)
						ri.Pixel(x,y)=ri.Pixel(ri.width-1,y);
				}

				//for (x=0; x<ri.width; x++)
					//if (x-shift>=0 && x-shift<ri.width && x-shift+1<ri.width)
						//ri.Pixel(x,y)=(1-fshift)*ris.Pixel(x-shift,y)+fshift*ris.Pixel(x-(shift+1),y);
			}
		}
	}
	
	if (ellipticalmask)
	{
		// mask out an ellipse
		//rx=55/ri.du;
		//ry=70/ri.dv;
		double rx=55/ri.du;
		double ry=70/ri.dv;
		double cx=ri.ou;
		double cy=ri.ov-20/ri.dv;
		int x,y;
		for (y=0; y<ri.height; y++)
		for (x=0; x<ri.width; x++)
		{
			if ((x-cx)*(x-cx)/(rx*rx)+(y-cy)*(y-cy)/(ry*ry)>1)
			{
				ri.Flag(x,y)=0;
				ri.Pixel(x,y)=0;
			}
		}
	}

	overallQuality=noseQuality*spikesQuality*holesQuality*noiseQuality*completenessQuality*registrationQuality;

	if (verbose)
	{
		fprintf(stderr,"noseQuality=%g\n",noseQuality);
		fprintf(stderr,"spikersQuality=%g\n",spikesQuality);
		fprintf(stderr,"holesQuality=%g\n",holesQuality);
		fprintf(stderr,"noiseQuality=%g\n",noiseQuality);
		fprintf(stderr,"completenessQuality=%g\n",completenessQuality);
		fprintf(stderr,"registrationQuality=%g\n",registrationQuality);
		fprintf(stderr,"overallQuality=%g\n",overallQuality);
	}
}

int Register::GetProfile(double maxd)
{
	profile.npoints=0;
	Plane3D yzplane(ri.u,ri.o*ri.u);
	Plane3D xzplane(ri.v,ri.o*ri.v);
	Plane3D xyplane(ri.u^ri.v,ri.o*(ri.u^ri.v));
	int i;
	for (i=0; i<ups.npoints; i++)
	{
		double d=yzplane.SignedDistance(ups[i]);
		double z=xyplane.SignedDistance(ups[i]);
		double y=xzplane.SignedDistance(ups[i]);

		if (fabs(d)<maxd)
			profile.AddPoint(Point3D(y,z,d));
	}

	profile.SortXYZ();

	return 1;
}

int Register::GetMaxProfile()
{
	if (profile.npoints==0)
		GetProfile();

	int i;

	// resample the profile with 1mm distance, taking the average
	int nmaxpoints=int(profile[profile.npoints-1].x-profile[0].x);

	double *ymean=new double[nmaxpoints];
	int *nmean=new int[nmaxpoints];
	for (i=0; i<nmaxpoints; i++)
	{
		ymean[i]=0.0;
		nmean[i]=0;
	}
	for (i=0; i<profile.npoints; i++)
	{
		Point3D p=profile[i];
		int j=int(p.x-profile[0].x+0.5);
		if (j<0 || j>=nmaxpoints)
			continue;
		ymean[j]+=p.y;
		nmean[j]++;
	}
	for (i=0; i<nmaxpoints; i++)
	{
		if (nmean[i]>0)
			ymean[i]/=nmean[i];
	}

	Point3D *maxpoints=new Point3D[nmaxpoints];
	for (i=0; i<nmaxpoints; i++)
		maxpoints[i]=Point3D(0,-100001,0);
	for (i=0; i<profile.npoints; i++)
	{
		Point3D p=profile[i];
		int j=int(p.x-profile[0].x+0.5);
		if (j<0 || j>=nmaxpoints)
			continue;
		if (p.y>maxpoints[j].y && fabs(p.y-ymean[j])<5)
			maxpoints[j]=p;
	}

	for (i=0; i<nmaxpoints; i++)
		if (maxpoints[i].y>-100000)
			maxprofile.AddPoint(maxpoints[i]);
	delete [] maxpoints;
	delete [] ymean;
	delete [] nmean;

	return 1;
}
}
