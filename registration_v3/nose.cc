#include "compat.h"
#include "math.h"
#include "nose.h"
#include "roi.h"

#include "spline.h"

namespace utw3dface {

extern int verbose;
extern int debug;

double nosemodel(Opt *opt)
{
	static int nrefs=0;
	double x=opt->pars[0];
	double y=opt->pars[1];
	double phi=opt->pars[2];

	UnorderedPointSet *ups=(UnorderedPointSet*)opt->data;
	double cosphi=cos(phi);
	double sinphi=sin(phi);
	
	Nose nose(x,y,phi);
	nose.SetSize(30,30,5,30,50); 

/*
	int i;
	double n=0;
	double sum=0;
	for (i=0; i<ups->npoints; i++)
	{
		Point3D p=(*ups)[i];
		if (!nose.PartOf(p))
			continue;
		sum+=nose.Dist(p);
		n++;
	}

	double val=(n<10) ? 1000 : sum/n;
*/
	double val=-nose.Corr(*ups);
//fprintf(stderr,"%d %g %g %g %g\n",nrefs++,x,y,phi,val);

	return val;
}

double fitnose2profile_old(UnorderedPointSet &profile,double noselength,double &nx,double &ny,double &nphi,double &nl,double &nh,double &nt)
{
	int i;

	printf("\n\n");
	for (i=0; i<profile.npoints; i++)
		fprintf(stdout,"%g %g %g\n",profile[i].x,profile[i].y,profile[i].z);
	printf("\n\n");

	// nose size search range and step
	const double NLMIN=20;
	//const double NLMAX=70;
	const double NLMAX=noselength;
	const double NHMIN=10;
	const double NHMAX=30;
	const double NLSTEP=2;
	const double NHSTEP=2;
	const double NTMIN=5;
	const double NTMAX=20;
	const double NTSTEP=2;

	const double NTIPLENGTH=5;
	const double NMOUTHLENGTH=10;
	const double NFHLENGTH=10;

	// nose position and orientation search range and step
	const double DFHN=noselength;
	const double DNMIN=-40;
	const double DNMAX=20;
	const double NSTEPSIZE=1;
	const int NSTEPS=int((DNMAX-DNMIN)/NSTEPSIZE);

	int iz;
	double *nz=new double[NSTEPS];
	double *fz=new double[NSTEPS];
	double *dtz=new double[NSTEPS];
	double *nnz=new double[NSTEPS];
	double *nfz=new double[NSTEPS];
	double *ndtz=new double[NSTEPS];

	for (iz=0; iz<NSTEPS; iz++)
	{
		nz[iz]=fz[iz]=dtz[iz]=0;
		nnz[iz]=nfz[iz]=ndtz[iz]=0;
	}

	// get forehead depths
	for (i=0; i<profile.npoints; i++)
	{
		Point3D p=profile[i];
		int iz=int(p.x-(DNMIN+DFHN)/NSTEPSIZE+0.5);
		if (iz<0 || iz>=NSTEPS)
			continue;
		fz[iz]+=p.y;
		nfz[iz]++;
	}
	for (iz=0; iz<NSTEPS; iz++)
		fz[iz]/=(nfz[iz]>0) ? nfz[iz]:1;

	// get nose base
	for (i=0; i<profile.npoints; i++)
	{
		Point3D p=profile[i];
		int iz=int(p.x-DNMIN/NSTEPSIZE+0.5);
		if (iz<0 || iz>=NSTEPS)
			continue;
		nz[iz]+=p.y;
		nnz[iz]++;
	}
	for (iz=0; iz<NSTEPS; iz++)
		nz[iz]/=(nnz[iz]>0) ? nnz[iz]:1;

	double best_mind=1000000;
	double best_ny=0;
	double best_nl=0;
	double best_nh=0;
	double best_nt=0;
	double best_nphi=0;
	double best_nz=0;

	if (debug>1)
		printf("\n\n");

	for (iz=0; iz<NSTEPS; iz++)
	{
		if (nz[iz]==0)
			continue;

		// get the best matching forehead value
		int izf;
		double d;
		for (izf=iz; izf>0; izf--)
		{
			d=sqrt((nz[iz]-fz[izf])*(nz[iz]-fz[izf])+(izf-iz+DFHN)*(izf-iz+DFHN));
			if (d<DFHN)
				break;
		}
		if (izf<iz)
			izf++;
		double nphi=atan2(fz[izf]-nz[iz],izf-iz+DFHN);
		double ny=iz*NSTEPSIZE+DNMIN;
nphi=0;
		// don't accept extreme angles  - added 11-12-07
		//if (nphi<M_PI/3 || nphi>M_PI/3)
			//continue;

		// this might help if we only want to use max values
		int i;
		int nmaxpoints=int(NLMAX+NMOUTHLENGTH+NFHLENGTH+20);
		double cosnphi=cos(nphi);
		double sinnphi=sin(nphi);

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
			Point3D q;
			q.x=cosnphi*(p.x-ny)+sinnphi*(p.y-nz[iz]);
			q.y=-sinnphi*(p.x-ny)+cosnphi*(p.y-nz[iz]);
			int j=int(q.x+NMOUTHLENGTH+10+0.5);
			if (j<0 || j>=nmaxpoints)
				continue;
			ymean[j]+=q.y;
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
			Point3D q;
			q.x=cosnphi*(p.x-ny)+sinnphi*(p.y-nz[iz]);
			q.y=-sinnphi*(p.x-ny)+cosnphi*(p.y-nz[iz]);
			q.z=p.z;
			int j=int(q.x+NMOUTHLENGTH+10+0.5);
			if (j<0 || j>=nmaxpoints)
				continue;
			if (q.y>maxpoints[j].y && fabs(q.y-ymean[j])<5)
				maxpoints[j]=q;
		}
		UnorderedPointSet maxprofile;
		for (i=0; i<nmaxpoints; i++)
			if (maxpoints[i].y>-100000)
				maxprofile.AddPoint(maxpoints[i]);
		double steepness=0;
		for (i=2; i<maxprofile.npoints-2; i++)
			if (maxprofile[i].x>0 && maxprofile[i-1].x<=0)
			{
				double steepness1=(maxprofile[i].y-maxprofile[i-1].y)/(maxprofile[i].x-maxprofile[i-1].x);
				double steepness2=(maxprofile[i+1].y-maxprofile[i-2].y)/(maxprofile[i+1].x-maxprofile[i-2].x);
				steepness=(steepness1+steepness2)/2;
				break;
			}
		//if (steepness>1)
			//continue;
		delete [] maxpoints;
		delete [] ymean;
		delete [] nmean;
		
		//Nose nose(ny,nz[iz],nphi);
		Nose nose(0,0,0);
		double t,l,h;
		double mind=1000000;
		double nl=0;
		double nh=0;
		double nt;
		for (l=NLMIN; l<NLMAX; l+=NLSTEP)
		for (h=NHMIN; h<NHMAX; h+=NHSTEP)
		for (t=NTMIN; t<NTMAX; t+=NTSTEP)
		{
			nose.SetSize(l,h,t,NMOUTHLENGTH,NFHLENGTH); 
			double d=nose.Dist(maxprofile);
			if (d<mind)
			{
				mind=d;
				nl=l;
				nh=h;
				nt=t;
			}
		}
		if (debug>1)
			printf("%g %g %g %g %g %g %g %d %d\n",ny,nz[iz],nl,nh,nt,nphi,mind,iz,izf);

		if (mind<best_mind)
		{
			best_mind=mind;
			best_nl=nl;
			best_nh=nh;
			best_nt=nt;
			best_nphi=nphi;
			best_ny=ny;
			best_nz=nz[iz];
		}
	}

	if (verbose)
		fprintf(stderr,"ny=%g nz=%g nphi=%g nl=%g nh=%g nt=%g mind=%g\n",best_ny,best_nz,best_nphi,best_nl,best_nh,best_nt,best_mind);

/*
fprintf(stdout,"\n\n");
int nmaxpoints=int(NLMAX+NMOUTHLENGTH+NFHLENGTH+20);
Point3D maxpoints[nmaxpoints];
for (i=0; i<nmaxpoints; i++)
	maxpoints[i]=Point3D(0,-100000,0);
double cosnphi=cos(best_nphi);
double sinnphi=sin(best_nphi);
for (i=0; i<profile.npoints; i++)
{
	Point3D p=profile[i];
	Point3D q;
	q.x=cosnphi*(p.x-best_ny)+sinnphi*(p.y-best_nz);
	q.y=-sinnphi*(p.x-best_ny)+cosnphi*(p.y-best_nz);
	q.z=p.z;
	int j=int(q.x+NMOUTHLENGTH+10+0.5);
	if (j<0 || j>=nmaxpoints)
		continue;
	if (q.y>maxpoints[j].y)
		maxpoints[j]=q;
}
UnorderedPointSet maxprofile;
for (i=0; i<nmaxpoints; i++)
	if (maxpoints[i].y>-100000)
		maxprofile.AddPoint(maxpoints[i]);
for (i=0; i<maxprofile.npoints; i++)
	printf("%g %g\n",maxprofile[i].x,maxprofile[i].y);

	fprintf(stdout,"\n\n");
	Nose nose(0,0,0);
	nose.SetSize(best_nl,best_nh,best_nt,NMOUTHLENGTH,NFHLENGTH); 
	int x;
	for (x=int(nose.Start()); x<int(nose.End()); x++)
	{
		Point3D p=nose.XY(x);
		fprintf(stdout,"%g %g\n",p.x,p.y);
	}
	fprintf(stdout,"\n\n");
*/
// now things could further be improved by determining a new angle by getting two points next to the nose and one on the forehead

	nx=best_ny;
	ny=best_nz;
	nphi=best_nphi;
	nl=best_nl;
	nh=best_nh;
	nt=best_nt;

	delete [] nz;
	delete [] fz;
	delete [] dtz;
	delete [] nnz;
	delete [] nfz;
	delete [] ndtz;

	return best_mind;
}

// calculates gamma_opt, nose_x,nose_y,nose_z, returns correlation
// LR determines left fit (1), rightfit (2) or leftright fit (0)
double fitnose_rough(RangeImage &rir,UnorderedPointSet &ups,int minflag,double theta_opt,double phi_opt,double dx_opt,double &gamma_opt,int &nose_x,int &nose_y,double &nose_z,int LR,int nearfrontal)
{
	double gamma_min=-0.2*M_PI;
	double gamma_max=0.2*M_PI;
	double gamma_step=0.025*M_PI;

	if (nearfrontal)
	{
		gamma_min*=0.5;
		gamma_max*=0.5;
	}

	double dx=rir.du;
	double dy=rir.dv;
	double cosphi=cos(phi_opt);
	double sinphi=sin(phi_opt);
	int width=rir.width;
	int height=rir.height;

	if (debug)
		printf("\n\n# rough nose fitting results\n");

// generate a nose model
	Nose nose(0,0,0);
	nose.SetSize(30,30,5,10,10); 
	int w=2*int((1.5*nose.height)/dx)+1;
	int h=int((nose.length+nose.mouthlength+nose.foreheadlength)/dy);
	int noy=int((nose.length+nose.foreheadlength)/dy)-1;
	int nox=w/2;
	RangeImage rinose(rir.o,Point3D(1,0,0),Point3D(0,-1,0),dx,dy,0.5*w,0.5*h,w,h);
	int x,y;
	for (y=0; y<h; y++)
	for (x=0; x<w; x++)
	{
		double nx=(noy-y)*dy;
		double z=nose.NY(nx);
		double d=fabs(dx*(x-nox));
		//rinose.Pixel(x,y)=(z-d*2>0) ? -z+2*d : 0;
		rinose.Pixel(x,y)=(z-1.5*d>0) ? -z+1.5*d : 0;
	}

	if (debug)
		rinose.WritePGM("nose.pgm",1);

	Point3D u=Point3D(cos(theta_opt),0,sin(theta_opt));
	Point3D v=Point3D(0,-1,0);
	rir.ou-=dx_opt;
	Point3D u1=cosphi*u+sinphi*v;
	Point3D v1=cosphi*v-sinphi*u;
	rir.u=u1;

	double gamma=0.0*M_PI;
	gamma_opt=0;
	double nmax=-1;
	nose_z=0;
	nose_y=0;
	nose_x=0;

	if (debug>1)
		printf("# rough nose detection\n");
	for (gamma=gamma_min; gamma<=gamma_max; gamma+=gamma_step)
	{
		rir.v=cos(gamma)*v1+sin(gamma)*(u1^v1);
		rir.calculateT();
		rir.AccumulateDepth(ups);
// Profile handling
//char s[20];
//sprintf(s,"im_%g_%g.pgm",theta_opt,gamma);
//rir.WritePGM(s,1);

		int y_opt=0;
		int x_opt=0;
		int n_opt=0;
		double z_opt=0;
		double *corry=new double[height];
		memset(corry,0,height*sizeof(double));
		int *x_opty=new int[height];
		double *z_opty=new double[height];
		int *n_opty=new int[height];

		// set x-matching range for left, right or both (if nose is occluded or missing)
		int xxmin=0, xxmax=w;
		switch (LR)
		{
			case 1: // left fit
				xxmin=0; xxmax=w/2; break;
			case 2: // right fit
				xxmin=w/2; xxmax=w; break;
			default: // left right fit
				xxmin=0; xxmax=w-1; break;
		}

		for (y=int(noy+30/dy); y<height-h+noy-30/dy; y++)
		{
			for (x=width/2-3; x<=width/2+3; x++)
			{
				int xx,yy;
				int n=0;
				double ncorr=0;
				double rirmean=0,rirvar=0;
				double nosemean=0,nosevar=0;
				for (yy=0; yy<h; yy++)
				for (xx=xxmin; xx<=xxmax; xx++)
				{
					if (rir.Flag(xx-nox+x,y+yy-noy)==0)
						continue;
					n++;
					double a=rinose.Pixel(xx,yy);
					nosemean+=a;
					nosevar+=a*a;
					double b=rir.Pixel(xx-nox+x,y+yy-noy);
					rirmean+=b;
					rirvar+=b*b;
					ncorr+=a*b;
				}
				if (n<0.5*w*h)
					continue;
				if (n<0.8*w*h)
				{
					//if left or right half is missing, still OK
					int nl=0,nr=0;
					for (yy=0; yy<h; yy++)
					{
						for (xx=0; xx<w/2; xx++)
							if (rir.Flag(xx-nox+x,y+yy-noy)>0)
								nl++;
						for (xx=w/2; xx<w; xx++)
							if (rir.Flag(xx-nox+x,y+yy-noy)>0)
								nr++;
					}
					if (nl<0.8*n && nr<0.8*n)
						continue;

				}

				nosemean/=n;
				rirmean/=n;
				nosevar=nosevar/n-nosemean*nosemean;
				rirvar=rirvar/n-rirmean*rirmean;
				ncorr=(ncorr/n-rirmean*nosemean)/sqrt(nosevar*rirvar);

				if (ncorr>corry[y])
				{
					corry[y]=ncorr;
					x_opty[y]=x;
					z_opty[y]=rirmean;
					n_opty[y]=n;
				}
			}
			if (debug>1)
				printf("%g %d %g %d\n",gamma,y,corry[y],n_opty[y]);
		}

		// get local maxima
		for (y=1; y<height; y++)
		{
			if (corry[y]>corry[y-1])
				corry[y-1]=0;
		}
		
		for (y=height-2; y>0; y--)
		{
			if (corry[y]>corry[y+1])
				corry[y+1]=0;
		}

		double ncorrmax=-1;
		for (y=0; y<height; y++)
		{
			if (corry[y]>ncorrmax)
			{
				ncorrmax=corry[y];
				y_opt=y;
				x_opt=x_opty[y];
				z_opt=z_opty[y];
				n_opt=n_opty[y];
			}
		}

		delete [] corry;
		delete [] x_opty;
		delete [] z_opty;
		delete [] n_opty;

		if (debug==1)
			printf("%g %g %d %d\n",gamma,ncorrmax,y_opt,n_opt);
		if (ncorrmax>nmax)
		{
			nmax=ncorrmax;
			gamma_opt=gamma;
			nose_y=y_opt;
			nose_x=x_opt;
			nose_z=z_opt;
		}
		if (debug>1)
			printf("\n\n");
	}

	return nmax;
}

double fitnose2profile(UnorderedPointSet &profile,double noselength,double &nx,double &ny,double &nphi,double &nl,double &nh,double &nt,int nosefitmethod=0)
{
	double quality=1;
	int i;

	// nose size search range and step
	const double NLMIN=20;
	const double NLMAX=80;
	const double NHMIN=10;
	const double NHMAX=30;
	const double NLSTEP=1;
	const double NHSTEP=1;
	const double NTMIN=5;
	const double NTMAX=20;
	const double NTSTEP=1;

	const double NMOUTHLENGTH=10;
	const double NFHLENGTH=0;

	// nose position and orientation search range and step
	const double DNMIN=-20;
	//const double DNMAX=20;
	const double DNMAX=30;
	const double DNSTEP=1;

	double best_mind=1000000;
	double best_ny=0;
	double best_nl=0;
	double best_nh=0;
	double best_nt=0;
	double best_nphi=0;
	double best_nz=0;

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
	UnorderedPointSet maxprofile;
	for (i=0; i<nmaxpoints; i++)
		if (maxpoints[i].y>-100000)
			maxprofile.AddPoint(maxpoints[i]);
	delete [] maxpoints;
	delete [] ymean;
	delete [] nmean;

	if (debug)
	{
		printf("\n\n# profile\n");
		for (i=0; i<profile.npoints; i++)
			printf("%g %g\n",profile[i].x,profile[i].y);
		printf("\n\n");
		for (i=0; i<maxprofile.npoints; i++)
			printf("%g %g\n",maxprofile[i].x,maxprofile[i].y);
		printf("\n\n");
	}
		
	// first check if the potential nose position has any big holes
	double maxholesize=0;
	int izmax=0;
	double zmax=-10000000;
	int ihole=0;
	for (i=1; i<maxprofile.npoints-1; i++)
	{
		if (maxprofile[i].x<DNMIN || maxprofile[i].x>DNMAX)
			continue;
		if (maxprofile[i].y>zmax)
		{
			zmax=maxprofile[i].y;
			izmax=i;
		}
		double dx=maxprofile[i+1].x-maxprofile[i].x;
		if (dx>maxholesize)
		{
			maxholesize=dx;
			ihole=i;
		}
		dx=maxprofile[i].x-maxprofile[i-1].x;
		if (dx>maxholesize)
		{
			maxholesize=dx;
			ihole=i-1;
		}
	}

	if (verbose)
		fprintf(stderr,"fitnose2profile: maxholesize=%g ihole=%d izmax=%d zmax=%g\n",maxholesize,ihole,izmax,zmax);

/* some code that might handle holes in noses; added 13/3/2013
	// if there is a hole within 5 mm below the max, probably the nosetip is in the hole
	if (ihole<izmax && maxprofile[izmax].x-maxprofile[ihole+1].x<5)
		izmax=ihole;
*/
	int imax=0;
	for (i=izmax; i<maxprofile.npoints-1; i++)
	{
		if (maxprofile[i].x-maxprofile[izmax].x>DNMAX)
		{
			imax=i;
			break;
		}
	}

	if (verbose)
		fprintf(stderr,"fitnose2profile: izmax=%d imax=%d\n",izmax,imax);

	if (debug)
		printf("# get line through nose bridge with max support\n");
		
	int j;
	int maxsupport=0;
	double mindtot=1000000;
	// proper initialisation of best_p and bestp_r; added 13-3-2013 
	Point3D best_p=maxprofile[izmax];
	Point3D best_r=Point3D(1,-tan(M_PI/6),0);
	for (i=izmax; i<imax; i++)
	for (j=izmax; j<imax; j++)
	{
		if (maxprofile[j].x<maxprofile[i].x+10)
			continue;
		// define a line through points i and j
		Point3D p=maxprofile[i];
		Point3D r=maxprofile[j]-maxprofile[i];
		r.y/=r.x;
		r.x=1;
		// get support
		int support=0;
		double dtot=0;
		int k;
		for (k=izmax; k<imax+30; k++)
		{
			if (k>=maxprofile.npoints)
				continue;
			double y=p.y+(maxprofile[k].x-p.x)*r.y;
			double d=fabs(y-maxprofile[k].y);
// 31-5-2013 this should be better and sometimes is, but not always
//d=fabs(p.x*r.y-p.y*r.x-maxprofile[k].x*r.y+maxprofile[k].y*r.x)/sqrt(1+r.y*r.y);
			if (d<2)
			{
				support++;	
				dtot+=d;
			}
		}
		dtot/=(support>0) ? support : 1;
		if (debug)
			printf("%g %g %g %d %g\n",maxprofile[i].x,maxprofile[j].x,r.y,support,dtot);
		if (support>maxsupport)
		//if (0.9*support>=maxsupport || support>=0.9*maxsupport && dtot<mindtot) /* yet another way to handle nose bridge, works better with little support */
		//if ((support>maxsupport) || (support>=maxsupport  && dtot<mindtot)) /* updated: 23-6-2010, more accurate nose bridge */
		{
			maxsupport=support;
			mindtot=dtot;
			best_r=r;
			best_p=p;
		}
	}

	if (verbose)
		fprintf(stderr,"fitnose2profile maxsupport=%d\n",maxsupport);

	if (debug)
	{
		printf("\n\n");
		double x;
		for (x=DNMIN-20; x<DNMAX+50; x++)
			printf("%g %g\n",x,best_p.y+(x-best_p.x)*best_r.y);
	}

	nphi=-M_PI/6-atan2(best_r.y,best_r.x);
	if (verbose)
		fprintf(stderr,"fitnose2profile: phi=%g\n",nphi);
	double cosphi=cos(nphi);
	double sinphi=sin(nphi);

/* just checking
{
double rx1=best_r.x*cosphi-best_r.y*sinphi;
double ry1=best_r.x*sinphi+best_r.y*cosphi;
fprintf(stderr,"****** %g %g %g %g\n",rx1,ry1,atan2(ry1,rx1),-M_PI/6);
}
*/

	// find the lowest point that still lies more or less on the line on the bridge of the nose
	for (i=izmax; i>0; i--)
	{
		if (maxprofile[i].x<maxprofile[izmax].x-20)
			break;
		double y=best_p.y+(maxprofile[i].x-best_p.x)*best_r.y;
		if (fabs(maxprofile[i].y-y)>2)
			break;
	}
	izmax=i;

	UnorderedPointSet maxrprofile;

{
	// resample the profile with 1mm distance, taking the average, but now rotated over phi
	UnorderedPointSet rprofile;
	for (i=0; i<profile.npoints; i++)
	{
		double x=cosphi*profile[i].x-sinphi*profile[i].y;
		double y=sinphi*profile[i].x+cosphi*profile[i].y;
		rprofile.AddPoint(Point3D(x,y,profile[i].z));
	}
	int nmaxpoints=int(rprofile[rprofile.npoints-1].x-rprofile[0].x);
	if (nmaxpoints<2)
	{
		if (verbose)
			fprintf(stderr,"fitnose2profile: nose detection completely failed\n");
		nmaxpoints=2;
		quality=0;
	}

	double *ymean=new double[nmaxpoints];
	int *nmean=new int[nmaxpoints];
	for (i=0; i<nmaxpoints; i++)
	{
		ymean[i]=0.0;
		nmean[i]=0;
	}
	for (i=0; i<rprofile.npoints; i++)
	{
		Point3D p=rprofile[i];
		int j=int(p.x-rprofile[0].x+0.5);
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
	for (i=0; i<rprofile.npoints; i++)
	{
		Point3D p=rprofile[i];
		int j=int(p.x-rprofile[0].x+0.5);
		if (j<0 || j>=nmaxpoints)
			continue;
		if (p.y>maxpoints[j].y && fabs(p.y-ymean[j])<5)
			maxpoints[j]=p;
	}
	for (i=0; i<nmaxpoints; i++)
		if (maxpoints[i].y>-100000)
			maxrprofile.AddPoint(maxpoints[i]);
	delete [] maxpoints;
	delete [] ymean;
	delete [] nmean;

	if (debug)
	{
		printf("\n\n");
		for (i=0; i<maxrprofile.npoints; i++)
			printf("%g %g\n",maxrprofile[i].x,maxrprofile[i].y);
		printf("\n\n");
	}
}

	// perform a median filter to remove small spikes
	UnorderedPointSet maxrprofile_f;
	for (i=1; i<maxrprofile.npoints-1; i++)
	{
		if (fabs(maxrprofile[i].x-maxrprofile[i-1].x)<5 && fabs(maxrprofile[i+1].x-maxrprofile[i].x)<5)
		{
			double low=maxrprofile[i-1].y;
			double med=maxrprofile[i].y;
			double high=maxrprofile[i+1].y;
			double tmp;
			if (low>med) { tmp=med; med=low; low=tmp; }
			if (high<med) { tmp=med; med=high; high=tmp; }
			if (low>med) { tmp=med; med=low; low=tmp; }
			if (high<med) { tmp=med; med=high; high=tmp; }
			Point3D p(maxrprofile[i]);
			p.y=med;
			maxrprofile_f.AddPoint(p);
		}
	}
	maxrprofile.CopyFrom(maxrprofile_f);
	
	Point3D tip=Point3D(0,-100000000,0);
	double xcenter=cosphi*maxprofile[izmax].x-sinphi*maxprofile[izmax].y;
	for (i=0; i<maxrprofile.npoints; i++)
	{
		if (fabs(maxrprofile[i].x-xcenter)>5)
			continue;
		if (maxrprofile[i].y>tip.y)
			tip=Point3D(maxrprofile[i].x,maxrprofile[i].y,0);
	}

	if (tip.y<=-99999999)
	{
		if (verbose)
			fprintf(stderr,"No tip of nose found near %g\n",maxprofile[izmax].x);
		quality=0;
		for (i=0; i<maxrprofile.npoints; i++)
		{
			if (maxrprofile[i].y>tip.y)
				tip=Point3D(maxrprofile[i].x,maxrprofile[i].y,0);
		}
		if (tip.y<=-99999999)
			tip=Point3D(0,0,0);
	}

	// set the tip of the nose to the intersection of the horizontal line
	// through the maximum point and the line along the nose bridge
	if (best_r.x!=0 && best_r.y != 0)
	{
		double px1=best_p.x*cosphi-best_p.y*sinphi;
		double py1=best_p.x*sinphi+best_p.y*cosphi;
		double rx1=best_r.x*cosphi-best_r.y*sinphi;
		double ry1=best_r.x*sinphi+best_r.y*cosphi;
		ry1/=(rx1!=0) ? rx1 : 0.0001;
		rx1=1;
		if (ry1==0)
			ry1=0.00001;
	
		double tipy1=tip.y;
		double tipx1=(tipy1-py1)/ry1+px1;
		
		if (verbose)
		{
			fprintf(stderr,"tipx=%g tipx1=%g\n",tip.x,tipx1);
			fprintf(stderr,"tipy=%g tipy1=%g\n",tip.y,tipy1);
			fprintf(stderr,"rx=%g rx1=%g\n",best_r.x,rx1);
			fprintf(stderr,"ry=%g ry1=%g\n",best_r.y,ry1);
		}
	
		if (debug)
		{
			printf("\n\n");
			double x;
			for (x=DNMIN-20; x<DNMAX+50; x++)
				printf("%g %g\n",x,py1+(x-px1)*ry1);
			printf("\n\n");
			for (i=0; i<maxprofile.npoints; i++)
			{
				double x=cosphi*maxprofile[i].x-sinphi*maxprofile[i].y;
				double y=sinphi*maxprofile[i].x+cosphi*maxprofile[i].y;
				printf("%g %g\n",x,y);
			}
			printf("\n\n");
		}
		tip.x=tipx1;
	}

// nosefitmethod: 
//	1 = dent in nose bridge as reference, re-estimate slope of nose bridge
//	2 = below nose point as reference, orientation determined by dent-below nose point
//	3 = dent in nose bridge as reference, orientation determined by dent and point on forehead
if (nosefitmethod>0)
{
	// some new experimental code to use dent in nose bridge of nose as reference
	// move upward to the forehead until the line along the bridge of the nose drops below the forehead
	FILE *f=NULL;
	if (debug)
		f=fopen("forehead.txt","w");

	double px1=best_p.x*cosphi-best_p.y*sinphi;
	double py1=best_p.x*sinphi+best_p.y*cosphi;
	double rx1=best_r.x*cosphi-best_r.y*sinphi;
	double ry1=best_r.x*sinphi+best_r.y*cosphi;
	ry1/=(rx1!=0) ? rx1 : 0.0001;
	rx1=1;
	if (ry1==0)
		ry1=0.00001;

	UnorderedPointSet maxrprofile_s;
	int hfilterwidth=5; // was 5
	for (i=0; i<maxrprofile.npoints; i++)
	{
		double x=maxrprofile[i].x;
		double y=maxrprofile[i].y;

		double meanx=x;
		double meany=y;
		int j,n=1;
		for (j=i-1; j>=0; j--)
		{
			if (maxrprofile[j].x<x-hfilterwidth)
				break;
			n++;
			meanx+=maxrprofile[j].x;
			meany+=maxrprofile[j].y;
		}
		for (j=i+1; j<maxrprofile.npoints; j++)
		{
			if (maxrprofile[j].x>x+hfilterwidth)
				break;
			n++;
			meanx+=maxrprofile[j].x;
			meany+=maxrprofile[j].y;
		}
		meanx/=n;
		meany/=n;
		
		if (debug)
			fprintf(f,"%d %g %g %g %g %g\n",i,x,y,meanx,meany,py1+(x-px1)*ry1);
		Point3D p(meanx,meany,0);
		maxrprofile_s.AddPoint(p);
	}

	/* added 20/3/13 to fill big holes in profile (bigger than 10mm: use 1st order extrapolation from top */
	if (maxholesize>10)
	{
		int n=maxrprofile.npoints;
		for (i=2; i<n-2; i++)
		{
			if (maxrprofile_s[i].x-maxrprofile_s[i-1].x>10)
			{
				int i2;
				for (i2=i+1; i2<n-1; i2++)
				{
					if (maxrprofile_s[i2].x-maxrprofile_s[i].x>2)
						break;
				}
				Point3D dq(maxrprofile_s[i2]-maxrprofile_s[i]);
				dq.y/=(dq.x<1) ? 1 : dq.x;	
				dq.x=1;

		 		Point3D p=maxrprofile_s[i-1];
		 		Point3D q=maxrprofile_s[i];
				for (;;)
				{
					q-=dq;
					maxrprofile_s.AddPoint(q);
					if (q.x-p.x<5)
						break;
				}
			}
		}
		maxrprofile_s.SortXYZ();

		if (debug)
		{
			fprintf(f,"\n\n");
			for (i=0; i<maxrprofile_s.npoints; i++)
				fprintf(f,"%d %g %g\n",i,maxrprofile_s[i].x,maxrprofile_s[i].y);
		}
	}

	for (i=1; i<maxrprofile_s.npoints; i++)
	{
		double dx=maxrprofile_s[i].x-maxrprofile_s[i-1].x;
		double dy=maxrprofile_s[i].y-maxrprofile_s[i-1].y;
		maxrprofile_s[i].z=maxrprofile_s[i-1].z+sqrt(dx*dx+dy*dy);
	}

	if (maxrprofile_s.npoints==0)
	{
		if (verbose)
			fprintf(stderr,"Found no nose profile at all!\n");
		return 0;
	}

// smooth the profile by fitting and resampling a spline
{
	Spline *spline=new Spline(2);
	for (i=1; i<maxrprofile_s.npoints; i++)
	{
		float p[2]={maxrprofile_s[i].x,maxrprofile_s[i].y};
		spline->AddPoint(p,i);
	}

	spline->Fit();
	spline->SetAbscissaGeoDist(10);
	spline->Fit();
	float t0=spline->abscissa->Get(0)[0];
	float t1=spline->abscissa->Get(spline->abscissa->npoints-1)[0];
	int n=int(t1-t0)-2;
	CoordList *clist=spline->InterpolateRange(t0+1,t1-1,n);
	delete spline;
	spline=new Spline(*clist,0);
	delete clist;
	spline->SetAbscissaGeoDist(10);
	spline->Fit();
	t0=spline->abscissa->Get(0)[0];
	t1=spline->abscissa->Get(spline->abscissa->npoints-1)[0];
	n=int(t1-t0)-2;
	clist=spline->InterpolateRange(t0+1,t1-1,n);
	delete spline;
	spline=new Spline(*clist,0);
	delete clist;
	spline->SetAbscissaGeoDist(10);
	spline->Fit();
	t0=spline->abscissa->Get(0)[0];
	t1=spline->abscissa->Get(spline->abscissa->npoints-1)[0];
	n=int(t1-t0)-2;
	clist=spline->InterpolateRange(t0+1,t1-1,n);

	UnorderedPointSet ups;
	for (i=0; i<clist->npoints; i++)
	{
		Point3D p(clist->Get(i)[0],clist->Get(i)[1],0);
		ups.AddPoint(p);
	}

	for (i=1; i<ups.npoints; i++)
	{
		double dx=ups[i].x-ups[i-1].x;
		double dy=ups[i].y-ups[i-1].y;
		ups[i].z=ups[i-1].z+sqrt(dx*dx+dy*dy);
	}

	delete clist;
	delete spline;
	maxrprofile_s.CopyFrom(ups);
}

	UnorderedPointSet maxrprofile_k;
	
	if (debug)
		fprintf(f,"\n\n");

	int kscale=10; //5;
	for (i=2*kscale; i<maxrprofile_s.npoints-2*kscale; i++)
	{
		double x=maxrprofile_s[i].x;
		double y=maxrprofile_s[i].y;
		double xs=(maxrprofile_s[i+kscale].x-maxrprofile_s[i-kscale].x)/(maxrprofile_s[i+kscale].z-maxrprofile_s[i-kscale].z);
		double ys=(maxrprofile_s[i+kscale].y-maxrprofile_s[i-kscale].y)/(maxrprofile_s[i+kscale].z-maxrprofile_s[i-kscale].z);
		double xsl=(maxrprofile_s[i].x-maxrprofile_s[i-2*kscale].x)/(maxrprofile_s[i].z-maxrprofile_s[i-2*kscale].z);
		double xsr=(maxrprofile_s[i+2*kscale].x-maxrprofile_s[i].x)/(maxrprofile_s[i+2*kscale].z-maxrprofile_s[i].z);
		double ysl=(maxrprofile_s[i].y-maxrprofile_s[i-2*kscale].y)/(maxrprofile_s[i].z-maxrprofile_s[i-2*kscale].z);
		double ysr=(maxrprofile_s[i+2*kscale].y-maxrprofile_s[i].y)/(maxrprofile_s[i+2*kscale].z-maxrprofile_s[i].z);
		double xss=(xsr-xsl)/(maxrprofile_s[i+kscale].z-maxrprofile_s[i-kscale].z);
		double yss=(ysr-ysl)/(maxrprofile_s[i+kscale].z-maxrprofile_s[i-kscale].z);
		
		double kdenom=pow(xs*xs+ys*ys,1.5);
		double k=(kdenom!=0) ? (xs*yss-ys*xss)/kdenom : 0;
		if (xs<=0 || xsl<=0 || xsr<=0)
			k=0;
		//if (ys>=0 || xs<=0 || xsl<=0 || xsr<=0)
		//	k=0;
		Point3D p(x,y,k);
		maxrprofile_k.AddPoint(p);

		if (debug)
			fprintf(f,"%d %g %g %g %g %g %g %g\n",i,maxrprofile_s[i].x,maxrprofile_s[i].y,xs,ys,xss,yss,k);
	}

	// remove spikes
	for (i=1; i<maxrprofile_k.npoints-1; i++)
	{
		if (maxrprofile_k[i].z>0 && (maxrprofile_k[i-1].z<0 || maxrprofile_k[i+1].z<0))
			maxrprofile_k[i].z=0;
	}
	// remove out of order points
	UnorderedPointSet copy;
	for (i=1; i<maxrprofile_k.npoints-1; i++)
	{
		if (maxrprofile_k[i].x>maxrprofile_k[i-1].x && maxrprofile_k[i].x<maxrprofile_k[i+1].x)
			copy.AddPoint(maxrprofile_k[i]);
	}
	maxrprofile_k.CopyFrom(copy);

	// handle case where there are no points in maxrprofile_k
	if (maxrprofile_k.npoints==0)
		 maxrprofile_k.AddPoint(Point3D(maxrprofile_s[0].x,maxrprofile_s[0].y,0));

	// find local maxima
	UnorderedPointSet lmax;
	lmax.CopyFrom(maxrprofile_k);
	double lmax_mindist=5;	// previously 10 mm
	for (i=1; i<maxrprofile_k.npoints-1; i++)
	{
		for (j=1; j<maxrprofile_k.npoints-1; j++)
		{
			if (i==j)
				continue;
			if (fabs(maxrprofile_k[i].x-maxrprofile_k[j].x)>lmax_mindist)
				continue;
			if (maxrprofile_k[i].z<maxrprofile_k[j].z)
			{
				lmax[i].z=0;
				break;
			}
		}
	}

	double mink=0.01,kmax=0;
	int imax=0;
	if (verbose)
		fprintf(stderr,"nose bridge candidates (picking the best upto 60 mm from tip with minimum curvature 0.01, if curvature larger than 0.03 take the first):\n");
	for (i=1; i<lmax.npoints-1; i++)
	{
		if (lmax[i].x<=tip.x)
			continue;
		if (lmax[i].x-lmax[i-1].x>5 || lmax[i+1].x-lmax[i].x>5)
			continue;
		if (lmax[i].z<mink)
			continue;
		if (lmax[i].z>0)
		{
			if (verbose)
				fprintf(stderr,"\t%d %g %g %g\n",i,lmax[i].x,lmax[i].y,lmax[i].z);
			if (lmax[i].x<tip.x+60 && kmax<0.03 && lmax[i].z>kmax)
			{
				kmax=lmax[i].z;
				imax=i;
			}
		}
	}
	if (verbose)
		fprintf(stderr,"\tpicked: %d %g %g %g\n",imax,lmax[imax].x,lmax[imax].y,lmax[imax].z);

	if (imax==0)
	{
		if (verbose)
			fprintf(stderr,"found no nose bridge setting to tip+40\n");
		for (i=0; i<maxrprofile_k.npoints; i++)
			if (maxrprofile_k[i].x>tip.x+40)
				break;
		imax=(i>=maxrprofile_k.npoints) ? maxrprofile_k.npoints-1 : i;
		quality=0;
	}

	double nbx=maxrprofile_k[imax].x;
	double nby=maxrprofile_k[imax].y;

	// find local maximum below nose tip
	int imax2=0;
	if (verbose)
		fprintf(stderr,"below nose candidates (picking the first):\n");
	for (i=lmax.npoints-2; i>=1; i--)
	{
		if (lmax[i].x>=tip.x)
			continue;
		if (lmax[i].x-lmax[i-1].x>5 || lmax[i+1].x-lmax[i].x>5)
			continue;
		if (lmax[i].z<mink)
			continue;
		if (lmax[i].z>0)
		{
			if (verbose)
				fprintf(stderr,"\t\%d %g %g %g\n",i,lmax[i].x,lmax[i].y,lmax[i].z);
			if (imax2==0)
				imax2=i;
		}
	}

	if (imax2==0)
	{
		if (verbose)
			fprintf(stderr,"found no below tip setting to tip-10\n");
		for (i=maxrprofile.npoints-1; i>=0; i--)
			if (maxrprofile[i].x<tip.x-10)
				break;
		imax2=(i<=0) ? 0 : (i>=maxrprofile_k.npoints) ? maxrprofile_k.npoints-1 : i;
		quality=0;
	}

	double bnx=maxrprofile_k[imax2].x;
	double bny=maxrprofile_k[imax2].y;

	// find local minima
	UnorderedPointSet lmin;
	lmin.CopyFrom(maxrprofile_k);
	for (i=1; i<maxrprofile_k.npoints-1; i++)
	{
		for (j=1; j<maxrprofile_k.npoints-1; j++)
		{
			if (i==j)
				continue;
			if (fabs(maxrprofile_k[i].x-maxrprofile_k[j].x)>10)
				continue;
			if (maxrprofile_k[i].z>=maxrprofile_k[j].z)
			{
				lmin[i].z=0;
				break;
			}
		}
	}

	double maxk=-0.01;
	int imin=0;
	if (verbose)
		fprintf(stderr,"forehead candidates (picking the first):\n");
	if (imax==0)
		imax=1;
	for (i=imax; i<lmin.npoints-1; i++)
	{
		if (lmin[i].x-lmin[i-1].x>5 || lmin[i+1].x-lmin[i].x>5)
			continue;
		if (lmin[i].z<0)
		{
			if (verbose)
				fprintf(stderr,"\t\%d %g %g %g\n",i,lmin[i].x,lmin[i].y,lmin[i].z);
			if (imin==0)
				imin=i;
		}
	}

	if (imin==0)
	{
		if (verbose)
			fprintf(stderr,"found no forehead setting to tip+70\n");
		for (i=0; i<maxrprofile.npoints; i++)
			if (maxrprofile[i].x>tip.x+70)
				break;
		imin=(i>=maxrprofile_k.npoints) ? maxrprofile_k.npoints-1 : i;
		quality=0;
	}

	double fhx=maxrprofile_k[imin].x;
	double fhy=maxrprofile_k[imin].y;

	int j;
	int maxsupport=0;
	double mindtot=1000000;
	Point3D best_p,best_r;
	for (i=0; i<maxrprofile.npoints; i++)
	for (j=0; j<maxrprofile.npoints; j++)
	{
		double x=maxrprofile[i].x;
		double y=maxrprofile[i].y;
		if (fabs(y-(py1+(x-px1)*ry1))<5)
			continue;
		if (x<tip.x-50 || x>tip.x+100)
			continue;
		x=maxrprofile[j].x;
		y=maxrprofile[j].y;
		if (fabs(y-(py1+(x-px1)*ry1))<5)
			continue;
		if (x<tip.x-50 || x>tip.x+100)
			continue;

		if (maxrprofile[j].x<tip.x+50 || maxrprofile[i].x>tip.x-15)
			continue;
		// define a line through points i and j
		Point3D p=maxrprofile[i];
		Point3D r=maxrprofile[j]-maxrprofile[i];
		if (r.x!=0)
			r.y/=r.x;
		r.x=1;
		// get support
		int support_l=0,support_r=0,support=0;
		double dtot=0;
		int k;
		for (k=0; k<maxrprofile.npoints; k++)
		{
			if (maxrprofile[k].x<tip.x-50 || maxrprofile[k].x>tip.x+100)
				continue;
			if (maxrprofile[k].x>tip.x-15 && maxrprofile[k].x<tip.x+50)
				continue;

			double y=p.y+(maxrprofile[k].x-p.x)*r.y;
			double d=fabs(y-maxrprofile[k].y);
			if (d<2)
			{
				if (maxrprofile[k].x<tip.x)
					support_l++;
				else
					support_r++;
				support++;	
				dtot+=d;
			}
		}
		dtot/=(support>0) ? support : 1;
		if (support>maxsupport && support_l>2 && support_r>2)
		{
			maxsupport=support;
			mindtot=dtot;
			best_r=r;
			best_p=p;
		}
	}

	//this sets the point on the forehead to the nosebridge + 20 mm
	if (nosefitmethod==3)
	{
		for (i=0; i<maxrprofile_k.npoints; i++)
		{
			if (maxrprofile_k[i].x>nbx)
				break;
		}
		for (; i<maxrprofile_k.npoints; i++)
		{
			double dx=maxrprofile_k[i].x-nbx;
			double dy=maxrprofile_k[i].y-nby;
			double d=sqrt(dx*dx+dy*dy);
			if (d>20)
				break;
		}

		if (i>=maxrprofile_k.npoints)
			i=maxrprofile_k.npoints-1;
		best_p=Point3D(nbx,nby,0);
		best_r=maxrprofile_k[i]-Point3D(nbx,nby,0);
		if (best_r.x!=0)
			best_r.y/=best_r.x;
		best_r.x=1;
	}

	// this uses nose bridge and forehead 
	if (nosefitmethod==1)
	{
		best_p=Point3D(nbx,nby,0);
		best_r=Point3D(fhx,fhy,0)-best_p;
		if (best_r.x!=0)
			best_r.y/=best_r.x;
		best_r.x=1;
	}

	// this uses nose bridge and point below nose tip
	if (nosefitmethod==2)
	{
		best_p=Point3D(bnx,bny,0);
		best_r=Point3D(nbx,nby,0)-best_p;
		if (best_r.x!=0)
			best_r.y/=best_r.x;
		best_r.x=1;

		//best_r.y=0;
	}

	double x,y;

	if (debug)
	{
		fprintf(f,"\n\n");

		for (x=-50; x<150; x++)
			fprintf(f,"%g %g\n",x,best_p.y+(x-best_p.x)*best_r.y);

		fclose(f);
	}

	x=(ry1-best_r.y!=0) ? (best_p.y-py1-best_p.x*best_r.y+px1*ry1)/(ry1-best_r.y) : 50;
	y=py1+(x-px1)*ry1;
	y=best_p.y+(x-best_p.x)*best_r.y;
	double phi=-atan2(best_r.y,best_r.x);

	if (nosefitmethod==1)	// refine orientation estimation based on nose bridge
	{
		phi=0;
		int ni=0;
		double nx=0,ny=0;
		for (i=0; i<maxrprofile_k.npoints-1; i++)
		{
			if (maxrprofile_k[i].x<tip.x+10 || maxrprofile_k[i].x>nbx-10)
				continue;
			nx+=maxrprofile_k[i+1].x-maxrprofile_k[i].x;
			ny+=maxrprofile_k[i+1].y-maxrprofile_k[i].y;
			ni++;
		}
		if (ni>0)
		{
			nx/=ni;
			ny/=ni;
			phi=-atan2(ny,nx);
			phi-=M_PI/6;
			if (verbose)
				fprintf(stderr,"==>>phi=%g\n",phi);
		}
	}

	// this derives nose tip from nose bridge dent
	if (nosefitmethod==1 || nosefitmethod==3)
	{
		tip.x=nbx*cos(phi)-nby*sin(phi) -40;
		tip.y=nbx*sin(phi)+nby*cos(phi);
	}

	// this derives nose tip from point below nose
	if (nosefitmethod==2)
	{
		tip.x=bnx*cos(phi)-bny*sin(phi)+5;
		tip.y=bnx*sin(phi)+bny*cos(phi);
	}

	nphi+=phi;

}
	cosphi=cos(nphi);
	sinphi=sin(nphi);


	if (debug)
	{
		FILE *f=fopen("dumpprof.txt","w");
		for (i=0; i<maxprofile.npoints; i++)
		{
			double x=cosphi*maxprofile[i].x-sinphi*maxprofile[i].y-tip.x;
			double y=sinphi*maxprofile[i].x+cosphi*maxprofile[i].y-tip.y;
			fprintf(f,"%g %g\n",x,y);
		}
		fclose(f);
	}

	if (debug)
	{
		printf("\n\n");
		for (i=0; i<maxprofile.npoints; i++)
		{
			double x=cosphi*maxprofile[i].x-sinphi*maxprofile[i].y-tip.x;
			double y=sinphi*maxprofile[i].x+cosphi*maxprofile[i].y-tip.y;
			printf("%g %g\n",x,y);
		}
		printf("\n\n");
	}
	if (debug)
	{
		printf("\n\n");
		for (i=0; i<maxrprofile.npoints; i++)
		{
			double x=maxrprofile[i].x-tip.x;
			double y=maxrprofile[i].y-tip.y;
			printf("%g %g\n",x,y);
		}
		printf("\n\n");
	}

	if (verbose)
		fprintf(stderr,"fitnose2profile: tip=%g %g\n",tip.x,tip.y);
		
	nphi=-nphi;
	nx=cos(nphi)*tip.x-sin(nphi)*tip.y;
	ny=sin(nphi)*tip.x+cos(nphi)*tip.y;
	nl=50;
	nh=30;
	nt=10;

	quality*=1-(maxholesize-1)/30;
	if (quality<0)
		quality=0;
	else if (quality>1)
		quality=1;

	return quality;
}
}
