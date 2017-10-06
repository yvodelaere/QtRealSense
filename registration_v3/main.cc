/*
	register 3D faces to common coordinate system using symmetry axis and nose geometry
	main program file
	Author Luuk Spreeuwers, University of Twente
	Date: 2006-2013

	26 March, 2013: added support for outputting multiple fine registrations as a single sfi image (-J option)
			currently 5x5 registrations with +/- 2 mm shift in origin and +/- 0.04 radians tilt variation
*/

#include <stdio.h>
#include <stdlib.h>
#ifdef linux
#include <unistd.h>
#include <malloc.h>
#endif
#include <string.h>
#include <math.h>
#include <time.h>

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

using namespace utw3dface;

namespace utw3dface {
extern int verbose;
extern int debug;
extern int Profile;
extern int nthreads;
extern double profile_theta;
}

const char *usage="register [OPTIONS]\n\
OPTIONS:\n\
\t-a file     output abs point file\n\
\t-b 0|1      remove background [1=yes]\n\
\t-c quality  compress/decompress result with quality 0..100 [90]\n\
\t-D          report in degrees instead of radians (together with -t)\n\
\t-F          assume near frontal images\n\
\t-g 0|1      remove glare (reflections) [0=no]\n\
\t-h 0|1      perform hole filling [1=yes]\n\
\t-i file     read abs point file or file with 1 xyz coords/line (pts) \n\
\t-I file     read ppm texture file\n\
\t-J          join number of registrations into 1 sfi file, currently 25 fine registrations\n\
\t-L          only use left part to find nose\n\
\t-m 0|1      apply elliptical mask [1=yes]\n\
\t-M float    max shift for motion correction (-S option) [7.5 mm]\n\
\t-n 0|1      normalise texture output image[1=yes]\n\
\t-N 0|1      choose nose fit model: 0=nose-tip+bridge slope,1=dent+bridge slope,2=below+dent-below slope, 3=dent+dent-forehead slope [0]\n\
\t-o file     write sfi range image file\n\
\t-O file     write sfi texture file (requires -I)\n\
\t-p file     write profile data\n\
\t-P theta    attempt profile processing; theta in through plane rotation in deg\n\
\t-q          print input file name and quality numbers\n\
\t-r float    resolution in xy plane in [0.67 mm]\n\
\t-R          only use right part to find nose\n\
\t-s 0|1|2+   perform spike removal 0=no, 1=yes, 2+=smoothing with n points [1=yes=25 points]\n\
\t-S float    symmetrize image (line based motion correction) number is average over rows [0 mm]\n\
\t-T float    threshold for motion correction, average detected shift must be > value [0 mm]\n\
\t-t          report rotation angles and translation\n\
\t-v          increase verbosity\n\
\t-X int      skip X: only use every Xth pixel from input abs file [1]\n\
\t-Y int      skip Y: only use every Yth pixel from input abs file [1]\n\
\t-z int      multithreaded [#threads=1]\n\
";

int main(int argc,char **argv)
{
	char *ifile=NULL;
	char *ofile=NULL;
	char *pfile=NULL;
	char *absfile=NULL;
	char *textureifile=NULL;
	char *textureofile=NULL;
	int spikeremoval=1;
	int holefilling=1;
	double resolution=0.67;
	int normalisetexture=1;
	int ellipticalmask=1;
	int reflectionremoval=0;
	int backgroundremoval=1;
	int printquality=0;
	int compressionquality=90;
	int compress=0;
	int reporttrafo=0;
	int degrees=0;
	int nosefitmethod=0;
	int LR=0;
	double symmetrize=0;
	double motion_threshold=0.0;
	int nearfrontal=0;
	double maxshift=7.5;
	int skipX=1;
	int skipY=1;
	int joinmultiple=0;

	int c;

	while ((c=getopt(argc,argv,"z:JX:Y:FM:S:P:LRN:Dtc:qb:m:n:I:O:a:dr:h:s:o:i:vp:g:T:"))!=-1)
	{
		switch (c)
		{
			case 'z': sscanf(optarg,"%d",&nthreads); break;
			case 'J': joinmultiple=1; break;
			case 'X': sscanf(optarg,"%d",&skipX); break;
			case 'Y': sscanf(optarg,"%d",&skipY); break;
			case 'L': LR=1; break;
			case 'R': LR=2; break;
			case 'N':
				sscanf(optarg,"%d",&nosefitmethod);
				break;
			case 'D':
				degrees=1;
				break;
			case 't':
				reporttrafo=1;
				break;
			case 'c':
				compress=1;
				sscanf(optarg,"%d",&compressionquality);
				break;
			case 'q':
				printquality=1;
				break;
			case 'n':
				sscanf(optarg,"%d",&normalisetexture);
				break;
			case 'd':
				debug++;
				break;
			case 'r':
				sscanf(optarg,"%lf",&resolution);
				break;
			case 's':
				sscanf(optarg,"%d",&spikeremoval);
				break;
			case 'b':
				sscanf(optarg,"%d",&backgroundremoval);
				break;
			case 'm':
				sscanf(optarg,"%d",&ellipticalmask);
				break;
			case 'h':
				sscanf(optarg,"%d",&holefilling);
				break;
			case 'g':
				sscanf(optarg,"%d",&reflectionremoval);
				break;
			case 'I':
				textureifile=optarg;
				break;
			case 'O':
				textureofile=optarg;
				break;
			case 'i':
				ifile=optarg;
				break;
			case 'a':
				absfile=optarg;
				break;
			case 'o':
				ofile=optarg;
				break;
			case 'v':
				verbose++;
				break;
			case 'p':
				pfile=optarg;
				break;
			case 'P':
				Profile=1;
				sscanf(optarg,"%lf",&profile_theta);
				profile_theta*=M_PI/180.0;
				break;
			case 'S':
				sscanf(optarg,"%lf",&symmetrize);
				break;
			case 'T':
				sscanf(optarg,"%lf",&motion_threshold);
				break;
			case 'M':
				sscanf(optarg,"%lf",&maxshift);
				break;
			case 'F':
				nearfrontal=1;
				break;
			default:
				fprintf(stderr,"%s",usage);
				exit(1);
		}
	}

	if (nthreads<0 || nthreads>128)
	{
		fprintf(stderr,"number of threads must be in range 1..128\n");
		exit(1);
	}

	clock_t t0=clock();

	OrderedPointSet ops;
	if (ifile==NULL)
	{
		fprintf(stderr,"-i is a required option\n");
		exit(1);
	}
	
	ops.Read(ifile);
	//ops.Resample(1);

	if (debug)
		ops.WriteVRML("foo.wrl");

	UnorderedPointSet ups;
	int x,y;
	for (y=0; y<ops.height; y++)
	for (x=0; x<ops.width; x++)
	{
		if (ops.ValidPoint(x,y) && x%skipX==0 && y%skipY==0)
			ups.AddPoint(ops.GetPoint(x,y));
	}

	if (verbose)
		fprintf(stderr,"Read %d points\n",ups.npoints);

	clock_t t1=clock();
        if (verbose)
                fprintf(stderr,"Reading and converting input data took %.6f seconds\n",double(t1-t0)/CLOCKS_PER_SEC);

	Register reg(ups,resolution,holefilling,spikeremoval,ellipticalmask,reflectionremoval,backgroundremoval,nosefitmethod,LR,symmetrize,motion_threshold,maxshift,nearfrontal);

	//reg.FitNose(2);
	//reg.PostProc();

	if (compress)
	{
		reg.ri.CompressDecompress(compressionquality);
		fprintf(stderr,"compressed size=%d\n",reg.ri.jpeg_buf_size);
	}

	if (printquality)
		printf("%s %g %g %g %g %g %g %g\n",ifile,reg.overallQuality,reg.spikesQuality,reg.holesQuality,reg.registrationQuality,reg.noiseQuality,reg.completenessQuality,reg.noseQuality);

	if (ofile!=NULL)
	{
		if (joinmultiple)
		{
			RangeImage rij(reg.ri.width*5,reg.ri.height*5);
			int x,y;
			for (y=0; y<reg.ri.width; y++)
			for (x=0; x<reg.ri.height; x++)
				rij[2*reg.ri.height+y][2*reg.ri.width+x]=reg.ri[y][x];

			verbose=0;
			Point3D o=reg.ri.o;
			Point3D v=reg.ri.v;
			Point3D u=reg.ri.u;
			double ou=reg.ri.ou;
			double ov=reg.ri.ov;
			int itilt=5,ishift=5;
			double tiltstep=0.02,tilt0=-0.04;
			double shiftstep=1,shift0=-2;
			for (itilt=0; itilt<5; itilt++)
			for (ishift=0; ishift<5; ishift++)
			{
				double tilt=itilt*tiltstep+tilt0;
				double shift=ishift*shiftstep+shift0;
				reg.ri.o=o-shift*v;
				reg.ri.v=cos(tilt)*v+sin(tilt)*(u^v);
				reg.ri.calculateT();
				reg.ri.AccumulateDepth(reg.ups);
				reg.PostProc();
				for (y=0; y<reg.ri.height; y++)
				for (x=0; x<reg.ri.width; x++)
					rij[itilt*reg.ri.height+y][ishift*reg.ri.width+x]=reg.ri[y][x];
			}

			rij.WriteSFI(ofile);
		}
		else
			reg.ri.WriteSFI(ofile);
	}

	
	if (reporttrafo)
	{
		float m[9],t[3];
		reg.GetRotationMatrix(m);
		reg.GetTranslationVector(t);
		//double phi=asin(m[2]);
		//double theta=asin(m[1]/cos(phi));
		//double gamma=-asin(m[7]/cos(phi));
		double phi=asin(m[1]);
		double theta=asin(m[2]/cos(phi));
		double gamma=asin(m[5]/cos(phi));

if (verbose || debug) // some extra code to check transformation and report origin in original coordinates
{
	float m00=m[0];
	float m10=m[1];
	float m20=m[2];
	float m01=m[3];
	float m11=m[4];
	float m21=m[5];
	float m02=m[6];
	float m12=m[7];
	float m22=m[8];

	float det = m00*m11*m22 + m01*m12*m20 + m02*m10*m21 - m00*m12*m21 - m01*m10*m22 - m02*m11*m20;

     	float i00 = (m11*m22 - m12*m21)/det;
     	float i01 = (m02*m21 - m01*m22)/det;
     	float i02 = (m01*m12 - m02*m11)/det;
     	float i10 = (m12*m20 - m10*m22)/det;
     	float i11 = (m00*m22 - m02*m20)/det;
     	float i12 = (m02*m10 - m00*m12)/det;
     	float i20 = (m10*m21 - m11*m20)/det;
     	float i21 = (m01*m20 - m00*m21)/det;
     	float i22 = (m00*m11 - m01*m10)/det;
	
	float p00=i00*m00+i10*m01+i20*m02;
	float p10=i00*m10+i10*m11+i20*m12;
	float p20=i00*m20+i10*m21+i20*m22;

	float p01=i01*m00+i11*m01+i21*m02;
	float p11=i01*m10+i11*m11+i21*m12;
	float p21=i01*m20+i11*m21+i21*m22;

	float p02=i02*m00+i12*m01+i22*m02;
	float p12=i02*m10+i12*m11+i22*m12;
	float p22=i02*m20+i12*m21+i22*m22;

	float ti0=-(i00*t[0]+i10*t[1]+i20*t[2]);
	float ti1=-(i01*t[0]+i11*t[1]+i21*t[2]);
	float ti2=-(i02*t[0]+i12*t[1]+i22*t[2]);

if (verbose>1)
{
	fprintf(stderr,"****************************************************\n");
	fprintf(stderr,"%4.3f %4.3f %4.3f\n",m00,m10,m20);
	fprintf(stderr,"%4.3f %4.3f %4.3f\n",m01,m11,m21);
	fprintf(stderr,"%4.3f %4.3f %4.3f\n",m02,m12,m22);
	fprintf(stderr,"\n");
	fprintf(stderr,"det=%4.3f\n\n",det);
	fprintf(stderr,"%4.3f %4.3f %4.3f\n",i00,i10,i20);
	fprintf(stderr,"%4.3f %4.3f %4.3f\n",i01,i11,i21);
	fprintf(stderr,"%4.3f %4.3f %4.3f\n",i02,i12,i22);
	fprintf(stderr,"\n");
	fprintf(stderr,"%4.3f %4.3f %4.3f\n",p00,p10,p20);
	fprintf(stderr,"%4.3f %4.3f %4.3f\n",p01,p11,p21);
	fprintf(stderr,"%4.3f %4.3f %4.3f\n",p02,p12,p22);
	fprintf(stderr,"\n");
	fprintf(stderr,"origin in original coordinates %4.3f %4.3f %4.3f\n",ti0,ti1,ti2);
	fprintf(stderr,"****************************************************\n");
}
else
	fprintf(stderr,"origin in original coordinates %4.3f %4.3f %4.3f\n",ti0,ti1,ti2);

if (debug)
{
double x,y,z;
for (x=-100; x<=100; x++)
	ups.AddPoint(Point3D(x+ti0,ti1,ti2));
for (y=-100; y<=100; y++)
	ups.AddPoint(Point3D(ti0,y+ti1,ti2));
for (z=-100; z<=100; z++)
	ups.AddPoint(Point3D(ti0,ti1,z+ti2));
ups.WriteVRML("origin.wrl");
}
}

		if (verbose)
		{
			fprintf(stderr,"phi=%g theta=%g gamma=%g (in deg: %g %g %g)\n",phi,theta,gamma,phi*180/M_PI,theta*180/M_PI,gamma*180/M_PI);
			fprintf(stderr,"t=(%g %g %g)\n",t[0],t[1],t[2]);
			fprintf(stderr,"quality=%g\n",reg.overallQuality);
		}
		else
		{
			if (degrees)
				printf("%g %g %g %g %g %g %g\n",phi*180/M_PI,theta*180/M_PI,gamma*180/M_PI,t[0],t[1],t[2],reg.overallQuality);
			else
				printf("%g %g %g %g %g %g %g\n",phi,theta,gamma,t[0],t[1],t[2],reg.overallQuality);
		}
	}

	if (absfile!=NULL)
	{
		RangeImage ri(reg.ri.o,reg.ri.u,reg.ri.v,1,1,110,170,220,260);
		ri.AccumulateDepth(ups);
		if (spikeremoval)
			filter(ri,ups,5);
		if (ellipticalmask)
		{
			// mask out an ellipse
			double rx=55/ri.du;
			double ry=70/ri.dv;
			double cx=ri.ou;
			double cy=ri.ov-20/ri.dv;
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

		RangeImage riO(reg.ri.o,reg.ri.u,-1*(reg.ri.v),1,1,0,0,220,260);
		int x,y;
		for (y=0; y<ops.height; y++)
		for (x=0; x<ops.width; x++)
		{
			int index=y*ops.width+x;
			if (ops.ValidPoint(x,y))
			{
				Point3D p=ops.GetPoint(x,y);
				Point3D q=ri.TransformPoint(p);
				double flag;
				double z=ri.GetDepth(q.x,q.y,flag);
				if (flag==0 || (spikeremoval && fabs(z-q.z)>5))
				{
					ops.flags[index]=0;
					ops.X[index]=0;
					ops.Y[index]=0;
					ops.Z[index]=0;
				}
				else
				{
					Point3D q=riO.TransformPoint(p);
					ops.X[index]=q.x;
					ops.Y[index]=q.y;
					ops.Z[index]=q.z;
				}
			}
		}

		ops.WriteAbs(absfile);
	}

	if (textureofile!=NULL)
	{
		if (textureifile==NULL)
		{
			fprintf(stderr,"Warning: texture output can only be generated if texture input is given\n");
			exit(0);
		}
		RangeImage texture;
		texture.ReadPNM(textureifile);

		if (debug)
			texture.WritePGM("textureraw.pgm",1);
		
		RangeImage texreg(reg.ri.o,reg.ri.u,reg.ri.v,1,1,55,85,110,130);
		texreg.Clear();

		int x,y;
		for (y=0; y<ops.height; y++)
		for (x=0; x<ops.width; x++)
		{
			int index=y*ops.width+x;
			if (ops.ValidPoint(x,y))
			{
				Point3D p=ops.GetPoint(x,y);
				Point3D q=reg.ri.TransformPoint(p);
				double flag;
				double z=reg.ri.GetDepth(q.x,q.y,flag);
				if (flag!=0 && fabs(z-q.z)<=5)
				{
					texreg.Pixel(int(q.x+0.5),int(q.y+0.5))=texture.Pixel(x,y);
					texreg.Flag(int(q.x+0.5),int(q.y+0.5))=1;
				}
			}
		}

		Nose nosef(-1000,-1000,0);
		fillholes(texreg,nosef);

		if (ellipticalmask)
		{
			double rx=55/reg.ri.du;
			double ry=70/reg.ri.dv;
			double cx=reg.ri.ou;
			double cy=reg.ri.ov-20/reg.ri.dv;
			for (y=0; y<reg.ri.height; y++)
			for (x=0; x<reg.ri.width; x++)
			{
				if ((x-cx)*(x-cx)/(rx*rx)+(y-cy)*(y-cy)/(ry*ry)>1)
				{
					texreg.Flag(x,y)=0;
					texreg.Pixel(x,y)=0;
				}
				else
					texreg.Flag(x,y)=1;
			}
		}

		if (normalisetexture)
		{
			int nbins=256;
			double min=100000,max=-100000;
			int n=0;
			for (y=0; y<texreg.height; y++)
			for (x=0; x<texreg.width; x++)
			{
				if (texreg.Flag(x,y)==0)
					continue;
				double l=texreg.Pixel(x,y);
				if (l<min)
					min=l;
				else if (l>max)
					max=l;
				n++;
			}
			
			double binsize=(max-min)/nbins;
			int bin;
			int hist[nbins];
			for (bin=0; bin<nbins; bin++)
				hist[bin]=0;

			for (y=0; y<texreg.height; y++)
			for (x=0; x<texreg.width; x++)
			{
				if (texreg.Flag(x,y)==0)
					continue;
				double l=texreg.Pixel(x,y);
				bin=int((l-min)/binsize);
				if (bin<0)
					bin=0;
				else if (bin>=nbins)
					bin=nbins-1;
				hist[bin]++;
			}

			double sum=0;
			double binstart[nbins],binend[nbins];
			for (bin=0; bin<nbins; bin++)
			{
				binstart[bin]=sum/n;
				sum+=hist[bin];
				binend[bin]=sum/n;
			}

			for (y=0; y<texreg.height; y++)
			for (x=0; x<texreg.width; x++)
			{
				if (texreg.Flag(x,y)==0)
				{
					texreg.Pixel(x,y)=0;
					continue;
				}
				double l=texreg.Pixel(x,y);
				double frac=(l-min)/binsize;
				bin=int(frac);
				frac-=bin;
				if (bin<0)
					bin=0;
				else if (bin>=nbins)
					bin=nbins-1;
				texreg.Pixel(x,y)=binstart[bin]*(1-frac)+binend[bin]*frac;
			}
		}

		texreg.WriteSFI(textureofile);
	}

	if (pfile!=NULL)
	{
		FILE *f=fopen(pfile,"w");
		if (f==NULL)
		{
			fprintf(stderr,"failed to open %s for output\n",pfile);
			exit(1);
		}

		reg.GetProfile();
		reg.GetMaxProfile();

		int i;
		for (i=0; i<reg.profile.npoints; i++)
			fprintf(f,"%g %g %g\n",reg.profile[i].x,reg.profile[i].y,reg.profile[i].z);

		fprintf(f,"\n\n");
		for (i=0; i<reg.maxprofile.npoints; i++)
			fprintf(f,"%g %g %g\n",reg.maxprofile[i].x,reg.maxprofile[i].y,reg.maxprofile[i].z);

		fprintf(f,"\n\n");
		Nose nose(0,0,0);
		nose.SetSize(reg.nl,reg.nh,reg.nt,10,20); 
		for (x=int(nose.Start()); x<int(nose.End()); x++)
		{
			Point3D p=nose.XY(x);
			fprintf(f,"%g %g\n",p.x,p.y);
		}

		fclose(f);
	}

}
