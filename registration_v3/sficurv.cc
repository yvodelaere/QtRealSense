#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include "rangeimage.h"

using namespace utw3dface;

const char *usage="sfidiff [OPTIONS]\n\
OPTIONS:\n\
\t-c float      clip curvature at +/- value [0.1]\n\
\t-d            debug\n\
\t-i file       input sfi file\n\
\t-t 0-5  type or curvature: 0=gaussian, 1=mean, 2=min, 3=max 4=|fxx,fyy| 5=fx,fy [0]\n\
\t-o file       output sfi file\n\
\t-r float      xy resolution in mm [1]\n\
\t-s sigma      first filter image with gaussian sigma in mm [0]\n\
\t-v            verbose\n\
";

int verbose=0;
int debug=0;

main(int argc,char **argv)
{
	char *ifile=NULL, *ofile=NULL;
	int curvtype=0;
	double xyres=1.0;
	double sigma=0;
	double clip=0.1;
	int c;

	while ((c=getopt(argc,argv,"di:t:o:s:vr:c:"))!=-1)
	{
		switch(c)
		{
			case 'i': ifile=optarg; break;
			case 't': sscanf(optarg,"%d",&curvtype); break;
			case 'o': ofile=optarg; break;
			case 'r': sscanf(optarg,"%lf",&xyres); break;
			case 's': sscanf(optarg,"%lf",&sigma); break;
			case 'c': sscanf(optarg,"%lf",&clip); break;
			case 'v': verbose++; break;
			case 'd': debug++; break;
			default: fprintf(stderr,usage); exit(1);
		}
	}

	if (ifile==NULL)
	{
		fprintf(stderr,"-i is a required option\n");
		exit(1);
	}

	RangeImage im;
	im.ReadSFI(ifile);

	// create a kernel
	sigma/=xyres;
	int nk=int(3*sigma)+1;

	double kernel[nk];
	double sum=0;
	int i;
	for (i=0; i<nk; i++)
	{
		kernel[i]=(sigma==0) ? 1 : exp(i*i/(2*sigma*sigma));
		sum+=kernel[i];
	}
	for (i=0; i<nk; i++)
		kernel[i]/=sum;

	RangeImage imf;
	imf.CopyDataFrom(im);

	// filter image
	int x,y;
	double norm;
	for (y=0; y<im.height; y++)
	for (x=0; x<im.width; x++)
	{
		sum=im.Pixel(x,y);
		for (i=1; i<nk; i++)
		{
			sum+=(x+i<im.width) ? im.Pixel(x+i,y)*kernel[i] : im.Pixel(im.width-1,y)*kernel[i];
			sum+=(x-i>=0) ? im.Pixel(x-i,y)*kernel[i] : im.Pixel(0,y)*kernel[i];
		}
		imf.Pixel(x,y)=sum;
	}

	im.CopyDataFrom(imf);
	for (y=0; y<im.height; y++)
	for (x=0; x<im.width; x++)
	{
		sum=im.Pixel(x,y);
		for (i=1; i<nk; i++)
		{
			sum+=(y+i<im.height) ? im.Pixel(x,y+i)*kernel[i] : im.Pixel(x,im.height)*kernel[i];
			sum+=(y-i>=0) ? im.Pixel(x,y-i)*kernel[i] : im.Pixel(x,0)*kernel[i];
		}
		imf.Pixel(x,y)=sum;
	}
	im.CopyDataFrom(imf);

	if (debug)
		im.WriteSFI("imf.sfi");

	// create derivatives
	// since f=z, fz=1, fzz=0, fzx=0, fzy=0
	double fx,fy,fz=1,fxx,fyy,fzz=0,fxy,fxz=0,fyz=0,fn;
	double g,m,det,kmax,kmin,curvature;

	for (y=2; y<im.height-2; y++)
	for (x=2; x<im.width-2; x++)
	{
		fx=(im.Pixel(x+1,y)-im.Pixel(x-1,y))/(2*xyres);
		fy=(im.Pixel(x,y+1)-im.Pixel(x,y-1))/(2*xyres);
		fxx=(im.Pixel(x+2,y)-2*im.Pixel(x,y)+im.Pixel(x-2,y))/(4*xyres);
		fyy=(im.Pixel(x,y+2)-2*im.Pixel(x,y)+im.Pixel(x,y-2))/(4*xyres);
		fxy=(im.Pixel(x+1,y+1)-im.Pixel(x-1,y+1)-im.Pixel(x+1,y-1)+im.Pixel(x-1,y-1))/(4*xyres);
		fn=(fx*fx+fy*fy+fz*fz);

		switch (curvtype)
		{
	        	case 0: // gaussian curvature
				curvature=(fx*fx*(fyy*fzz-fyz*fyz) + fy*fy*(fxx*fzz-fxz*fxz) + fz*fz*(fxx*fyy-fxy*fxy) + 2*(fx*fy*(fxz*fyz-fxy*fzz) + fy*fz*(fxy*fxz-fyz*fxx) + fx*fz*(fxy*fyz-fxz*fyy)))/(fn*fn); 
				break; 
			case 1: //mean curvature 
				curvature=((fyy+fzz)*fx*fx + (fxx+fzz)*fy*fy + (fxx+fyy)*fz*fz - 2*fx*fy*fxy - 2*fx*fz*fxz- 2*fy*fz*fyz)/(fn*sqrt(fn)); 
				break; 
			case 2: // min curvature 
				g=(fx*fx*(fyy*fzz-fyz*fyz) + fy*fy*(fxx*fzz-fxz*fxz) + fz*fz*(fxx*fyy-fxy*fxy) + 2*(fx*fy*(fxz*fyz-fxy*fzz) + fy*fz*(fxy*fxz-fyz*fxx) + fx*fz*(fxy*fyz-fxz*fyy)))/(fn*fn); 
				m=((fyy+fzz)*fx*fx + (fxx+fzz)*fy*fy + (fxx+fyy)*fz*fz - 2*fx*fy*fxy - 2*fx*fz*fxz - 2*fy* fz*fyz)/(fn*sqrt(fn));
				det=m*m-g; 
				if (det>=0) 
				{ 
					det=sqrt(det); 
					kmin=m-det; 
					kmax=m+det; 
					curvature=fabs(kmin)<fabs(kmax) ? kmin : kmax; 
				} 
				else 
					curvature=0; 
				break;
			case 3: // max curvature
				g=(fx*fx*(fyy*fzz-fyz*fyz) + fy*fy*(fxx*fzz-fxz*fxz) + fz*fz*(fxx*fyy-fxy*fxy) + 2*(fx*fy*(fxz*fyz-fxy*fzz) + fy*fz*(fxy*fxz-fyz*fxx) + fx*fz*(fxy*fyz-fxz*fyy)))/(fn*fn); 
				m=((fyy+fzz)*fx*fx + (fxx+fzz)*fy*fy + (fxx+fyy)*fz*fz - 2*fx*fy*fxy - 2*fx*fz*fxz - 2*fy* fz*fyz)/(fn*sqrt(fn));
				det=m*m-g; 
				if (det>=0) 
				{ 
					det=sqrt(det); 
					kmin=m-det; 
					kmax=m+det; 
					curvature=fabs(kmin)>fabs(kmax) ? kmin : kmax; 
				} 
				else 
					curvature=0; 
				break;
			case 4: // magnitude of fxx,fyy
				curvature=sqrt(fxx*fxx+fyy*fyy);
				break;
			case 5: // alternating fx and fy
				curvature=(x&0x01) ? fx : fy;
				break;
		}

		imf.Pixel(x,y)=(curvature>clip || curvature<-clip) ? 0 : curvature;
	}
	
	for (y=0; y<im.height; y++)
		imf.Pixel(0,y)=imf.Pixel(1,y)=imf.Pixel(imf.width-1,y)=imf.Pixel(imf.width-2,y)=0;
	for (x=0; x<im.width; x++)
		imf.Pixel(x,0)=imf.Pixel(x,1)=imf.Pixel(x,imf.height-1)=imf.Pixel(x,imf.height-2)=0;
	for (x=0; x<im.width; x++)
		imf.Pixel(x,imf.height-3)=imf.Pixel(x,imf.height-4)=0;

	if (ofile)
		imf.WriteSFI(ofile);
}
