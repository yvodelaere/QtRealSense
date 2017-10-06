#include <stdio.h>
#ifdef linux
#include <values.h>
#endif
#include <math.h>
#include "compat.h"

#include "roi.h"
#include "filter.h"

namespace utw3dface {

extern int verbose;
extern int debug;

double filter(RangeImage &ri,UnorderedPointSet &ups,double maxd,int lesssmoothing)
{
	int minn=(lesssmoothing) ? lesssmoothing : 25;

	double meanflag=0,maxflag=0,n=0;
	int x,y;
	for (y=0; y<ri.height; y++)
	for (x=0; x<ri.width; x++)
	{
		double flag=ri.Flag(x,y);
		if (flag)
		{
			meanflag+=flag;
			n++;
		}
		if (flag>maxflag)
			maxflag=flag;
	}
	if (n>0)
		meanflag/=n;

	if (debug)
		fprintf(stderr,"meanflag=%g maxflag=%g\n",meanflag,maxflag);

	double *f=new double[ri.width*ri.height];
	for (y=0; y<ri.height; y++)
	for (x=0; x<ri.width; x++)
	{
		double sum=ri.Flag(x,y)*ri.Pixel(x,y);
		double n=ri.Flag(x,y);
		int r;
		for (r=1; r<ri.width; r++)
		{
			int xx,yy;
			for (xx=-r+1; xx<=r-1; xx++)
			{
				if (x+xx<0 || x+xx>ri.width)
					continue;
				if (y-r>0)
				{
					n+=ri.Flag(x+xx,y-r);
					sum+=ri.Flag(x+xx,y-r)*ri.Pixel(x+xx,y-r);
				}
				if (y+r<ri.height)
				{
					n+=ri.Flag(x+xx,y+r);
					sum+=ri.Flag(x+xx,y+r)*ri.Pixel(x+xx,y+r);
				}
			}
			for (yy=-r+1; yy<=r-1; yy++)
			{
				if (y+yy<0 || y+yy>=ri.height)
					continue;
				if (x-r>0)
				{
					n+=ri.Flag(x-r,y+yy);
					sum+=ri.Flag(x-r,y+yy)*ri.Pixel(x-r,y+yy);
				}
				if (x+r<ri.width)
				{
					n+=ri.Flag(x+r,y+yy);
					sum+=ri.Flag(x+r,y+yy)*ri.Pixel(x+r,y+yy);
				}
			}
			if (n>meanflag*minn)
				break;
		}

		f[y*ri.width+x]=sum/n;
	}

	memcpy(ri.pixel,f,sizeof(double)*ri.width*ri.height);

	delete [] f;

	if (debug)
		ri.WritePGM("filtered.pgm",1);
	
	UnorderedPointSet upsf;
	int i;
	int upsnpoints=0;
	for (i=0; i<ups.npoints; i++)
	{
		Point3D q=ri.TransformPoint(ups[i]);
		double flag;
		if (fabs(ri.GetDepth(q.x,q.y,flag)-q.z)<maxd)
			upsf.AddPoint(ups[i]);
		if (flag)
			upsnpoints++;
	}


	ups.npoints=0;
	for (i=0; i<upsf.npoints; i++)
		ups.AddPoint(upsf[i]);

	if (debug)
		fprintf(stderr,"%d points after filtering\n",upsf.npoints);

	ri.AccumulateDepth(ups);
	if (debug)
		ri.WritePGM("filtered2.pgm",1);
	if (debug)
		ri.WriteFlagPGM("filtered2flags.pgm");

	return (upsf.npoints<=0) ? 0 : double(upsf.npoints)/upsnpoints;
}

void trace_boundary(RangeImage &edges,RangeImage &xedges,RangeImage &yedges,int startx,int starty,double mingrad=1)
{
	int w=xedges.width;
	int h=xedges.height;
	int x=startx;
	int y=starty;
	
	// directions: 0=w,1=nw,2=n,3=ne,4=e,5=se,6=s,7=se
	// starting direction e is the correct setting for the left most pixel
	double gradx=xedges.Pixel(x,y);
	double grady=yedges.Pixel(x,y);
	int initdir=(fabs(gradx)>fabs(grady)) ? 2 : 4;
	int grad=(initdir==2) ? (gradx>0) ? 1 : -1 : (grady>0) ? 1 : -1;

	edges.Pixel(x,y)=100;

	int dir=initdir;
	for (;;)
	{
		int newdir,xx,yy;
		double maxg=0;
		int maxgdir=dir;
		int maxg_xx=x,maxg_yy=y;
		int d;
		for (d=dir-3; d<dir+4; d++)
		{
			newdir=(d+8)%8;
			switch (newdir)
			{
				case 0: xx=x-1; yy=y; gradx=0; grady=-1; break;
				case 1: xx=x-1; yy=y-1; gradx=1; grady=-1; break;
				case 2: xx=x;   yy=y-1; gradx=1; grady=0; break;
				case 3: xx=x+1; yy=y-1; gradx=1; grady=1; break;
				case 4: xx=x+1; yy=y; gradx=0; grady=1; break;
				case 5: xx=x+1; yy=y+1; gradx=-1; grady=1; break;
				case 6: xx=x;   yy=y+1; gradx=-1; grady=0; break;
				case 7: xx=x-1; yy=y+1; gradx=-1; grady=-1; break;
			}
			
			if (xx>=0&&yy>=0&&xx<w&&yy<h)
			{
				double gx=gradx*grad*xedges.Pixel(xx,yy);
				double gy=grady*grad*yedges.Pixel(xx,yy);
				double g=(gx>gy) ? gx : gy;

				if (g>maxg)
				{
					maxg=g;
					maxgdir=newdir;
					maxg_xx=xx;
					maxg_yy=yy;
				}
			}
		}

		if (maxg>mingrad)
		{
			newdir=maxgdir;
			xx=maxg_xx;
			yy=maxg_yy;
		}

		// stop if we're back at the beginning or if we don't move
		if ((xx==startx && yy==starty) || (xx==x && yy==y))
			break;

		dir=newdir;
		x=xx; y=yy;

		if (edges.Pixel(x,y)==100)
			break;

		edges.Pixel(x,y)=100;
	}

	dir=(initdir==2) ? 6 : 0;
	x=startx;
	y=starty;
	for (;;)
	{
		int newdir,xx,yy;
		double maxg=0;
		int maxgdir=dir;
		int maxg_xx=x,maxg_yy=y;
		int d;
		for (d=dir-3; d<dir+4; d++)
		{
			newdir=(d+8)%8;
			switch (newdir)
			{
				case 0: xx=x-1; yy=y; gradx=0; grady=1; break;
				case 1: xx=x-1; yy=y-1; gradx=-1; grady=1; break;
				case 2: xx=x;   yy=y-1; gradx=-1; grady=0; break;
				case 3: xx=x+1; yy=y-1; gradx=-1; grady=-1; break;
				case 4: xx=x+1; yy=y; gradx=0; grady=-1; break;
				case 5: xx=x+1; yy=y+1; gradx=1; grady=-1; break;
				case 6: xx=x;   yy=y+1; gradx=1; grady=0; break;
				case 7: xx=x-1; yy=y+1; gradx=1; grady=1; break;
			}
			if (xx>=0&&yy>=0&&xx<w&&yy<h)
			{
				double gx=gradx*grad*xedges.Pixel(xx,yy);
				double gy=grady*grad*yedges.Pixel(xx,yy);
				double g=(gx>gy) ? gx : gy;

				if (g>maxg)
				{
					maxg=g;
					maxgdir=newdir;
					maxg_xx=xx;
					maxg_yy=yy;
				}
			}
		}

		if (maxg>mingrad)
		{
			newdir=maxgdir;
			xx=maxg_xx;
			yy=maxg_yy;
		}

		// stop if we're back at the beginning or if we don't move
		if ((xx==startx && yy==starty) || (xx==x && yy==y))
			break;

		dir=newdir;
		x=xx; y=yy;

		if (edges.Pixel(x,y)==100)
			break;

		edges.Pixel(x,y)=100;
	}
}

double reflectionfilter(RangeImage &ri,UnorderedPointSet &ups,double maxd,double mind)
{
	// we assume the nose is in ri.ou,ri.ov
	int y0=1;
	int y1=int(ri.ov-30/ri.dv);

	double meanflag=0,maxflag=0,n=0;
	int x,y;
	for (y=0; y<ri.height; y++)
	for (x=0; x<ri.width; x++)
	{
		double flag=ri.Flag(x,y);
		if (flag)
		{
			meanflag+=flag;
			n++;
		}
		if (flag>maxflag)
			maxflag=flag;
	}
	if (n>0)
		meanflag/=n;

	if (debug)
		fprintf(stderr,"meanflag=%g maxflag=%g\n",meanflag,maxflag);

	double *f=new double[ri.width*ri.height];
	for (y=0; y<ri.height; y++)
	for (x=0; x<ri.width; x++)
	{
		double sum=ri.Flag(x,y)*ri.Pixel(x,y);
		double n=ri.Flag(x,y);
		int r;
		for (r=1; r<ri.width; r++)
		{
			int xx,yy;
			for (xx=-r+1; xx<=r-1; xx++)
			{
				if (x+xx<0 || x+xx>ri.width)
					continue;
				if (y-r>0)
				{
					n+=ri.Flag(x+xx,y-r);
					sum+=ri.Flag(x+xx,y-r)*ri.Pixel(x+xx,y-r);
				}
				if (y+r<ri.height)
				{
					n+=ri.Flag(x+xx,y+r);
					sum+=ri.Flag(x+xx,y+r)*ri.Pixel(x+xx,y+r);
				}
			}
			for (yy=-r+1; yy<=r-1; yy++)
			{
				if (y+yy<0 || y+yy>=ri.height)
					continue;
				if (x-r>0)
				{
					n+=ri.Flag(x-r,y+yy);
					sum+=ri.Flag(x-r,y+yy)*ri.Pixel(x-r,y+yy);
				}
				if (x+r<ri.width)
				{
					n+=ri.Flag(x+r,y+yy);
					sum+=ri.Flag(x+r,y+yy)*ri.Pixel(x+r,y+yy);
				}
			}
			if (n>meanflag/2)
				break;
		}

		f[y*ri.width+x]=sum/n;
	}

	memcpy(ri.pixel,f,sizeof(double)*ri.width*ri.height);
	memset(f,0,sizeof(double)*ri.width*ri.height);

	if (debug)
		ri.WriteSFI("discont.sfi");

	// this part is new

	int xn=int(ri.ou+0.5);
	int yn=int(ri.ov+0.5);
	int wn=int(13.5/ri.du);
	y0=int(yn-50/ri.dv);
	if (y0<1)
		y0=1;

	RangeImage candidates;
	candidates.CopyDataFrom(ri);
	RangeImage xedges(ri.width,ri.height);
	RangeImage yedges(ri.width,ri.height);

	for (y=1; y<yn; y++)
	for (x=1; x<ri.width; x++)
	{
		double dx=ri.Pixel(x,y)-ri.Pixel(x-1,y);
		double dy=ri.Pixel(x,y)-ri.Pixel(x,y-1);

		xedges.Pixel(x,y)=dx;
		yedges.Pixel(x,y)=dy;
	}

	if (debug)
	{
		xedges.WriteSFI("xedges.sfi");
		yedges.WriteSFI("yedges.sfi");
	}
/*
{
	RangeImage xedges2,yedges2;
	xedges2.CopyDataFrom(xedges);
	yedges2.CopyDataFrom(yedges);
	for (y=1; y<yn; y++)
	for (x=1; x<ri.width; x++)
	{
		if (xedges2.Pixel(x,y)>1)
			xedges2.Pixel(x,y)=1;
		else if (xedges2.Pixel(x,y)<-1)
                        xedges2.Pixel(x,y)=-1;
		if (yedges2.Pixel(x,y)>1)
			yedges2.Pixel(x,y)=1;
		else if (yedges2.Pixel(x,y)<-1)
                        yedges2.Pixel(x,y)=-1;
	}
	xedges2.WriteSFI("xedges2.sfi");
	yedges2.WriteSFI("yedges2.sfi");
}
*/

	// detect large jumps in brightness and trace connected boundaries
	for (y=1; y<yn; y++)
	for (x=1; x<ri.width; x++)
	{
		if (fabs(xedges.Pixel(x,y))>2 || fabs(yedges.Pixel(x,y))>2)
			trace_boundary(candidates,xedges,yedges,x,y,1);
	}

	if (debug)
		candidates.WriteSFI("reflection_edges.sfi");

	// fit parabola to horizontal intersections of the nose
	double *a=new double[ri.height];
	for (y=yn; y>1; y--)
	{
		int x,xi;
		double z0=ri.Pixel(xn,y);
		int nmax=0;
		double dmin=1000;
		double a_opt,b_opt,c_opt;
		for (xi=int(ri.ou)-wn;xi<=int(ri.ou)+wn; xi++)
		{
			if (xi==xn)
				continue;
			double zi=ri.Pixel(xi,y);
			double a=(zi-z0)/((xi-xn)*(xi-xn));
			if (a<=0)
				continue;
			double b=-2*a*xn;
			double c=a*xn*xn+z0;
			int n=0;
			double dtot=0;
			int yi;
			for (yi=y-1; yi<=y+1; yi++)
			for (x=xn-wn;x<=xn+wn; x++)
			{
				double z=ri.Pixel(x,yi);
				double d=fabs(z-(a*x*x+b*x+c));
				if (d<1)
				{
					n++;
					dtot+=d;
				}
			}
			if (n>nmax)
			{
				nmax=n;
				a_opt=a;
				b_opt=b;
				c_opt=c;
				dmin=dtot;
			}
			else if (n==nmax && dtot<dmin)
			{
				dmin=dtot;
				a_opt=a;
				b_opt=b;
				c_opt=c;
			}
		}
		a[y]=a_opt;

		//fprintf(stderr,"%d %g %g %g %d %g\n",y,a_opt,b_opt,c_opt,nmax,dmin);
	}

	// filter a, i.e. the steepness of the parabola in y-direction
	double *af=new double[ri.height];
	for (y=yn; y>1; y--)
	{
		int dy,n=0;
		double sum=0;
		for (dy=-10; dy<=10; dy++)
		{
			if (y+dy<1 || y+dy>int(ri.ov))
				continue;
			sum+=a[y+dy];
			n++;
		}
		af[y]=(n>0) ? sum/n : 0;
	}
	
	for (y=yn; y>y0; y--)
	{
		double z0=ri.Pixel(xn,y);
		double a_opt=af[y];
		double b_opt=-2*a_opt*xn;
		double c_opt=a_opt*xn*xn+z0;

		//fprintf(stderr,"%d %g %g %g %g\n",y,a_opt,b_opt,c_opt,af[y]);

		int x0;
		for (x0=0; x0<xn; x0++)
		{
			if (candidates.Pixel(x0,y)!=100 || xedges.Pixel(x0,y)>-1)
				continue;

			int x1;
			for (x1=xn-1; x1>x0; x1--)
				if (candidates.Pixel(x1,y)==100 && xedges.Pixel(x1,y)>1)
					break;
			if (x1==x0)
			{
				double z1=ri.Pixel(x0+1,y);
				double D=b_opt*b_opt-4*a_opt*(c_opt-z1);
				x1=(D<0) ? xn-1 : int((-b_opt-sqrt(D))/(2*a_opt));
				if (x1<x0)
					x1=x0;
			}
			for (x=x0; x<x1; x++)
				candidates.Pixel(x,y)=100;
		}
		for (x0=ri.width-1; x0>xn+1; x0--)
		{
			if (candidates.Pixel(x0,y)!=100 || xedges.Pixel(x0,y)<1)
				continue;
			int x1;
			for (x1=xn+1; x1<x0; x1++)
				if (candidates.Pixel(x1,y)==100 && xedges.Pixel(x1,y)<-1)
					break;
			if (x1==x0)
			{
				double z1=ri.Pixel(x0-1,y);
				double D=b_opt*b_opt-4*a_opt*(c_opt-z1);
				x1=(D<=0) ? xn-1 : int((-b_opt+sqrt(D))/(2*a_opt));
				if (x1>x0)
					x1=x0;
			}
			for (x=x0; x>x1; x--)
				candidates.Pixel(x,y)=100;
		}
	}

	delete [] a;
	delete [] af;

	RangeImage can2;
	can2.CopyDataFrom(candidates);

	for (y=yn-1; y>y0+1; y--)
	for (x=1; x<candidates.width-1; x++)
	{
		if (can2.Pixel(x,y)==100)
		{
			int dx,dy;
			for (dy=-1; dy<=1; dy++)
			for (dx=-1; dx<=1; dx++)
				candidates.Pixel(x+dx,y+dy)=100;
		}

	}

	if (debug)
		candidates.WriteSFI("reflections.sfi");

//#define DEVELOPMENT 1
/* development stuff */
#ifdef DEVELOPMENT
{
	RangeImage locmax;
	locmax.CopyDataFrom(ri);
	for (x=0; x<locmax.width; x++)
	{
		for (y=0; y<locmax.height-1; y++)
		{
			if (locmax.Pixel(x,y)<locmax.Pixel(x,y+1))
				locmax.Pixel(x,y)=0;
		}
		for (y=locmax.height-1; y>0; y--)
		{
			if (locmax.Pixel(x,y)<locmax.Pixel(x,y-1))
				locmax.Pixel(x,y)=0;
		}
	}
	locmax.WriteSFI("locmax.sfi");

	FILE *f=fopen("discontline_x.txt","w");

	for (y=0; y<ri.height-1; y+=1)
	{
		double line[ri.width];
		for (x=0; x<ri.width; x++)
			line[x]=ri.Pixel(x,y);
		double sline[ri.width];
		for (x=3; x<ri.width-3; x++)
		{
			sline[x]=0;
			int dx;
			for (dx=-3; dx<=3; dx++)
				sline[x]+=line[x+dx];
			sline[x]/=7;
		}
		for (x=2; x<ri.width-2; x++)
		{
			double dzdxc=(ri.Pixel(x+1,y)-ri.Pixel(x-1,y))/2;
			double dzdxl=(ri.Pixel(x,y)-ri.Pixel(x-1,y))/1;
			double dzdxr=(-ri.Pixel(x,y)+ri.Pixel(x+1,y))/1;
			double ddzdx=dzdxr-dzdxl;
			
			fprintf(f,"%d %d %g %g %g\n",x,y,line[x],dzdxl,ddzdx);
		}
		fprintf(f,"\n\n");
	}
	fclose(f);

	f=fopen("discontline_y.txt","w");

	for (x=0; x<ri.width; x+=1)
	{
		for (y=2; y<ri.height-2; y++)
		{
			double dzdy=ri.Pixel(x,y+1)+ri.Pixel(x,y+2)-ri.Pixel(x,y-1)-ri.Pixel(x,y-2);
			double dzdyl=ri.Pixel(x,y)-ri.Pixel(x,y-2);
			double dzdyr=ri.Pixel(x,y+2)-ri.Pixel(x,y);
			fprintf(f,"%d %d %g %g %g\n",x,y,ri.Pixel(x,y),dzdy,dzdyr-dzdyl);
		}
		fprintf(f,"\n\n");
	}
	fclose(f);

	f=fopen("histz.txt","w");
	y=89;
	int histz[500];
	int z;
	for (z=0; z<500; z++)
		histz[z]=0;
	double r=200;
	for (x=0; x<ri.width; x++)
	{
		double dx=0.5*ri.width-x;
		double z=ri.Pixel(x,y)+sqrt(r*r-dx*dx);
		if (z>=0 && z<500)
			histz[int(z+0.5)]++;
	}

	for (z=0; z<500; z++)
		fprintf(f,"%d %d\n",z,histz[z]);

	
	fclose(f);

	// trace nose boundaries
	double min_dzdx=0,max_dzdx=0;
	int xmin=int(ri.ou),xmax=int(ri.ou);
	y=int(ri.ov);
	for (x=ri.width/8; x<7*ri.width/8; x++)
	{
		double dzdx=ri.Pixel(x+1,y)+ri.Pixel(x+2,y)-ri.Pixel(x-1,y)-ri.Pixel(x-2,y);
		if (dzdx<min_dzdx)
		{
			min_dzdx=dzdx;
			xmin=x;
		}
		else if (dzdx>max_dzdx)
		{
			max_dzdx=dzdx;
			xmax=x;
		}
	}
	ri.Pixel(xmin,y)=0;
	ri.Pixel(xmax,y)=0;

	for (y=int(ri.ov)-1; y>=0; y--)
	{
		min_dzdx=max_dzdx=0;
		for (x=xmin-1; x<=xmin+3; x++)
		{
			if (x+2>=ri.width || x-2<0)
				continue;
			double dzdx=ri.Pixel(x+1,y)+ri.Pixel(x+2,y)-ri.Pixel(x-1,y)-ri.Pixel(x-2,y);
			if (dzdx<min_dzdx)
			{
				min_dzdx=dzdx;
				xmin=x;
			}
		}
		min_dzdx=max_dzdx=0;
		for (x=xmax-3; x<=xmax+1; x++)
		{
			if (x+2>=ri.width || x-2<0)
				continue;
			double dzdx=ri.Pixel(x+1,y)+ri.Pixel(x+2,y)-ri.Pixel(x-1,y)-ri.Pixel(x-2,y);
			if (dzdx>max_dzdx)
			{
				max_dzdx=dzdx;
				xmax=x;
			}
		}
		ri.Pixel(xmin,y)=0;
		ri.Pixel(xmax,y)=0;
	}
	ri.WritePGM("noseboundaries.pgm");
}
#endif /*DEVELOPMENT*/

	UnorderedPointSet upsf;
	int i;
	for (i=0; i<ups.npoints; i++)
	{
		Point3D q=ri.TransformPoint(ups[i]);
		if (candidates.GetDepth(q.x,q.y)<100)
			upsf.AddPoint(ups[i]);
	}

	int upsnpoints=ups.npoints;

	ups.npoints=0;
	for (i=0; i<upsf.npoints; i++)
		ups.AddPoint(upsf[i]);

	if (debug)
		fprintf(stderr,"%d points after filtering\n",upsf.npoints);

	ri.AccumulateDepth(ups);
	if (debug)
		ri.WritePGM("reflfiltered2.pgm",1);
	if (debug)
		ri.WriteFlagPGM("reflfiltered2flags.pgm");

	return (upsf.npoints<=0) ? 0 : double(upsf.npoints)/upsnpoints;

// old poor method
#ifdef OLD_METHOD
	for (y=y0; y<y1; y++)
	{
		double *pf=f+y*ri.width;
		for (x=1; x<ri.width-1; x++)
			pf[x]=ri.Pixel(x-1,y)-ri.Pixel(x+1,y);
		for (x=1; x<ri.width-1; x++)
			if (pf[x]>0 && pf[x+1]>pf[x])
				pf[x]=0;
		for (x=ri.width-2; x>0; x--)
			if (pf[x]>0 && pf[x-1]>pf[x])
				pf[x]=0;
		for (x=1; x<ri.width-1; x++)
			if (pf[x]<0 && pf[x+1]<pf[x])
				pf[x]=0;
		for (x=ri.width-2; x>0; x--)
			if (pf[x]<0 && pf[x-1]<pf[x])
				pf[x]=0;
		for (x=1; x<ri.width-1; x++)
			if (pf[x]>maxd)
			{
				int dx;
				for (dx=1; dx<0.25*ri.width; dx++)
				{
					if (x+dx<ri.width-1 && pf[x+dx]<-mind)
					{
						if (x<ri.width/2 && x+dx>ri.width/2)
							break;
						int xx;
						for (xx=x; xx<=x+dx; xx++)
							pf[xx]=maxd;
						break;
					}
				}
			}
			if (pf[x]<-maxd)
			{
				int dx;
				for (dx=1; dx<0.25*ri.width; dx++)
				{
					if (x-dx>0 && pf[x-dx]>mind)
					{
						if (x>ri.width/2 && x-dx<ri.width/2)
							break;
						int xx;
						for (xx=x; xx>=x-dx; xx--)
							pf[xx]=maxd;
						break;
					}
					
				}
			}
	}

	memcpy(ri.pixel,f,sizeof(double)*ri.width*ri.height);

	for (y=1; y<ri.height-1; y++)
	for (x=1; x<ri.width-1; x++)
	{
		if (f[y*ri.width+x]<maxd || f[(y-1)*ri.width+x]<maxd || f[(y+1)*ri.width+x]<maxd)
			ri.Pixel(x,y)=0;
	}

	if (debug)
		ri.WritePGM("discontborders.pgm",1);

	delete [] f;

	UnorderedPointSet upsf;
	int i;
	for (i=0; i<ups.npoints; i++)
	{
		Point3D q=ri.TransformPoint(ups[i]);
		if (ri.GetDepth(q.x,q.y)<maxd)
			upsf.AddPoint(ups[i]);
	}

	int upsnpoints=ups.npoints;

	ups.npoints=0;
	for (i=0; i<upsf.npoints; i++)
		ups.AddPoint(upsf[i]);

	if (debug)
		fprintf(stderr,"%d points after filtering\n",upsf.npoints);

	ri.AccumulateDepth(ups);
	if (debug)
		ri.WritePGM("reflfiltered2.pgm",1);
	if (debug)
		ri.WriteFlagPGM("reflfiltered2flags.pgm");

	return (upsf.npoints<=0) ? 0 : double(upsf.npoints)/upsnpoints;
#endif /*OLD_METHOD*/
}

double backgroundfilter(RangeImage &ri,UnorderedPointSet &ups,double maxd)
{
	// we assume the nose is in ri.ou,ri.ov
	int y0=1;
	int y1=ri.height-1;

	double meanflag=0,maxflag=0,n=0;
	int x,y;
	for (y=0; y<ri.height; y++)
	for (x=0; x<ri.width; x++)
	{
		double flag=ri.Flag(x,y);
		if (flag)
		{
			meanflag+=flag;
			n++;
		}
		if (flag>maxflag)
			maxflag=flag;
	}
	if (n>0)
		meanflag/=n;

	if (debug)
		fprintf(stderr,"meanflag=%g maxflag=%g\n",meanflag,maxflag);

	double *f=new double[ri.width*ri.height];
	for (y=0; y<ri.height; y++)
	for (x=0; x<ri.width; x++)
	{
		double sum=ri.Flag(x,y)*ri.Pixel(x,y);
		double n=ri.Flag(x,y);
		int r;
		for (r=1; r<ri.width; r++)
		{
			int xx,yy;
			for (xx=-r+1; xx<=r-1; xx++)
			{
				if (x+xx<0 || x+xx>ri.width)
					continue;
				if (y-r>0)
				{
					n+=ri.Flag(x+xx,y-r);
					sum+=ri.Flag(x+xx,y-r)*ri.Pixel(x+xx,y-r);
				}
				if (y+r<ri.height)
				{
					n+=ri.Flag(x+xx,y+r);
					sum+=ri.Flag(x+xx,y+r)*ri.Pixel(x+xx,y+r);
				}
			}
			for (yy=-r+1; yy<=r-1; yy++)
			{
				if (y+yy<0 || y+yy>=ri.height)
					continue;
				if (x-r>0)
				{
					n+=ri.Flag(x-r,y+yy);
					sum+=ri.Flag(x-r,y+yy)*ri.Pixel(x-r,y+yy);
				}
				if (x+r<ri.width)
				{
					n+=ri.Flag(x+r,y+yy);
					sum+=ri.Flag(x+r,y+yy)*ri.Pixel(x+r,y+yy);
				}
			}
			if (n>meanflag/2)
				break;
		}

		f[y*ri.width+x]=sum/n;
	}

	memcpy(ri.pixel,f,sizeof(double)*ri.width*ri.height);
	memset(f,0,sizeof(double)*ri.width*ri.height);

	if (debug)
		ri.WritePGM("background.pgm");

	for (y=y0; y<y1; y++)
	{
		double *pf=f+y*ri.width;
		for (x=1; x<ri.width-1; x++)
			pf[x]=ri.Pixel(x-1,y)-ri.Pixel(x+1,y);
		for (x=1; x<ri.width-1; x++)
			if (pf[x]>0 && pf[x+1]>pf[x])
				pf[x]=0;
		for (x=ri.width-2; x>0; x--)
			if (pf[x]>0 && pf[x-1]>pf[x])
				pf[x]=0;
		for (x=1; x<ri.width-1; x++)
			if (pf[x]<0 && pf[x+1]<pf[x])
				pf[x]=0;
		for (x=ri.width-2; x>0; x--)
			if (pf[x]<0 && pf[x-1]<pf[x])
				pf[x]=0;
		for (x=1; x<ri.width-1; x++)
			if (pf[x]<-maxd)
			{
				int xx;
				for (xx=x; xx<ri.width; xx++)
					pf[xx]=maxd;
			}
			if (pf[x]>maxd)
			{
				int xx;
				for (xx=x; xx>=0; xx--)
					pf[xx]=maxd;
			}
	}

	memcpy(ri.pixel,f,sizeof(double)*ri.width*ri.height);

	if (debug)
		ri.WritePGM("background2.pgm",1);

	delete [] f;

	UnorderedPointSet upsf;
	int i;
	for (i=0; i<ups.npoints; i++)
	{
		Point3D q=ri.TransformPoint(ups[i]);
		if (ri.GetDepth(q.x,q.y)<maxd)
			upsf.AddPoint(ups[i]);
	}

	int upsnpoints=ups.npoints;
	ups.CopyFrom(upsf);

	if (debug)
		fprintf(stderr,"%d points after filtering\n",upsf.npoints);

	ri.AccumulateDepth(ups);
	if (debug)
		ri.WritePGM("bgfiltered2.pgm",1);
	if (debug)
		ri.WriteFlagPGM("bgfiltered2flags.pgm");

	if (verbose)
		fprintf(stderr,"removed %d background points\n",upsnpoints-ups.npoints);

	return (upsf.npoints<=0) ? 0 : double(upsf.npoints)/upsnpoints;
}

int filtershoulders(UnorderedPointSet &ups)
{
	#define MAXHEADWIDTH	300
	#define MAXHEADHEIGHT	350
	#define MINHEADHEIGHT	200
	#define NSTRIPES 	40
	#define STRIPEHEIGHT	50

	int i;
	ups.GetStats();

	// check if there is significantly more than just the head
	double width=ups.maxX-ups.minX;
	double height=ups.maxY-ups.minY;
	if (width>MAXHEADWIDTH || height>MAXHEADHEIGHT)
	{
		if (verbose)
			fprintf(stderr,"image contains more than just the head attempting to find the head %g[mm]x%g[mm]\n",width,height);

		// check if most of the extra width is in the bottom of the image
		double left[NSTRIPES],right[NSTRIPES];
		for (i=0; i<NSTRIPES; i++)
			left[i]=right[i]=ups.minX+width/2;
		for (i=0; i<ups.npoints; i++)
		{
			int stripe=int((ups[i].y-ups.minY)/STRIPEHEIGHT);
			if (stripe<0 || stripe>=NSTRIPES)
				continue;
			if (ups[i].x<left[stripe])
				left[stripe]=ups[i].x;
			else if (ups[i].x>right[stripe])
				right[stripe]=ups[i].x;
			
		}

		for (i=NSTRIPES-1; i>=0; i--)
		{
			if (right[i]-left[i]>MAXHEADWIDTH)
				break;
		}
		double shoulder=ups.minY+i*STRIPEHEIGHT;
		if (ups.maxY-(shoulder)>MINHEADHEIGHT)
		{
			if (verbose)
				fprintf(stderr,"found shoulders at y=%g; dropping points below shoulder\n",shoulder);
			UnorderedPointSet ups2;
			for (i=0; i<ups.npoints; i++)
			{
				if (ups[i].y>shoulder)
					ups2.AddPoint(ups[i]);
			}
			ups.CopyFrom(ups2);
		}
		else if (verbose)
			fprintf(stderr,"could not find shoulder, just continuing\n");
	}

	return 1;
}

double filterpickfrontal(RangeImage &ri,UnorderedPointSet &ups,double maxvar,double maxd)
{
	if (verbose)
		fprintf(stderr,"-->>filterpickfrontal\n");

	int upsnpoints=ups.npoints;

	UnorderedPointSet upsf;
	UnorderedPointSet ups_check;

	int i;
	for (i=0; i<ups.npoints; i++)
	{
		Point3D q=ri.TransformPoint(ups[i]);
		if (ri.GetStddev(q.x,q.y)<maxvar)
			upsf.AddPoint(ups[i]);
		else
		{
			int x=int(q.x+0.5);
			int y=int(q.y+0.5);
			if (x>=0 && y>=0 && x<ri.width && y<ri.height)
				ups_check.AddPoint(ups[i]);
		}
	}
	if (verbose)
		fprintf(stderr,"%d occluding points out of %d\n",ups_check.npoints,upsnpoints);

	int x,y;
	for (y=0; y<ri.height; y++)
	for (x=0; x<ri.width; x++)
	{
		if (ri.Stddev(x,y)>maxvar)
		{
			ri.Pixel(x,y)=+MAXFLOAT;
			ri.Flag(x,y)=0;
		}
	}

	for (i=0; i<ups_check.npoints; i++)
	{
		Point3D q=ri.TransformPoint(ups_check[i]);
		int x=int(q.x+0.5);
		int y=int(q.y+0.5);
		if (q.z<ri.Pixel(x,y))
			ri.Pixel(x,y)=q.z;
	}
	for (i=0; i<ups_check.npoints; i++)
	{
		Point3D q=ri.TransformPoint(ups_check[i]);
		int x=int(q.x+0.5);
		int y=int(q.y+0.5);
		if (fabs(q.z-ri.Pixel(x,y))<maxd)
		{
			upsf.AddPoint(ups_check[i]);
			ri.Flag(x,y)++;
		}
	}
	if (verbose)
		fprintf(stderr,"dropping %d occluded points\n",upsnpoints-upsf.npoints);

	ups.CopyFrom(upsf);

	return (upsf.npoints<=0) ? 0 : double(upsf.npoints)/upsnpoints;
}

double filterhair(RangeImage &ri,UnorderedPointSet &ups)
{
	int x,y;

	int y0=64;
	for (x=0; x<ri.width; x++)
	{
		if (ri.Flag(x,y0))
			printf("%d %g\n",x,ri.Pixel(x,y0));
	}
	printf("\n\n");

	RangeImage ri2(ri.o,ri.u,ri.v,ri.du,ri.dv,ri.ou,ri.ov,ri.width,ri.height);

	for (y=1; y<ri.height-1; y++)
	for (x=1; x<ri.width-1; x++)
	{
		ri2.Pixel(x,y)=-ri.Pixel(x-1,y-1)+ri.Pixel(x+1,y-1)-ri.Pixel(x-1,y)+ri.Pixel(x+1,y)-ri.Pixel(x-1,y+1)+ri.Pixel(x+1,y+1);
		if (x>ri.width/2)
			ri2.Pixel(x,y)*=-1;
	}

	if (debug)
		ri2.WriteSFI("ri2.sfi");

	int xx;
	for (y=1; y<ri.height-1; y++)
	{
		for (x=int(0.2*ri.width); x>0; x--)
		{
			if (ri2.Pixel(x,y)>5)
			{
				// pick 3 points and fit a cylinder (x-x0)^2/a+(z-z0)^2/b=1
				double x0=ri.width/2;
				double x1=x+1-x0,x2=x+1+7.5/ri.du-x0,x3=x+1+15/ri.du-x0;
				double z1=ri.Pixel(int(x1+x0),y);
				double z2=ri.Pixel(int(x2+x0),y);
				double z3=ri.Pixel(int(x3+x0),y);
				double a=(0.6*ri.width)*(0.6*ri.width);
				double k1=1-x3*x3/a;
				double k2=-(1-x1*x1/a);
				double z0=((k1*z1+k2*z3)+sqrt((k1*z1+k2*z3)*(k1*z1+k2*z3)-(k1+k2)*(k1*z1*z1+k2*z3*z3)))/(k1+k2);
				double b=(z1-z0)*(z1-z0)/(1-x1*x1/a);
if (y==y0)
{
	fprintf(stderr,"x1=%g x2=%g x3=%g z1=%g z2=%g z3=%g\n",x1,x2,x3,z1,z2,z3);
	fprintf(stderr,"z0=%g; a=%g; b=%g\n",z0,a,b);
}
				int nfit=0,nnotfit=0;
				for (xx=int(x1+x0); xx<=int(x3+x0); xx++)
				{
					if (fabs(z0-sqrt(b*(1-(xx-x0)*(xx-x0)/a))-ri.Pixel(xx,y))>1)
						nnotfit++;
					else
						nfit++;
				}
				fprintf(stderr,"y=%d nfit=%d nnotfit=%d\n",y,nfit,nnotfit);
				if (nfit==0 || nnotfit/nfit>0.95)
				{
					for (xx=0; xx<=x+1; xx++)
					{
						ri.Pixel(xx,y)=0;
						ri.Flag(xx,y)=0;
					}
				}
				else
				{
					for (xx=0; xx<=x+1; xx++)
					{
						ri.Pixel(xx,y)=z0-sqrt(b*(1-(xx-x0)*(xx-x0)/a));
						ri.Flag(xx,y)=1;
					}
				}
				break;
			}
		}
		for (x=int(0.8*ri.width); x<ri.width-1; x++)
		{
			if (ri2.Pixel(x,y)>5)
			{
				// pick 3 points and fit a cylinder (x-x0)^2/a+(z-z0)^2/b=1
				double x0=ri.width/2;
				double x1=x-1-x0,x2=x-1-7.5/ri.du-x0,x3=x-1-15/ri.du-x0;
				double z1=ri.Pixel(int(x1+x0),y);
				double z2=ri.Pixel(int(x2+x0),y);
				double z3=ri.Pixel(int(x3+x0),y);
				double a=(0.6*ri.width)*(0.6*ri.width);
				double k1=1-x3*x3/a;
				double k2=-(1-x1*x1/a);
				double z0=((k1*z1+k2*z3)+sqrt((k1*z1+k2*z3)*(k1*z1+k2*z3)-(k1+k2)*(k1*z1*z1+k2*z3*z3)))/(k1+k2);
				double b=(z1-z0)*(z1-z0)/(1-x1*x1/a);
if (y==y0)
{
	fprintf(stderr,"x1=%g x2=%g x3=%g z1=%g z2=%g z3=%g\n",x1,x2,x3,z1,z2,z3);
	fprintf(stderr,"z0=%g; a=%g; b=%g\n",z0,a,b);
}
				int nfit=0,nnotfit=0;
				for (xx=int(x3+x0); xx<=int(x1+x0); xx++)
				{
					if (fabs(z0-sqrt(b*(1-(xx-x0)*(xx-x0)/a))-ri.Pixel(xx,y))>1)
						nnotfit++;
					else
						nfit++;
				}
				fprintf(stderr,"y=%d nfit=%d nnotfit=%d\n",y,nfit,nnotfit);
printf("%d %g %g %g\n",y,z0,a,b);
				if (nfit==0 || nnotfit/nfit>0.95)
				{
					for (xx=x-1; xx<ri.width; xx++)
					{
						ri.Pixel(xx,y)=0;
						ri.Flag(xx,y)=0;
					}
				}
				else
				{
					for (xx=x-1; xx<ri.width; xx++)
					{
						ri.Pixel(xx,y)=z0-sqrt(b*(1-(xx-x0)*(xx-x0)/a));
						ri.Flag(xx,y)=1;
					}
				}
				break;
			}
		}
	}
	printf("\n\n");

	for (x=0; x<ri.width; x++)
	{
		if (ri.Flag(x,y0))
			printf("%d %g\n",x,ri.Pixel(x,y0));
	}
	printf("\n\n");

	if (debug)
		ri.WriteSFI("ri.sfi");

	return 1.0;
}

}
