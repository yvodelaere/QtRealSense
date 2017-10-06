#include <stdio.h>
#include "holes.h"

namespace utw3dface {

extern int verbose;
extern int debug;

static double ThresholdMeanflagFraction=0.25;

double fillholes(RangeImage &ri,Nose &nose)
{
	double meanflag=0;
	int n=0;
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
	}
	if (n>0)
		meanflag/=n;

	for (y=0; y<ri.height; y++)
	for (x=0; x<ri.width; x++)
	{
		if (ri.Flag(x,y)>meanflag+1)
			ri.Flag(x,y)=int(meanflag+1);
	}
	
	int threshold=(ThresholdMeanflagFraction*meanflag<1) ? 1 : int(ThresholdMeanflagFraction*meanflag);

	int nholes=0;
	for (y=0; y<ri.height; y++)
	for (x=0; x<ri.width; x++)
	{
		if (ri.Flag(x,y)>=threshold)
			continue;
		nholes++;
	}

	for (y=1; y<ri.height-1; y++)
	for (x=1; x<ri.width-1; x++)
	{
		if (ri.Flag(x,y)>0 || ri.Flag(x+1,y)>0 || ri.Flag(x-1,y)>0 || ri.Flag(x,y-1)>0 || ri.Flag(x,y+1)>0)
			continue;
		if (nose.PartOf(ri.ov-y))
		{
			double z=nose.NY(ri.ov-y);
			if (z-1*fabs(x-ri.ou)<0)
				continue;
			ri.Pixel(x,y)=1*fabs(x-ri.ou)-z;
			ri.Flag(x,y)=(threshold<2) ? 1 : int(0.5*threshold);
		}
	}

	double badness=0;
	for (y=0; y<ri.height; y++)
	for (x=0; x<ri.width; x++)
	{
		if (ri.Flag(x,y)>=threshold)
			continue;
		double sum=0;
		double n=0;
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
					n+=ri.Flag(x+xx,y-r)/double(r);
					sum+=ri.Flag(x+xx,y-r)*ri.Pixel(x+xx,y-r)/r;
				}
				if (y+r<ri.height)
				{
					n+=ri.Flag(x+xx,y+r)/double(r);
					sum+=ri.Flag(x+xx,y+r)*ri.Pixel(x+xx,y+r)/r;
				}
			}
			for (yy=-r+1; yy<=r-1; yy++)
			{
				if (y+yy<0 || y+yy>=ri.height)
					continue;
				if (x-r>0)
				{
					n+=ri.Flag(x-r,y+yy)/double(r);
					sum+=ri.Flag(x-r,y+yy)*ri.Pixel(x-r,y+yy)/r;
				}
				if (x+r<ri.width)
				{
					n+=ri.Flag(x+r,y+yy)/double(r);
					sum+=ri.Flag(x+r,y+yy)*ri.Pixel(x+r,y+yy)/r;
				}
			}
			if (n>=threshold && n>4)
				break;
		}
		badness+=(r*r/9>1) ? 1 : (r*r/16);

		ri.Pixel(x,y)=sum/n;
	}

	
	return double(badness)/(ri.width*ri.height);
}

double fill_big_holes_using_symmetry(RangeImage &ri)
{
#define MIND 1
	double meanflag=0;
	int n=0;
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
	}
	if (n>0)
		meanflag/=n;

	for (y=0; y<ri.height; y++)
	for (x=0; x<ri.width; x++)
	{
		if (ri.Flag(x,y)>meanflag+1)
			ri.Flag(x,y)=int(meanflag+1);
	}
	
	int threshold=(ThresholdMeanflagFraction*meanflag<1) ? 1 : int(ThresholdMeanflagFraction*meanflag);

	int nholes=0;
	for (y=MIND; y<ri.height-MIND; y++)
	for (x=MIND; x<ri.width-MIND; x++)
	{
		if (ri.Flag(x,y)>=threshold)
			continue;
		int neighbours=0;
		int dx,dy;
		for (dy=-MIND; dy<=MIND; dy++)
		for (dx=-MIND; dx<=MIND; dx++)
		{
			if (ri.Flag(x+dx,y+dy)>=threshold)
				neighbours++;
		}
		if (neighbours==0)
		{
			// try to get the value from the other side of the symmetry axis
			if (ri.Flag(ri.width-1-x,y)>=threshold)
			{
				ri.Flag(x,y)=-ri.Flag(ri.width-1-x,y);
				ri.Pixel(x,y)=ri.Pixel(ri.width-1-x,y);
			}
			nholes++;
		}
	}

	for (y=1; y<ri.height-1; y++)
	for (x=1; x<ri.width-1; x++)
	{
		if (ri.Flag(x,y)<0)
			ri.Flag(x,y)=-ri.Flag(x,y);
	}

	return 2*double(nholes)/(ri.width*ri.height);
}
}
