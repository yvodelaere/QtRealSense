/*
 * 11-3-2013: changed max dz in symmetryfunction from 10 to 5
 */

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#include "unorderedpointset.h"
#include "rangeimage.h"
#include "nose.h"
#include "symmetry.h"

#include "opt.h"

namespace utw3dface {

extern int nthreads;
extern int verbose;
extern int debug;

double symmetryfunction(Opt *opt)
{
	static int nrefs=0;
	double ox=opt->pars[0];
	double phi=opt->pars[1];
	RangeImage *ri=(RangeImage *)opt->data;
	UnorderedPointSet *ups=ri->ups;

	Point3D n(cos(phi),sin(phi),0);
	Point3D o(ox,ri->ov,0);
	double k=n*o;
	double sum=0;
	int nhits=0;
	double flag;
	int i;
	int di=ups->npoints/1000;
	if (di<1)
		di=1;
	for (i=0; i<ups->npoints; i+=di)
	{
		Point3D p=(*ups)[i];
		Point3D q=ri->TransformPoint(p);
		Point3D q1=q+2*(k-q*n)*n;
		q1.z=ri->GetDepth(q1.x,q1.y,flag);
		double dz=fabs(q1.z-q.z);
		//if (flag && dz<10 && dz>0.000000001)
		if (flag && dz<5 && dz>0.000000001)
		{
			sum+=dz;
			nhits++;
		}
	}
//fprintf(stderr,"symmetryfunction %d: %g %g %g %d\n",nrefs++,ox,phi,sum/nhits/nhits,nhits);

	return sum/nhits/nhits;
}

double symmetryfunction2(Opt *opt2)
{
	static int nrefs=0;
	double theta=opt2->pars[0];
	Opt *opt=(Opt *)opt2->data;
	RangeImage *ri=(RangeImage *)opt->data;

	Point3D u=ri->u;
	Point3D v=ri->v;

	ri->u=cos(theta)*u+sin(theta)*(u^v);
	ri->calculateT();
	ri->AccumulateDepth();

	opt->SetPar(0,0.5*(opt->min[0]+opt->max[0]),opt->min[0],opt->max[0],opt->acc[0]);
	opt->SetPar(1,0.5*(opt->min[1]+opt->max[1]),opt->min[1],opt->max[1],opt->acc[1]);
	opt->ParabolicFit(10,3);

	ri->u=u;
	ri->v=v;

	return opt->GetMinimum();
}

class RoughSymmetryThreadContainer
{
public:
	int itheta;
	double theta_min,theta_step;
	int nphi;
	double phi_min,phi_step;
	int dx_min,dx_max,dx_step;
	int symwidth;
	double range;
	RangeImage ri;
	ArrayOfSymmetryScores *scores;
	UnorderedPointSet *ups;

	void Init(int itheta,double theta_min,double theta_step,int nphi,double phi_min,double phi_step,int dx_min,int dx_max,int dx_step,int symwidth,double range, RangeImage &ri,ArrayOfSymmetryScores &scores,UnorderedPointSet &ups)
	{
		this->itheta=itheta;
		this->theta_min=theta_min;
		this->theta_step=theta_step;
		this->nphi=nphi;
		this->phi_min=phi_min;
		this->phi_step=phi_step;
		this->dx_min=dx_min;
		this->dx_max=dx_max;
		this->dx_step=dx_step;
		this->symwidth=symwidth;
		this->range=range;
		this->ri.CopyDataFrom(ri);
		this->scores=&scores;
		this->ups=&ups;
	}
};

void *threadtestroughsymmetry(void *p)
{
	RoughSymmetryThreadContainer *rs=(RoughSymmetryThreadContainer*)p;

	double theta=rs->theta_min+rs->itheta*rs->theta_step;
	double min=rs->range,phisum=0;
	int dxsum=0;
	int maxn=0,nsum=0;
	rs->ri.u=Point3D(cos(theta),0,sin(theta));
	rs->ri.calculateT();
	rs->ri.AccumulateDepth(*(rs->ups));

	int width=int(sqrt((double)rs->ri.width*rs->ri.width+rs->ri.height*rs->ri.height));
	int height=width;
	double ox=0.5*width;
	double oy=0.5*height;
	RangeImage rir(rs->ri.o,Point3D(1,0,0),Point3D(0,-1,0),rs->ri.du,rs->ri.dv,ox,oy,width,height);

	int iphi;
	for (iphi=0; iphi<rs->nphi; iphi++)
	{
		double phi=rs->phi_min+iphi*rs->phi_step;
		double cosphi=cos(phi);
		double sinphi=sin(phi);
		double tx=ox-ox*cosphi+oy*sinphi-0.5*(rir.width-rs->ri.width);
		double ty=oy-ox*sinphi-oy*cosphi-0.5*(rir.height-rs->ri.height);

		int x,y;
		for (y=0; y<rir.height; y++)
		for (x=0; x<rir.width; x++)
		{
			double x1=x*cosphi-y*sinphi+tx;
			double y1=x*sinphi+y*cosphi+ty;
			if (x1<0 || x1>=rs->ri.width-1 || y1<0 || y1>=rs->ri.height-1)
			{
				rir.Flag(x,y)=0;
				rir.Pixel(x,y)=0;
				continue;
			}
			int ix1=int(x1);
			int iy1=int(y1);
			double dx1=x1-ix1;
			double dy1=y1-iy1;
			double sum=0;
			double normaliser=0;
			if (rs->ri.Flag(ix1,iy1))
			{
				sum+=rs->ri.Pixel(ix1,iy1)*(1-dx1)*(1-dy1);
				normaliser+=(1-dx1)*(1-dy1);
			}
			if (rs->ri.Flag(ix1+1,iy1))
			{
				sum+=rs->ri.Pixel(ix1+1,iy1)*(dx1)*(1-dy1);
				normaliser+=(dx1)*(1-dy1);
			}
			if (rs->ri.Flag(ix1+1,iy1+1))
			{
				sum+=rs->ri.Pixel(ix1+1,iy1+1)*dx1*dy1;
				normaliser+=(dx1)*(dy1);
			}
			if (rs->ri.Flag(ix1,iy1+1))
			{
				sum+=rs->ri.Pixel(ix1,iy1+1)*(1-dx1)*dy1;
				normaliser+=(1-dx1)*(dy1);
			}
			if (normaliser==0)
			{
				rir.Flag(x,y)=0;
				rir.Pixel(x,y)=0;
			}
			else
			{
				rir.Flag(x,y)=1;
				rir.Pixel(x,y)=sum/normaliser;
			}
		}
		int dx;
		for (dx=rs->dx_min; dx<=rs->dx_max; dx+=rs->dx_step)
		{
			int n=0;
			double sum=0;

			int xmirror=int(ox+dx);
			int xmin=(xmirror<width/2) ? 0 : xmirror-(width-xmirror);
			if (xmirror-xmin>rs->symwidth/rs->ri.du)
				xmin=xmirror-int(rs->symwidth/rs->ri.du);
			int x,y;
			for (y=0; y<height; y++)
			for (x=xmin; x<xmirror; x++)
			{
				int xm=xmirror+(xmirror-x);
				if (rir.Flag(x,y)==0 || rir.Flag(xm,y)==0)
					continue;
				double dz=fabs(rir.Pixel(x,y)-rir.Pixel(xm,y));
				if (dz<rs->range && dz>0.0000001)
				{
					sum+=dz;
					n++;
				}
			}
			if (n<maxn/2 || n==0)
				continue;
			if (n>maxn*1.5)
				min=rs->range;
			if (n>maxn)
				maxn=n;
			sum/=(n>0) ? n*n : 1;
			if (sum<min && n>0)
			{
				min=sum;
				nsum=n;
				phisum=phi;
				dxsum=dx;
			}

		}

		(*(rs->scores))[rs->itheta*rs->nphi+iphi]=SymmetryScore(theta,phisum,dxsum,min);
	}

	if (debug)
		printf("%g %g %d %d %g %d\n",theta,min,maxn,nsum,phisum,dxsum);
}

int rough_symmetry(RangeImage &ri,UnorderedPointSet &ups,ArrayOfSymmetryScores &scores,int nearfrontal)
{
	double range=10;	// max distance to use with comparison
	double symwidth=50;	// 0.5*width to determine symmetry in
	double theta_min=-0.25*M_PI;
	double theta_max=0.25*M_PI;
//theta_min=theta_max=-0.39;
	double theta_step=0.025*M_PI;
	double phi_min=-0.25*M_PI;
	double phi_max=0.25*M_PI;
//phi_max=phi_min=-0.31;
//phi_max=phi_min=-0.39;
	double phi_step=0.025*M_PI;
	int dx_min=-3*ri.width/4;
	int dx_max=3*ri.width/4;
	int dx_step=1;

	if (nearfrontal)
	{
		theta_min=0.5*theta_min;
		theta_max=0.5*theta_max;
		phi_min=0.5*phi_min;
		phi_max=0.5*phi_max;
	}

	int ntheta=int((theta_max-theta_min)/theta_step)+1;
	int nphi=int((phi_max-phi_min)/phi_step)+1;
	int ndx=int((dx_max-dx_min)/dx_step)+1;

	scores.Init(ntheta*nphi);
	int i;
	for (i=0; i<ntheta*nphi; i++)
		scores[i].score=range;

	// first find a rough estimate for the through plane rotation
	int itheta,iphi;
	RoughSymmetryThreadContainer rs[ntheta];
	pthread_t thread[ntheta];
	for (itheta=0; itheta<ntheta; itheta++)
	{
		rs[itheta].Init(itheta,theta_min,theta_step,nphi,phi_min,phi_step,dx_min,dx_max,dx_step,symwidth,range,ri,scores,ups);
		if (nthreads==1)
			threadtestroughsymmetry((void *)&rs[itheta]);
		else
			pthread_create(&thread[itheta],NULL,threadtestroughsymmetry,(void*)&rs[itheta]);
#ifdef NOTHREADS
		double theta=theta_min+itheta*theta_step;
		double min=range,phisum=0;
		int dxsum=0;
		int maxn=0,nsum=0;
		ri.u=Point3D(cos(theta),0,sin(theta));
		ri.calculateT();
		ri.AccumulateDepth(ups);

		int width=int(sqrt((double)ri.width*ri.width+ri.height*ri.height));
		int height=width;
		double ox=0.5*width;
		double oy=0.5*height;
		RangeImage rir(ri.o,Point3D(1,0,0),Point3D(0,-1,0),ri.du,ri.dv,ox,oy,width,height);
	
		for (iphi=0; iphi<nphi; iphi++)
		{
			double phi=phi_min+iphi*phi_step;
/*
        Point3D u=Point3D(cos(theta),0,sin(theta));
        Point3D v=Point3D(0,-1,0);
        Point3D u1=cos(phi)*u+sin(phi)*v;
        Point3D v1=cos(phi)*v-sin(phi)*u;
        rir.u=u1;
        rir.v=v1;
        rir.calculateT();
        rir.AccumulateDepth(ups,1);
*/
			double cosphi=cos(phi);
			double sinphi=sin(phi);
			double tx=ox-ox*cosphi+oy*sinphi-0.5*(rir.width-ri.width);
			double ty=oy-ox*sinphi-oy*cosphi-0.5*(rir.height-ri.height);

			int x,y;
			for (y=0; y<rir.height; y++)
			for (x=0; x<rir.width; x++)
			{
				//int x1=int(x*cosphi-y*sinphi+tx);
				//int y1=int(x*sinphi+y*cosphi+ty);
				double x1=x*cosphi-y*sinphi+tx;
				double y1=x*sinphi+y*cosphi+ty;
				if (x1<0 || x1>=ri.width-1 || y1<0 || y1>=ri.height-1)
				{
					rir.Flag(x,y)=0;
					rir.Pixel(x,y)=0;
					continue;
				}
				int ix1=int(x1);
				int iy1=int(y1);
				double dx1=x1-ix1;
				double dy1=y1-iy1;
				double sum=0;
				double normaliser=0;
				if (ri.Flag(ix1,iy1))
				{
					sum+=ri.Pixel(ix1,iy1)*(1-dx1)*(1-dy1);
					normaliser+=(1-dx1)*(1-dy1);
				}
				if (ri.Flag(ix1+1,iy1))
				{
					sum+=ri.Pixel(ix1+1,iy1)*(dx1)*(1-dy1);
					normaliser+=(dx1)*(1-dy1);
				}
				if (ri.Flag(ix1+1,iy1+1))
				{
					sum+=ri.Pixel(ix1+1,iy1+1)*dx1*dy1;
					normaliser+=(dx1)*(dy1);
				}
				if (ri.Flag(ix1,iy1+1))
				{
					sum+=ri.Pixel(ix1,iy1+1)*(1-dx1)*dy1;
					normaliser+=(1-dx1)*(dy1);
				}
				if (normaliser==0)
				{
					rir.Flag(x,y)=0;
					rir.Pixel(x,y)=0;
				}
				else
				{
					rir.Flag(x,y)=1;
					rir.Pixel(x,y)=sum/normaliser;
				}
				//rir.Pixel(x,y)=ri.Pixel(x1,y1);
				//rir.Flag(x,y)=ri.Flag(x1,y1);
			}
//rir.WritePGM("rir.pgm",0);
//ri.WritePGM("ri.pgm",0);
			int dx;
			for (dx=dx_min; dx<=dx_max; dx+=dx_step)
			{
				int n=0;
				double sum=0;

				int xmirror=int(ox+dx);
				int xmin=(xmirror<width/2) ? 0 : xmirror-(width-xmirror);
				if (xmirror-xmin>symwidth/ri.du)
					xmin=xmirror-int(symwidth/ri.du);
				int x,y;
				for (y=0; y<height; y++)
				for (x=xmin; x<xmirror; x++)
				{
					int xm=xmirror+(xmirror-x);
					if (rir.Flag(x,y)==0 || rir.Flag(xm,y)==0)
						continue;
					double dz=fabs(rir.Pixel(x,y)-rir.Pixel(xm,y));
					if (dz<range && dz>0.0000001)
					{
						sum+=dz;
						n++;
					}
				}
				if (n<maxn/2 || n==0)
					continue;
				if (n>maxn*1.5)
					min=range;
				if (n>maxn)
					maxn=n;
				sum/=(n>0) ? n*n : 1;
				if (sum<min && n>0)
				{
					min=sum;
					nsum=n;
					phisum=phi;
					dxsum=dx;
				}
//printf("%g %g %d %d %g %g\n",theta,phi,dx,n,sum,min);

			}

			scores[itheta*nphi+iphi]=SymmetryScore(theta,phisum,dxsum,min);
		}
//printf("\n\n");

		if (debug)
		{
			printf("%g %g %d %d %g %d\n",theta,min,maxn,nsum,phisum,dxsum);
			//printf("\n\n");
		}
#endif
	}

	if (nthreads>1)
	{
		for (i=0; i<itheta; i++)
			pthread_join(thread[i],NULL);
	}

	if (debug)
		printf("\n\n");

/*
	int i;
	for (i=1; i<ntheta; i++)
		if (scores[i].score>0 && scores[i-1].score>scores[i].score)
			scores[i-1].score=0;

	for (i=ntheta-1; i>0; i--)
		if (scores[i-1].score>0 && scores[i].score>scores[i-1].score)
			scores[i].score=0;
*/

	for (i=0; i<ntheta*nphi; i++)
		scores[i].flag=1;
	for (itheta=0; itheta<ntheta; itheta++)
	for (iphi=1; iphi<nphi; iphi++)
		if (scores[itheta*nphi+iphi].score<scores[itheta*nphi+iphi-1].score)
			scores[itheta*nphi+iphi-1].flag=0;
	for (itheta=0; itheta<ntheta; itheta++)
	for (iphi=nphi-1; iphi>0; iphi--)
		if (scores[itheta*nphi+iphi-1].score<=scores[itheta*nphi+iphi].score)
			scores[itheta*nphi+iphi].flag=0;
	for (iphi=0; iphi<nphi; iphi++)
	for (itheta=1; itheta<ntheta; itheta++)
		if (scores[itheta*nphi+iphi].score<scores[(itheta-1)*nphi+iphi].score)
			scores[(itheta-1)*nphi+iphi].flag=0;
	for (iphi=0; iphi<nphi; iphi++)
	for (itheta=ntheta-1; itheta>0; itheta--)
		if (scores[(itheta-1)*nphi+iphi].score<=scores[itheta*nphi+iphi].score)
			scores[itheta*nphi+iphi].flag=0;

/*
	for (itheta=1; itheta<ntheta-1; itheta++)
	for (iphi=1; iphi<nphi-1; iphi++)
	{
		double c=scores[itheta*nphi+iphi].score;
		scores[itheta*nphi+iphi].flag=1;
		int di,dj;
		for (di=-1; di<=1; di++)
		for (dj=-1; dj<=1; dj++)
		{
			if (scores[(itheta+di)*nphi+iphi+dj].score<c)
			{
				scores[itheta*nphi+iphi].flag=0;
				break;
			}
		}
	}
*/

/*
{
	FILE *ff=fopen("sym.plot","w");
	for (itheta=0; itheta<ntheta-0; itheta++)
	for (iphi=0; iphi<nphi-0; iphi++)
	{
		double c=scores[itheta*nphi+iphi].score;
		int f=scores[itheta*nphi+iphi].flag;
		fprintf(ff,"%g %g %g %d\n",theta_min+itheta*theta_step,phi_min+iphi*phi_step,c,f);
	}
	fclose(ff);
}
*/

	double min=range;	// this is the maximum score we can get
	int imin=0;
	for (i=0; i<ntheta*nphi; i++)
	{
		if (scores[i].score>0 && scores[i].flag)
		{
			if (scores[i].score<min)
			{
				imin=i;
				min=scores[i].score;
			}
		}
	}

	scores.Sort();

	int n=0;
	int nsignificant=0;
	for (i=0; i<ntheta*nphi; i++)
	{
		if (scores[i].score>0 && scores[i].flag)
		{
			if (debug)
				fprintf(stderr,"minimum %d %g %g\n",n,scores[i].score,scores[i].score/min);
			if (scores[i].score/min<2.5)
			{
				scores[nsignificant]=scores[i];
				nsignificant++;
			}
			n++;
		}
	}

	if (verbose)
		fprintf(stderr,"rough_symmetry: found %d local minima of which %d significant\n",n,nsignificant);

//scores[0]=
	return nsignificant;
}

}
