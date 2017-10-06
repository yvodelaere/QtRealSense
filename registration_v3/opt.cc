#include <stdlib.h>
#include <stdio.h>
#ifdef linux
#include <malloc.h>
#endif
#include <limits.h>
#include <math.h>
#define HAVE_FLOAT
#ifdef HAVE_FLOAT
#include <float.h>
#endif

#include "opt.h"

namespace utw3dface {

Opt::Opt(int npars,double (*func)(Opt *opt))
{
	pars=min=max=acc=NULL;
	this->npars=npars;

	pars=new double[npars];
	min=new double[npars];
	max=new double[npars];
	acc=new double[npars];

	if (pars==NULL || min==NULL || max==NULL || acc==NULL)
	{
		fprintf(stderr,"Opt::Opt: failed to allocate memory\n");
		exit(1);
	}

	data=NULL;
	this->func=func;
	minimum=FLT_MAX;
}

Opt::~Opt()
{
	if (pars != NULL)
		delete pars;
	if (min != NULL)
		delete min;
	if (max != NULL)
		delete max;
	if (acc != NULL)
		delete acc;
}

int Opt::SetPar(int index,double value,double minval,double maxval,double acc)
{
	if (index<0 || index>=npars)
		return 0;

	pars[index]=value;
	min[index]=minval;
	max[index]=maxval;
	this->acc[index]=acc;

	return 1;
}

int Opt::SetParValue(int index,double value)
{
	if (index<0 || index>=npars)
		return 0;

	pars[index]=value;

	return 1;
}

int Opt::SetFunc(double (*func)(Opt *opt))
{
	this->func=func;

	return 1;
}

int Opt::SetData(void *data)
{
	this->data=data;

	return 1;
}

double Opt::GetParValue(int index)
{
	if (index<0 || index>=npars)
		return 0.0;

	return pars[index];
}

double Opt::GetFuncValue(double *pars)
{
	double *tpars;
	double value;
	int i;

	if (pars==NULL)
		return FLT_MAX;

	tpars=new double[npars];
	if (tpars==NULL)
	{
		fprintf(stderr,"Opt::getFuncValue - failed to allocate memory\n");
		exit(1);
	}

	for (i=0; i<npars; i++)
	{
		tpars[i]=pars[i];
		this->pars[i]=pars[i];
	}

	value=func(this);

	// restore pars
	for (i=0; i<npars; i++)
		this->pars[i]=tpars[i];
	
	delete tpars;

	return value;
}

double Opt::GetMinimum()
{
	return minimum;
}

int Opt::ParabolicFit1D(int par,int maxiter)
{
	double fa,fb,fc,fd;
	double a,b,c,d;
	int iter;
	double eps=1e-10;
	double acc;

	a=(min[par]<max[par]) ? min[par] : max[par];
	c=(min[par]<max[par]) ? max[par] : min[par];
	b=(pars[par]>a && pars[par]<c) ? pars[par] : 0.5*(a+c);
	acc=this->acc[par];

	pars[par]=a; fa=func(this);
	pars[par]=b; fb=func(this);
	pars[par]=c; fc=func(this);

	// find a b between a and c for which either f(b)<f(c) or f(b)<f(a)
	if (fb>=fa && fb>=fc)
	{
		// check the gradient in b
		d=b+acc;
		pars[par]=d; fd=func(this);

		if (fd-fb>0)
		{
			c=d;
			fc=fd;
		}
		else
		{
			a=d;
			fa=fd;
		}
	}			

	for (iter=0; iter<maxiter; iter++)
	{
		double det;
		double k0,k1;
		double d,fd;
		double a2,b2,c2;

		if (c-a<acc)
			break;

		// fit a parabola y=k0*x^2+k1*x+k2 through (a,fa) (b,fb) and (c,fc)
		// this is done by:
		//	|fa|   | a^2 a 1 ||k0|       |k0|   | a^2 a 1 |-1 |fa|
		//	|fb| = | b^2 b 1 ||k1|	=>   |k1| = | b^2 b 1 |   |fb|
		//	|fc|   | c^2 c 1 ||k2|       |k2|   | c^2 c 1 |   |fc|
	
		a2=a*a;
		b2=b*b;
		c2=c*c;
		det=a2*b+a*c2+b2*c-a2*c-a*b2-b*c2;
	
		if (fabs(det)<eps)
			break;
		
		k0=((b-c)*fa+(c-a)*fb+(a-b)*fc)/det;
		k1=((c2-b2)*fa+(a2-c2)*fb+(b2-a2)*fc)/det;
		//k2=((b2*c-b*c2)*fa+(a*c2-a2*c)*fb+(a2*b-a*b2)*fc)/det;
	
		if (k0==0)
			break;
	
		d=-0.5*k1/k0;
		
		if (d<a)
		{
		// bc is much steeper than ab
			c=b;
			fc=fb;
			d=0.5*(a+c);
		}

		if (d>c)
		{
		// ab is much steeper than bc
			a=b;
			fa=fb;
			d=0.5*(a+c);
		}

		pars[par]=d; fd=func(this);
		
		if (fd<fb)
		{
			if (d<b)
			{
				c=b;
				fc=fb;
			}
			else
			{
				a=b;
				fa=fb;
			}
			b=d;
			fb=fd;
		}
		else
		{
			if (d>b)
			{
				c=d;
				fc=fd;
			}
			else
			{
				a=d;
				fa=fd;
			}
		}
	}

	pars[par]=b;

	minimum=fb;

	return 1;
}
	
int Opt::ParabolicFit(int maxpariter,int niter)
{
	int i,par;

	for (i=0; i<npars; i++)
		if (acc[i]<=0)
			return 0;
	
	minimum=func(this);

	for (i=0; i<niter; i++)
	{
		for (par=0; par<npars; par++)
		{
			if (min[par]==max[par])
				continue;
        		ParabolicFit1D(par,maxpariter);
		}
	}

	return 1;
}
}
