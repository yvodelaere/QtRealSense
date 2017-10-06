#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef linux
#include <values.h>
#include <malloc.h>
#endif
#include <math.h>
#include "compat.h"

#include "spline.h"

namespace utw3dface {

#define bool int

static double interpolate(int n,double *t,double *y,double *z,double x,double tension,bool periodic);
static void fit(int n,double *t,double *y,double *z,double k,double tension,bool periodic);

Spline::Spline(int dim)
{
	if (dim<1)
	{
		fprintf(stderr,"Cannot create spline with dim<1\n");
		exit(1);
	}

	points= new CoordList(dim);
	abscissa= new CoordList(1);

	closed=0;
	tension=0.0;
	boundary_condition=1.0;
	y=NULL;
	z=NULL;
	t=NULL;
}

Spline::~Spline()
{
	delete points;
	delete abscissa;

	if (y != NULL)
		free(y);
	if (z != NULL)
		free(z);
	if (z != NULL)
		free(t);
}

int Spline::AddPoint(float *coords,float t)
{
	points->Add(coords);
	abscissa->Add(&t);

	return 1;
}

int Spline::SetAbscissaGeoDist(int npoints)
{
	int i,n;
	CoordList t(1);
	double dist;

	if (npoints<0)
		return 0;
	npoints +=2;

	n=points->npoints;

	t.Resize(n);

	if (abscissa->npoints != n)
		abscissa->Resize(n);

	// first just fit a spline 
	Fit();
	
	dist=0;
	t[0][0]=0;
	for (i=1; i<n; i++)
	{
		CoordList *cl=InterpolateRange((float)this->t[i-1],(float)this->t[i],npoints);
		int j;
		
		for (j=0; j<cl->npoints-1; j++)
		{
			double d=0; 
			int dim;

			for (dim=0; dim<cl->dim; dim++)
			{
				float ddim=cl->Get(j+1)[dim]-cl->Get(j)[dim];
				d += ddim*ddim;
			}
			dist += sqrt(d);
		}
		t.coords[i]=(float)dist;

		delete cl;
	}

	for (i=0; i<n; i++)
		abscissa->coords[i]=t.coords[i];

	return 1;
}

int Spline::SetAbscissa(float *t)
{
	int i,n;

	n=points->npoints;

	if (abscissa->npoints != n)
		abscissa->Resize(n);

	if (t==NULL) // auto, range [0,1]
	{
		for (i=0; i<n; i++)
			abscissa->coords[i]=(float)i/(float)(n-1);
	}
	else
	{
		for (i=0; i<n; i++)
		{
			if (i>0 && t[i]<=t[i-1])
				return 0;
			abscissa->coords[i]=t[i];
		}
	}

	return 1;
}

int Spline::Fit()
{
	int dim,n,i,d;

	n=points->npoints;
	dim=points->dim;

	if (n<2+(closed!=0))
		return 0;

	if (closed)
	{
		float *first=points->coords;
		float *last=points->coords+(n-1)*dim;
		float t;

		this->boundary_condition=0.0;

		// ensure periodicity
		for (d=0; d<dim; d++)
		{
			if (first[d] != last[d])
				last[d]=first[d];
		}

		// add pseudo-point
		points->Add(first+dim);
		if (abscissa->npoints==n)
		{
			t=abscissa->coords[n-1] + abscissa->coords[1]- abscissa->coords[0];
			abscissa->Add(&t);
		}
		n++;
	}

	if (abscissa->npoints != n)
		SetAbscissa(NULL);
	
	if (this->y!=NULL)
		free(this->y);
	this->y=(double*)malloc(sizeof(double)*n*dim);
	if (this->z!=NULL)
		free(this->z);
	this->z=(double*)malloc(sizeof(double)*n*dim);
	if (this->t!=NULL)
		free(this->t);
	this->t=(double*)malloc(sizeof(double)*n);

	for (d=0; d<dim; d++)
	for (i=0; i<n; i++)
		this->y[d*n+i]=this->points->coords[i*dim+d];
	for (i=0; i<n; i++)
		this->t[i]=this->abscissa->coords[i];

	// do the actual fit 
	for (d=0; d<dim; d++)
		fit(n-1,this->t,this->y+d*n,this->z+d*n,this->boundary_condition,this->tension,this->closed);

	if (this->closed)
	{
		// remove the pseudo-point
		points->Remove(n-1);
		abscissa->Remove(n-1);
	}

	return 1;
}

Spline::Spline(CoordList &points,int closed)
{
	this->points=new CoordList(points);
	this->abscissa= new CoordList(1);

	tension=0.0;
	boundary_condition=1.0;
	y=NULL;
	z=NULL;
	t=NULL;
	this->closed=closed;

	if (closed)
	{
		int i;
		float *first=points.coords;
		float *last=points.coords+(points.npoints-1)*points.dim;
		float d=0;

		for (i=0; i<points.dim; i++)
			d += (first[i]-last[i])*(first[i]-last[i]);
		if (d>1e-10)	// take into account rounding errors
			this->points->Add(first);
	}

	SetAbscissa(NULL);
	Fit();
}

int Spline::Interpolate(float t,float *coords)
{
	int dim,n,d,nc;

	if (this->y==NULL || this->z==NULL || this->t==NULL)
		return 0;

	n=this->points->npoints;
	dim=this->points->dim;

	// correct for pseudo point in closed curves
	nc = (this->closed) ? n+1 : n;

	if (t<this->t[0] || t>this->t[n-1])
		return 0;

	for (d=0; d<this->points->dim; d++)
		coords[d]=(float)interpolate(n-1,this->t,this->y+d*nc,this->z+d*nc,t,this->tension,this->closed);

	return 1;
}

CoordList *Spline::InterpolateRange(float t0,float t1,int n)
{
	CoordList *cl;
	int dim,i;
	float dt;

	if (t1<=t0 || n<1)
		return NULL;

	if (this->y==NULL || this->z==NULL || this->t==NULL)
		return NULL;

	if (t0<this->t[0] || t1>this->t[this->points->npoints-1])
		return NULL;
	dim=this->points->dim;

	cl=new CoordList(dim);
	cl->Resize(n);

	dt = (t1-t0)/(n-1);

	for (i=0; i<n; i++)
		Interpolate(t0+i*dt,cl->coords+i*dim);

	return cl;
}

int Spline::Normal2D(float t,float dt,float *n)
{
	// returns normal for 2-D spline according to central differential scheme
	float left[2],right[2],norm;
	float tl,tr;

	if (this->points->dim != 2)
		return 0;

	tr=t+dt;
	tl=t-dt;
	if (this->closed)
	{
		float tmin=(float)this->t[0];
		float tmax=(float)this->t[this->abscissa->npoints-1];
		float trange=tmax-tmin;
		while (tl<tmin) tl += trange;
		while (tr<tmin) tr += trange;
		while (tl>tmax) tl -= trange;
		while (tr>tmax) tr -= trange;
	}

	if (!(Interpolate(tr,right)  && Interpolate(tl,left)))
		return 0;

	n[1]=(left[0]-right[0])/(2*dt);
	n[0]=(right[1]-left[1])/(2*dt);
	norm=(float)sqrt(n[0]*n[0]+n[1]*n[1]);
	n[1]/=norm;
	n[0]/=norm;

	return 1;
}

float Spline::Curvature2D(float t,float dt)
{
	// returns curvature for 2-D spline according to central difference scheme
	float left[2],right[2],middle[2];
	double xt,yt,xtt,ytt,norm;
	float tl,tr;

	if (this->points->dim != 2)
		return 0;

	tr=t+dt;
	tl=t-dt;
	if (this->closed)
	{
		float tmin=(float)this->t[0];
		float tmax=(float)this->t[this->abscissa->npoints-1];
		float trange=tmax-tmin;
		while (tl<tmin) tl += trange;
		while (tr<tmin) tr += trange;
		while (tl>tmax) tl -= trange;
		while (tr>tmax) tr -= trange;
		while (t>tmax) t -= trange;
		while (t<tmin) t += trange;
	}

	if (!(Interpolate(tr,right) && Interpolate(t,middle) && Interpolate(tl,left)))
		return 0;
	
	xt=(right[0]-left[0])/(2*dt);
	yt=(right[1]-left[1])/(2*dt);
	xtt=(right[0]-2*middle[0]+left[0])/(dt*dt);
	ytt=(right[1]-2*middle[1]+left[1])/(dt*dt);
	norm=sqrt(xt*xt+yt*yt);

	return (float)((xt*ytt - yt*xtt)/(norm*norm*norm));
}

//////////////////////////////////////////////////////////////////////////////

// this part is from the Gnu spline program, part of the libplot library. 

// some macro's 

#define xmalloc 	malloc
#define xrealloc 	realloc
#define bool		int
#define true		1
#define false		0
#define EXIT_FAILURE	1
#define _HAVE_PROTOS

static double quotient_sin_func (double x,double y);
static double quotient_sinh_func (double x,double y);
static double sin_func (double x);
static double sinh_func (double x);
static double tan_func (double x);
static double tanh_func (double x);

/* This program, spline, interpolates scalar or vector-valued input data
   using splines with tension, including piecewise cubic (zero-tension)
   splines.  When acting as a real-time filter, it uses cubic Bessel
   interpolation instead.  Written by Robert S. Maier
   <rsm@math.arizona.edu>, based on earlier work by Rich Murphey.
   Copyright (C) 1989-1999 Free Software Foundation, Inc.

   References:

   D. Kincaid and [E.] W. Cheney, Numerical Analysis, Brooks/Cole,
   2nd. ed., 1996, Section 6.4.

   C. de Boor, A Practical Guide to Splines, Springer-Verlag, 1978, 
   Chapter 4.

   A. K. Cline, "Scalar and Planar-Valued Curve Fitting Using Splines under
   Tension", Communications of the ACM 17 (1974), 218-223.

   The tension in a spline is set with the -T (i.e., --tension) option.  By
   definition, a one-dimensional spline with tension satisfies the
   differential equation y''''=sgn(tension)*(tension**2)y''.  The default
   value for the tension is zero.  If tension=0 then a spline with tension
   reduces to a conventional piecewise cubic spline.  In the limits
   tension->+infinity and tension->-infinity, a spline with tension reduces
   to a piecewise linear (`broken line') interpolation.

   To oversimplify a bit, 1.0/tension is the maximum abscissa range over
   which the spline likes to curve, at least when tension>0.  So increasing
   the tension far above zero tends to make the spline contain short curved
   sections, separated by sections that are almost straight.  The curved
   sections will be centered on the user-specified data points.  The
   behavior of the spline when tension<0 is altogether different: it will
   tend to oscillate, though as tension->-infinity the oscillations are
   damped out.

   Tension is a `dimensionful' quantity.  If tension=0 (the cubic spline
   case), then the computation of the spline is scale-invariant.  But if
   the tension is nonzero, then when the abscissa values are multiplied by
   some common positive factor, the tension should be divided by the same
   factor to obtain a scaled version of the original spline.

   The algorithms of Kincaid and Cheney have been extended to include
   support for periodicity.  To obtain a periodic spline, with or without
   tension, the user uses the -p (i.e., --periodic) option and supplies
   input data satisfying y[n]=y[0].  Also, in the non-periodic case the
   algorithms have been extended to include support for a parameter k,
   which appears in the two boundary conditions y''[0]=ky''[1] and
   y''[n]=ky''[n-1].  The default value of k is 1.0.  The parameter k,
   which is specified with the -k (i.e. --boundary-condition) option, is
   ignored for periodic splines (using the -k option with the -p option
   will elicit a warning).

   If the -f option is specified, then an altogether different (real-time)
   algorithm for generating interpolating points will be used, so that this
   program can be used as a real-time filter.  If -f is specified then the
   -t option, otherwise optional, must also be used.  (I.e., the minimum
   and maximum abscissa values for the interpolating points must be
   specified, and optionally the spacing between them as well.  If the
   spacing is not specified on the command line, then the interval
   [tmin,tmax] will be subdivided into a default number of intervals [100],
   unless the default number of intervals is overridden with the -n option.

   The real-time algorithm that is used when the -f option is specified is
   cubic Bessel interpolation.  (The -T, -p, and -k options are ignored
   when -f is specified; using them will elicit a warning.)  Interpolation
   in this case is piecewise cubic, and the slopes at either end of each
   sub-interval are found by fitting a parabola through each successive
   triple of points.  That is, the slope at t=t_n is found by fitting a
   parabola through the points at t_(n-1), t_n, and t_(n+1).  This
   interpolation scheme yields a spline that is only once, rather than
   twice, continuously differentiable.  However, it has the feature that
   all computations are local rather than global, so it is suitable for
   real-time work.

   Since the above was written, the -d option has been added, to permit the
   splining of multidimensional data.  All components of a d-dimensional
   data set (a d-dimensional vector y is specified at each t) are splined
   in the same way, as if they were one-dimensional functions of t.  All
   options that apply to 1-dimensional datasets, such as -T, -p, -k, -f,
   etc., apply to d-dimensional ones also. */

/* Minimum value for magnitude of x, for such functions as x-sinh(x),
   x-tanh(x), x-sin(x), and x-tan(x) to have acceptable accuracy.  If the
   magnitude of x is smaller than this value, these functions of x will be
   computed via power series to accuracy O(x**6). */
#define TRIG_ARG_MIN 0.001

/* Maximum value for magnitude of x, beyond which we approximate
   x/sinh(x) and x/tanh(x) by |x|exp(-|x|). */
#define TRIG_ARG_MAX 50.0

const char *progname = "spline"; /* name of this program */

/* fit() computes the array z[] of second derivatives at the knots, i.e.,
   internal data points.  The abscissa array t[] and the ordinate array y[]
   are specified.  On entry, have n+1 >= 2 points in the t, y, z arrays,
   numbered 0..n.  The knots are numbered 1..n-1 as in Kincaid and Cheney.
   In the periodic case, the final knot, i.e., (t[n-1],y[n-1]), has the
   property that y[n-1]=y[0]; moreover, y[n]=y[1].  The number of points
   supplied by the user was n+1 in the non-periodic case, and n in the
   periodic case.  When this function is called, n>=1 in the non-periodic
   case, and n>=2 in the periodic case. */

/* Algorithm: the n-1 by n-1 tridiagonal matrix equation for the vector of
   2nd derivatives at the knots is reduced to upper diagonal form.  At that
   point the diagonal entries (pivots) of the upper diagonal matrix are in
   the vector u[], and the vector on the right-hand side is v[].  That is,
   the equation is of the form Ay'' = v, where a_(ii) = u[i], and a_(i,i+1)
   = alpha[i].  Here i=1..n-1 indexes the set of knots.  The matrix
   equation is solved by back-substitution for y''[], i.e., for z[]. */

void
#ifdef _HAVE_PROTOS
fit (int n, double *t, double *y, double *z, double k, double tension,
     bool periodic)
#else
fit (n, t, y, z, k, tension, periodic)
     int n;
     double *t, *y, *z;
     double k;			/* y''_1 = k y''_0, etc. */
     double tension;
     bool periodic;
#endif
{
  double *h, *b, *u, *v, *alpha, *beta;
  double *uu = NULL, *vv = NULL, *s = NULL;
  int i;

  if (n == 1)			/* exactly 2 points, use straight line */
    {
      z[0] = z[1] = 0.0;
      return;
    }

  h = (double *)xmalloc (sizeof(double) * n);
  b = (double *)xmalloc (sizeof(double) * n);
  u = (double *)xmalloc (sizeof(double) * n);
  v = (double *)xmalloc (sizeof(double) * n);
  alpha = (double *)xmalloc (sizeof(double) * n);
  beta = (double *)xmalloc (sizeof(double) * n);

  if (periodic)
    {
      s = (double *)xmalloc (sizeof(double) * n); 
      uu = (double *)xmalloc (sizeof(double) * n); 
      vv = (double *)xmalloc (sizeof(double) * n); 
    }

  for (i = 0; i <= n - 1 ; ++i)
    {
      h[i] = t[i + 1] - t[i];
      b[i] = 6.0 * (y[i + 1] - y[i]) / h[i]; /* for computing RHS */
    }

  if (tension < 0.0)		/* must rule out sin(tension * h[i]) = 0 */
    {
      for (i = 0; i <= n - 1 ; ++i)
	if (sin (tension * h[i]) == 0.0)
	  {
	    fprintf (stderr, "%s: error: specified negative tension value is singular\n", progname);
	    exit (EXIT_FAILURE);
	  }
    }
  if (tension == 0.0)
    {
      for (i = 0; i <= n - 1 ; ++i)
	{
	  alpha[i] = h[i];	/* off-diagonal = alpha[i] to right */
	  beta[i] = 2.0 * h[i];	/* diagonal = beta[i-1] + beta[i] */
	}
    }
  else
    if (tension > 0.0)
      /* `positive' (really real) tension, use hyperbolic trig funcs */
      {
	for (i = 0; i <= n - 1 ; ++i)
	  {
	    double x = tension * h[i];
	    double xabs = (x < 0.0 ? -x : x);

	    if (xabs < TRIG_ARG_MIN)
	      /* hand-compute (6/x^2)(1-x/sinh(x)) and (3/x^2)(x/tanh(x)-1)
                 to improve accuracy; here `x' is tension * h[i] */
	      {
		alpha[i] = h[i] * sinh_func(x);
		beta[i] = 2.0 * h[i] * tanh_func(x);
	      }
	    else if (xabs > TRIG_ARG_MAX)
	      /* in (6/x^2)(1-x/sinh(x)) and (3/x^2)(x/tanh(x)-1),
		 approximate x/sinh(x) and x/tanh(x) by 2|x|exp(-|x|)
		 and |x|, respectively */
	      {
		int sign = (x < 0.0 ? -1 : 1);

		alpha[i] = ((6.0 / (tension * tension))
			   * ((1.0 / h[i]) - tension * 2 * sign * exp(-xabs)));
		beta[i] = ((6.0 / (tension * tension))
			   * (tension - (1.0 / h[i])));
	      }
	    else
	      {
		alpha[i] = ((6.0 / (tension * tension))
			    * ((1.0 / h[i]) - tension / sinh(x)));
		beta[i] = ((6.0 / (tension * tension))
			   * (tension / tanh(x) - (1.0 / h[i])));
	      }
	  }
      }
    else				/* tension < 0 */
      /* `negative' (really imaginary) tension,  use circular trig funcs */
      {
	for (i = 0; i <= n - 1 ; ++i)
	  {
	    double x = tension * h[i];
	    double xabs = (x < 0.0 ? -x : x);

	    if (xabs < TRIG_ARG_MIN)
	      /* hand-compute (6/x^2)(1-x/sin(x)) and (3/x^2)(x/tan(x)-1)
                 to improve accuracy; here `x' is tension * h[i] */
	      {
		alpha[i] = h[i] * sin_func(x);
		beta[i] = 2.0 * h[i] * tan_func(x);
	      }
	    else
	      {
		alpha[i] = ((6.0 / (tension * tension))
		           * ((1.0 / h[i]) - tension / sin(x)));
		beta[i] = ((6.0 / (tension * tension))
			   * (tension / tan(x) - (1.0 / h[i])));
	      }
	  }
      }
  
  if (!periodic && n == 2)
      u[1] = beta[0] + beta[1] + 2 * k * alpha[0];
  else
    u[1] = beta[0] + beta[1] + k * alpha[0];

  v[1] = b[1] - b[0];
  
  if (u[1] == 0.0)
    {
      fprintf (stderr, 
	       "%s: error: as posed, problem of computing spline is singular\n",
	       progname);
      exit (EXIT_FAILURE);
    }

  if (periodic)
    {
      s[1] = alpha[0];
      uu[1] = 0.0;
      vv[1] = 0.0;
    }

  for (i = 2; i <= n - 1 ; ++i)
    {
      u[i] = (beta[i] + beta[i - 1]
	      - alpha[i - 1] * alpha[i - 1] / u[i - 1]
	      + (i == n - 1 ? k * alpha[n - 1] : 0.0));

      if (u[i] == 0.0)
	{
	  fprintf (stderr, 
		   "%s: error: as posed, problem of computing spline is singular\n",
		   progname);
	  exit (EXIT_FAILURE);
	}


      v[i] = b[i] - b[i - 1] - alpha[i - 1] * v[i - 1] / u[i - 1];

      if (periodic)
	{
	  s[i] = - s[i-1] * alpha[i-1] / u[i-1];
	  uu[i] = uu[i-1] - s[i-1] * s[i-1] / u[i-1];
	  vv[i] = vv[i-1] - v[i-1] * s[i-1] / u[i-1];
	}
    }
      
  if (!periodic)
    {
      /* fill in 2nd derivative array */
      z[n] = 0.0;
      for (i = n - 1; i >= 1; --i)
	z[i] = (v[i] - alpha[i] * z[i + 1]) / u[i];
      z[0] = 0.0;
      
      /* modify to include boundary condition */
      z[0] = k * z[1];
      z[n] = k * z[n - 1];
    }
  else		/* periodic */
    {
      z[n-1] = (v[n-1] + vv[n-1]) / (u[n-1] + uu[n-1] + 2 * s[n-1]);
      for (i = n - 2; i >= 1; --i)
	z[i] = ((v[i] - alpha[i] * z[i + 1]) - s[i] * z[n-1]) / u[i];

      z[0] = z[n-1];
      z[n] = z[1];
    }

  if (periodic)
    {
      free (vv);
      free (uu);
      free (s);
    }
  free (beta);
  free (alpha);
  free (v);
  free (u);
  free (b);
  free (h);
}

/* interpolate() computes an approximate ordinate value for a given
   abscissa value, given an array of data points (stored in t[] and y[],
   containing abscissa and ordinate values respectively), and z[], the
   array of 2nd derivatives at the knots (i.e. internal data points).
   
   On entry, have n+1 >= 2 points in the t, y, z arrays, numbered 0..n.
   The number of knots (i.e. internal data points) is n-1; they are
   numbered 1..n-1 as in Kincaid and Cheney.  In the periodic case, the
   final knot, i.e., (t[n-1],y[n-1]), has the property that y[n-1]=y[0];
   also, y[n]=y[1].  The number of data points supplied by the user was n+1
   in the non-periodic case, and n in the periodic case.  When this
   function is called, n>=1 in the non-periodic case, and n>=2 in the
   periodic case. */

double
#ifdef _HAVE_PROTOS
interpolate (int n, double *t, double *y, double *z, double x, 
	     double tension, bool periodic)
#else
interpolate (n, t, y, z, x, tension, periodic)
     int n;
     double *t, *y, *z, x;
     double tension;
     bool periodic;
#endif
{
  double diff, updiff, reldiff, relupdiff, h;
  double value;
  int is_ascending = (t[n-1] < t[n]);
  int i = 0, k;

  /* in periodic case, map x to t[0] <= x < t[n] */
  if (periodic && (x - t[0]) * (x - t[n]) > 0.0)
    x -= ((int)(floor( (x - t[0]) / (t[n] - t[0]) )) * (t[n] - t[0]));

  /* do binary search to find interval */
  for (k = n - i; k > 1;)
    {
      if (is_ascending ? x >= t[i + (k>>1)] : x <= t[i + (k>>1)])
	{
	  i = i + (k>>1);
	  k = k - (k>>1);
	}
      else
	k = k>>1;
    }

  /* at this point, x is between t[i] and t[i+1] */
  h = t[i + 1] - t[i];
  diff = x - t[i];
  updiff = t[i+1] - x;
  reldiff = diff / h;
  relupdiff = updiff / h;

  if (tension == 0.0)
  /* evaluate cubic polynomial in nested form */
    value =  y[i] 
      + diff
	* ((y[i + 1] - y[i]) / h - h * (z[i + 1] + z[i] * 2.0) / 6.0
	   + diff * (0.5 * z[i] + diff * (z[i + 1] - z[i]) / (6.0 * h)));
  
  else if (tension > 0.0)
    /* `positive' (really real) tension, use sinh's */
    {
      if (fabs(tension * h) < TRIG_ARG_MIN)
	/* hand-compute (6/y^2)(sinh(xy)/sinh(y) - x) to improve accuracy;
	   here `x' means reldiff or relupdiff and `y' means tension*h */
	value = (y[i] * relupdiff + y[i+1] * reldiff
		 + ((z[i] * h * h / 6.0) 
		    * quotient_sinh_func (relupdiff, tension * h))
		 + ((z[i+1] * h * h / 6.0) 
		    * quotient_sinh_func (reldiff, tension * h)));
      else if (fabs(tension * h) > TRIG_ARG_MAX)
	/* approximate 1/sinh(y) by 2 sgn(y) exp(-|y|) */
	{
	  int sign = (h < 0.0 ? -1 : 1);

	  value = (((z[i] * (exp (tension * updiff - sign * tension * h) 
			     + exp (-tension * updiff - sign * tension * h))
		     + z[i + 1] * (exp (tension * diff - sign * tension * h) 
				   + exp (-tension * diff - sign * tension*h)))
		    * (sign / (tension * tension)))
		   + (y[i] - z[i] / (tension * tension)) * (updiff / h)
		   + (y[i + 1] - z[i + 1] / (tension * tension)) * (diff / h));
	}
      else
	value = (((z[i] * sinh (tension * updiff) 
		   + z[i + 1] * sinh (tension * diff))
		  / (tension * tension * sinh (tension * h)))
		 + (y[i] - z[i] / (tension * tension)) * (updiff / h)
		 + (y[i + 1] - z[i + 1] / (tension * tension)) * (diff / h));
    }
  else
    /* `negative' (really imaginary) tension, use sin's */
    {
      if (fabs(tension * h) < TRIG_ARG_MIN)
	/* hand-compute (6/y^2)(sin(xy)/sin(y) - x) to improve accuracy;
	   here `x' means reldiff or relupdiff and `y' means tension*h */
	value = (y[i] * relupdiff + y[i+1] * reldiff
		 + ((z[i] * h * h / 6.0) 
		    * quotient_sin_func (relupdiff, tension * h))
		 + ((z[i+1] * h * h / 6.0) 
		    * quotient_sin_func (reldiff, tension * h)));
      else
	value = (((z[i] * sin (tension * updiff) 
		   + z[i + 1] * sin (tension * diff))
		  / (tension * tension * sin (tension * h)))
		 + (y[i] - z[i] / (tension * tension)) * (updiff / h)
		 + (y[i + 1] - z[i + 1] / (tension * tension)) * (diff / h));
    }
  
  return value;
}

/* Following four functions compute (6/x^2)(1-x/sinh(x)),
   (3/x^2)(x/tanh(x)-1), (6/x^2)(1-x/sin(x)), and (3/x^2)(x/tan(x)-1) via
   the first three terms of the appropriate power series.  They are used
   when |x|<TRIG_ARG_MIN, to avoid loss of significance.  Errors are
   O(x**6). */
double
#ifdef _HAVE_PROTOS
sinh_func (double x) 
#else
sinh_func (x)
     double x;
#endif
{
  /* use 1-(7/60)x**2+(31/2520)x**4 */
  return 1.0 - (7.0/60.0)*x*x + (31.0/2520.0)*x*x*x*x;
}

double
#ifdef _HAVE_PROTOS
tanh_func (double x) 
#else
tanh_func (x)
     double x;
#endif
{
  /* use 1-(1/15)x**2+(2/315)x**4 */
  return 1.0 - (1.0/15.0)*x*x + (2.0/315.0)*x*x*x*x;
}

double
#ifdef _HAVE_PROTOS
sin_func (double x) 
#else
sin_func (x)
     double x;
#endif
{
  /* use -1-(7/60)x**2-(31/2520)x**4 */
  return -1.0 - (7.0/60.0)*x*x - (31.0/2520.0)*x*x*x*x;
}

double
#ifdef _HAVE_PROTOS
tan_func (double x) 
#else
tan_func (x)
     double x;
#endif
{
  /* use -1-(1/15)x**2-(2/315)x**4 */
  return -1.0 - (1.0/15.0)*x*x - (2.0/315.0)*x*x*x*x;
}


/* Following two functions compute (6/y^2)(sinh(xy)/sinh(y)-x) and
   (6/y^2)(sin(xy)/sin(y)-x), via the first three terms of the appropriate
   power series in y.  They are used when |y|<TRIG_ARG_MIN, to avoid loss
   of significance.  Errors are O(y**6). */
double
#ifdef _HAVE_PROTOS
quotient_sinh_func (double x, double y) 
#else
quotient_sinh_func (x, y)
     double x, y;
#endif
{
  return ((x*x*x-x) + (x*x*x*x*x/20.0 - x*x*x/6.0 + 7.0*x/60.0)*(y*y)
	  + (x*x*x*x*x*x*x/840.0 - x*x*x*x*x/120.0 + 7.0*x*x*x/360.0
	     -31.0*x/2520.0)*(y*y*y*y));
}

double
#ifdef _HAVE_PROTOS
quotient_sin_func (double x, double y) 
#else
quotient_sin_func (x, y)
     double x, y;
#endif
{
  return (- (x*x*x-x) + (x*x*x*x*x/20.0 - x*x*x/6.0 + 7.0*x/60.0)*(y*y)
	  - (x*x*x*x*x*x*x/840.0 - x*x*x*x*x/120.0 + 7.0*x*x*x/360.0
	     -31.0*x/2520.0)*(y*y*y*y));
}

}
