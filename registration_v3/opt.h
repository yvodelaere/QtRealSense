////////////////////////////////////////////////////////////////////
//
// A class for parabolic parameter optimisation in any number of dimensions
//
// First the function to optimise must be set. This takes a parameter of 
// type Opt, i.e. of the class itself. In this way all class variables are
// available to the function (most important of course the pars variable).
// Next the parameter search ranges, initial values, required accuracies are set.
// Optional data may be passed to the structure (like e.g. samples of a signal).
// Then the optimiser: ParabolicFit() is called. This performs a number (niter) of 
// 1D optimisations for each parameter. The maximum number of iterations for each 1D 
// optimisations is maxpariter (if the accuracy, i.e. a change smaller than the set
// accuracy is reached before maxpariter, the 1D optimisation is finished before maxpariter
// iterations). 1D optimisations for a certain parameter can be performed separately
// if necessary using ParabolicFit1D().
//
// Author: Luuk Spreeuwers, University of Twente
// Date: 2005-2009
// Last update: 30-07-2009
// version: 2.0
//
////////////////////////////////////////////////////////////////////////

#ifndef OPT_INCLUDED
#define OPT_INCLUDED

#include <math.h>

namespace utw3dface {

class Opt
{
public:
	double *pars;
	double *min;
	double *max;
	double *acc;

	int npars;

	double minimum;

	double (*func)(Opt *opt);

	/// data that may be used by the function func()
	void *data;

	Opt(int npars,double (*func)(Opt *opt));
	~Opt();

	int SetPar(int index,double value,double minval,double maxval,double acc);
	int SetParValue(int index,double value);
	int SetFunc(double (*func)(Opt *opt));
	int SetData(void *data);
	double GetParValue(int index);
	double GetFuncValue(double *pars);
	double GetMinimum();

	int ParabolicFit1D(int par,int maxiter);
	int ParabolicFit(int maxpariter,int niter);
};

}
#endif // OPT_INCLUDED
