#ifndef SYMMETRY_INCLUDED
#define SYMMETRY_INCLUDED

#include <stdio.h>
#include "opt.h"

namespace utw3dface {

class SymmetryScore
{
public:
	double theta;
	double phi;
	double dx;
	double score;
	int flag;

	SymmetryScore():theta(0),phi(0),dx(0),score(0) {}
	SymmetryScore(double theta,double phi,double dx,double score):theta(theta),phi(phi),dx(dx),score(score){}
};

class ArrayOfSymmetryScores
{
public:
	int nscores;
	SymmetryScore *scores;
	
	int Init(int n)
	{
		if (scores!=NULL)
			delete [] scores;

		scores=new SymmetryScore[n];
		nscores=n;
		if (scores==NULL)
		{
			fprintf(stderr,"symmetry::ArrayOfSymmetryScores::constructor - memory allocation failed\n");
			exit(1);
		}

		return 1;
	}

	ArrayOfSymmetryScores()
	{
		nscores=0;
		scores=NULL;
	}

	~ArrayOfSymmetryScores()
	{
		if (scores!=NULL)
			delete [] scores;
	}

	SymmetryScore &operator[](int index)
	{
		if (index<0 || index>=nscores)
		{
			fprintf(stderr,"symmetry::ArrayOfSymmetryScores::operator[]() - index out of range\n");
			exit(1);
		}

		return scores[index];
	}

	void Sort()
	{
		int swapped;
		do
		{
			swapped=0;
			int i;
			for (i=1; i<nscores; i++)
			{
				if (scores[i].score<scores[i-1].score)
				{
					SymmetryScore dummy=scores[i];
					scores[i]=scores[i-1];
					scores[i-1]=dummy;
					swapped++;
				}
			}
		} while (swapped);
	}
};

double symmetryfunction(Opt *opt);
double symmetryfunction2(Opt *opt2);
int rough_symmetry(RangeImage &ri,UnorderedPointSet &ups,ArrayOfSymmetryScores &scores,int nearfrontal=0);
}

#endif // SYMMETRY_INCLUDED
