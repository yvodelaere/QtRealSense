#ifndef FILTER_INCLUDED
#define FILTER_INCLUDED

#include "rangeimage.h"
#include "unorderedpointset.h"

namespace utw3dface {
double filter(RangeImage &ri,UnorderedPointSet &ups,double maxd=50,int lesssmoothing=0);
double reflectionfilter(RangeImage &ri,UnorderedPointSet &ups,double maxd=5,double mind=2);
double backgroundfilter(RangeImage &ri,UnorderedPointSet &ups,double maxd=100);
int filtershoulders(UnorderedPointSet &ups);
double filterpickfrontal(RangeImage &ri,UnorderedPointSet &ups,double maxvar=5,double maxd=5);
double filterhair(RangeImage &ri,UnorderedPointSet &ups);
}

#endif // FILTER_INCLUDED
