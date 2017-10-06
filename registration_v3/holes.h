#ifndef HOLES_INCLUDED
#define HOLES_INCLUDED

#include "rangeimage.h"
#include "nose.h"

namespace utw3dface {
double fillholes(RangeImage &ri,Nose &nose);
double fill_big_holes_using_symmetry(RangeImage &ri);
}

#endif // HOLES_INCLUDED
