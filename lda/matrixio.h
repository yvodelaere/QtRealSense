#include <cv.h>

#include "asciimatrixio.h"
#include "beematrixio.h"

#ifndef MATRIXIO_DEFINED
#define MATRIXIO_DEFINED

namespace utw3dface {
int write_matrix(const char *filename,CvMat *matrix);
CvMat *read_matrix(const char *filename, int transpose=0);
}

#endif // MATRIXIO_DEFINED
