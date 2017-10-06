#include <cv.h>

#ifndef ASCIIMATRIXIO_DEFINED
#define ASCIIMATRIXIO_DEFINED

namespace utw3dface {
int write_asciimatrix(const char *filename,CvMat *matrix);
CvMat *read_asciimatrix(const char *filename, int transpose=0);
}

#endif // ASCIIMATRIXIO_DEFINED
