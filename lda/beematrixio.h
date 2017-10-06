#ifndef BEEMATRIXIO_DEFINED
#define BEEMATRIXIO_DEFINED

#include "cv.h"
namespace utw3dface {
void write_beematrix(const char *filename,CvMat *mat,char *probepath,char *gallerypath);

// note: probepath and gallerypath should be allocated by caller and have length>=1024
CvMat *read_beematrix(const char *filename,char *probepath,char *gallerypath);
}

#endif //BEEMATRIXIO_DEFINED
