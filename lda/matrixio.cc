#include <stdio.h>
#ifdef linux
#include <unistd.h>
#endif
#include <stdlib.h>
#include <string.h>

#include "matrixio.h"

namespace utw3dface {

int write_matrix(const char *filename,CvMat *matrix)
{
	return write_asciimatrix(filename,matrix);
}

CvMat *read_matrix(const char *filename, int transpose)
{
	int readbee=0;

	FILE *f=fopen(filename,"r");
	if (f==NULL)
	{
		fprintf(stderr,"matrixio: failed to open %s for reading\n",filename);
		exit(1);
	}

	double a;
	if (fscanf(f,"%lf",&a)!=1)
		readbee=1;
	fclose(f);

	if (readbee)
	{
		char ppath[1024],gpath[1024];
		return read_beematrix(filename,ppath,gpath);
	}

	return read_asciimatrix(filename,transpose);
}
}
