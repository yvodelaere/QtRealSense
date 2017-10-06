#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <malloc.h>
#include <string.h>
#include <values.h>
#include <math.h>
#include <time.h>

#include "orderedpointset.h"

namespace utw3dface
{
	int verbose;
	int debug;
}	

char *usage="absscale [OPTIONS]\n\
OPTIONS:\n\
\t-i file     read abs point file\n\
\t-o file     write abs file\n\
\t-s float    scaling factor\n\
\t-v          increase verbosity\n\
";

using namespace utw3dface;

int main(int argc,char **argv)
{
	char *ifile=NULL;
	char *ofile=NULL;
	double factor=1;

	int c;

	while ((c=getopt(argc,argv,"o:i:vs:"))!=-1)
	{
		switch (c)
		{
			case 'i':
				ifile=optarg;
				break;
			case 'o':
				ofile=optarg;
				break;
			case 'v':
				verbose++;
				break;
			case 's':
				sscanf(optarg,"%lf",&factor);
				break;
			default:
				fprintf(stderr,usage);
				exit(1);
		}
	}

	OrderedPointSet ops;
	if (ifile==NULL)
	{
		fprintf(stderr,"-i is a required option\n");
		exit(1);
	}
	
	ops.ReadAbs(ifile);
	int x,y;
	for (y=0; y<ops.height; y++)
	for (x=0; x<ops.width; x++)
	{
		ops.X[y*ops.width+x]*=factor;
		ops.Y[y*ops.width+x]*=factor;
		ops.Z[y*ops.width+x]*=factor;

	}

	if (ofile!=NULL)
	{
		ops.WriteAbs(ofile);
	}
}
