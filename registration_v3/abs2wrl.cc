#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#ifdef linux
#include <malloc.h>
#include <values.h>
#endif
#include <string.h>
#include <math.h>
#include <time.h>

#include "orderedpointset.h"

namespace utw3dface {
int verbose=0;
int debug=0;
}

char *usage="abs2wrl [OPTIONS]\n\
OPTIONS:\n\
\t-i file     read abs point file\n\
\t-o file     write wrl file\n\
\t-r          add 3d points by interpolation\n\
\t-s          write surface instead of point cloud\n\
\t-v          increase verbosity\n\
";

using namespace utw3dface;

int main(int argc,char **argv)
{
	char *ifile=NULL;
	char *ofile=NULL;
	int surface=0;
	int resample=0;

	int c;

	while ((c=getopt(argc,argv,"o:i:vsr"))!=-1)
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
				surface=1;
				break;
			case 'r':
				resample=1;
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
	
	ops.Read(ifile);

	if (ofile!=NULL)
	{
		if (surface)
			ops.WriteVRMLsurface(ofile);
		else
			ops.WriteVRML(ofile);
	}
}
