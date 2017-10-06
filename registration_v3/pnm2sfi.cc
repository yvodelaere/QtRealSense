#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include "rangeimage.h"

using namespace utw3dface;

const char *usage="pbm2sfi [OPTIONS]\n\
OPTIONS:\n\
\t-i file       input pnm file\n\
\t-o file       output sfi file\n\
\t-v            verbose\n\
";

int verbose=0;
int debug=0;

main(int argc,char **argv)
{
	char *ifile1=NULL,*ofile=NULL;
	int c;

	while ((c=getopt(argc,argv,"i:o:v"))!=-1)
	{
		switch(c)
		{
			case 'i': ifile1=optarg; break;
			case 'o': ofile=optarg; break;
			case 'v': verbose++; break;
			default: fprintf(stderr,usage); exit(1);
		}
	}

	if (ifile1==NULL)
	{
		fprintf(stderr,"-i is a required options\n");
		exit(1);
	}

	RangeImage rim1;
	rim1.ReadPNM(ifile1);

	if (ofile)
		rim1.WriteSFI(ofile);
}
