#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include "rangeimage.h"

using namespace utw3dface;

const char *usage="sfidiff [OPTIONS]\n\
OPTIONS:\n\
\t-i file       input sfi file\n\
\t-I file       input sfi file\n\
\t-o file       output sfi file\n\
\t-v            verbose\n\
";

int verbose=0;
int debug=0;

main(int argc,char **argv)
{
	char *ifile1=NULL,*ifile2=NULL,*ofile=NULL;
	int c;

	while ((c=getopt(argc,argv,"i:I:o:v"))!=-1)
	{
		switch(c)
		{
			case 'i': ifile1=optarg; break;
			case 'I': ifile2=optarg; break;
			case 'o': ofile=optarg; break;
			case 'v': verbose++; break;
			default: fprintf(stderr,usage); exit(1);
		}
	}

	if (ifile1==NULL || ifile2==NULL)
	{
		fprintf(stderr,"-i and -I are required options\n");
		exit(1);
	}

	RangeImage rim1,rim2;
	rim1.ReadSFI(ifile1);
	rim2.ReadSFI(ifile2);

	int x,y;
	for (y=0; y<rim1.height; y++)
	for (x=0; x<rim1.width; x++)
		rim1.Pixel(x,y)-=rim2.Pixel(x,y);

	if (ofile)
		rim1.WriteSFI(ofile);
}
