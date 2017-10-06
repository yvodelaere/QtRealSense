#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>

//#define LITTLE_ENDIAN 1

typedef struct {char a,b,c,d;} EL;
typedef union {float f; EL elem; } NUM;

main(int argc, char **argv)
{
	if (argc!=2)
	{
		fprintf(stderr,"usage: sfi2pgm file.sfi\n");
		exit(1);
	}

    	int i,j,c,width,height,channels;
    	FILE *f;

	char *fname=argv[1];

    	NUM junk;
    	char firstline[1000];
    	char ftype[1000];

    	f = fopen( fname, "rb" );
    	if ( !f ) { printf("Can't open %s\n", fname); exit(1); }

    	fgets(firstline,1000,f);
    	sscanf(firstline,"%s %d %d %d", ftype,&width, &height, &channels);

    	if( !(strcmp(ftype,"CSU_SFI") == 0) )
	{
        	fprintf(stderr,"Wrong filetype: %s\n", ftype);
        	exit(1);
	}

	if (channels!=1)
	{
		fprintf(stderr,"can only convert 1 channel images\n");
		exit(1);
	}

	fprintf(stderr,"converting %d x %d image %s to pgm\n",width,height,fname);

	float *im=new float[width*height];

    	for(j=0;j<height;j++)
        for(i=0;i<width;i++)
	{
                /* read in the correct byte order for floating point values */
                if (BYTE_ORDER == LITTLE_ENDIAN)
		{
                    fread(&(junk.elem.d),1,1,f);
                    fread(&(junk.elem.c),1,1,f);
                    fread(&(junk.elem.b),1,1,f);
                    fread(&(junk.elem.a),1,1,f);
                }
                else{
                    fread(&(junk.elem.a),1,1,f);
                    fread(&(junk.elem.b),1,1,f);
                    fread(&(junk.elem.c),1,1,f);
                    fread(&(junk.elem.d),1,1,f);
                }

		im[j*width+i]=junk.f;
	}

	float min=im[0],max=im[0];
	for (j=0; j<height; j++)
	for (i=0; i<width; i++)
	{
		if (im[j*width+i]<min)
			min=im[j*width+i];
		else if (im[j*width+i]>max)
			max=im[j*width+i];
	}

	fprintf(stderr,"min=%g max=%g\n",min,max);

	printf("P2\n%d %d\n%d\n",width,height,255);

	for (j=0; j<height; j++)
	for (i=0; i<width; i++)
		printf("%d\n",int((im[j*width+i]-min)*255/(max-min)));
}

