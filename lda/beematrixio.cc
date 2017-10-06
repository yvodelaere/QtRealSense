#include <stdio.h>
#include <stdlib.h>
#ifdef linux
#include <unistd.h>
#endif
#include <string.h>

#include "cv.h"

namespace utw3dface {

extern int verbose;

static int fgetint32(FILE *f)
{
	int i;
	fread(&i,sizeof(int),1,f);
/*
	unsigned char *pi=&i;
	fgetc(pi[0],f);
	fgetc(pi[1],f);
	fgetc(pi[2],f);
	fgetc(pi[3],f);
*/

	return i;
}

static float fgetfloat(FILE *f)
{
	float fl;
	fread(&fl,sizeof(float),1,f);
/*
	unsigned char *pfl=&fl;
	fgetc(pfl[0],f);
	fgetc(pfl[1],f);
	fgetc(pfl[2],f);
	fgetc(pfl[3],f);
*/

	return fl;
}

static void fputint32(int i,FILE *f)
{
	fwrite(&i,sizeof(int),1,f);
/*
	unsigned char *pi=&i;
	fputc(pi[0],f);
	fputc(pi[1],f);
	fputc(pi[2],f);
	fputc(pi[3],f);
*/
}

static void fputfloat(float fl,FILE *f)
{
	fwrite(&fl,sizeof(float),1,f);
/*
	unsigned char *pfl=&fl;
	fputc(pfl[0],f);
	fputc(pfl[1],f);
	fputc(pfl[2],f);
	fputc(pfl[3],f);
*/
}

void write_beematrix(const char *filename,CvMat *mat,char *probepath,char *gallerypath)
{
	FILE *f;

	f=fopen(filename,"w");
	if (f==NULL)
	{
		fprintf(stderr,"write_beematrix: failed to open %s for writing\n",filename);
		exit(1);
	}

	// write header stuff
	fprintf(f,"%s\n",probepath);
	fprintf(f,"%s\n",gallerypath);

	// xml part closed with '\0'
	fputc(0,f);
		
	int storageID=4;
	fputint32(storageID,f);

	int version=1;
	fputint32(version,f);

	fputc('M',f);	// this is a matrix
	fputc('F',f);	// of floats

	fputc(' ',f);

	fprintf(f,"%d ",mat->rows);
	fprintf(f,"%d\n",mat->cols);

	//int delimiter=0;
	//fputc(delimiter,f);

	// write data as floats
	int i,j;
	for (i=0; i<mat->rows; i++)
	for (j=0; j<mat->cols; j++)
		fputfloat(float(cvGetReal2D(mat,i,j)),f);

	fclose(f);

	if (verbose)
		fprintf(stderr,"wrote %d x %d matrix to %s\n",mat->rows,mat->cols,filename);
}


void skipauditheader(FILE *fin)
{
	int c=fgetc(fin);
	ungetc(c,fin);
	if (c!='<')
		return;
		
	if (verbose)
		fprintf(stderr,"beematrixio: skipping auditheader\n");

	char s[1024];
	if (fgets(s,1024,fin)==NULL)
	{
		fprintf(stderr,"beematrixio::skipauditheader: file seems truncated\n");
		exit(1);
	}

	if (s[0]!='<' || s[1]!='?' || s[2]!='x' || s[3]!='m' || s[4]!='l')
	{
		fprintf(stderr,"beematrixio::skipauditheader: no xml header found\n");
		exit(1);
	}

	for (;;)
	{
		if (fgets(s,1024,fin)==NULL)
		{
			fprintf(stderr,"beematrixio::skipauditheader: seems truncated\n");
			exit(1);
		}

		if (strstr(s,"END AUDIT TRAIL") != NULL)
			break;
	}

	// now skip any '0' characters
	while ((c=fgetc(fin))=='\0');
	ungetc(c,fin);
}

CvMat *read_beematrix(const char *filename,char *probepath,char *gallerypath)
{
	FILE *f;

	f=fopen(filename,"r");
	if (f==NULL)
	{
		fprintf(stderr,"read_beematrix: failed to open %s for reading\n",filename);
		exit(1);
	}

	skipauditheader(f);

	// read header stuff
	fgets(probepath,1023,f);
	if (probepath[strlen(probepath)-1]=='\n')
		probepath[strlen(probepath)-1]='\0';
	if (verbose>1)
		fprintf(stderr,"read_beematrix: probepath=%s\n",probepath);
	fgets(gallerypath,1023,f);
	if (gallerypath[strlen(gallerypath)-1]=='\n')
		gallerypath[strlen(gallerypath)-1]='\0';
	if (verbose>1)
		fprintf(stderr,"read_beematrix: gallerypath=%s\n",gallerypath);

	// xml part closed with '\0'
	while (fgetc(f) != '\0');
	if (verbose>1)
		fprintf(stderr,"read_beematrix: skipped xml part\n");

	int storageID=fgetint32(f);
	if (verbose>1)
		fprintf(stderr,"read_beematrix: storageID=%d\n",storageID);
	if (storageID!=4 && storageID!=0x33)
	{
		fprintf(stderr,"read_beematrix: storageID != 4 or 0x33, cannot read this matrix\n");
		exit(1);
	}

	int version=fgetint32(f);
	if (version!=1)
	{
		fprintf(stderr,"read_beematrix: version != 1, cannot read this matrix\n");
		exit(1);
	}
	if (verbose>1)
		fprintf(stderr,"read_beematrix: version=%d\n",version);

	if (fgetc(f)!='M')
	{
		fprintf(stderr,"read_beematrix: not a matrix\n");
		exit(1);
	}
	int type=fgetc(f);
	if (type!='F' && type!='B')
	{
		fprintf(stderr,"read_beematrix: not a matrix of floats/bytes\n");
		exit(1);
	}

	char rowscols[1024];
	fgets(rowscols,1023,f);	// read until end of line

	int rows,cols;
	sscanf(rowscols,"%d",&rows);
	sscanf(rowscols,"%d",&cols);

	if (verbose>1)
		fprintf(stderr,"read_beematrix: rows=%d cols=%d\n",rows,cols);

	// create matrix and read data
	CvMat *mat;
	int i,j;
	switch (type)
	{
		case 'B': 
			mat=cvCreateMat(rows,cols,CV_8U); 
			for (i=0; i<mat->rows; i++)
			for (j=0; j<mat->cols; j++)
				cvSetReal2D(mat,i,j,(double)fgetc(f));
			break;
		case 'F': 
			mat=cvCreateMat(rows,cols,CV_32F); 
			for (i=0; i<mat->rows; i++)
			for (j=0; j<mat->cols; j++)
				cvSetReal2D(mat,i,j,fgetfloat(f));
			break;
	}

	fclose(f);

	if (verbose)
		fprintf(stderr,"read %d x %d matrix from %s\n",mat->rows,mat->cols,filename);

	return mat;
}
}

