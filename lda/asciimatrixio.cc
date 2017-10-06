#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#ifdef linux
#include <unistd.h>
#endif
#include <time.h>

#include "cv.h"

namespace utw3dface {
#define MAXLINELENGTH	1000000
#define MAXCOLS		100000
#define MAXROWS		100000

extern int verbose;

int write_asciimatrix(const char *filename,CvMat *matrix)
{
	FILE *f;
	
	f=fopen(filename,"w");
	if (f==NULL)
	{
		fprintf(stderr,"could not open %s for output\n",filename);
		exit(1);
	}

	int i,j;

	for (i=0; i<matrix->rows; i++)
	{
		for (j=0; j<matrix->cols; j++)
			fprintf(f,"%g ",cvGetReal2D(matrix,i,j));
		fprintf(f,"\n");
	}

	fclose(f);

	if (verbose)
		fprintf(stderr,"wrote %dx%d matrix to %s\n",matrix->rows,matrix->cols,filename);
	return 1;
}

CvMat *read_asciimatrix(const char *filename, int transpose=0)
{
	FILE *f;
	char *line=new char[MAXLINELENGTH];
	double *row=new double[MAXCOLS];
	double **mat=new double*[MAXROWS];
	int i,j;

	f=fopen(filename,"r");
	if (f==NULL)
	{
		fprintf(stderr,"could not open %s for input\n",filename);
		exit(1);
	}

	int rows=0,cols=-1;
	while (fgets(line,MAXLINELENGTH-1,f) != NULL)
	{
		char *s;
		int i=0;
		for (s=strtok(line," \t\n"); s!=NULL; s=strtok(NULL," \t\n"))
		{
			row[i++]=atof(s);
			if (i>=MAXCOLS)
			{
				fprintf(stderr,"%s: too many columns\n",filename);
				exit(1);
			}
		}
		if (cols==-1)
			cols=i;
		else if (i!=cols)
		{
			fprintf(stderr,"%s: lines with different numbers of colums!\n",filename);
			exit(1);
		}
		mat[rows]=new double[cols];
		for (i=0; i<cols; i++)
			mat[rows][i]=row[i];
		rows++;
		if (rows>=MAXROWS)
		{
			fprintf(stderr,"%s: too many rows\n",filename);
			exit(1);
		}
	}

	fclose(f);

	CvMat *matrix;
	if (transpose)
	{
		matrix=cvCreateMat(cols,rows,CV_64F);
		for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
			cvSetReal2D(matrix,j,i,mat[i][j]);
	}
	else
	{
		matrix=cvCreateMat(rows,cols,CV_64F);
		for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
			cvSetReal2D(matrix,i,j,mat[i][j]);
	}

	for (i=0; i<rows; i++)
		delete [] mat[i];

	if (verbose)
		fprintf(stderr,"read %dx%d matrix from %s\n",matrix->rows,matrix->cols,filename);

	delete [] line;
	delete [] mat;
	delete [] row;

	return matrix;
}
}
