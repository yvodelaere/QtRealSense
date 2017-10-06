#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "point3d.h"
#include "orderedpointset.h"
#include "rangeimage.h"

const char *usage="transformabs [OPTIONS]\n\
OPTIONS\n\
\t-i file     input abs point file\n\
\t-o file     output abs point file\n\
\t-s float    scale \n\
\t-x float    translation in x-direction\n\
\t-y float    translation in y-direction\n\
\t-z float    translation in z-direction\n\
\t-X float    rotation around x-axis (gamma)\n\
\t-Y float    rotation around y-axis (phi)\n\
\t-Z float    rotation around z-axis (theta)\n\
\t-v          increase verbosity\n\
";

using namespace utw3dface;

namespace utw3dface {
int verbose=0;
int debug=0;
}

int main(int argc,char **argv)
{
	int c;
	char *ifile=NULL;
	char *ofile=NULL;
	double tx=0,ty=0,tz=0,phi=0,theta=0,gamma=0,scale=1;

	while ((c=getopt(argc,argv,"vs:i:o:x:y:z:X:Y:Z:"))!=-1)
	{
		switch (c)
		{
			case 'v':
				verbose++;
				break;
			case 'i':
				ifile=optarg;
				break;
			case 'o':
				ofile=optarg;
				break;
			case 's':
				sscanf(optarg,"%lf",&scale);
				break;
			case 'x':
				sscanf(optarg,"%lf",&tx);
				break;
			case 'y':
				sscanf(optarg,"%lf",&ty);
				break;
			case 'z':
				sscanf(optarg,"%lf",&tz);
				break;
			case 'X':
				sscanf(optarg,"%lf",&gamma);
				break;
			case 'Y':
				sscanf(optarg,"%lf",&phi);
				break;
			case 'Z':
				sscanf(optarg,"%lf",&theta);
				break;
			default:
				fprintf(stderr,"%s\n",usage);
				exit(1);
		}
	}

	Point3D u,v,t;
	phi=-phi;
	theta=-theta;
	u.x=cos(phi)*cos(theta);
	u.y=cos(phi)*sin(theta);
	u.z=sin(phi);
	v.x=sin(phi)*cos(theta)*sin(gamma)-sin(theta)*cos(gamma);
	v.y=sin(phi)*sin(theta)*sin(gamma)+cos(theta)*cos(gamma);
	v.z=-cos(phi)*sin(gamma);
	t=Point3D(tx,ty,tz);

	RangeImage ri(t,u,v,1,1,0,0,220,260);
	
	OrderedPointSet ops;
	if (ifile==NULL)
	{
		fprintf(stderr,"-i is a required option\n");
		exit(1);
	}
	
	ops.Read(ifile);

	int x,y;
	for (y=0; y<ops.height; y++)
	for (x=0; x<ops.width; x++)
	{
		int index=y*ops.width+x;
		if (ops.ValidPoint(x,y))
		{
			Point3D p=ops.GetPoint(x,y);
			Point3D q=ri.TransformPoint(p);
			ops.X[index]=q.x*scale;
			ops.Y[index]=q.y*scale;
			ops.Z[index]=q.z*scale;
		}
	}

	if (ofile!=NULL)
		ops.WriteAbs(ofile);
}

