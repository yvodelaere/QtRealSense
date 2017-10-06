#include <stdio.h>
#include <spline.h>

using namespace utw3dface;

main()
{
	CoordList points(2);

	float a[]={1,1};
	float b[]={10,30};
	float c[]={70,40};
	float d[]={70,0};
	points.Add(a);
	points.Add(b);
	points.Add(c);
	points.Add(d);

	Spline spline(points,1);

	CoordList *interp=spline.InterpolateRange(0,spline.abscissa->coords[spline.abscissa->npoints-1],500);

	{
		float t;
		for (t=0; t<1; t+=0.1)
		{
			float c=spline.Curvature2D(t,0.01);
			float coords[2];
			
			spline.Interpolate(t,coords);
			printf("t=%6.3g x=%6.3g y=%6.3g c=%6.3g\n",t,coords[0],coords[1],c);
		}
	}

	interp->Write("interp");
	spline.abscissa->Write("abscissa");
	spline.points->Write("pointsout");

	spline.SetAbscissaGeoDist(0);
	spline.abscissa->Write("abscissa0");
	spline.SetAbscissaGeoDist(1);
	spline.abscissa->Write("abscissa1");
	spline.SetAbscissaGeoDist(5);
	spline.abscissa->Write("abscissa5");
	spline.SetAbscissaGeoDist(10);
	spline.abscissa->Write("abscissa10");
	spline.abscissa->Write("abscissa100");
	spline.SetAbscissaGeoDist(100);

	delete interp;
}
