#include <stdio.h>
#include <jpeglib.h>
#include <setjmp.h>
#include <X11/Xarch.h>

#include "rangeimage.h"

namespace utw3dface {

enum ImageType {PBM,PGM,PPM};

void RangeImage::calculateT()
{
	Point3D n=u^v;
	n.Normalise();
	T[0][0]=u.x/du; T[0][1]=u.y/du; T[0][2]=u.z/du; T[0][3]=-(u*o)/du+ou;
	T[1][0]=v.x/dv; T[1][1]=v.y/dv; T[1][2]=v.z/dv; T[1][3]=-(v*o)/dv+ov;
	T[2][0]=n.x; T[2][1]=n.y; T[2][2]=n.z; T[2][3]=-(n*o);
}

void RangeImage::Clear()
{
	memset(front,0,sizeof(double)*width*height);
	memset(back,0,sizeof(double)*width*height);
	memset(pixel,0,sizeof(double)*width*height);
	memset(stddev,0,sizeof(double)*width*height);
	memset(flag,0,sizeof(double)*width*height);
	memset(frontflag,0,sizeof(double)*width*height);
	memset(backflag,0,sizeof(double)*width*height);
}

void RangeImage::Free()
{
	if (pixel!=NULL)
		delete [] pixel;
	if (stddev!=NULL)
		delete [] stddev;
	if (flag!=NULL)
		delete [] flag;
	if (frontflag!=NULL)
		delete [] frontflag;
	if (backflag!=NULL)
		delete [] backflag;
	if (front!=NULL)
		delete [] front;
	if (back!=NULL)
		delete [] back;
	if (jpeg_buf!=NULL)
		delete [] jpeg_buf;
	pixel=NULL;
	flag=NULL;
	stddev=NULL;
	front=NULL;
	back=NULL;
	frontflag=NULL;
	backflag=NULL;
}

void RangeImage::Init(const Point3D &o,const Point3D &u,const Point3D &v,double du,double dv,double ou,double ov,int width,int height)
{
	Free();
	this->o=o;
	this->u=u;
	this->v=v;
	this->du=du;
	this->dv=dv;
	this->ou=ou;
	this->ov=ov;
	this->width=width;
	this->height=height;

	this->u.Normalise();
	this->v.Normalise();
	calculateT();

	pixel=new double[width*height];
	stddev=new double[width*height];
	flag=new double[width*height];
	front=new double[width*height];
	back=new double[width*height];
	frontflag=new double[width*height];
	backflag=new double[width*height];
	jpeg_buf=new unsigned char[width*height];
	Clear();
	ups=NULL;
}

RangeImage::RangeImage(const Point3D &o,const Point3D &u,const Point3D &v,double du,double dv,double ou,double ov,int width,int height)
{
	pixel=NULL;
	stddev=NULL;
	flag=NULL;
	front=NULL;
	back=NULL;
	frontflag=NULL;
	backflag=NULL;
	jpeg_buf=NULL;
	Init(o,u,v,du,dv,ou,ov,width,height);
}

RangeImage::RangeImage()    
{
	o=Point3D(0,0,0);
	u=Point3D(1,0,0);
	v=Point3D(0,-1,0);
	du=dv=1;
	ou=ov=0;
	width=height=0;
	pixel=NULL;
	stddev=NULL;
	flag=NULL;
	front=NULL;
	back=NULL;
	frontflag=NULL;
	backflag=NULL;
	ups=NULL;
	jpeg_buf=NULL;
}

RangeImage::RangeImage(int width,int height)
{
	o=Point3D(0,0,0);
	u=Point3D(1,0,0);
	v=Point3D(0,-1,0);
	du=dv=1;
	ou=ov=0;
	pixel=NULL;
	stddev=NULL;
	flag=NULL;
	front=NULL;
	back=NULL;
	frontflag=NULL;
	backflag=NULL;
	ups=NULL;
	jpeg_buf=NULL;
	Init(o,u,v,du,dv,ou,ov,width,height);
}

void RangeImage::CopyDataFrom(RangeImage &ri)
{
	if (ri.width!=width || ri.height!=height)
		Init(ri.o,ri.u,ri.v,ri.du,ri.dv,ri.ou,ri.ov,ri.width,ri.height);

	memcpy(pixel,ri.pixel,sizeof(double)*ri.width*ri.height);
}

RangeImage::~RangeImage()   
{
	Free();
}

Point3D RangeImage::TransformPoint(const Point3D &p)
{
	Point3D uvd;
	uvd.x=T[0][0]*p.x+T[0][1]*p.y+T[0][2]*p.z+T[0][3];
	uvd.y=T[1][0]*p.x+T[1][1]*p.y+T[1][2]*p.z+T[1][3];
	uvd.z=T[2][0]*p.x+T[2][1]*p.y+T[2][2]*p.z+T[2][3];

	return uvd;
}

void RangeImage::SetDepth(const Point3D &p)
{
	Point3D q=TransformPoint(p);
	int x=int(q.x+0.5);
	int y=int(q.y+0.5);
	if (x<0 || x>=width || y<0 || y>=height)
		return;
	pixel[y*width+x]=q.z;
	flag[y*width+x]=255;
}

void RangeImage::AccumulateDepth(const Point3D &p)
{
	Point3D q=TransformPoint(p);
/*      this is some code to distribute the contribution of a single point over 4 neighbouring gridpoints
it turned out that this did not give significant improvement over simple rounding and it is much slower
the code is left here for reference. Note that it is older than the code below
	int x=int(q.x);
	int y=int(q.y);
	if (x<0 || x>=width || y<0 || y>=height)
		return;
	double dx=q.x-x;
	double dy=q.y-y;
	int index;
	double fraction;

	index=y*width+x;
	fraction=(1-dx)*(1-dy);
	pixel[index]+=q.z*fraction;
	stddev[index]+=q.z*q.z*fraction;
	flag[index]+=fraction;

	if (x+1<width)
	{
		index=y*width+x+1;
		fraction=dx*(1-dy);
		pixel[index]+=q.z*fraction;
		stddev[index]+=q.z*q.z*fraction;
		flag[index]+=fraction;
	}

	if (y+1<height)
	{
		index=(y+1)*width+x;
		fraction=(1-dx)*dy;
		pixel[index]+=q.z*fraction;
		stddev[index]+=q.z*q.z*fraction;
		flag[index]+=fraction;
	}
	
	if (x+1<width && y+1<height)
	{
		index=(y+1)*width+x+1;
		fraction=dx*dy;
		pixel[index]+=q.z*fraction;
		stddev[index]+=q.z*q.z*fraction;
		flag[index]+=fraction;
	}
*/
	int x=int(q.x+0.5);
	int y=int(q.y+0.5);
	if (x<0 || x>=width || y<0 || y>=height)
		return;
	int index=y*width+x;
	double n=flag[index];
	if (n>0)
	{
		double d=q.z-pixel[index];
		if (d>50)
		{
			back[index]=(back[index]*backflag[index]+q.z)/(backflag[index]+1);
			backflag[index]++;
		}
		else if (d<-50)
		{
			back[index]=(back[index]*backflag[index]+pixel[index]*n)/(n+backflag[index]+1);
			backflag[index]+=n;
			pixel[index]=q.z;
			stddev[index]=q.z*q.z;
			flag[index]=1;
		}
		else
		{
			pixel[index]=(pixel[index]*n+q.z)/(n+1);
			stddev[index]=(stddev[index]*n+q.z*q.z)/(n+1);
			flag[index]++;
		}
	}
	else
	{
		pixel[index]=(pixel[index]*n+q.z)/(n+1);
		stddev[index]=(stddev[index]*n+q.z*q.z)/(n+1);
		flag[index]++;
	}
}

void RangeImage::AccumulateDepth(double minflag)
{
	if (ups==NULL)
		return;
	Clear();
	int i;
	for (i=0; i<ups->npoints; i++)
		AccumulateDepth((*ups)[i]);
	if (minflag<=0)
		minflag=0.001;
	for (i=0; i<height*width; i++)
	{
		if (flag[i]>=minflag)
			stddev[i]=sqrt(stddev[i]-pixel[i]*pixel[i]);
		else
		{
			if (backflag[i]>=minflag)
			{
				pixel[i]=back[i];
				stddev[i]=0;
				flag[i]=backflag[i];
			}
			else
			{
				pixel[i]=0;
				stddev[i]=0;
				flag[i]=0;
			}
		}
	}
}


void RangeImage::SetDepth()
{
	if (ups==NULL)
		return;
	Clear();
	int i;
	for (i=0; i<ups->npoints; i++)
		SetDepth((*ups)[i]);
}

void RangeImage::SetDepth(UnorderedPointSet &ups)
{
	this->ups=&ups;
	SetDepth();
}

double RangeImage::GetDepth(double dx,double dy,double &flag)
{
	int x=int(dx+0.5);
	int y=int(dy+0.5);
	if (x<0 || x>=width || y<0 || y>=height)
	{
		flag=0;
		return 0;
	}
	flag=this->flag[y*width+x];
	return pixel[y*width+x];
}

double RangeImage::GetStddev(double dx,double dy)
{
	int x=int(dx+0.5);
	int y=int(dy+0.5);
	if (x<0 || x>=width || y<0 || y>=height)
		return 0;

	return stddev[y*width+x];
}

double RangeImage::GetDepth(double dx,double dy)
{
	double flag;
	return GetDepth(dx,dy,flag);
}

void RangeImage::AccumulateDepth(UnorderedPointSet &ups,double minflag)
{
	this->ups=&ups;
	AccumulateDepth(minflag);
}

int RangeImage::WritePGM(const char *filename,int stretch)
{
	FILE *f;
	f=fopen(filename,"w");
	if (f==NULL)
	{
		fprintf(stderr,"RangeImage::WritePGM - failed to open %s\n",filename);
		exit(1);
	}

	double min=pixel[0];
	double max=pixel[0];
	int i;
	for (i=0; i<width*height; i++)
	{
		if (pixel[i]>max)
			max=pixel[i];
		else if (pixel[i]<min)
			min=pixel[i];
	}

	if (!stretch)
	{
		max=min+255;
	}
	//fprintf(stderr,"min=%g max=%g\n",min,max);

	double offset=min;
	double scale=(255.0/(max-min));
	fprintf(f,"P2\n%d %d\n255\n",width,height);
	for (i=0; i<width*height; i++)
	{
		int l=int((pixel[i]-offset)*scale);
		if (l<0) l=0; else if (l>255) l=255;
		fprintf(f,"%d ",l);
		if (i%10==9)
			fprintf(f,"\n");
	}

	fclose(f);

	return 1;
}

int RangeImage::ReadPNM(const char *filename)
{
	FILE *f;
	f=fopen(filename,"r");
	if (f==NULL)
	{
		fprintf(stderr,"RangeImage::ReadPGM: could not open %s for input\n",filename);
		exit(1);
	}

	int type,raw;

	// read header
	char magic[64];
	fscanf(f,"%s",magic);
	if (magic[0]!='P' || magic[2]!='\0')
	{
		fprintf(stderr,"RangeImage::ReadPGM: %s is not a ppm or pgm file\n", filename);
		exit(1);
	}

	switch (magic[1])
	{
		case '1': type=PBM; raw=0; break;
		case '2': type=PGM; raw=0; break;
		case '3': type=PPM; raw=0; break;
		case '4': type=PBM; raw=1; break;
		case '5': type=PGM; raw=1; break;
		case '6': type=PPM; raw=1; break;
		default:
			fprintf(stderr,"RangeImage::ReadPGM: %s is not a ppm or pgm file\n",filename);
			exit(1);
	}

	int maxval=1;
	char s[1024];
	for (;;)
	{ 
		if (fscanf(f,"%s",s)!=1)
		{
			fprintf(stderr,"RangeImage::ReadPGM: error reading %s\n",filename);
			exit(1);
		}
		if (s[0]=='#') 
			fgets(s,1024,f); 
		else
			break;
	} 
	sscanf(s,"%d",&width);
	
	for (;;)
	{ 
		if (fscanf(f,"%s",s)!=1)
		{
			fprintf(stderr,"RangeImage::ReadPGM: error reading %s\n",filename);
			exit(1);
		}
		if (s[0]=='#') 
			fgets(s,1024,f); 
		else
			break;
	} 
	sscanf(s,"%d",&height);

	if (type!=PBM)
	{
		for (;;)
		{ 
			if (fscanf(f,"%s",s)!=1)
			{
				fprintf(stderr,"RangeImage::ReadPGM: error reading %s\n",filename);
				exit(1);
			}
			if (s[0]=='#') 
				fgets(s,1024,f); 
			else
				break;
		} 
		sscanf(s,"%d",&maxval);
	}

	// skip to the end of the line
	fgets(s,1024,f); 

	// allocate the image
	Init(Point3D(0,0,0),Point3D(1,0,0),Point3D(0,-1,0),1,1,0.5*width,0.5*height,width,height);

	int x,y,val,r,g,b;
	switch (type)
	{
		case PBM:
			if (raw)
			{
				for (y=0; y<height; y++)
				for (x=0; x<width; x+=8)
				{
					val=fgetc(f);
					for (b=0; b<8; b++)
					{
						if (x+b<width)
							Pixel(x+b,y)=(val & 1<<(7-b)) ? 1 : 0;
					}
				}
			}
			else
			{
				for (y=0; y<height; y++)
				for (x=0; x<width; x++)
				{
					fscanf(f,"%d",&val);
						Pixel(x,y)=val;
				}
			}
			break;
		case PGM:
			if (raw)
			{
				for (y=0; y<height; y++)
				for (x=0; x<width; x++)
				{
					val=fgetc(f);
					Pixel(x,y)=val;
				}
			}
			else
			{
				for (y=0; y<height; y++)
				for (x=0; x<width; x++)
				{
					fscanf(f,"%d",&val);
					Pixel(x,y)=val;
				}
			}
			break;
		case PPM:
			if (raw)
			{
				for (y=0; y<height; y++)
				for (x=0; x<width; x++)
				{
					r=fgetc(f);
					g=fgetc(f);
					b=fgetc(f);
					Pixel(x,y)=(r+g+b)/3;
				}
			}
			else
			{
				for (y=0; y<height; y++)
				for (x=0; x<width; x++)
				{
					fscanf(f,"%d",&r);
					fscanf(f,"%d",&g);
					fscanf(f,"%d",&b);
					Pixel(x,y)=(r+g+b)/3;
				}
			}
			break;
	}

	return 1;
}

typedef struct {char a,b,c,d;} EL;
typedef union {float f; EL elem; } NUM;

int RangeImage::WriteSFI(const char *filename)
{
	FILE *f;
	f=fopen(filename,"wb");
	if (f==NULL)
	{
		fprintf(stderr,"RangeImage::WriteSFI - failed to open %s\n",filename);
		exit(1);
	}

	fprintf(f,"%s %d %d %d\n","CSU_SFI",width,height,1);

	NUM junk;
	int x,y;
	for(y=0;y<height;y++)
	for(x=0;x<width;x++)
	{
		junk.f=(float)Pixel(x,y);

		if (BYTE_ORDER == LITTLE_ENDIAN)
		{
			fwrite(&(junk.elem.d),1,1,f);
			fwrite(&(junk.elem.c),1,1,f);
			fwrite(&(junk.elem.b),1,1,f);
			fwrite(&(junk.elem.a),1,1,f);
		}
		else
		{
			fwrite(&(junk.elem.a),1,1,f);
			fwrite(&(junk.elem.b),1,1,f);
			fwrite(&(junk.elem.c),1,1,f);
			fwrite(&(junk.elem.d),1,1,f);
		}
	}

	fclose(f);

	return 1;
}

int RangeImage::ReadSFI(const char *fname)
{
	int i,j,channels;
	FILE *f;

	NUM junk;
	char firstline[1000];
	char ftype[1000];

	f = fopen( fname, "rb" );
	if ( !f ) { fprintf(stderr,"Can't open %s\n", fname); exit(1); }

	fgets(firstline,1000,f);
	sscanf(firstline,"%s %d %d %d", ftype,&width, &height, &channels);

	if( !(strcmp(ftype,"CSU_SFI") == 0) )
	{
		fprintf(stderr,"Wrong filetype: %s\n", ftype);
		exit(1);
	}

	if (channels!=1)
	{
		fprintf(stderr,"can only display 1 channel images\n");
		exit(1);
	}

	Init(Point3D(0,0,0),Point3D(1,0,0),Point3D(0,-1,0),1,1,0.5*width,0.5*height,width,height);

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

		Pixel(i,j)=junk.f;
	}

	fclose(f);

	return 1;
}

int RangeImage::WriteFlagPGM(const char *filename)
{
	FILE *f;
	f=fopen(filename,"w");
	if (f==NULL)
	{
		fprintf(stderr,"RangeImage::WriteFlagsPGM - failed to open %s\n",filename);
		exit(1);
	}

	fprintf(f,"P2\n%d %d\n255\n",width,height);
	int i;
	for (i=0; i<width*height; i++)
		fprintf(f,"%d ",int(flag[i]));

	fclose(f);

	return 1;
}

int RangeImage::WriteStddevPGM(const char *filename,int stretch)
{
	FILE *f;
	f=fopen(filename,"w");
	if (f==NULL)
	{
		fprintf(stderr,"RangeImage::WriteStddevPGM - failed to open %s\n",filename);
		exit(1);
	}

	double min=0;
	double max=stddev[0];
	int i;
	for (i=0; i<width*height; i++)
	{
		if (stddev[i]>max)
			max=stddev[i];
	}

	if (!stretch)
		max=min+255;

	double offset=min;
	double scale=(255.0/(max-min));
	fprintf(f,"P2\n%d %d\n255\n",width,height);
	for (i=0; i<width*height; i++)
	{
		int l=int((stddev[i]-offset)*scale);
		if (l<0) l=0; else if (l>255) l=255;
		fprintf(f,"%d ",l);
		if (i%10==9)
			fprintf(f,"\n");
	}

	fclose(f);

	return 1;
}

static void memdest_init(j_compress_ptr cinfo)
{
}

static boolean memdest_empty(j_compress_ptr cinfo)
{
	return TRUE;
}

static void memdest_term(j_compress_ptr cinfo)
{
}

static void memsrc_init(j_decompress_ptr cinfo)
{
}

static boolean memsrc_fill(j_decompress_ptr cinfo)
{
	return TRUE;
}

static void memsrc_skip(j_decompress_ptr cinfo, long num_bytes)
{
	cinfo->src->next_input_byte+=num_bytes;
	cinfo->src->bytes_in_buffer-=num_bytes;
}

static void memsrc_term(j_decompress_ptr cinfo)
{
}

void RangeImage::CompressDecompress(int quality)
{
	jpeg_buf_size=width*height;
	unsigned char *line_buf=new unsigned char[width];

	struct jpeg_destination_mgr memdest={jpeg_buf,jpeg_buf_size,memdest_init,memdest_empty,memdest_term};

	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;

	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);

	cinfo.dest=&memdest;

	cinfo.image_width = width;
	cinfo.image_height = height;
	cinfo.input_components = 1;	
	cinfo.in_color_space = JCS_GRAYSCALE;
	jpeg_set_defaults(&cinfo);
	jpeg_set_quality(&cinfo, quality,TRUE);

	jpeg_start_compress(&cinfo, TRUE);

	while (cinfo.next_scanline < cinfo.image_height) 
	{
		int i;
		for (i=0; i<width; i++)
		{
			double l=3*Pixel(i,cinfo.next_scanline);
			line_buf[i]=(l<0) ? 0 : (l>255) ? 255 : int(l);
		}
		jpeg_write_scanlines(&cinfo, &line_buf, 1);
	}

	jpeg_finish_compress(&cinfo);
	jpeg_buf_size=width*height-cinfo.dest->free_in_buffer;

	jpeg_destroy_compress(&cinfo);

	struct jpeg_decompress_struct dinfo;

  	dinfo.err = jpeg_std_error(&jerr);
  	jpeg_create_decompress(&dinfo);

	struct jpeg_source_mgr memsrc={jpeg_buf,jpeg_buf_size,memsrc_init,memsrc_fill,memsrc_skip,jpeg_resync_to_restart,memsrc_term};

	dinfo.src=&memsrc;

  	jpeg_read_header(&dinfo, TRUE);
  	jpeg_start_decompress(&dinfo);

	if (dinfo.output_width!=width || dinfo.output_height!=height || dinfo.output_components!=1)
	{
		fprintf(stderr,"rangeimage: error in compressdecompress\n");
		exit(1);
	}

  	while (dinfo.output_scanline < dinfo.output_height) 
	{
		int y=dinfo.output_scanline;
		jpeg_read_scanlines(&dinfo, &line_buf, 1);
		int x;
		for (x=0; x<width; x++)
			Pixel(x,y)=line_buf[x]/3.0;
	}	

  	jpeg_finish_decompress(&dinfo);
  	jpeg_destroy_decompress(&dinfo);

	delete [] line_buf;
}
}
