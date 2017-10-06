/*!
* \file TDFmodules_utw.cc
* \version 0.1.1
* \author Luuk Spreeuwers
*         <a href="mailto:l.j.spreeuwers@utwente.nl">
*         l.j.spreeuwers@utwente.nl</a>,
*         University of Twente, Netherlands
* \date 2007/10/22
* \brief Application programming interface for 3D Face project.
*
* This file describes the application programming interface within the
* 3D face project. This header file is distributed with the corrosponding
* dynamic link library. All API documentation is generated out of this
* source file with the open source tool DOXYGEN. Please do not edit any
* automatically generated HTML (and CHM), LaTeX, PDF, XML or RTF file.
* API changes should <b>only</b> be done in the header file. Any code must
* be documented in DOXYGEN style.
*

* @{
*/
///////////////////////////////////////////////////////////////////////////////

/*! @} */
/*
namespace utw3dface {
int verbose=0;
}
*/

#include "TDFExtract_utw.h"
#include "lda.h"
#include "matrixio.h"
#include "register.h"
#include "rangeimage.h"
#include "holes.h"
#include "filter.h"
#include "nose.h"

using namespace utw3dface;

/*!
* \brief Initialise the Feature Extraction Module
*
* Initialise the Feature Extraction Module. The function reads the necessary
* parameters in a file and allocates a buffer for future parameters communication.
* The parameter file can be 3D Face Partner dependent. This buffer can be
* pre-filled by the function if necessary.
*
*	\param[in] szParamFileName - is the name of the parameter file.
*	             This will not be evaluated for the UTW module, a description for the
*                parameter file can therefore not be given.
*
*	\param[in] pParamBuffer - is the pointer to the memory allocated by the function.
*                This will not be evaluated for the UTW module, a description for pParamBuffer
*                can therefore not be given.
*
* \return TDF_ SUCCESS, TDF_ ERROR, TDF_FAIL_FILE, TDF_FAIL_ALLOC
*/
TDF_EXTRACT_UTW_API int tdfInitExtractModule_utw(const char* szParamFileName,void** pParamBuffer)
{
	LDA *lda=new LDA;
	lda->load(szParamFileName);
	
	*pParamBuffer=(void*)lda;

	return TDF_SUCCESS;
}

/*!
* \brief Returns the  module  version number of the module
*
* Returns the the  module  major and minor version number of the extraction  module.
* A function like this should be implemented for each module / partner
*
* \param[out] majorNumber - is the major part of the version number (represented as integer value)
* \param[out] minorNumber - is the minor part of the version number (represented as integer value)
* \param[out] revision    - is the revision number (represented as integer value)
*
\* return TDF_ SUCCESS, TDF_ ERROR,TDF_FAIL_ALLOC.
*/
TDF_EXTRACT_UTW_API int tdfGetExtractModuleVersion_utw(int *majorNumber, int *minorNumber, int* revision )
{
	*majorNumber=0;
	*minorNumber=1;
	*revision=1;

	return TDF_SUCCESS;
}

/*!
* brief Performs the deletion of  allocated resources
*
* Performs the deletion of previously allocated resources, in particular pParamBuffer.
*
* \param[in] pParamBuffer - is the pointer to the memory allocated by the init function.
*							This pointer is deleted by  tdfTerminateExtractModule.
*							
* \return TDF_ SUCCESS, TDF_ ERROR, TDF_FAIL_FILE, TDF_FAIL_ALLOC.
*/
TDF_EXTRACT_UTW_API int tdfTerminateExtractModule_utw(void** pParamBuffer)
{
	LDA *lda=(LDA*)(*pParamBuffer);
	delete lda;

	return TDF_SUCCESS;
}

//3D FEATURE EXTRACTION

/*!
* \brief Computes and returns the 3D reference.
*
* Computes and returns the 3D reference. The reference is an one dimensional array of float
* values. The quality of the reference is returned along with the reference.
*
* REMARK: HR Texture parameters have been removed as HR Texture should be accessible through dataset
*	    structure in future (discussed on telco 20070612).
*
* \param[in] data3d - represents the normalized 3D data (and attached 2D data) for which the
*					  reference will be computed.
* \param[out] p3DReference - is the pointer to the reference. It must be allocated by the caller.
* \param[out] pQuality - is the pointer to the Quality structure which is
*	               allocated and initialized in the function, must
*	               be deallocated using \see tdfFreeQualityScore.
* \param[in] pParamBuffer - is the pointer to the memory containing some function parameters.
*
* \return TDF_ SUCCESS, TDF_ ERROR, TDF_FAIL_ALLOC.
*/
TDF_EXTRACT_UTW_API int tdfExtract3DReference_utw(const TDFDataset3* data3d, unsigned char* p3DReference, TDFQualityScore** pQuality, void* pParamBuffer )
{
	LDA *lda=(LDA*)pParamBuffer;
	double *features=(double*)p3DReference;

	UnorderedPointSet ups;
	int i;
	for (i=0; i<data3d->numberOfPoints; i++)
		ups.AddPoint(Point3D(data3d->pPoints[i].x,data3d->pPoints[i].y,data3d->pPoints[i].z));

	RangeImage *pri;	
	int x,y;

#ifdef REGISTER
	Register reg(ups,1,1);
	pri=&reg.ri;
#else
	RangeImage ri(Point3D(0,0,0),Point3D(1,0,0),Point3D(0,-1,0),1,1,55,85,110,130);
	ri.AccumulateDepth(ups);
	filter(ri,ups,5);
	Nose nose(10000,10000,0);
	fillholes(ri,nose);

	double rx=55/ri.du;
	double ry=70/ri.dv;
	double cx=ri.ou;
	double cy=ri.ov-20/ri.dv;
	for (y=0; y<ri.height; y++)
	for (x=0; x<ri.width; x++)
	{
		if ((x-cx)*(x-cx)/(rx*rx)+(y-cy)*(y-cy)/(ry*ry)>1)
		{
			ri.Flag(x,y)=0;
			ri.Pixel(x,y)=0;
		}
	}
	pri=&ri;
#endif
	
	double lowerpart=0.25;
	double upperpart=0.75;
	int width=pri->width;
	int height=pri->height;
	int y_min=int(height*lowerpart);
	int y_max=int(height*upperpart);
	CvMat *v=cvCreateMat((y_max-y_min)*width,1,CV_64F);

	for (y=y_min; y<y_max; y++)
	for (x=0; x<width; x++)
		cvSetReal1D(v,(y-y_min)*width+x,pri->Pixel(x,y));

	CvMat *f=lda->features(v);

	for (x=0; x<lda->featurelength(); x++)
		features[x]=cvGetReal1D(f,x);

	cvReleaseMat(&f);
	cvReleaseMat(&v);
	
	return TDF_SUCCESS;

}

/*!
* \brief Retrieves the size of the 3D reference.
*
* Retrieves the size of the 3D reference.
*
*
* \param[out] nSize - is the size of a 3D reference in bytes.
* \param[in] pParamBuffer - is the pointer to the memory containing some function parameters.
*
* \return  TDF_ SUCCESS, TDF_ ERROR.
*/
TDF_EXTRACT_UTW_API int tdfGet3DReferenceSize_utw(int *nSize, void* pParamBuffer)
{
	LDA *lda=(LDA*)pParamBuffer;

	*nSize=lda->featurelength()*sizeof(double);

	return TDF_SUCCESS;
}
