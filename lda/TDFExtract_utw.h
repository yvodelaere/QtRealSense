/*!
* \file TDFmodules_utw.h
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

#ifndef THREE_D_FACE_EXTRACT_MODULE_UTW_H
#define THREE_D_FACE_EXTRACT_MODULE_UTW_H

#include "TDFdef.h"

#ifdef _WIN32
  #ifdef  TDF_EXTRACT_UTW_EXPORTS
    #define TDF_EXTRACT_UTW_API __declspec(dllexport)
  #else
    #define TDF_EXTRACT_UTW_API __declspec(dllimport)
  #endif
#endif




///////////////////////////////////////////////////////////////////////////////
/*!
* \addtogroup FEATURE_EXTRACTION
*
* Feature extraction uses 3D data to produce the references used
* during the comparison procedure (authentication). Feature extraction module
* is mainly a black box, however it has to propose some interfaces to allow
* the system to get the references.
*
* PROPOSAL: As the feature extraction and matching (authentication) functios are
*	         strongly related and cannot be mixed for different version of the modules
*			 the feature extraction and the authent1 modules are integrated into and
*           provided as a single module.
*			 E.g. only one shared library (dll) containing all functions is provided.
* @{
*/
///////////////////////////////////////////////////////////////////////////////

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
TDF_EXTRACT_UTW_API int tdfInitExtractModule_utw(const char* szParamFileName,void** pParamBuffer);

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
TDF_EXTRACT_UTW_API int tdfGetExtractModuleVersion_utw(int *majorNumber, int *minorNumber, int* revision );
/*!
* brief Performs the deletion of  allocated resources
*
* Performs the deletion of previously allocated resources, in particular pParamBuffer.
*
* \param[in] pParamBuffer - is the pointer to the memory allocated by the init function.
*							This pointer is deleted by  tdfTerminateExtractModule.
*							Will not be evaluated for the UTW module.
*
* \return TDF_ SUCCESS, TDF_ ERROR, TDF_FAIL_FILE, TDF_FAIL_ALLOC.
*/
TDF_EXTRACT_UTW_API int tdfTerminateExtractModule_utw(void** pParamBuffer);

//3D FEATURE EXTRACTION

/*!
* \brief Computes and returns the 3D reference.
*
* Computes and returns the 3D reference. The reference is an one dimensional array of float
* values. The quality of the reference is returned along with the reference.
*
* REMARK: HR Texture parameters have been removed as HR Texture should be accessible through dataset
*            structure in future (discussed on telco 20070612).
*
* \param[in] data3d - represents the normalized 3D data (and attached 2D data) for which the
*					  reference will be computed.
* \param[out] p3DReference - is the pointer to the reference. It must be allocated by the caller.
* \param[out] pQuality - is the pointer to the Quality structure which is
*                       allocated and initialized in the function, must
*                       be deallocated using \see tdfFreeQualityScore.
* \param[in] pParamBuffer - is the pointer to the memory containing some function parameters.
*
* \return TDF_ SUCCESS, TDF_ ERROR, TDF_FAIL_ALLOC.
*/
TDF_EXTRACT_UTW_API int tdfExtract3DReference_utw(const TDFDataset3* data3d, unsigned char* p3DReference, TDFQualityScore** pQuality, void* pParamBuffer );

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
TDF_EXTRACT_UTW_API int tdfGet3DReferenceSize_utw(int *nSize, void* pParamBuffer);


/*! @} */


#endif // THREE_D_FACE_EXTRACT_MODULE_UTW_H
