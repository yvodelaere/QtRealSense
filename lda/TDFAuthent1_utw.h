/*!
* \file TDFmodules_utw.h
* \version 0.1.1
* \author Luuk Spreeuwers,
*         <a href="mailto:l.j.spreeuwers@utwente.nl">
*         l.j.spreeuwers@utwente.nl</a>,
*         University of Twente, Netherlands
* \date 2007/02/28
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

#ifndef THREE_D_FACE_AUTH1_MODULE_UTW_H
#define THREE_D_FACE_AUTH1_MODULE_UTW_H

#include "TDFdef.h"



#ifdef _WIN32
  #ifdef  TDF_AUTH1_UTW_EXPORTS
    #define TDF_AUTH1_UTW_API __declspec(dllexport)
  #else
    #define TDF_AUTH1_UTW_API __declspec(dllimport)
  #endif
#endif



///////////////////////////////////////////////////////////////////////////////
/*!
* \addtogroup AUTHENTICATION_MODULES
*
* This is the comparison module. It takes references, compares them and
* retrieves some authentication decision or score. A second step can be added,
* this step performs the fusion of the formers comparison scores or decisions.
* @{
*/
///////////////////////////////////////////////////////////////////////////////

/*!
{
* \brief Initialises the module
*
* Initialises the Authent1 module. The function reads the necessary parameters in a file and allocates a buffer
* for future parameters communication. The parameter file can be 3D Face Partner dependent. This buffer can be
* pre-filled by the function if necessary.
*
*	\param[in] szParamFileName - is the name of the parameter file.
*	             This will not be evaluated for the UTW module, a description for the parameter file can
*                therefore not be given.
*
*	\param[in] pParamBuffer - is the pointer to the memory allocated by the function.
*                This will not be evaluated for the UTW module, a description for pParamBuffer
*                can therefore not be given.
*
* \return TDF_ SUCCESS, TDF_ ERROR, TDF_FAIL_FILE, TDF_FAIL_ALLOC.
*/
TDF_AUTH1_UTW_API int tdfInitAuthent1Module_utw(const char* szParamFileName, void** pParamBuffer);



/*!
*  Returns the  module  version number of the module
*
* Returns the the  module  major and minor version number of the (combined extraction and authentication) module.
* A function like this should be implemented for each Module / Partner
*
* \param[out] majorNumber - is the major part of the version number (represented as integer value)
* \param[out] minorNumber - is the minor part of the version number (represented as integer value)
*
\* return TDF_ SUCCESS, TDF_ ERROR,TDF_FAIL_ALLOC.
*/
TDF_AUTH1_UTW_API int tdfGetAuthent1ModuleVersion_utw(int *majorNumber, int *minorNumber);

/*!
* \brief Performs the deletion of unnecessary memory
*
* Performs the deletion of unnecessary memory, in particular pParamBuffer.
*
* \param[in] pParamBuffer - is the pointer to the memory allocated by the init
*							function. This pointer is deleted by
*							tdfTerminateAuthent1Module.
*                           This will not be evaluated for the UTW module.
*
* \return TDF_ SUCCESS, TDF_ ERROR, TDF_FAIL_FILE, TDF_FAIL_ALLOC.
*
*/

TDF_AUTH1_UTW_API int tdfTerminateAuthent1Module_utw(void** pParamBuffer);

/*!
* brief Authentication based on two 3D references.
*
*         REMARK: Singnature has been changed to use const references as these should not be mofified within
*         this function.
*         Authentication based on two 3D references.
*
* \param[in] p3DReference_1 - is the probe 3D reference used by authentication process.
* \param[in] p3DReference_2 - is the gallery 3D reference used by authentication process.
* \param[out] pScore - is the comparison score. Scores range from zero to one.
* \param[out] pQuality - is the quality of the computed score, which depends from the given references.
* \param[in] pParamBuffer - is the pointer to the memory containing some function parameters.
*
* \return TDF_ SUCCESS, TDF_ ERROR, TDF_FAIL_ALLOC, TDF_FAIL_SIZE, TDF_FAIL_AUTHENT.
*/
TDF_AUTH1_UTW_API int tdf3D3DAuthent1_utw(const unsigned char* p3DReference_1, const unsigned char* p3DReference_2, float* pScore, TDFQualityScore** pQuality, void* pParamBuffer);

/*! @} */

#endif // THREE_D_FACE_AUTH1_MODULE_UTW_H
