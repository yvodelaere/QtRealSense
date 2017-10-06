////////////////////////////////////////////////////////////////////
//
// This header defines some constants/functions etc. to ensure 
// compatibility with MS Visual C++ and other systems
//
// Author: Luuk Spreeuwers, University of Twente
// Date: 2005-2009
// Last update: 29-07-2009
// version: 2.0
//
////////////////////////////////////////////////////////////////////////


#ifndef COMPAT_INCLUDED
#define COMPAT_INCLUDED

#ifdef WIN32
#ifndef __LITTLE_ENDIAN
#define __LITTLE_ENDIAN	8421
#endif
#ifndef __BYTE_ORDER
#define __BYTE_ORDER	__LITTLE_ENDIAN
#endif
#endif

// Visual C++ does not know about times() and the tms structure
#include <time.h>
#ifdef WIN32
struct tms
{
	clock_t tms_utime;
	clock_t tms_stime;
	clock_t tms_cutime;
	clock_t tms_cstime;
};
#define times(buf) ((buf)->tms_utime=clock(),(buf)->tms_stime=(buf)->tms_stime=(buf)->tms_cstime=0,(buf)->tms_utime)
#else
#include <sys/times.h>
#endif

// CLK_TCK is not available on some systems
#ifndef CLK_TCK
#define CLK_TCK CLOCKS_PER_SEC	
#endif

// braindead Visual C++ needs __declspec(dllexport) to export symbols to dll's
#ifndef EXPORT
#ifdef WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT /* nothing */
#endif
#endif
#ifndef IMPORT
#ifdef WIN32
#define IMPORT __declspec(dllimport)
#else
#define IMPORT /* nothing */
#endif
#endif

// to use math constants visual C++ needs _USE_MATH_DEFINES
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

// on win32 systems use float.h to get FLT_MAX etc.
#if WIN32 || linux
#define HAVE_FLOAT
#endif

#ifdef WIN32
#ifndef MAXFLOAT
#define MAXFLOAT      3.402823466e+38F
#endif
#endif

// for Visual C 
#ifdef WIN32
#define isnan(number) _isnan(number)
#endif 

// on win32 systems there is no drand48
#ifdef WIN32
#define drand48() ((double)rand()/(double)RAND_MAX)
#endif

#endif // COMPAT_INCLUDED
