//===========================================================================//
// This file is part of the DIP Solver Framework.                            //
//                                                                           //
// DIP is distributed under the Eclipse Public License as part of the        //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Authors: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)      //
//          Ted Ralphs, Lehigh University (ted@lehigh.edu)                   //
//          Jiadong Wang, Lehigh University (jiw408@lehigh.edu)              //
//                                                                           //
// Copyright (C) 2002-2019, Lehigh University, Matthew Galati, Ted Ralphs    //
// All Rights Reserved.                                                      //
//===========================================================================//

/* include the COIN-OR-wide system specific configure header */
#include "configall_system.h"

/* this needs to come before the include of config_dip_default.h */
#ifndef DIPLIB_EXPORT
#if defined(_WIN32) && defined(DLL_EXPORT)
#define DIPLIB_EXPORT __declspec(dllexport)
#else
#define DIPLIB_EXPORT
#endif
#endif

/* include the public project specific macros */
#include "config_dip_default.h"

/***************************************************************************/
/*             HERE DEFINE THE PROJECT SPECIFIC MACROS                     */
/*    These are only in effect in a setting that doesn't use configure     */
/***************************************************************************/

/***************************************************************************/
/*             HERE DEFINE MS PRAGMAS TO DISABLE SOME WARNINGS             */
/***************************************************************************/

#include "CoinPragma.hpp"

#if defined(_MSC_VER)
// warning C4290: C++ exception specification ignored except to indicate
// a function is not __declspec(nothrow)
# pragma warning(disable:4290)
//warning C4996: 'std::xxx' was declared deprecated
# pragma warning(disable:4996)
#endif

/***************************************************************************/
/*             HERE DEFINE THE CONFIGURATION SPECIFIC MACROS               */
/***************************************************************************/

/* Define to 1 if the Cgl package is used */
#define DIP_HAS_CGL 1

/* Define to 1 if the CoinUtils package is used */
#define DIP_HAS_COINUTILS 1

/* Define to 1 if the Osi package is available */
#define DIP_HAS_OSI 1

/* Define to the debug sanity check level (0 is no test) */
//#define COIN_DECOMP_CHECKLEVEL 0

/* Define to the debug verbosity level (0 is no output) */
//#define COIN_DECOMP_VERBOSITY 0


