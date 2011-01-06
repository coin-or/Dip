//===========================================================================//
// This file is part of the DIP Solver Framework.                            //
//                                                                           //
// DIP is distributed under the Eclipse Public License as part of the        //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)       //
//                                                                           //
// Conceptual Design: Matthew Galati, SAS Institute Inc.                     //
//                    Ted Ralphs, Lehigh University                          //
//                                                                           //
// Copyright (C) 2002-2009, Lehigh University, Matthew Galati, Ted Ralphs    //
// All Rights Reserved.                                                      //
//===========================================================================//

/*
 * Include file for the configuration of Decomp.
 *
 * On systems where the code is configured with the configure script
 * (i.e., compilation is always done with HAVE_CONFIG_H defined), this
 * header file includes the automatically generated header file, and
 * undefines macros that might configure with other Config.h files.
 *
 * On systems that are compiled in other ways (e.g., with the
 * Developer Studio), a header files is included to define those
 * macros that depend on the operating system and the compiler.  The
 * macros that define the configuration of the particular user setting
 * (e.g., presence of other COIN packages or third party code) are set
 * here.  The project maintainer needs to remember to update this file
 * and choose reasonable defines.  A user can modify the default
 * setting by editing this file here.
 *
 */

#ifndef __DECOMPCONFIG_H__

#ifdef HAVE_CONFIG_H
#include "config_dip.h"

/* undefine macros that could conflict with those in other config.h
   files */
#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION

#else /* HAVE_CONFIG_H */

/* include the COIN-wide system specific configure header */
#include "configall_system.h"

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

/* Define to the debug sanity check level (0 is no test) */
//#define COIN_DECOMP_CHECKLEVEL 0

/* Define to the debug verbosity level (0 is no output) */
//#define COIN_DECOMP_VERBOSITY 0

#define __DECOMP_LP_CLP__
#define __DECOMP_IP_CBC__


#endif /* HAVE_CONFIG_H */

#ifndef DIP_VERSION
#define DIP_VERSION "0.81.2"
#endif

#endif /*__DECOMPCONFIG_H__ */
