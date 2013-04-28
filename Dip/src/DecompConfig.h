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
// Copyright (C) 2002-2013, Lehigh University, Matthew Galati, Ted Ralphs    //
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

#ifndef __DIPCONFIG_H__

#ifdef HAVE_CONFIG_H
#ifdef DIP_BUILD
#include "config.h"
#else
#include "config_dip.h"
#endif

#else /* HAVE_CONFIG_H */

#ifdef DIP_BUILD
#include "config_default.h"
#else
#include "config_dip_default.h"
#endif

#endif /* HAVE_CONFIG_H */

#endif /*__DIPCONFIG_H__*/
