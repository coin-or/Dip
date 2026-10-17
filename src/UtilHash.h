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


#ifndef UTIL_HASH_INCLUDED
#define UTIL_HASH_INCLUDED

#include "DecompConfig.h"

#include <string>

DIPLIB_EXPORT std::string UtilCreateStringHash(const int      len,
                                               const double* els,
                                               const int      precision = 6);

DIPLIB_EXPORT std::string UtilCreateStringHash(const int      len,
                                               const int*     ind,
                                               const double* els,
                                               const int      precision = 6);
DIPLIB_EXPORT std::string UtilCreateStringHash(const int      len,
                                               const int*     ind,
                                               const double   els,
                                               const int      precision = 6);
DIPLIB_EXPORT std::string UtilCreateStringHash(const int      len,
                                               const int*     ind,
                                               const double* els,
                                               const char     sense,
                                               const double   rhs,
              				 const double   infinity,
                                               const int      precision = 6);

#endif
