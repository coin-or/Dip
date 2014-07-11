//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)       //
//                                                                           //
// Copyright (C) 2002-2007, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//


#ifndef UTIL_HASH_INCLUDED
#define UTIL_HASH_INCLUDED

#include <string>
using namespace std;

string UtilCreateStringHash(const int      len,
                            const double* els,
                            const int      precision = 6);

string UtilCreateStringHash(const int      len,
                            const int*     ind,
                            const double* els,
                            const int      precision = 6);

string UtilCreateStringHash(const int      len,
                            const int*     ind,
                            const double* els,
                            const char     sense,
                            const double   rhs,
                            const int      precision = 6);

#endif
