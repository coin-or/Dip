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


#ifndef UTIL_HASH_INCLUDED
#define UTIL_HASH_INCLUDED

#include <string>

std::string UtilCreateStringHash(const int      len,
                                 const double* els,
                                 const int      precision = 6);

std::string UtilCreateStringHash(const int      len,
                                 const int*     ind,
                                 const double* els,
                                 const int      precision = 6);
std::string UtilCreateStringHash(const int      len,
                                 const int*     ind,
                                 const double   els,
                                 const int      precision = 6);
std::string UtilCreateStringHash(const int      len,
                                 const int*     ind,
                                 const double* els,
                                 const char     sense,
                                 const double   rhs,
                                 const int      precision = 6);

#endif
