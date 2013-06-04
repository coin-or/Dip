//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, Lehigh University                                 //
//                                                                           //
// Copyright (C) 2002-2007, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

// --------------------------------------------------------------------- //
#include <sstream>
#include <string>
#include <iomanip>
using namespace std;

#include "UtilMacros.h"

//http://burtleburtle.net/bob/hash/evahash.html or just map?
//if this is not really doing hashing, then move this function into
//util macros

//---
//--- NOTE:
//---  There is a memory leak in stringstream with MSVS-2005.
//---  MSVS-2005/SP1 fixes the issue.
//---

// --------------------------------------------------------------------- //
string UtilCreateStringHash(const int      len,
                            const double* els,
                            const int      precision)
{
   stringstream ss;
   ss << setprecision(precision);

   for (int i = 0; i < len; i++) {
      if (!UtilIsZero(els[i])) {
         ss << i << "_" << els[i] << "_";
      }
   }

   return ss.str();
}

// --------------------------------------------------------------------- //
string UtilCreateStringHash(const int      len,
                            const int*     ind,
                            const double* els,
                            const int      precision)
{
   stringstream ss;
   ss << setprecision(precision);

   for (int i = 0; i < len; i++) {
      if (!UtilIsZero(els[i])) {
         ss << ind[i] << "_" << els[i] << "_";
      }
   }

   return ss.str();
}

// --------------------------------------------------------------------- //
string UtilCreateStringHash(const int      len,
                            const int*     ind,
                            const double* els,
                            const char     sense,
                            const double   rhs,
                            const int      precision)
{
   stringstream ss;
   ss << setprecision(precision);
   ss << rhs << "_" << sense << "_";
   ss << UtilCreateStringHash(len, ind, els, precision);
   return ss.str();
}


