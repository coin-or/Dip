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

//===========================================================================//
#ifndef DecompMemPool_h_
#define DecompMemPool_h_

#include "CoinError.hpp"

// --------------------------------------------------------------------- //
class DecompMemPool {
public:
   double* dblArrNCoreCols;
   double* dblArrNCoreRows;

public:
   void allocateMemory(const int nCoreCols,
                       const int nCoreRows) {
      if (nCoreCols > 0) {
         dblArrNCoreCols = new double[nCoreCols];
         CoinAssertHint(dblArrNCoreCols, "Error: Out of Memory");
      }

      if (nCoreRows > 0) {
         dblArrNCoreRows = new double[nCoreRows];
         CoinAssertHint(dblArrNCoreRows, "Error: Out of Memory");
      }
   }

public:
   DecompMemPool() :
      dblArrNCoreCols(0),
      dblArrNCoreRows(0) {
   }
   ~DecompMemPool() {
      UTIL_DELARR(dblArrNCoreCols);
      UTIL_DELARR(dblArrNCoreRows);
   }
};

#endif
