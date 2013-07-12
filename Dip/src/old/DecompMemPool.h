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

#ifndef DECOMP_MEMPOOL_INCLUDED
#define DECOMP_MEMPOOL_INCLUDED

#include "CoinError.hpp"

// --------------------------------------------------------------------- //
class DecompMemPool {
public:
   double * dblArrNCoreCols;
   double * dblArrNCoreRows;
   
public:
   void allocateMemory(const int nCoreCols,
                       const int nCoreRows) throw(CoinError) {
      if(nCoreCols > 0){
         dblArrNCoreCols = new double[nCoreCols];
         CoinAssertHint(dblArrNCoreCols, "Error: Out of Memory");
      }
      if(nCoreRows > 0){
         dblArrNCoreRows = new double[nCoreRows];
         CoinAssertHint(dblArrNCoreRows, "Error: Out of Memory");
      }
   }
   
public:
   DecompMemPool() :
      dblArrNCoreCols(0),
      dblArrNCoreRows(0)
   {
   }
   ~DecompMemPool()
   {
      UTIL_DELARR(dblArrNCoreCols);
      UTIL_DELARR(dblArrNCoreRows);
   }
};

#endif
