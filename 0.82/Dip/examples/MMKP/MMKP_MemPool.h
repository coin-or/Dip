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

#ifndef MMKP_MEMPOOL_INCLUDED
#define MMKP_MEMPOOL_INCLUDED

#include "CoinError.hpp"

// --------------------------------------------------------------------- //
class MMKP_MemPool {
 public:
   double            * dblArrNCoreCols;
   
 public:
   void allocateMemory(const int nCoreCols) {

      if(nCoreCols > 0){
         //intArrNCoreCols     = new int[nCoreCols];
         dblArrNCoreCols     = new double[nCoreCols];
         //pIntDblArrNCoreCols = new pair<int,double>[nCoreCols];
         //CoinAssertHint(intArrNCoreCols &&
         //               dblArrNCoreCols &&
         //               pIntDblArrNCoreCols,
         //               "Error: Out of Memory");
      }
   }
   
 public:
   MMKP_MemPool() :      
      dblArrNCoreCols(0)
      {}
   ~MMKP_MemPool()
      {
         UTIL_DELARR(dblArrNCoreCols);       
      }
};

#endif
