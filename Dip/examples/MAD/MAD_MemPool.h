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

#ifndef MAD_MEMPOOL_INCLUDED
#define MAD_MEMPOOL_INCLUDED

#include "CoinError.hpp"

// --------------------------------------------------------------------- //
class MAD_MemPool {
 public:
   int               * intArrNOrigRows;
   int               * intArrNCoreCols;
   double            * dblArrNCoreCols;
   pair<int, double> * pIntDblArrNCoreCols;
   double            * dblArrNBlocks;
   
 public:
   void allocateMemory(const int nOrigRows,
                       const int nCoreCols,
                       const int nBlocks    ) throw(CoinError) {
      if(nOrigRows > 0){
         intArrNOrigRows     = new int[nOrigRows];
         CoinAssertHint(intArrNOrigRows, "Error: Out of Memory");
      }  
      if(nCoreCols > 0){
         intArrNCoreCols     = new int[nCoreCols];
         dblArrNCoreCols     = new double[nCoreCols];
         pIntDblArrNCoreCols = new pair<int,double>[nCoreCols];
         CoinAssertHint(intArrNCoreCols &&
                        dblArrNCoreCols &&
                        pIntDblArrNCoreCols,
                        "Error: Out of Memory");
      }
      if(nBlocks > 0){
         dblArrNBlocks = new double[nBlocks];
         CoinAssertHint(dblArrNBlocks, "Error: Out of Memory");
      }
   }
   
public:
   MAD_MemPool() :
      intArrNOrigRows(0),
      intArrNCoreCols(0),
      dblArrNCoreCols(0),
      pIntDblArrNCoreCols(0),
      dblArrNBlocks(0)
      {
      }
   ~MAD_MemPool()
      {
         UTIL_DELARR(intArrNOrigRows);
         UTIL_DELARR(intArrNCoreCols);
         UTIL_DELARR(dblArrNCoreCols);
         UTIL_DELARR(pIntDblArrNCoreCols);
         UTIL_DELARR(dblArrNBlocks);
      }
};

#endif
