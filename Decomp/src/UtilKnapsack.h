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

#ifndef UTIL_KNAPSACK_INCLUDED
#define UTIL_KNAPSACK_INCLUDED


/*==========================================================================*/
/*                           SOR_IntDblArr                                  */
/*==========================================================================*/

typedef struct SOR_IntDblT{
  int    i;
  double x;
} SOR_IntDbl;

typedef struct SOR_IntDblArrT{
  SOR_IntDbl * arr;
  int          len;
  int          size;
} SOR_IntDblArr;

typedef SOR_IntDblArr   SOR_IntDblArr;
typedef SOR_IntDblArr * SOR_IntDblArrPtr;

SOR_IntDblArrPtr SOR_IntDblArrNew (int    size,
                                   int  * pstatus);
void SOR_IntDblSwap(SOR_IntDbl * A,
                    SOR_IntDbl * B);
void SOR_IntDblArrPrint(const SOR_IntDblArr * A);                        
void SOR_IntDblArrFree(SOR_IntDblArrPtr * A);


void KnapsackSortRatioOut(const int       n,
                          const double  * p,
                          const double  * w,
                          double        * psort,
                          double        * wsort,
                          SOR_IntDbl    * ratio);
int KnapsackOptimizeHS(const int      n,
                       const double   c,
                       double       * p,
                       double       * w,
                       int          * x,
                       double       * z,
                       int          * pstatus);



#endif
