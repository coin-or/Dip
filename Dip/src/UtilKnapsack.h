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

#ifndef UTIL_KNAPSACK_INCLUDED
#define UTIL_KNAPSACK_INCLUDED


/*==========================================================================*/
/*                           SOR_IntDblArr                                  */
/*==========================================================================*/

typedef struct SOR_IntDblT {
   int    i;
   double x;
} SOR_IntDbl;

typedef struct SOR_IntDblArrT {
   SOR_IntDbl* arr;
   int          len;
   int          size;
} SOR_IntDblArr;

typedef SOR_IntDblArr   SOR_IntDblArr;
typedef SOR_IntDblArr* SOR_IntDblArrPtr;

SOR_IntDblArrPtr SOR_IntDblArrNew (int    size,
                                   int*   pstatus);
void SOR_IntDblSwap(SOR_IntDbl* A,
                    SOR_IntDbl* B);
void SOR_IntDblArrPrint(const SOR_IntDblArr* A);
void SOR_IntDblArrFree(SOR_IntDblArrPtr* A);


void KnapsackSortRatioOut(const int       n,
                          const double*   p,
                          const double*   w,
                          double*         psort,
                          double*         wsort,
                          SOR_IntDbl*     ratio);
int KnapsackOptimizeHS(const int      n,
                       const double   c,
                       double*        p,
                       double*        w,
                       int*           x,
                       double*        z,
                       int*           pstatus);



#endif
