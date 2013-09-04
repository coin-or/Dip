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


//===========================================================================
//---
//--- Knapsack Problem:
//--- max  sum{j = 0..n-1} p[j] x[j]
//--- s.t. sum{j = 0..n-1} w[j] x[j] <= c
//---      x binary
//---
//===========================================================================
#include "UtilMacros.h"
#include "UtilKnapsack.h"

//===========================================================================
//TODO: change this to a class
#define OK             0
#define ERR_NO_MEMORY  1

#define KP_OPTIMAL     0
#define KP_INFEASIBLE  1
#define KP_ERROR       2
#define KP_ITERLIMIT   3

#define EPSILON        1.0e-6
#define HS_MAXITER     1000000
#define BIGM           1.0e+17
#define IS_SMALL( x )  ( fabs( x ) < EPSILON )
#define IS_TINY( x )   ( fabs( x ) < 1.0e-15 )

//===========================================================================
//#define DEBUG_KNAP(x)
#define DEBUG_KNAP(x) x


/*==========================================================================*/
/*                           SOR_IntDblArr                                  */
/*==========================================================================*/

/*==========================================================================*/
SOR_IntDblArrPtr SOR_IntDblArrNew (int    size,
                                   int*   pstatus)
{
   SOR_IntDblArrPtr A = NULL;
   *pstatus = OK;
   A = (SOR_IntDblArrPtr) malloc(sizeof(SOR_IntDblArr));

   if (A == NULL) {
      *pstatus = ERR_NO_MEMORY;
      return A;
   }

   A->len          = 0;
   A->size         = size;
   A->arr          = (SOR_IntDbl*) malloc(size * sizeof(SOR_IntDbl));

   if (A->arr == NULL) {
      *pstatus = ERR_NO_MEMORY;
      return A;
   }

   return A;
}

/*==========================================================================*/
void SOR_IntDblSwap(SOR_IntDbl* A,
                    SOR_IntDbl* B)
{
   SOR_IntDbl temp;
   temp.i = B->i;
   temp.x = B->x;
   B->i   = A->i;
   B->x   = A->x;
   A->i   = temp.i;
   A->x   = temp.x;
}

/*==========================================================================*/
void SOR_IntDblArrPrint(const SOR_IntDblArr*    A)
{
   int i;

   for (i = 0; i < A->len; i++) {
      printf("index: %d, i: %d, x: %g\n",
             i, A->arr[i].i, A->arr[i].x);
   }
}

/*==========================================================================*/
void SOR_IntDblArrFree(SOR_IntDblArrPtr* A)
{
   if (A == NULL) {
      return;
   }

   if ((*A) == NULL) {
      return;
   }

   UTIL_DELARR((*A)->arr);
   UTIL_DELPTR(*A);
}


/*==========================================================================*/
static void SOR_QSortIntDblDec(SOR_IntDbl* item, int lft, int rght)
{
   int i  = lft;
   int j  = rght;
   int lr = (lft + rght) / 2;
   SOR_IntDbl x;
   x.i = item[lr].i;
   x.x = item[lr].x;

   do {
      while (item[i].x > x.x && i < rght) {
         i++;
      }

      while (item[j].x < x.x && j > lft) {
         j--;
      }

      if (i <= j) {
         SOR_IntDblSwap(&item[i], &item[j]);
         i++;
         j--;
      }
   } while (i <= j);

   if (lft < j) {
      SOR_QSortIntDblDec(item, lft, j);
   }

   if (i < rght) {
      SOR_QSortIntDblDec(item, i, rght);
   }
}





//===========================================================================
int KnapsackSortRatio(const int      n,
                      const double* p,
                      const double* w,
                      double*        psort,
                      double*        wsort)
{
   //---
   //--- Sort non-increasing order of p[j]/w[j]. Return the new
   //--- order in psort, wsort which will be of size n+1 to conform
   //--- to the routine SOR_KnapsackOptimizeHS.
   //---
   int i;
   int status = 0;
   SOR_IntDblArrPtr ratio = NULL;
   ratio = SOR_IntDblArrNew(n, &status);

   if (status != OK) {
      return status;
   }

   for (i = 0; i < n; i++) {
      assert(w[i] >= 0);
      assert(p[i] >= 0);
      assert(!IS_TINY(w[i]));
      ratio->arr[i].x = p[i] / w[i];
      ratio->arr[i].i = i;
   }

   ratio->len = n;

   if (n > 1) {
      SOR_QSortIntDblDec(ratio->arr, 0, n - 1);
   }

   for (i = 0; i < n; i++) {
      psort[i] = p[ratio->arr[i].i];
      wsort[i] = w[ratio->arr[i].i];
      printf("i:%d j:%d p:%g w:%g\n",
             i, ratio->arr[i].i, psort[i], wsort[i]);
   }

   psort[n] = 0;
   wsort[n] = BIGM;
   SOR_IntDblArrFree(&ratio);
   return status;
}

//===========================================================================
void KnapsackSortRatioOut(const int       n,
                          const double*   p,
                          const double*   w,
                          double*         psort,
                          double*         wsort,
                          SOR_IntDbl*     ratio)
{
   //---
   //--- Sort non-increasing order of p[j]/w[j]. Return the new
   //--- order in psort, wsort which will be of size n+1 to conform
   //--- to the routine SOR_KnapsackOptimizeHS.
   //---
   //--- This version accepts a pointer to ratio and returns it after
   //--- the sort. We need this version if we need to construct the
   //--- solution at the end, since the KP solver's return is not in
   //--- the original order.
   //---
   int i;

   for (i = 0; i < n; i++) {
      assert(!IS_SMALL(w[i]));
      ratio[i].x = p[i] / w[i];
      ratio[i].i = i;
   }

   if (n > 1) {
      SOR_QSortIntDblDec(ratio, 0, n - 1);
   }

   for (i = 0; i < n; i++) {
      psort[i] = p[ratio[i].i];
      wsort[i] = w[ratio[i].i];
      printf("i:%d j:%d p:%g w:%g\n",
             i, ratio[i].i, psort[i], wsort[i]);
   }

   psort[n] = 0;
   wsort[n] = BIGM;
}

//===========================================================================
int KnapsackOptimizeHS(const int      n,
                       const double   c,
                       double*        p,
                       double*        w,
                       int*           x,
                       double*        z,
                       int*           pstatus)
{
   //---
   //--- Knapsack Problems (Algorithms and Computer Implementations)
   //--- Horowitz-Sanhi DFS Branch and Bound (Martello/Toth p30-31)
   //---
   //--- INPUT:
   //---  p and w should be sent in as size n+1
   //---  the input is ordered nonincreasing in p[j]/w[j]
   //---  x should be sent in calloc'd
   //---
   //--- OUTPUT:
   //---  x is the best solution (**in permuted order, NOT original order**)
   //---  z is the best solution value = sum{j = 0..n-1} p[j] x[j]
   //---
   //TODO: better return codes
   int i, j, r, k, iter;
   double zhat; // curr solution value = sum{j = 0..n-1} p[j] xhat[j]
   double chat; // curr residual capacity = c - sum{j = 0..n-1} w[j] xhat[j]
   double wSum, pSum, u;
   int* xhat = NULL;
   *pstatus = OK;
   DEBUG_KNAP(
      printf("\n" "//---in HS");

   for (i = 0; i < n; i++) {
   printf("\n" "p[%d]: %12.10f, w[%d]: %12.10f",
          i, p[i], i, w[i]);
   }
   printf("\n" "CAP: %12.10f", c);
   fflush(stdout);
   printf("\n" "//---out HS");
   );

   //---
   //--- we assume a[j] > 0, for all j, and x in {0,1}, so if b < 0, INF
   //---
   if (c < 0) {
      return KP_INFEASIBLE;
   }

   memset(x, 0, n * sizeof(int));

   if (n == 0) {
      return KP_OPTIMAL;
   }

   assert(n >= 0);

   //---
   //--- deal with some trivial cases here
   //---
   if (n == 1) {
      //---
      //--- if w[j] <= b, y = 1, psi = coverEl
      //--- else        , y = 0, psi = 0
      //---
      if ((w[0] - c ) > EPSILON) {
         *z   = 0;
         x[0] = 0;
      } else {
         *z   = p[0];
         x[1] = 1;
      }

      return KP_OPTIMAL;
   } else if (n == 2) {
      //---
      //--- Possible solutions: (0,0), (0,1), (1,0), (1,1).
      //--- Check each one, keep the best feasible one.
      //---
      //--- We have ??
      //--- already checked above that w[i] < cap, for each i.
      //---
      *z   = 0.0;
      x[0] = 0;
      x[1] = 0;

      if (((w[1] - c) <= EPSILON) && (p[1] > *z)) {
         *z   = p[1];
         x[0] = 0;
         x[1] = 1;
      }

      if (((w[0] - c) <= EPSILON) && (p[0] > *z)) {
         *z   = p[0];
         x[0] = 1;
         x[1] = 0;
      }

      if (((w[0] + w[1] - c) <= EPSILON) && ((p[0] + p[1]) > *z)) {
         *z = p[0] + p[1];
         x[0] = 1;
         x[1] = 1;
      }

      return KP_OPTIMAL;
   }

   xhat = (int*) calloc(n, sizeof(int));

   if (xhat == NULL) {
      *pstatus = ERR_NO_MEMORY;
      return KP_ERROR;
   }

   //---
   //--- (1) initialize
   //---
   *z    = 0;
   zhat  = 0;
   chat  = c;
   p[n]  = 0;
   w[n]  = BIGM;
   j     = 0;
   iter  = 0;

   while (1) {
      iter++;

      if (iter > HS_MAXITER) {
         return KP_ITERLIMIT;
      }

      //---
      //--- (2) compute the upper bound U1
      //---
      //---  r = arg min{i : sum{k = j..i} w[k] > chat}
      //---  u = sum{k = j..r-1} p[k]
      //---      + floor((chat - sum{k = j..r-1} w[k]) p[r]/w[r])
      //---  if(z >= zhat + u) then go to 5
      r    = j;
      wSum = w[r];
      pSum = p[r];

      //while((wSum < (chat - EPSILON)) && r <= n){
      //if wSum <= chat, continue
      //if wSum  > chat, stop
      while ((wSum <= chat) && r <= n) {
         r++;
         wSum += w[r];
         pSum += p[r];
         printf("r=%d, wSum=%12.10f, chat=%12.10f\n",
                r, wSum, chat);
      }

      assert(r >= 0 && r <= n);
      wSum -= w[r];
      pSum -= p[r];
      u = pSum + floor((chat - wSum) * p[r] / w[r]);
      DEBUG_KNAP(printf("\n" "z, %g, zhat: %g, u: %g r: %d ws: %g, ps: %g",
                        *z, zhat, u, r, wSum, pSum);
                );

      //if(*z + EPSILON <= zhat + u){
      if (*z <= zhat + u) {
         //if(!(*z >= zhat + u)){
         //---
         //--- (3) perform a forward step
         //---
         do {
            DEBUG_KNAP(printf("\n" "w[%d]: %9.6f, caphat: %9.6f",
                              j, w[j], chat);
                      );

            //while(w[j] - chat <= EPSILON){
            while (w[j] <= chat) {
               chat -= w[j];
               zhat += p[j];
               xhat[j] = 1;
               j++;
            }

            if (j < n) {
               xhat[j] = 0;
               j++;
            }
         } while (j == n - 1);

         //---
         //--- if j < n-1, then goto 2
         //---
         if (j < n - 1) {
            continue;
         }

         //---
         //--- (4) update the best solution so far
         //---
         DEBUG_KNAP(printf("\n" "update best sol zhat %g", zhat);
                   );

         //if(zhat > *z){
         if (zhat >= *z) {
            *z = zhat;

            for (k = 0; k < n; k++) {
               x[k] = xhat[k];
            }
         }

         j = n - 1;

         if (xhat[n - 1] == 1) {
            chat += w[n - 1];
            zhat -= p[n - 1];
            xhat[n - 1] = 0;
         }
      }

      DEBUG_KNAP(printf("\n" "backtrack");
                );
      //---
      //--- (5) backtrack (find i = max arg{k < j : x[k] = 1}
      //---
      i = std::max(0, j - 1); // need max, in case i = 0

      while ((xhat[i] != 1) && i > 0) {
         i--;
         DEBUG_KNAP(printf("\n" "i : %d xhat[i]: %d", i, xhat[i]);
                   );
      }

      if (i == 0) {
         break;
      }

      chat += w[i];
      zhat -= p[i];
      xhat[i] = 0;
      j = i + 1;
   }

   DEBUG_KNAP(
      wSum = 0.0;
      pSum = 0.0;

   for (k = 0; k < n; k++) {
   if (x[k] == 1) {
         wSum += w[k];
         pSum += p[k];
         printf("\n" "x[%d]=%d, w: %12.4f, p: %12.4f, ",
                k, x[k], w[k], p[k]);
         printf("wSum: %12.4f, pSum: %12.4f",
                wSum, pSum);
      }
   }
   );
   UTIL_DELPTR(xhat);
   return KP_OPTIMAL;
}

