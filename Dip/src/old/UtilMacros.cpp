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

//===========================================================================//
#include "UtilMacros.h"
#include "DecompApp.h"
#define SPACES " \t\r\n"
//===========================================================================//

// =========================================================================
// Graph Macros
// =========================================================================
// ------------------------------------------------------------------------- //
pair<int, int> UtilBothEndsU(const int index)
{
   int i = (int)(floor(sqrt(1 + 8.0 * index) / 2 + .500000001));
   return make_pair( i, index - (i * (i - 1) / 2) );
}

// =========================================================================
// COIN Macros
// =========================================================================
// ------------------------------------------------------------------------- //
CoinPackedVector* UtilPackedVectorFromDense(const int      len,
      const double* dense,
      const double   etol)
{
   //TODO: test for dup? - efficiency vs debug
   //TODO: insert is slow - better to use setVector?
   int i;
   CoinPackedVector* v = new CoinPackedVector();

   for (i = 0; i < len; i++) {
      if (fabs(dense[i]) > etol) {
         v->insert(i, dense[i]);
      }
   }

   return v;
}

// ------------------------------------------------------------------------- //
void UtilPackedVectorFromDense(const int          len,
                               const double*      dense,
                               const double       etol,
                               CoinPackedVector& v)
{
   //TODO: test for dup? - efficiency vs debug
   //TODO: insert is slow - better to use setVector?
   int i;

   for (i = 0; i < len; i++) {
      if (fabs(dense[i]) > etol) {
         v.insert(i, dense[i]);
      }
   }
}

// ------------------------------------------------------------------------- //
void UtilPrintPackedVector(const CoinPackedVector& v,
                           ostream*                 os,
                           DecompApp*               app)
{
   (*os).precision(2);
   const int*     inds  = v.getIndices();
   const double* elems = v.getElements();
   const int      len   = v.getNumElements();

   for (int i = 0; i < len; i++) {
      if (!app) {
         (*os) << elems[i] << " x[" << inds[i] << "]  ";
      } else {
         (*os) << elems[i] << " ";
         app->printOriginalColumn(inds[i], os);
         (*os) << "  ";
      }

      if ((i + 1) % 5 == 0) {
         (*os) << "\n";
      }
   }

   (*os) << endl;
}


// =========================================================================
// Other Macros
// =========================================================================

//TODO: look into CUT_scaleRow
//combine these two!
/*==========================================================================*/
int UtilScaleDblToIntArr(const int      arrLen,
                         const double* arrDbl,
                         int*           arrInt,
                         const double   oneDbl,
                         int*           oneInt,
                         const double   epstol)
{
   //---
   //--- A very simple function to scale an array of doubles to integers.
   //---    Note: epstol denotes the preferred accuracy,
   //---    so, we will scale by 1.0/epstol, unless something smaller works.
   //--- It constructs the scaled array and returns the scale factor.
   //---
   //--- It can also scale oneDbl to oneInt wrt to the array (e.g..,
   //--- the rhs of row). If oneInt == NULL, then this part is skipped.
   //---
   int      i, scaleFactor = 1, n_aFrac = 0, factorToBig = 0;
   double* arrFrac = NULL;
   double   fractionalPart;
   double   oneOverEps = 1.0 / epstol;
   //TODO: pass in arrFrac?
   arrFrac = new double[arrLen + 1];
   CoinAssertHint(arrFrac, "Error: Out of Memory");

   for (i = 0; i < arrLen; i++) {
      fractionalPart = UtilFracPart(arrDbl[i]);

      if (!UtilIsZero(fractionalPart)) {
         fractionalPart     *= oneOverEps;
         arrFrac[n_aFrac++]  = (int)fractionalPart * (double)epstol;
      }
   }

   if (oneInt) {
      fractionalPart = UtilFracPart(oneDbl);

      if (!UtilIsZero(fractionalPart)) {
         fractionalPart     *= oneOverEps;
         arrFrac[n_aFrac++]  = (int)fractionalPart * (double)epstol;
      }
   }

   for (i = 0; i < n_aFrac; i++) {
      CoinAssertDebug(arrFrac[i] < (INT_MAX / scaleFactor));
      arrFrac[i] *= scaleFactor;

      while (!UtilIsZero(UtilFracPart(arrFrac[i]))) {
         scaleFactor *= 10;

         if (scaleFactor >= oneOverEps) {
            factorToBig = 1;
            break;
         }

         CoinAssertDebug(arrFrac[i] < (INT_MAX / 10));
         arrFrac[i]    *= 10;
         CoinAssertDebug(arrFrac[i] >= 0);
      }

      if (factorToBig) {
         break;
      }
   }

   for (i = 0; i < arrLen; i++) {
      arrInt[i] = (int)(arrDbl[i] * scaleFactor);
   }

   if (oneInt) {
      *oneInt = (int)(oneDbl * scaleFactor);
   }

   UTIL_DELARR(arrFrac);
   return scaleFactor;
}


/*==========================================================================*/
int UtilScaleDblToIntArr(const int      arrLen,
                         const double* arrDbl,
                         int*           arrInt,
                         const double   epstol)
{
   //---
   //--- A very simple function to scale an array of doubles to integers.
   //---    Note: epstol denotes the preferred accuracy,
   //---    so, we will scale by 1.0/epstol, unless something smaller works.
   //--- It constructs the scaled array and returns the scale factor.
   //---
   //--- It can also scale oneDbl to oneInt wrt to the array (e.g..,
   //--- the rhs of row). If oneInt == NULL, then this part is skipped.
   //---
   int      i, scaleFactor = 1, n_aFrac = 0, factorToBig = 0;
   double* arrFrac = NULL;
   double   fractionalPart;
   double   oneOverEps = 1.0 / epstol;
   //TODO: pass in arrFrac?
   arrFrac = new double[arrLen];
   CoinAssertHint(arrFrac, "Error: Out of Memory");

   for (i = 0; i < arrLen; i++) {
      fractionalPart = UtilFracPart(arrDbl[i]);

      if (!UtilIsZero(fractionalPart)) {
         fractionalPart     *= oneOverEps;
         arrFrac[n_aFrac++]  = (int)fractionalPart * (double)epstol;
      }
   }

   for (i = 0; i < n_aFrac; i++) {
      CoinAssertDebug(arrFrac[i] < (INT_MAX / scaleFactor));
      arrFrac[i] *= scaleFactor;

      while (!UtilIsZero(UtilFracPart(arrFrac[i]))) {
         scaleFactor *= 10;

         if (scaleFactor >= oneOverEps) {
            factorToBig = 1;
            break;
         }

         CoinAssertDebug(arrFrac[i] < (INT_MAX / 10));
         arrFrac[i]    *= 10;
         CoinAssertDebug(arrFrac[i] >= 0);
      }

      if (factorToBig) {
         break;
      }
   }

   for (i = 0; i < arrLen; i++) {
      arrInt[i] = (int)(arrDbl[i] * scaleFactor);
   }

   UTIL_DELARR(arrFrac);
   return scaleFactor;
}

