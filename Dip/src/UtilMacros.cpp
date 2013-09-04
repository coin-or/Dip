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
#include "UtilMacros.h"
#define SPACES " \t\r\n"
//===========================================================================//

// =========================================================================
// Graph Macros
// =========================================================================
// ------------------------------------------------------------------------- //
std::pair<int, int> UtilBothEndsU(const int index)
{
   int i = (int)(floor(sqrt(1 + 8.0 * index) / 2 + .500000001));
   return std::make_pair( i, index - (i * (i - 1) / 2) );
}

// =========================================================================
// Other Macros
// =========================================================================

// ------------------------------------------------------------------------- //
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
   int      i, scaleFactor = 1, n_aFrac = 0, factorTooBig = 0;
   double* arrFrac = NULL;
   double   fractionalPart;
   double   oneOverEps = 1.0 / epstol;
   arrFrac = new double[arrLen + 1];
   assert(arrFrac);

   for (i = 0; i < arrLen; i++) {
      fractionalPart = UtilFracPart(arrDbl[i]);

      if (!UtilIsZero(fractionalPart)) {
         fractionalPart     *= oneOverEps;
         arrFrac[n_aFrac++]  = static_cast<int>(round(fractionalPart))
                               * (double)epstol;
      }
   }

   if (oneInt) {
      fractionalPart = UtilFracPart(oneDbl);

      if (!UtilIsZero(fractionalPart)) {
         fractionalPart     *= oneOverEps;
         arrFrac[n_aFrac++]  = static_cast<int>(round(fractionalPart))
                               * (double)epstol;
      }
   }

   for (i = 0; i < n_aFrac; i++) {
      assert(arrFrac[i] < (INT_MAX / scaleFactor));
      arrFrac[i] *= scaleFactor;

      while (!UtilIsZero(UtilFracPart(arrFrac[i]))) {
         scaleFactor *= 10;

         if (scaleFactor >= oneOverEps) {
            factorTooBig = 1;
            break;
         }

         assert(arrFrac[i] < (INT_MAX / 10));
         arrFrac[i]    *= 10;
         assert(arrFrac[i] >= 0);
      }

      if (factorTooBig) {
         break;
      }
   }

   for (i = 0; i < arrLen; i++) {
      arrInt[i] = static_cast<int>(round(arrDbl[i] * scaleFactor));
   }

   if (oneInt) {
      *oneInt = static_cast<int>(round(oneDbl * scaleFactor));
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
   int      i, scaleFactor = 1, n_aFrac = 0, factorTooBig = 0;
   double* arrFrac = NULL;
   double   fractionalPart;
   double   oneOverEps = 1.0 / epstol;
   //TODO: pass in arrFrac?
   arrFrac = new double[arrLen];
   assert(arrFrac);

   for (i = 0; i < arrLen; i++) {
      fractionalPart = UtilFracPart(arrDbl[i]);

      //printf("arrDbl[%d]=%10.5f fracPart=%6.5f\n",
      //     i, arrDbl[i], fractionalPart);
      if (!UtilIsZero(fractionalPart)) {
         fractionalPart     *= oneOverEps;
         //printf("fracPart is not zero oneOverEps= %10.5f fracPart= %10.5f\n",
         //oneOverEps, fractionalPart);
         arrFrac[n_aFrac++]  = static_cast<int>(round(fractionalPart))
                               * (double)epstol;
         //printf("arrFrac[%d] = %10.5f\n", (n_aFrac-1), arrFrac[n_aFrac-1]);
      }
   }

   for (i = 0; i < n_aFrac; i++) {
      assert(arrFrac[i] < (INT_MAX / scaleFactor));
      arrFrac[i] *= scaleFactor;

      while (!UtilIsZero(UtilFracPart(arrFrac[i]))) {
         scaleFactor *= 10;

         if (scaleFactor >= oneOverEps) {
            factorTooBig = 1;
            break;
         }

         assert(arrFrac[i] < (INT_MAX / 10));
         arrFrac[i]    *= 10;
         assert(arrFrac[i] >= 0);
      }

      if (factorTooBig) {
         break;
      }
   }

   //---
   //--- must be careful not to trunc here
   //---   so, we want to round
   //---
   for (i = 0; i < arrLen; i++) {
      arrInt[i] = static_cast<int>(round(arrDbl[i] * scaleFactor));
   }

   UTIL_DELARR(arrFrac);
   return scaleFactor;
}

#if 0
/*==========================================================================*/
//---
//--- taken from Concorde's integerize_vector in LOCALCUT/first.c
//---
/*==========================================================================*/
#define INTEGERIZE_MUL (16*9*5*7*11*13*17)
#define CC_SWAP(a,b,t) (((t)=(a)),((a)=(b)),((b)=(t)))

/*==========================================================================*/
static int CCutil_our_gcd (int a, int b)
{
   int c;

   if (a < 0) {
      a = -a;
   }

   if (b < 0) {
      b = -b;
   }

   if (a > b) {
      CC_SWAP (a, b, c);
   }

   while (a) {
      c = b % a;
      b = a;
      a = c;
   }

   return b;
}
/*==========================================================================*/
static void integerize_check (int count, double* dvec, int* ivec)
{
   double error;
   double scale;
   int i;
   scale = 1.0;

   for (i = 0; i < count; i++) {
      if (ivec[i]) {
         scale = fabs(dvec[i] / ivec[i]);
         break;
      }
   }

   error = 0.0;

   for (i = 0; i < count; i++) {
      error += fabs(dvec[i] / scale - ivec[i]);
   }

   if (error > 0.001) {
      printf ("bad integerize error %f\n", error);
      printf ("from:");

      for (i = 0; i < count; i++) {
         printf (" %f", dvec[i]);
      }

      printf ("\n");
      printf ("to:");

      for (i = 0; i < count; i++) {
         printf (" %d", ivec[i]);
      }

      printf ("\n");
      fflush (stdout);
   }
}

/*==========================================================================*/
static void cfrac_approx (double d, int* p_num, int* p_den)
{
   int sgn = 1;
   int a;
   int g2 = 0;
   int g1 = 1;
   int g = 0;
   int h2 = 1;
   int h1 = 0;
   int h = 1;
   int bestg = 0;
   int besth = 1;
   double besterr = 1.0;
   int lim;

   if (d < 0.0) {
      d = -d;
      sgn = -1;
   }

   for (;;) {
      a = (int) d;

      if (a) {
         lim = INT_MAX / a;
      } else {
         lim = INT_MAX;
      }

      if (g1 > lim || h1 > lim ||
            g2 > INT_MAX - g1 * a || h2 > INT_MAX - h1 * a) {
         break;
      }

      g = a * g1 + g2;
      h = a * h1 + h2;
      g2 = g1;
      h2 = h1;
      g1 = g;
      h1 = h;

      if (d - (double) a < besterr) {
         bestg = g;
         besth = h;
         besterr = d - (double) a;

         if (besterr < 0.0001) {
            break;
         }
      }

      d = 1 / (d - a);
   }

   *p_num = bestg * sgn;
   *p_den = besth;
}

/*==========================================================================*/
static void integerize_vector2 (int count, double* dvec, int* ivec)
{
   int i;
   double scale;
   int divv;
#if 0
   printf ("integerize(2):");

   for (i = 0; i < count; i++) {
      printf (" %f", dvec[i]);
   }

   printf ("\n");
#endif
   scale = fabs(dvec[0]);

   for (i = 1; i < count; i++) {
      if (fabs(dvec[i]) > scale) {
         scale = fabs(dvec[i]);
      }
   }

   if (scale == 0.0) {
      for (i = 0; i < count; i++) {
         ivec[i] = 0;
      }

      return;
   }

   for (i = 0; i < count; i++) {
      if (dvec[i] < 0.0) {
         ivec[i] = -((int) (-dvec[i] * INTEGERIZE_MUL / scale + 0.5));
      } else {
         ivec[i] = ((int) (dvec[i] * INTEGERIZE_MUL / scale + 0.5));
      }
   }

   divv = ivec[0];

   for (i = 1; i < count; i++) {
      divv = CCutil_our_gcd (divv, ivec[i]);
   }

   if (divv > 1) {
      for (i = 0; i < count; i++) {
         ivec[i] /= divv;
      }
   }

#if 0
   printf ("integerized");

   for (i = 0; i < count; i++) {
      printf (" %d", ivec[i]);
   }

   printf ("\n");
   integerize_check (count, dvec, ivec);
#endif
}

/*==========================================================================*/
static void integerize_vector (int count, double* dvec, int* ivec)
{
   int i;
   int num;
   int den;
   int g;
   int scale = 1;
   double dmax;
#if 0
   printf ("integerize:");

   for (i = 0; i < count; i++) {
      printf (" %f", dvec[i]);
   }

   printf ("\n");
#endif
   dmax = 0.0;

   for (i = 0; i < count; i++) {
      if (fabs(dvec[i]) > dmax) {
         dmax = fabs(dvec[i]);
      }
   }

   dmax = 1.0 / dmax;

   for (i = 0; i < count; i++) {
      cfrac_approx (dvec[i]*dmax, &num, &den);
      g = CCutil_our_gcd (den, scale);
      scale /= g;

      if (INT_MAX / abs(scale) < abs(den)) {
         goto FAILURE;
      }

      scale *= den;
   }

   for (i = 0; i < count; i++) {
      cfrac_approx (dvec[i]*dmax, &num, &den);
      den = scale / den;

      if (INTEGERIZE_MUL / abs(den) < abs(num)) {
         goto FAILURE;
      }

      ivec[i] = num * den;
   }

#if 0
   printf ("integerized");

   for (i = 0; i < count; i++) {
      printf (" %d", ivec[i]);
   }

   printf ("\n");
   integerize_check (count, dvec, ivec);
#endif
   return;
FAILURE:
   printf ("Overflow in integerize_vector.  Calling integerize_vector2\n");
   integerize_vector2 (count, dvec, ivec);
}



/*==========================================================================*/
//---
//--- taken from concorde integerize_vector in LOCALCUT/first.c
//---
void UtilScaleDblToIntArrCC(const int      arrLen,
                            double* arrDbl,
                            int*           arrInt)
{
   integerize_vector(arrLen, arrDbl, arrInt);
}
#endif
