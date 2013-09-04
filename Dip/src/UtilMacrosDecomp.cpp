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
#include "UtilMacrosDecomp.h"
#include "DecompApp.h"
#define SPACES " \t\r\n"
//===========================================================================//

using namespace std;

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

// ------------------------------------------------------------------------- //
void UtilPrintPackedVector(const CoinPackedVector& v,
                           ostream*                 os,
                           const vector<string>&    colNames,
                           const double*            value)
{
   (*os).precision(2);
   const int*     inds  = v.getIndices();
   const double* elems = v.getElements();
   const int      len   = v.getNumElements();
   int    namesSize = static_cast<int>(colNames.size());
   double sum       = 0.0;

   for (int i = 0; i < len; i++) {
      if (!namesSize)
         (*os) << setw(10) << UtilDblToStr(elems[i], 4)
               << " x[" << setw(6) <<  inds[i] << "]  ";
      else {
         (*os) << setw(10) << UtilDblToStr(elems[i], 4)
               << " " << setw(10) << colNames[inds[i]] << " ";
      }

      if (value) {
         sum += elems[i] * value[inds[i]];
         (*os) << " --> " << setw(10) << UtilDblToStr(value[inds[i]], 4);
         (*os) << " --> " << setw(10) << UtilDblToStr(sum, 4);
      }

      (*os) << "\n";
   }

   if (value) {
      (*os) << "dot product = " << UtilDblToStr(sum, 4) << endl;
   }

   (*os) << endl;
}


