//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)       //
//                                                                           //
// Copyright (C) 2002-2007, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

#include "DecompVar.h"

// --------------------------------------------------------------------- //
void DecompVar::fillDenseArr(int      len,
                             double* arr)
{
   CoinFillN(arr, len, 0.0);
   const int      sz    = m_s.getNumElements();
   const int*     inds  = m_s.getIndices();
   const double* elems = m_s.getElements();

   for (int i = 0; i < sz; ++i) {
      arr[inds[i]] = elems[i];
   }
}

// --------------------------------------------------------------------- //
void
DecompVar::print(ostream*    os,
                 DecompApp* app) const
{
   double lb = getLowerBound();
   double ub = getUpperBound();
   (*os) << "\nVAR c: " << m_origCost
         << " rc: "     << m_redCost
         << " eff: "    << m_effCnt;

   if (lb > -DecompInf) {
      (*os) << " lb:  "    << getLowerBound();
   } else {
      (*os) << " lb: -INF";
   }

   if (ub < DecompInf) {
      (*os) << " ub:  "    << getUpperBound();
   } else {
      (*os) << " ub:  INF";
   }

   (*os) << "\n";
   UtilPrintPackedVector(m_s, os, app);
}
