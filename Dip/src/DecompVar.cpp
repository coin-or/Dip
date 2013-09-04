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

#include "DecompVar.h"
#include "DecompModel.h"

using namespace std;

// --------------------------------------------------------------------- //
//Design question - we use a check over the entire space to simplify
// storage of branching bounds. But, in reality, for node 1, for example,
// these has only been one branch - so we only need to check one number.
//This could be much faster if we are willing to do more accounting related
// to branching.
bool DecompVar::doesSatisfyBounds(int                     denseLen,
                                  double*                 denseArr,
                                  const DecompAlgoModel& model,
                                  const double*           lbs,
                                  const double*           ubs)
{
   int            j;
   double         xj;//, lb, ub;
   vector<int> ::const_iterator it;
   map<int, int>::const_iterator mcit;
   DecompConstraintSet* modelRelax    = model.getModel();
   const vector<int>&   activeColumns  = modelRelax->getActiveColumns();
   //---
   //--- activeColumns are in original space
   //---    denseArr, lbs, ubs are all in original space
   //---
   fillDenseArr(denseLen, denseArr);//TODO: expensive...

   for (it = activeColumns.begin(); it != activeColumns.end(); it++) {
      j  = *it;
      xj = denseArr[j];

      if (xj < (lbs[j] - DecompEpsilon) ||
            xj > (ubs[j] + DecompEpsilon)) {
         return false;
      }
   }

   return true;
}

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
   (*os) << "\nVAR c: "   << m_origCost
         << " rc: "       << m_redCost
         << " eff: "      << m_effCnt
         << " block: "    << m_blockId
         << " colIndex: " << m_colMasterIndex;

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

// --------------------------------------------------------------------- //
void
DecompVar::print(ostream*               os,
                 const vector<string>& colNames,
                 const double*          value) const
{
   double lb = getLowerBound();
   double ub = getUpperBound();
   (*os) << "\nVAR c: "   << m_origCost
         << " rc: "       << m_redCost
         << " eff: "      << m_effCnt
         << " block: "    << m_blockId
         << " colIndex: " << m_colMasterIndex;

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
   UtilPrintPackedVector(m_s, os, colNames, value);
}
