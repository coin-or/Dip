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
// Copyright (C) 2002-2009, Lehigh University, Matthew Galati, Ted Ralphs    //
// All Rights Reserved.                                                      //
//===========================================================================//

#include "DecompVar.h"

// --------------------------------------------------------------------- //
//this is not good enough - this only checks the nonzero elements
//   what if we branch to fix a component to 1, need to make sure that
//   is in this set - will have to compare vs dense array
bool DecompVar::doesSatisfyBounds(int            denseLen,
				  double       * denseArr,
				  const double * lbs,
				  const double * ubs){

   int            j;
   double         xj;
   fillDenseArr(denseLen, denseArr);
   for (j = 0; j < denseLen; ++j){
      xj = denseArr[j];
      //printf("j:%d xj:%g lb:%g ub:%g --> ",
      //     j, xj, lbs[j], ubs[j]);
      if(xj < (lbs[j] - DecompEpsilon) ||
	 xj > (ubs[j] + DecompEpsilon)){
	//printf(" false\n");
	return false;	 
      }
      else{
	//printf("\n");
      }
   }
   return true;
}

// --------------------------------------------------------------------- //
void DecompVar::fillDenseArr(int      len,
                             double * arr){
   CoinFillN(arr, len, 0.0);
   const int      sz    = m_s.getNumElements();
   const int    * inds  = m_s.getIndices();
   const double * elems = m_s.getElements();
   for (int i = 0; i < sz; ++i)
      arr[inds[i]] = elems[i];
}

// --------------------------------------------------------------------- //
void 
DecompVar::print(ostream   * os,
                 DecompApp * app) const {

   double lb = getLowerBound();
   double ub = getUpperBound();
   (*os) << "\nVAR c: "   << m_origCost 
         << " rc: "       << m_redCost 
         << " eff: "      << m_effCnt
	 << " block: "    << m_blockId
         << " colIndex: " << m_colMasterIndex;
   if(lb > -DecompInf)
      (*os) << " lb:  "    << getLowerBound();
   else
      (*os) << " lb: -INF";
   if(ub < DecompInf)
      (*os) << " ub:  "    << getUpperBound();
   else
      (*os) << " ub:  INF";
   (*os) << "\n";   
   UtilPrintPackedVector(m_s, os, app);
}

// --------------------------------------------------------------------- //
void 
DecompVar::print(ostream              * os,
		 const vector<string> & colNames,
		 const double         * value) const {

   double lb = getLowerBound();
   double ub = getUpperBound();
   (*os) << "\nVAR c: "   << m_origCost 
         << " rc: "       << m_redCost 
         << " eff: "      << m_effCnt
	 << " block: "    << m_blockId
         << " colIndex: " << m_colMasterIndex;
   if(lb > -DecompInf)
      (*os) << " lb:  "    << getLowerBound();
   else
      (*os) << " lb: -INF";
   if(ub < DecompInf)
      (*os) << " ub:  "    << getUpperBound();
   else
      (*os) << " ub:  INF";
   (*os) << "\n";   
   UtilPrintPackedVector(m_s, os, colNames, value);
}
