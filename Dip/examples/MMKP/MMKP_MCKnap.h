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

#ifndef MMKP_MCKNAP_INCLUDED
#define MMKP_MCKNAP_INCLUDED

//---
//--- The MMKP MCKnap -- interface class for
//---     Pisinger's Multi-Choice Knapsack Solver
//---       http://www.diku.dk/~pisinger/codes.html
//---

//---
//--- Multi-Choice Knapsack Problem
//---  min  sum{i in 1..n, j in 1..l[i]}  c[i,j] x[i,j]
//---  s.t. sum{i in 1..n, j in 1..l[i]}  w[i,j] x[i,j] <= b
//---       sum{j in 1..l[i]}                      x[i,j]  = 1   , i in 1..n
//---


// --------------------------------------------------------------------- //
extern "C"{
#include "mcknap.h"
}

//extern solstruct optsol;

// --------------------------------------------------------------------- //
#include "UtilMacros.h"
using namespace std;

// --------------------------------------------------------------------- //
class MMKP_MCKnap{
 private:

   //---
   //--- typedef int           itype;   // item profits and weights 
   //--- typedef long          stype;   // sum of pofit or weight
   //--- typedef unsigned long vtype;   // solution vector
   //---
   //---
   //--- typedef struct {
   //---    itemset  *fset;
   //---    itemset  *lset;
   //---    ntype    size;
   //--- } isetset;
   //---
   isetset        * m_setset;        // set of itemsetsets
   double         * m_costDbl;       // temp storage
   int            * m_cost;          // cost vector (integer)
   int            * m_weight;        // weight vector (integer)
   int              m_nCols;         // size of mc-knap problem
   int              m_nGroupRows;
   int              m_nGroupCols;
   int              m_capacity;      // capacity (integer)
   int              m_cscale;        // cost vector scale factor
   int              m_wscale;        // weight vector scale factor

 public:
   //separate function just to reset cost...from one call to next
   void solveTrivialMaxSum(const double * redCost,
			   const double * origCost,
			   vector<int>  & solInd,
			   double       & varRedCost,
			   double       & varOrigCost);   

   inline const int getIndexIJ(const int i,
                               const int j) const{
      return (i * m_nGroupCols) + j;
   }

   inline pair<int,int> getIndexInv(const int index) const {      
      return make_pair(index / m_nGroupCols, index % m_nGroupCols);
   }

   void setMCKnapData(const double   capacity,
		      const double * weight);

   void solveMCKnap(const double   * redCost,
		    const double   * origCost,
		    vector<int>    & solInd,
		    vector<double> & solEls,
		    double         & varRedCost,
		    double         & varOrigCost);

 public:
   MMKP_MCKnap(const int nGroupRows,
	       const int nGroupCols) :
      m_setset   (NULL),
      m_costDbl  (NULL),
      m_cost     (NULL),
      m_weight   (NULL),
      m_nCols    (nGroupRows * nGroupCols),
      m_nGroupRows(nGroupRows),
      m_nGroupCols(nGroupCols),
      m_capacity (0),
      m_cscale   (1),
      m_wscale   (1)      
      {
	 m_costDbl = new double[m_nCols];
	 m_cost    = new int[m_nCols];
	 m_weight  = new int[m_nCols];
	 assert(m_costDbl && m_cost && m_weight);		
      }
   ~MMKP_MCKnap() {      
      if(m_setset){         
         itemset * setPtr = m_setset->fset;
         if(setPtr){
            int       i;
            for(i = 0; i < m_setset->size; i++){
               if(setPtr->fset)
                  free(setPtr->fset);
               setPtr++;
            }
            free(m_setset->fset);
         }
         free(m_setset);
      }
      UTIL_DELARR(m_costDbl);
      UTIL_DELARR(m_cost);
      UTIL_DELARR(m_weight);
   }
};


#endif
