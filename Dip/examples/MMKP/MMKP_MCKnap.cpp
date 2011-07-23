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
#include "MMKP_MCKnap.h"
#include "UtilMacrosDecomp.h"
extern "C"{
#include "mcknap.h"
}

// --------------------------------------------------------------------- //
//#define MCKP_EPSILON      1.0e-4 //still causing overflow
#define MCKP_EPSILON      1.0e-3
//#define MMKP_MCKNAP_DEBUG

// --------------------------------------------------------------------- //
#include "UtilMacros.h"

// --------------------------------------------------------------------- //
void MMKP_MCKnap::solveTrivialMaxSum(const double * redCost,
				     const double * origCost,
				     vector<int>  & solInd,
				     double       & varRedCost,
				     double       & varOrigCost){
   
   double minRedCost;
   int    i, j, minRedCostInd, ijIndex, minWeight, totalWeight;

   //---
   //--- Pisinger's code breaks on this trivial case.
   //--- 
   //---  In the case where maxwsum <= c, then we can trivially
   //---  pick the max profit / min cost element from each group
   //---   
   //--- We have to be careful here because there might be ties to deal
   //---   with. The max profit might be p* for two different choices
   //---   but only one of those gave a weight that was under capacity.
   //--- So, we have to find the alternative choice with lowest weight.
   //---
   //--- NOTE: this algorithm was added because mcknap algorithm seems
   //---  to crash on this trivial case. It would be better if we could
   //---  get this case fixed. TODO: send test case to Pisinger.
   //---
   varRedCost  = 0.0;
   varOrigCost = 0.0;
   solInd.reserve(m_nGroupRows);
   
   ijIndex     = 0;
   totalWeight = 0;
   for(i = 0; i < m_nGroupRows; i++){
      ijIndex       = getIndexIJ(i, 0);
      minRedCostInd = ijIndex;
      minRedCost    = redCost[ijIndex];
      minWeight     = m_weight[ijIndex];//need if we have ties
#ifdef MMKP_MCKNAP_DEBUG
      printf("i:%d j:%d redCost:%g minRedCost:%g wt:%d minWeight:%d\n",
           i, 0, redCost[ijIndex], minRedCost, m_weight[ijIndex],
           minWeight);
#endif

      ijIndex++;
      for(j = 1; j < m_nGroupCols; j++){
	if((redCost[ijIndex] - minRedCost) < -MCKP_EPSILON){
	  minRedCost    = redCost[ijIndex];
	  minWeight     = m_weight[ijIndex];
	  minRedCostInd = ijIndex;
	}else if( UtilIsZero(redCost[ijIndex] - minRedCost, MCKP_EPSILON) ){
	  //---
	  //--- break ties with element with least weight
	  //---
	  if(minWeight > m_weight[ijIndex]){
	    minWeight     = m_weight[ijIndex];
	    minRedCostInd = ijIndex;
	  }
	}
#ifdef MMKP_MCKNAP_DEBUG
	printf("i:%d j:%d redCost:%g minRedCost:%g wt:%d minWeight:%d\n",
	     i, j, redCost[ijIndex], minRedCost, m_weight[ijIndex],
	     minWeight);
#endif

	ijIndex++;
      }      
      assert((minRedCostInd >= 0) && (minRedCostInd < m_nCols));      
      totalWeight += m_weight[minRedCostInd];
      solInd.push_back(minRedCostInd);
#ifdef MMKP_MCKNAP_DEBUG
      printf("i:%d totalWeight:%d cap:%d mincostInd:%d redCost:%g origCost:%g\n",
             i, totalWeight, m_capacity, minRedCostInd, minRedCost, 
           origCost[minRedCostInd]);fflush(stdout);
#endif
      assert(totalWeight <= m_capacity);
      varRedCost  += minRedCost;
      varOrigCost += origCost[minRedCostInd];
   }
}

// --------------------------------------------------------------------- //
void MMKP_MCKnap::setMCKnapData(const double   capacity,
                                const double * weight){

   int       i, j, colIndex;      
   itemset * setPtr    = NULL;
   itemrec * recPtr    = NULL;
      
   //TODO: allow pass in arrFrac to this util function, for speed
   //TODO: really need to use malloc here?
      
   //---
   //--- Pisinger's code assume cost/weight is integer type.
   //---   So, we need to scale the cost/weight to nearest integer first.
   //---
   m_wscale = UtilScaleDblToIntArr(m_nCols,
                                   weight, m_weight,
                                   capacity, &m_capacity, MCKP_EPSILON);
   
   m_setset = (isetset*) malloc(sizeof(isetset));
   assert(m_setset);
      
   m_setset->size  = m_nGroupRows;
   m_setset->fset = (itemset*) malloc(m_setset->size * sizeof(itemset));
   assert(m_setset->fset);
   
   setPtr = m_setset->fset;
   colIndex = 0;
   for(i = 0; i < m_setset->size; i++){
      setPtr->size = m_nGroupCols;
      setPtr->fset = (itemrec*) malloc(setPtr->size * sizeof(itemrec));
      assert(m_setset->fset);

      recPtr       = setPtr->fset;
      for(j = 0; j < setPtr->size; j++){
	 recPtr->i    = i;
	 recPtr->j    = j;
         recPtr->wsum = m_weight[colIndex];
         recPtr->psum = m_cost[colIndex];
         recPtr++;
         colIndex++;
      }
      setPtr->lset = setPtr->fset + setPtr->size - 1;
      setPtr++;
   }
   m_setset->lset = m_setset->fset + m_setset->size - 1;
}

// --------------------------------------------------------------------- //
void MMKP_MCKnap::solveMCKnap(const double   * redCost,
			      const double   * origCost,
			      vector<int>    & solInd,
			      vector<double> & solEls,
			      double         & varRedCost,
			      double         & varOrigCost){
   
   int i, j, colIndex;

   //---
   //--- Pisinger's code (mcknap) solves a max problem. And, it assumes 
   //--- all positive costs/profits and weights.
   //---   
   memcpy(m_costDbl, redCost, m_nCols * sizeof(double));
#ifdef MMKP_MCKNAP_DEBUG
   for(i = 0; i < m_nCols; i++){
      pair<int,int> ij = getIndexInv(i);
      printf("\ncostDbl[%d: %d, %d]: %g",
	     i, ij.first, ij.second, m_costDbl[i]);
   }
#endif

   
   //---
   //--- flip reduced costs (max c == min -c)
   //---
   UtilNegateArr(m_nCols, m_costDbl);

   //---
   //--- add a constant so that all vertex weights are positive, inc alpha
   //---
   double offset = 0.0;
   double minrc  = *min_element(m_costDbl, m_costDbl + m_nCols);
#ifdef MMKP_MCKNAP_DEBUG
   printf("\nminrc   = %g", minrc);
#endif   
   if(minrc <= 0){
      offset = -minrc + 1;
      UtilAddOffsetArr(m_nCols, offset, m_costDbl);
   }

   //---
   //--- now scale the double array to an integer array
   //---
   //TODO: magic number - have to be careful of overflow... 
   m_cscale = UtilScaleDblToIntArr(m_nCols, m_costDbl, m_cost, MCKP_EPSILON);

#ifdef MMKP_MCKNAP_DEBUG
   double diff;
   printf("\noffset   = %g", offset);
   printf("\nm_cscale = %d", m_cscale);
   printf("\nm_wscale = %d", m_wscale);
   printf("\ncapacity = %d", m_capacity);
   for(i = 0; i < m_nCols; i++){
      pair<int,int> ij = getIndexInv(i);
      diff = fabs((m_costDbl[i]*m_cscale) - m_cost[i]);
      printf("\n[%d: %d, %d]: dbl-> %12.5f int-> %8d diff-> %12.5f",
	     i, ij.first, ij.second, m_costDbl[i], m_cost[i], diff);
      assert( diff < 0.99 );
   }
#endif

   //---
   //--- sanity check
   //---  if any cost value becomes negative that
   //---  denotes an overflow happened
   //--- TODO: not sure how to do this scaling safely and
   //---  accurately
   //---
   for(i = 0; i < m_nCols; i++){
      if(m_cost[i] < 0){
         throw UtilException("negative cost value",
                             "solveMCKnap", "MMKP_MCKnap");
      }
   }
      
   //---
   //--- setup the data structures for mcknap
   //---
   itemset * setPtr    = m_setset->fset;
   itemrec * recPtr    = NULL;
   itemrec * recSolPtr = NULL;
   
   //THINK: reset - assume memory is still there
   m_setset->size  = m_nGroupRows;
   setPtr = m_setset->fset;
   for(i = 0; i < m_setset->size; i++){
      setPtr->size = m_nGroupCols;
      recPtr       = setPtr->fset;    
      setPtr->lset = setPtr->fset + setPtr->size - 1;
      setPtr++;
   }
   m_setset->lset = m_setset->fset + m_setset->size - 1;
   

   colIndex            = 0;
   setPtr = m_setset->fset;
   for(i = 0; i < m_setset->size; i++){
      recPtr    = setPtr->fset;
      for(j = 0; j < setPtr->size; j++){
	 recPtr->i = i;
	 recPtr->j = j;
         recPtr->psum = m_cost[colIndex];
         recPtr->wsum = m_weight[colIndex];
#ifdef MMKP_MCKNAP_DEBUG
         printf("\ncolIndex: %d i: %d, j: %d, p: %d, w: %d",
		colIndex, i, j, recPtr->psum, recPtr->wsum);
#endif
         recPtr++;
         colIndex++;
      }
      setPtr++;
   }

   itemrec * solRec = new itemrec[m_setset->size];
   
   double minObj = 99999;
//   long minObj = 99999;
   int status  = minmcknapSolve(m_capacity, m_setset, solRec, &minObj);

   solInd.reserve(m_nGroupRows);
   solEls.reserve(m_nGroupRows);
   UtilFillN<double>(solEls, m_nGroupRows, 1.0);
   varRedCost  = 0.0;
   varOrigCost = 0.0;

   //---
   //--- TODO:
   //--- this is painful to get optimal assignments
   //---    wrote Dr. Pisinger for help (7/4/07)
   //--- NOTE: optsol is NOT reentrant
   //---
   //CoinAssert(optsol.size == 1); //TODO
#ifdef MMKP_MCKNAP_DEBUG
   printf("\nstatus=%d minObj=%g\n", status, minObj);
#endif

   switch(status){
   case MCKNAP_RC_INF:
      assert(status != MCKNAP_RC_INF);
      break;
   case MCKNAP_RC_OK:
      {
	 //---
	 //--- need to unravel:
	 //---    s * ((-x) + offset)
	 //---
	  
	 double solObj = minObj / static_cast<double>(m_cscale);
	 solObj -= (offset * m_setset->size);
	 solObj  = -solObj;
	  
#ifdef MMKP_MCKNAP_DEBUG
	 printf("\nminObj = %g, solObj = %g", minObj, solObj);
#endif
	  
	  
	 int c, i, j, g;
	 for(g = 0; g < m_setset->size; g++){
	    recSolPtr = &solRec[g];
	    i         = recSolPtr->i;//was missing - STOP
	    j         = recSolPtr->j;//was missing
	    c         = getIndexIJ(i, j);
	    solInd.push_back(c);
	    varRedCost  += redCost[c];
	    varOrigCost += origCost[c];
#ifdef MMKP_MCKNAP_DEBUG	     
	    printf("\nc: %d = (%d,%d) redCost: %g cost: %d wt: %d",
		   c, i, j, 
		   redCost[c], m_cost[c], m_weight[c]);
#endif
	    /*for(c = 0; c < m_nCols; c++){
	    //but could have more than one equal in p and w!
	    //this is NOT the right way to get back optimal
	    if((m_cost[c]   == recSolPtr->psum) &&
	    (m_weight[c] == recSolPtr->wsum)){
	    //TODO: why should the user have to calc origCost?
	    //framework should probably do that
	    solInd.push_back(c);
	    varRedCost  += redCost[c];
	    varOrigCost += origCost[c];
	    #ifdef MMKP_MCKNAP_DEBUG
	    pair<int,int> ij = getIndexInv(c);
	    printf("\nc: %d = (%d,%d) redCost: %g cost: %d wt: %d",
	    c,
	    ij.first, ij.second,
	    redCost[c], m_cost[c], m_weight[c]);
	    #endif
	    break;
	    }
	    }*/
	 }
	 //UtilPrintVector<int>(solInd);
	 break;
      }
   case MCKNAP_RC_TRIVIAL_MAXSUM:
     fflush(stdout);
      //maxwsum <= c, so just pick the max elements from each group
      solveTrivialMaxSum(redCost, origCost, solInd, 
			 varRedCost, varOrigCost);
#ifdef MMKP_MCKNAP_DEBUG	     
      printf("trivial sum varRedCost=%g varOrigCost=%g\n",
	     varRedCost, varOrigCost);
#endif
      break;
   default:
      assert(0);
   }
   UTIL_DELARR(solRec);
}


