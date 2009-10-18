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

#ifndef GAP_KNAP_PISINGER_INCLUDED
#define GAP_KNAP_PISINGER_INCLUDED

// --------------------------------------------------------------------- //
/** 
 * \class GAP_KnapPisinger
 * 
 * An interface class for Pisinger's Combo Algorithm for 
 * 0-1 Knapsack Problem. See \url{http://www.diku.dk/~pisinger/codes.html}.
 *
 * The 0-1 Knapsack Problem:
 *   max  sum{j = 0..n-1} p[j] x[j]
 *   s.t. sum{j = 0..n-1} w[j] x[j] <= b
 *        x binary
 */

/**
 * Types and protos for combo algorithm API:
 *
 * typedef int           boolean; // logical variable         
 * typedef int           ntype;   // number of states/items   
 * typedef long          itype;   // item profits and weights 
 * typedef long          stype;   // sum of profit or weight  
 *
 * typedef struct {
 *   itype   p;              // profit                  
 *   itype   w;              // weight                  
 *   boolean x;              // solution variable       
 * } item;
 *
 *
 * extern stype combo(item  * f, 
 *                    item  * l, 
 *                    stype   c, 
 *                    stype   lb, 
 *                    stype   ub,
 *                    boolean def, 
 *                    boolean relx);
 * f,l : first, last item                                               
 * c   : capacity of knapsack                                           
 * lb  : lower bound. Solution vector is only updated if better z found 
 * ub  : upper bound. When upper bound is reached terminate immediately 
 * def : should solution vector be defined or just find objective value 
 * relx: relaxed problem is solved (no more relaxations will be made)   
 * returns the objective value of the problem                           
 *
 */

// --------------------------------------------------------------------- //
extern "C"{
#include "combo.h"
}

// --------------------------------------------------------------------- //
#include "UtilMacros.h"

// --------------------------------------------------------------------- //
class GAP_KnapPisinger{

private:
   int      m_nItems;   //< number of items (== nTasks)
   item   * m_items;    //< array of items for Pisinger solver   
   int      m_capacity; //< capacity of knapsack  
   int    * m_weight;   //< weight of items
   int    * m_profit;   //< profit of items 
   double * m_workDbl;  //< temp work space 

private:
   inline const int getIndexIJ(const int i,
                               const int j) const {
      return (i * m_nItems) + j;
   }
   
   //1.0e-6 seems to cause overflow
   int calcScaleFactor(const double * redCost,
		       const double   epsTol = 1.0e-4){
      
      //---
      //--- A very simple function to scale an array of doubles to integers.
      //---    Note: epstol denotes the preferred accuracy,
      //---    so, we will scale by 1.0/epstol, unless something smaller works.
      //--- It returns the scale factor.
      //---
      //--- Since the knapsack solver solves the max problem, we will
      //--- flip the sign of the redCost first (max p == min -p).
      //---
      //--- This also means that we can look just at redCost[i] <= 0.
      //---      
      int      i, scaleFactor = 1, n_aFrac = 0, factorToBig = 0;
      double   fractionalPart, rcFlipped;
      double   oneOverEps = 1.0 / epsTol;      
      for(i = 0; i < m_nItems; i++){
	 if(redCost[i] >= 0)
	    continue;
	 rcFlipped = -redCost[i];
	 fractionalPart = UtilFracPart(rcFlipped);
	 if(!UtilIsZero(fractionalPart)){
	    fractionalPart     *= oneOverEps;
	    m_workDbl[n_aFrac++]  = (int)fractionalPart * (double)epsTol;
	 }
      }
      for(i = 0; i < n_aFrac; i++){
	 CoinAssertDebug(m_workDbl[i] < (COIN_INT_MAX / scaleFactor));
	 m_workDbl[i] *= scaleFactor;
	 while(!UtilIsZero(UtilFracPart(m_workDbl[i]))){
	    scaleFactor *= 10;
	    if(scaleFactor >= oneOverEps){
	       factorToBig = 1;
	       break;
	    }
	    CoinAssertDebug(m_workDbl[i] < (COIN_INT_MAX / 10));
	    m_workDbl[i]    *= 10;
	    CoinAssertDebug(m_workDbl[i] >= 0);
	 }
	 if(factorToBig)
	    break;
      }
      return scaleFactor;
   }

public:
   void solve(const int        blockB,
	      const double   * redCost,
	      const double   * origCost,
	      vector<int>    & solInd,
	      vector<double> & solEls,
	      double         & varRedCost,
	      double         & varOrigCost){


      //---
      //--- Since we want to find min redCost, and the constraint is a 
      //--- simple knapsack constraint, we know that if redCost[i] >= 0
      //--- then we can fix x[i]=0. It makes no sense to select.
      //---
      //--- The knapsack solver solves maximization problem, so we also 
      //--- need to fip the sign of the redCost.
      //--- 
      //--- The knapsack solver also expects integers for weight and profit, 
      //--- the weights are already integer, but since the profits are coming
      //--- from reduced cost, we need to scale them to integers
      //---
      int i, j;
      int nItems      = 0;
      int scaleFactor = calcScaleFactor(redCost);            

      double      rcScaled;
      int         weightSum = 0;
      vector<int> shrunkToOrig;
      for(i = 0; i < m_nItems; i++){
	 if(redCost[i] < 0.0) {	    
	    rcScaled    = (-redCost[i]) * scaleFactor;	    
	    long rcLong = static_cast<long>(floor(rcScaled + 0.5));
	    CoinAssert(rcScaled < COIN_INT_MAX);
	    
	    m_items[nItems].p = rcLong;
	    CoinAssert(m_items[nItems].p >= 0);
	    
	    m_items[nItems].w = m_weight[i];
	    m_items[nItems].i = nItems;
	    m_items[nItems].x = 0;

	    weightSum += m_weight[i];
	    
	    shrunkToOrig.push_back(i);
	    nItems++;
	 }
      }

      //---
      //--- print the current problem
      //---
      /*
      printf("COMBO Knap Max (nItems = %d, cap=%d)\n", nItems, m_capacity);
      for(i = 0; i < nItems; i++){
	 j = shrunkToOrig[i];
	 printf("[%d->%d] pOrig:%g, p:%ld, w:%ld\n",
		i, j, -redCost[j], m_items[i].p, m_items[i].w);
      }
      */

      //---
      //--- solve the knapsack problem
      //---   first deal with some trivial cases
      if(nItems == 1){
	 CoinAssert(m_items[0].w <= m_capacity);
	 m_items[0].x = 1;
      }
      else if(weightSum <= m_capacity){
	 //---
	 //--- then put all items in the knapsack (to max profit)
	 //---
	 for(i = 0; i < nItems; i++)
	    m_items[i].x = 1;
      }
      else if(nItems == 2){
	 CoinAssert(m_items[0].w <= m_capacity);
	 CoinAssert(m_items[1].w <= m_capacity);
	 int w0 = m_items[0].w;
	 int w1 = m_items[1].w;
	 int p0 = m_items[0].p;
	 int p1 = m_items[1].p;
	 double bestObj = 0;
	 if( (w1 <= m_capacity) && (p1 > bestObj) ){
	    bestObj      = p1;
	    m_items[0].x = 0;
	    m_items[1].x = 1;
	 }
	 if( (w0 <= m_capacity) && (p0 > bestObj) ){
	    bestObj      = p0;
	    m_items[0].x = 1;
	    m_items[1].x = 0;
	 }
	 if( (w0 + w1 <= m_capacity) && (p0 + p1 > bestObj) ){
	    bestObj      = p0 + p1;
	    m_items[0].x = 1;
	    m_items[1].x = 1;
	 }	 
      }
      else if(nItems > 2){
	 item * fItem  = m_items;
	 item * lItem  = m_items + (nItems-1);
	 long   obj    = 0;
	 obj = combo(fItem,           //first item
		     lItem,           //last item
		     m_capacity,      //capacity of knapsack
		     -COIN_INT_MAX,   //obj lower bound
		     COIN_INT_MAX,    //obj upper bound
		     1,               //define sol vector
		     0);              //do not solve relaxed prob
	 //printf("combo obj = %ld\n", obj);
      }
      
      int origColIndex;
      /*int wtSum = 0;
	for(i = 0; i < nItems; i++)
	printf("m_items[%d].x=%d i=%d p=%ld w=%ld\n", 
	i, 
	m_items[i].i,
	m_items[i].x,
	m_items[i].p,
	m_items[i].w);*/
      
      for(i = 0; i < nItems; i++){
	 if(m_items[i].x == 1){
	    j            = shrunkToOrig[m_items[i].i];
	    origColIndex = getIndexIJ(blockB, j);
	    varRedCost  += redCost[j];
	    varOrigCost += origCost[j];
	    solInd.push_back(origColIndex);
	    solEls.push_back(1.0);
	 }
      }
   }
   
public:
   GAP_KnapPisinger(const int   nItems,
		    const int   capacity,
		    const int * weight,
		    const int * profit) :
      m_nItems   (nItems),
      m_items    (0),
      m_capacity (capacity),
      m_weight   (0),
      m_profit   (0),
      m_workDbl  (0)
   {
      m_items    = new item[nItems];
      m_weight   = new int[nItems];
      m_profit   = new int[nItems];
      m_workDbl  = new double[nItems];
      CoinAssertHint(m_items && m_weight && m_profit && m_workDbl, 
		     "Error: Out of Memory");
      
      memcpy(m_weight, weight, nItems * sizeof(int));
      memcpy(m_profit, profit, nItems * sizeof(int));
   }
   ~GAP_KnapPisinger() {
      UTIL_DELARR(m_items)
	 UTIL_DELARR(m_weight);
      UTIL_DELARR(m_profit);
      UTIL_DELARR(m_workDbl);
   }
};


#endif
