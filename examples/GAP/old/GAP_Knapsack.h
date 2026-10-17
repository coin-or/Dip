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

#ifndef GAP_KNAPSACK_INCLUDED
#define GAP_KNAPSACK_INCLUDED

//TODO: consider one container object rather than K and then 
//  load constraint and profit each time?

// --------------------------------------------------------------------- //
/*!
 * \class GAP_Knapsack
 * A class to store an instance of the 
 *     Knapsack Problem
 *
 *   max  sum{j = 0..n-1} p[j] x[j]
 *   s.t. sum{j = 0..n-1} w[j] x[j] <= b
 *        x binary
 */

// --------------------------------------------------------------------- //
#include "UtilKnapsack.h"

// --------------------------------------------------------------------- //
class GAP_Knapsack {

 private:
   /** GAP_Knapsack problem instance data */   
   int             m_length;     //n (use j index)   
   double          m_capacity;   //b
   double        * m_profit;     //p[j]
   double        * m_weight;     //w[j]
   double        * m_profitSort; //temp storage for KP HS algo
   double        * m_weightSort; //temp storage for KP HS algo
   int           * m_solution;   //storage for solution
   double          m_optObj;     //storage for optimal objective
   SOR_IntDblArr * m_ratio;      //temp storage for ratios

 public:
   /** @name Access get methods. */
   inline const int       getLength     ()  const {return m_length;     }
   inline const double    getCapacity   ()  const {return m_capacity;   }
   inline const double *  getProfit     ()  const {return m_profit;     }
   inline const double *  getWeight     ()  const {return m_weight;     }
   inline const double *  getProfitSort ()  const {return m_profitSort; }
   inline const double *  getWeightSort ()  const {return m_weightSort; }
   inline const int    *  getSolution   ()  const {return m_solution;   }
   inline const double    getOptObj     ()  const {return m_optObj;     }
   inline const SOR_IntDblArr * getRatio() const  {return m_ratio;      }
   inline double       *  getMutableProfit ()     {return m_profit;     }
   inline double       *  getMutableWeight ()     {return m_weight;     }
   
   /** @name Set methods. */
   void setLength(const int length) {m_length = length;}

 public:
   /** @name Helper Methods. */
   void populateAndSortRatio(const double * profit = 0,
                             const double * weight = 0){
      const double * p = profit;
      const double * w = weight;
      if(!p) 
         p = m_profit;
      if(!w) 
         w = m_weight;
      KnapsackSortRatioOut(m_length, p, w, 
                           m_profitSort, 
                           m_weightSort, 
                           m_ratio->arr);
   }

   int solveKnapsack(){
      int    statusKP = 0;
      int    status   = 0;
      statusKP = 
         KnapsackOptimizeHS(m_length, m_capacity,
                            m_profitSort, m_weightSort, 
                            m_solution, &m_optObj, &status);
      CoinAssert(!statusKP && !status);
      //printf("Knapsack OptObj = %g\n", m_optObj);
      return status;
   }

 private:
   /** @name Copy Constructors
    *
    * Disable the default copy constructors.
    *
    */
   GAP_Knapsack(const GAP_Knapsack &);
   GAP_Knapsack & operator=(const GAP_Knapsack &);

 public:
   /** @name Constructor and Destructor */
   
   /** Default constructor. Takes an instance of UtilParameters */
   GAP_Knapsack(const int      length,
                const double   capacity,
                const double * profit,
                const double * weight) :
      m_length    (length),
      m_capacity  (capacity),
      m_profit    (0),
      m_weight    (0),
      m_profitSort(0),
      m_weightSort(0),
      m_solution  (0),
      m_ratio     (0)
      {
         int status;
         CoinAssert(length >= 0);
         m_profit     = new double[length];
         m_weight     = new double[length];
         m_profitSort = new double[length+1];
         m_weightSort = new double[length+1];
         m_solution   = new int[length];
         CoinAssertHint(m_profit     && m_weight && 
                        m_profitSort && m_weightSort && m_solution, 
                        "Error: Out of Memory");
         memcpy(m_profit, profit, length * sizeof(double));
         memcpy(m_weight, weight, length * sizeof(double));

         m_ratio = SOR_IntDblArrNew(length, &status);
         CoinAssert(!status);
      }
   
   ~GAP_Knapsack() {      
      UTIL_DELARR(m_profit);
      UTIL_DELARR(m_weight);
      UTIL_DELARR(m_profitSort);
      UTIL_DELARR(m_weightSort);
      UTIL_DELARR(m_solution);
      SOR_IntDblArrFree(&m_ratio);
   };
};

#endif
