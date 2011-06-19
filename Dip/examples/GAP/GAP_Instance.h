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

#ifndef GAP_INSTANCE_INCLUDED
#define GAP_INSTANCE_INCLUDED

//===========================================================================//
#include "UtilMacros.h"

using namespace std;

//===========================================================================//
/*!
 * \class GAP_Instance
 * A class to store an instance of the 
 *     Generalized Assignment Problem (GAP).
 * 
 * Find the maximum profit assignment of n tasks to m machines
 * such that each task is assinged to precisely one machine subject
 * to capacity restrictions of the machine.
 *
 *    max  sum{i in 1..m, j in 1..n} p[i,j] x[i,j]
 *    s.t. sum{           j in 1..n} w[i,j] x[i,j] <= b[i], i in 1..m
 *         sum{i in 1..m           }        x[i,j]  = 1   , j in 1..n
 *         x[i,j] in {0,1}, i in 1..m, j in 1..n
 *
 *    x[i,j]=1 means assign task j to agent i
 *
 * Note: DIP does min, so, we solve for min sum{ij} -p[i,j] x[i,j].
 *
 */

//===========================================================================//
class GAP_Instance {

private:
   /** GAP_Instance problem instance data */   
   int    m_nTasks;       //n (use j index)
   int    m_nMachines;    //m (use i index)
   int *  m_capacity;     //b[i = 1..m]
   int *  m_profit;       //p[i,j]
   int *  m_weight;       //w[i,j]
   
   /** GAP_Instance best known LB/UB */
   bool      m_isProvenOptimal;
   double    m_bestKnownLB;
   double    m_bestKnownUB;


public:
   /** @name Access methods. */
   inline const int    getNTasks   ()  const {return m_nTasks;   }
   inline const int    getNMachines()  const {return m_nMachines;}
   inline const int *  getCapacity ()  const {return m_capacity; }
   inline const int *  getProfit   ()  const {return m_profit;   }
   inline const int *  getWeight   ()  const {return m_weight;   }

public:
   /** @name Helper Methods. */
   void readInstance(string & filename);
   void readBestKnown(string & fileName,
                      string & instanceName);

   inline void initMembers(){
      m_nTasks    = 0;
      m_nMachines = 0;
      m_capacity  = NULL;
      m_profit    = NULL;
      m_weight    = NULL;
      m_isProvenOptimal =  false;
      m_bestKnownLB     = -1.e20;
      m_bestKnownUB     =  1.e20;
   }

   inline const int getIndexIJ(const int i,
                               const int j) const{
      return (i * m_nTasks) + j;
   }
   
   inline pair<int,int> getIndexInv(const int index) const {      
      return make_pair(index / m_nTasks, index % m_nTasks);
   }

   inline const double getBestKnownLB() const {return m_bestKnownLB;}
   inline const double getBestKnownUB() const {return m_bestKnownUB;}
   
public:
   /** @name Constructor and Destructor */

   /** Default constructor. */
   GAP_Instance(){
      initMembers();
   };
   
   /** Default constructor. Takes an instance of UtilParameters */
   GAP_Instance(string & fileName) {
      initMembers();
      readInstance(fileName);
   }

   ~GAP_Instance() {
      UTIL_DELARR(m_capacity);
      UTIL_DELARR(m_profit);
      UTIL_DELARR(m_weight);
   };
};

#endif
