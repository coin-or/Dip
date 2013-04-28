//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, Lehigh University                                 //
//                                                                           //
// Copyright (C) 2002-2013, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

#ifndef MCF_INSTANCE_INCLUDED
#define MCF_INSTANCE_INCLUDED

//===========================================================================//
#include "UtilMacros.h"
//===========================================================================//
class MCF_Param;
using namespace std;
//===========================================================================//

//===========================================================================//
/*!
 * \class MCF_Instance
 * A class to store an instance of the 
 *     (Integer) Multi-Commodity Flow Problem (MCF).
 *
 *     min  sum{i in 1..n} v[i,j]   x[i,j]
 *     s.t. sum{i in 1..n, j in 1..l[i]} r[k,i,j] x[i,j] <= b[k], k in 1..m
 *          sum{j in 1..l[i]}                     x[i,j]  = 1   , i in 1..n
 *          x[i,j] in {0,1}, i in 1..n, j in 1..l[i]
 *
 */

//===========================================================================//
class MCF_Instance {
public:
   /** MCF_Instance problem instance data */   
   struct arc {
      int    tail;
      int    head;
      int    lb;
      int    ub;
      double weight;
   };
   struct commodity {
      int    source;
      int    sink;
      int    demand;
   };   
   string            m_problemName;
   arc             * m_arcs;
   commodity       * m_commodities;
   int               m_numNodes;
   int               m_numArcs;
   int               m_numCommodities;

public:
   /** @name Helper Methods. */
   int readInstance(string & fileName,
                    bool     addDummyArcs = true);

   inline void initMembers(){
      m_problemName     = "";
      m_arcs            = NULL;
      m_commodities     = NULL;
      m_numNodes        = 0;
      m_numArcs         = 0;
      m_numCommodities  = 0;
   }

public:
   /** @name Constructor and Destructor */

   /** Default constructor. */
   MCF_Instance(){
      initMembers();
   };
   
   /** Default constructor. Takes an instance of UtilParameters */
   MCF_Instance(string & fileName) {
      initMembers();
      readInstance(fileName);
   }
   
   /** Default destructor. */
   ~MCF_Instance() {
      UTIL_DELARR(m_arcs);
      UTIL_DELARR(m_commodities);
   };
};

#endif
