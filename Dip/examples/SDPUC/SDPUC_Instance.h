//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, Lehigh University                                 //
//                                                                           //
// Copyright (C) 2002-2009, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

#ifndef SDPUC_INSTANCE_INCLUDED
#define SDPUC_INSTANCE_INCLUDED

//===========================================================================//
#include "UtilMacros.h"
//===========================================================================//
class SDPUC_Param;
using namespace std;
//===========================================================================//

//===========================================================================//
/*!
 * \class SILCEP_Instance
 * A class to store an instance of the 
 *     Switch Investment and Line Capacity Expansion Problem.
 *
 *     min  sum{i in 1..n} v[i,j]   x[i,j]
 *     s.t. sum{i in 1..n, j in 1..l[i]} r[k,i,j] x[i,j] <= b[k], k in 1..m
 *          sum{j in 1..l[i]}                     x[i,j]  = 1   , i in 1..n
 *          x[i,j] in {0,1}, i in 1..n, j in 1..l[i]
 *
 */

//===========================================================================//
class SDPUC_Instance {
public:
   /** SDPUC_Instance problem instance data */   
   struct arc {
      int    tail;
      int    head;
      double lb;
      double ub;
      double weight; //e.g. reactance
	  double mcost;  //marginal cost
	  double fcost1; //arc installation cost
	  double fcost2; //switch installation cost
  	  int	 tscap;
	  int	 tscost;
	  int	 acline;  // arc is included in kirchoffs constraints (set to 0 for supply and hvdc-arcs)
	  int	 switchable;  // arc is switchable
   };
   struct node {
      int	   id;
      double   demand;
	  int	   tsdemand;
   };   
   struct timeseries {
	  int	 id;
	  double * values;
   };
   string       m_problemName;
   arc          * m_arcs;
   node			* m_nodes;
   timeseries   * m_timeseries;
   int          m_numNodes;
   int          m_numArcs;
   int		    m_numTimeseries;
   int		    m_numTimeperiods;
   int			m_numSwitchings;    // max. no. of simultaneously employed switches

public:
   /** @name Helper Methods. */
   int readInstance(string & fileName,
                    bool     addDummyArcs = true);

   inline void initMembers(){
      m_problemName     = "";
      m_arcs            = NULL;
	  m_nodes			= NULL;
      m_numNodes        = 0;
      m_numArcs         = 0;
	  m_numTimeseries	= 0;
	  m_numTimeperiods	= 0;
	  m_numSwitchings	= 0;
   }

public:
   /** @name Constructor and Destructor */

   /** Default constructor. */
   SDPUC_Instance(){
      initMembers();
   };
   
   /** Default constructor. Takes an instance of UtilParameters */
   SDPUC_Instance(string & fileName) {
      initMembers();
      readInstance(fileName);
   }
   
   /** Default destructor. */
   ~SDPUC_Instance() {
      UTIL_DELARR(m_arcs);
      UTIL_DELARR(m_nodes);
	  UTIL_DELARR(m_timeseries);
   };
};

#endif
