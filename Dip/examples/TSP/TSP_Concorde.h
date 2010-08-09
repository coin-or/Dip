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

#ifndef TSP_CONCORDE_INCLUDED
#define TSP_CONCORDE_INCLUDED

//---
//--- Interface class to Concorde (TSP solver).
//---

// --------------------------------------------------------------------- //
extern "C"{
#include "TSP_ConcordeI.h"
}

// --------------------------------------------------------------------- //
#include <vector>
using namespace std;

// --------------------------------------------------------------------- //
#include "UtilMacros.h"

// --------------------------------------------------------------------- //
class ConcordeSubtourCut {
public:
   vector<int>  S;
   vector<bool> inS;
public:
   void print(){
      cout << "ConcordeSubtour: ";
      for(int i = 0; i < static_cast<int>(inS.size()); i++){
	 if(inS[i]){
	    cout << i << " ";
	 }
      }
      cout << endl;
   }
   ConcordeSubtourCut(const int nVerts) :
      S(),
      inS(nVerts, false) {}
   
   ~ConcordeSubtourCut() {}
};

// --------------------------------------------------------------------- //
class TSP_Concorde {
public:
   //define concorde subgraph
   int            m_nVerts;
   vector<int>    m_edgeList;
   vector<double> m_edgeValue;
   
public:
   /** @name Helper Functions */
   void clearSubGraph(){
      m_nVerts = 0;
      m_edgeList.clear();
      m_edgeValue.clear();
   }

   void buildSubGraph(const int      nVerts,
		      const int      nEdges,
                      const double * edgeWeight,
                      const double   tol = 1.0e-6){      

      int i;
      clearSubGraph();
      m_nVerts = nVerts;
      m_edgeList.reserve(2*nEdges);
      m_edgeValue.reserve(nEdges);      
      for(i = 0; i < nEdges; i++){
         if(edgeWeight[i] <= tol)
            continue;
         pair<int,int> uv = UtilBothEndsU(i);
         m_edgeValue.push_back(edgeWeight[i]);
         m_edgeList.push_back(uv.first);
         m_edgeList.push_back(uv.second);
      }
   }

   int generateCutsSubtour(vector<ConcordeSubtourCut> & subtourCuts){
      CCtsp_lpcut_in * tsp_cuts     = NULL;
      CCtsp_lpcut_in * tsp_cut      = NULL;
      CCtsp_lpcut_in * tsp_cut_next = NULL;
      
      int i, c, tmp;
      int n_subtour  = 0;
      int n_edgeList = static_cast<int>(m_edgeValue.size());
      int retCode    = CCtsp_exact_subtours(&tsp_cuts, 
					    &n_subtour, 
					    m_nVerts,
					    n_edgeList,
					    &m_edgeList[0], 
					    &m_edgeValue[0]);
      assert(!retCode);
      
      //cout << "CONCORDE found n_subtours : " << n_subtour << endl;

      tsp_cut = tsp_cuts;      
      for(c = 0; c < n_subtour; c++){
	 ConcordeSubtourCut subtourCut(m_nVerts);	 
	 assert(tsp_cut);
	 CC_FOREACH_NODE_IN_CLIQUE(i, tsp_cut->cliques[0], tmp){
	    subtourCut.S.push_back(i);
	    subtourCut.inS[i] = true;
	 }
	 //subtourCut.print();
	 subtourCuts.push_back(subtourCut);
	 tsp_cut = tsp_cut->next;
      }
      
      for (tsp_cut = tsp_cuts; tsp_cut; tsp_cut = tsp_cut_next){
	 tsp_cut_next = tsp_cut->next;
	 CCtsp_free_lpcut_in(tsp_cut);
	 free(tsp_cut); //need UtilFree? for C style alloc/free?
      }
      tsp_cuts = NULL;
      return static_cast<int>(subtourCuts.size());
   }
   
public:
   TSP_Concorde() :
      m_nVerts   (0),
      m_edgeList (),
      m_edgeValue()
   {}
   ~TSP_Concorde() 
   {}
};

#endif
