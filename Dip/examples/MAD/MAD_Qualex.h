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

#ifndef MAD_QUALEX_INCLUDED
#define MAD_QUALEX_INCLUDED

//---
//--- The MAD Qualex -- interface class for Qualex-MS solver
//---      http://www.busygin.dp.ua/npc.html
//---

// --------------------------------------------------------------------- //
#include <algorithm>
using namespace std;

// --------------------------------------------------------------------- //
#include "qualex.h"
#include "bool_vector.h"
#include "preproc_clique.h"
#include "greedy_clique.h"

// --------------------------------------------------------------------- //
class MAD_Qualex{
 public:
   
   //private:
   MaxCliqueInfo  * m_info;
   Graph          * m_graph;         //temp storage for compl graph
   Graph          * m_graphOrig;     //temp storage for compl graph (orig)
   double         * m_ijMatrix;      //temp storage for qualex-ms

   vector<int>      m_residual;
   list<int>        m_preselected;
   list<int>        m_clique;
   double           m_preselected_weight;
   double           m_cliqueWeight;

   


 public:  
   inline void addEdge(const int u,
                       const int v){
      //---
      //--- add an edge to the original graph
      //---   NOTE: qualex will set edge u-v and v-u
      //---
      m_graphOrig->add_edge(u,v);
   }

   inline void removeEdge(const int u,
                          const int v){
      //---
      //--- remove an edge to the original graph
      //---
      m_graphOrig->remove_edge(u,v);
   }

   void removeNonNegVertices(const double * weights){
      //---
      //--- "remove" vertices with non-negative weights
      //--- since the algorithm will look for cliques, we can 
      //--- just remove all the edges incident to the vertex
      //--- therefore, it can't be selected as part of a clique
      //---
      int i, j;
      for(i = 0; i < m_graph->n; i++){
         if(weights[i] >= 0){
	    //printf("\nZERO OUT i: %d, weight: %g", i, weights[i]);

            vector<bool_vector>::iterator bi = m_graph->mates.begin() + i;
            for(j = 0; j < m_graph->n; j++){
               if(bi->at(j)){
                  m_graph->remove_edge(i,j);
               }
            }
            //m_graph->mates[i].zero();
         }
      }
   }
   
   void setUpForSolve(const double * vertexWeight){
      
      //---
      //--- cleanup old copies, if exist
      //---
      UTIL_DELPTR(m_graph);
      UTIL_DELPTR(m_info);
      m_residual.clear();
      m_preselected.clear();
      m_clique.clear();
      m_cliqueWeight = 0.0;
      
      //---
      //--- allocate new graph
      //---
      m_graph = new Graph(m_graphOrig->n);
      CoinAssertHint(m_graph, "Error: Out of Memory");

      //---
      //--- copy orig graph adjMatrix into active graph
      //---  NOTE: this is neccessary because preprocessor
      //---        can change the structure of active graph
      //---
      u_long * dataSrc  = NULL;
      u_long * dataDst  = NULL;
      unsigned int i;
      int      dataSize = 0;
      for(i = 0; i < m_graphOrig->mates.size(); i++){
         dataSrc  = m_graphOrig->mates[i].data;
         dataDst  = m_graph->mates[i].data;         
         dataSize = m_graph->mates[i].data_size;
         memcpy(dataDst, dataSrc, dataSize * sizeof(u_long));         
      }

      //---
      //--- copy vertex weights
      //---
      copy(vertexWeight, 
           vertexWeight + m_graph->n, 
           m_graph->weights.begin());
   }
   
   void findMaxIndSetGreedy(const double * vertexWeight,
                            bool           forClique    = true){

      //printf("==== START findMaxIndSetGreedy\n");
      //fflush(stdout);
      
      //---
      //--- setup
      //---
      //setUpForSolve(vertexWeight);

      //---
      //--- preprocess the current graph
      //---
#if 0
      //bug in preproc_clique?
      m_preselected_weight = preproc_clique(*m_graph, 
                                            m_residual,
                                            m_preselected,
                                            m_cliqueWeight,
                                            m_clique);
#else
      int i;
      m_residual.resize(m_graph->n);
      for(i = 0; i < m_graph->n; i++)
	 m_residual[i] = i;
#endif

#if 0
      printf("\npreselected   = %d", m_preselected.size());
      printf("\nresidual size = %d", m_residual.size());
      UtilPrintVector<int>(m_residual);
#endif

      if(!m_residual.empty()){
         //---
         //--- setup MaxCliqueInfo object
         //---
         m_info = new MaxCliqueInfo(*m_graph, forClique);
         CoinAssertHint(m_info, "Error: Out of Memory");
         
         //---
         //--- call greedy heuristic
         //---
         meta_greedy_clique(*m_info);  

	 //---
	 //--- if greedy improved upon ?, save the clique
	 //---
	 list<int>::iterator it;
	 if(m_info->lower_clique_bound > m_cliqueWeight) {
	    m_cliqueWeight = m_info->lower_clique_bound;
	    m_clique.erase(m_clique.begin(), m_clique.end());
	    for(it = m_info->clique.begin(); it != m_info->clique.end(); it++)
	       m_clique.push_back(m_residual[*it]);
	 }
      }
      
      //---
      //--- join preselected vertices to maximum clique
      //---
      m_clique.splice(m_clique.begin(), m_preselected);
      m_cliqueWeight += m_preselected_weight;

#if 0
      printf("\npreselected   = %d", m_preselected.size());
      printf("\nresidual size = %d", m_residual.size());
      printf("\nm_clique  size = %d, wt = %g",
	     m_clique.size(), m_cliqueWeight);
      UtilPrintVector<int>(m_residual);
      UtilPrintList<int>(m_clique);
#endif

      //printf("\n==== END findMaxIndSetGreedy\n");
      //fflush(stdout);
   }
   
   void findMaxIndSetQualexMS(){

      if(!m_residual.empty()){
         printf("==== START findMaxIndSetQualexMS\n");
         fflush(stdout);
         
         //this should not be called before Greedy
         int      i,j;
         int      n     = m_graph->n;
         double * sqrtw = m_info->sqrtw;
         
         memset(m_ijMatrix, 0, sizeof(double) * n * n);
         for(i = 0; i < n; i++) {
            m_ijMatrix[i*(n+1)] = m_graph->weights[i] - m_info->w_min;
            bit_iterator bi(m_graph->mates[i]);
            while((j = bi.next()) > -1) {
               if(j > i) 
                  break;
               m_ijMatrix[i*n+j] = m_ijMatrix[j*n+i] = sqrtw[i] * sqrtw[j];
            }
         }      
         qualex_ms(*m_info, m_ijMatrix);

         //---
         //--- if qualex improved upon greedy, save the clique
         //---
         list<int>::iterator it;
         if(m_info->lower_clique_bound > m_cliqueWeight) {
            m_cliqueWeight = m_info->lower_clique_bound;
            m_clique.erase(m_clique.begin(), m_clique.end());
            for(it = m_info->clique.begin(); it != m_info->clique.end(); it++)
               m_clique.push_back(m_residual[*it]);
         }

         //---
         //--- join preselected vertices to maximum clique
         //---
         m_clique.splice(m_clique.begin(), m_preselected);
         m_cliqueWeight += m_preselected_weight;

#if 0
	 printf("\npreselected   = %d", m_preselected.size());
	 printf("\nresidual size = %d", m_residual.size());
	 printf("\nm_clique  size = %d, wt = %g",
		m_clique.size(), m_cliqueWeight);
	 UtilPrintVector<int>(m_residual);
	 UtilPrintList<int>(m_clique);
#endif
         
         printf("==== END findMaxIndSetQualexMS\n");
         fflush(stdout);
      }
   }


   void printGraph(Graph * G){
      int                   i, j;
      int                   n     = G->n;
      vector<bool_vector> & mates = G->mates;
      for(i = 0; i < n; i ++){
         printf("%2d -> ", i);
         bool_vector & matesI = mates[i];
         assert(matesI.size == n);
         for(j = 0; j < matesI.size; j++){
            if(matesI.at(j))
               printf(" %d", j);
         }
		 printf("\n");
      }
   }
   
   
 public:
   MAD_Qualex(const int nVertices) :
      m_info              (NULL),
      m_graph             (NULL),
      m_graphOrig         (NULL),
      m_ijMatrix          (NULL), 
      m_residual          (),
      m_preselected       (),
      m_clique            (),
      m_preselected_weight(0.0),
      m_cliqueWeight      (0.0)
      
   {
      CoinAssert(nVertices > 0);

      m_ijMatrix = new double[nVertices * nVertices];
      CoinAssertHint(m_ijMatrix, "Error: Out of Memory");
      
      m_graphOrig = new Graph(nVertices);
      CoinAssertHint(m_graphOrig, "Error: Out of Memory");
   }
   ~MAD_Qualex() {
      UTIL_DELPTR(m_info);
      UTIL_DELPTR(m_graph);
      UTIL_DELPTR(m_graphOrig);
      UTIL_DELARR(m_ijMatrix);
   }
};


#endif
