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

#ifndef TSP_BOOST_INCLUDED
#define TSP_BOOST_INCLUDED

//---
//--- Interface class to Boost/Graph (an open source graph library).
//---

//THINK design:
//design wise there should be a class that just wraps boost algos
//then one that carries boost data specific to your application

// --------------------------------------------------------------------- //
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
using namespace boost;

typedef property<edge_weight_t, double,
                 property<edge_index_t, int> > edge_prop;
typedef adjacency_list<vecS, vecS, undirectedS,
                       no_property, edge_prop> Graph;

// --------------------------------------------------------------------- //
class TSP_Boost {
public:
   Graph m_sg;
   Graph m_cgV;

public:
   TSP_Boost() :
      m_sg (),
      m_cgV()  {}
   ~TSP_Boost() {}
   
public:
   /** @name Helper Functions */
   void buildSubGraph(const int      len,
                      const double * val,
                      const double   tol = 1.0e-6){      

      property_map<Graph, edge_weight_t>::type e_weight 
	 = get(edge_weight, m_sg);
      property_map<Graph, edge_index_t>::type  e_index  
	 = get(edge_index,  m_sg);
      graph_traits<Graph>::edge_descriptor ed; bool inserted;      
      for (int c = 0; c < len; ++c) {
	 if(val[c] <= tol)
	    continue;
	 pair<int,int> uv  = UtilBothEndsU(c);
	 tie(ed, inserted) = add_edge(uv.first, uv.second, m_sg);
	 e_weight[ed]      = val[c];
	 e_index[ed]       = c;
      }      
   }

   void buildCompleteGraphMinusVert(const int vert,
				    const int nVerts){
      
      property_map<Graph, edge_index_t>::type  e_index  
	 = get(edge_index,  m_cgV);
      graph_traits<Graph>::edge_descriptor ed; bool inserted;            
      int u, v;
      int index = 0;
      for(u = 1; u < nVerts; u++){
	 for(v = 0; v < u; v++){
	    if(u != vert && v != vert){
	       tie(ed, inserted) = add_edge(u, v, m_cgV);
	       e_index[ed]       = index;
	    }
	    index++;    
	 }
      }
   }

   inline void clearSubGraph(){
      m_sg.clear();
   }

   inline int findConnectedComponents(vector<int> & component){
      return connected_components(m_sg, &component[0]);
   }

   inline int getDegree(const int nodeIndex){
      return static_cast<int>(degree(nodeIndex, m_sg));
   }
   
};

#endif
