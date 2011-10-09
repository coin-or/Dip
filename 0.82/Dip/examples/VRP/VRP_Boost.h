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

#ifndef VRP_BOOST_INCLUDED
#define VRP_BOOST_INCLUDED

//---
//--- Interface class to Boost/Graph (an open source graph library).
//---

//THINK design:
//design wise there should be a class that just wraps boost algos
//then one that carries boost data specific to your application

// --------------------------------------------------------------------- //
#include <boost/graph/graphviz.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
using namespace boost;
using namespace std;

typedef property<edge_weight_t, double,
                 property<edge_index_t, int> > edge_prop;
typedef adjacency_list<vecS, vecS, undirectedS,
                       no_property, edge_prop> Graph;

// --------------------------------------------------------------------- //
//---
//--- utility classes for writing to graphviz format (edges)
//---   TODO: convert double to closest fraction?
//---
template < class EdgeWeight >
class sg_label_writer {
public:
   sg_label_writer(EdgeWeight e_) : weight(e_) {}
   template <class VertexOrEdge>
   void operator()(ostream & out, const VertexOrEdge & e) const {   
      if(weight[e] >= 0.9999999)
         out << "[style=solid]";
    else
       out << "[label=\"" << weight[e] << "\", stype=dashed]";
   }
private:
   EdgeWeight weight;
};
template <class EdgeWeight >
inline sg_label_writer<EdgeWeight> sg_make_label_writer(EdgeWeight e) {
  return sg_label_writer<EdgeWeight>(e);
}


// --------------------------------------------------------------------- //
//---
//--- utility classes for writing to graphviz format (vertices)
//---   TODO: convert double to closest fraction?
//---
class sg_label_writer_vertex {
public:
   sg_label_writer_vertex(const int * vertex_wt) 
      : demand(vertex_wt) {}
   template <class VertexOrEdge>
   void operator()(ostream & out, const VertexOrEdge & v) const {
      if(v == 0)
         out << " [style=filled, fillcolor=black]";
      else{
         out << " [label=\"" << v << " : " 
             << demand[v] << "\" width=\"0.5\" height=\"0.5\"]";
      }
   }
private:
   const int * demand;
};
inline sg_label_writer_vertex 
sg_make_label_writer_vertex(const int * vertex_wt) {
   return sg_label_writer_vertex(vertex_wt);
}

// --------------------------------------------------------------------- //
class VRP_Boost {
public:
   Graph  m_sg;          //for storage of support graph
   Graph  m_sg0;         //for storage of support graph minus depot
   double m_depotDegree; //degree at depot (since we can have single-customers)

public:
   VRP_Boost() :
      m_sg         (),
      m_sg0        (),
      m_depotDegree(0.0)
   {}
   ~VRP_Boost() {}
   
public:
   /** @name Helper Functions */
   void buildSubGraph(const int      len,
                      const double * val,
                      const double   tol = 1.0e-6){      

      //---
      //--- build the boost subgraph as the support graph for x=val
      //---  at the same time, calculate the depot degree since we can 
      //---  have single-customer routes with x=2.0
      //---
      property_map<Graph, edge_weight_t>::type e_weight 
	 = get(edge_weight, m_sg);
      property_map<Graph, edge_index_t>::type  e_index  
	 = get(edge_index,  m_sg);
      graph_traits<Graph>::edge_descriptor ed; bool inserted;      

      int c;      
      m_depotDegree = 0.0;
      for (c = 0; c < len; ++c) {
	 if(val[c] <= tol)
	    continue;
	 pair<int,int> uv  = UtilBothEndsU(c);
	 
	 if(uv.first == 0 || uv.second == 0){
	    m_depotDegree += val[c];
	 }

	 tie(ed, inserted) = add_edge(uv.first, uv.second, m_sg);
	 e_weight[ed]      = val[c];
	 e_index[ed]       = c;
      }      
   }

   /*
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
	       e_index[ed]    = index;
	    }
	    index++;    
	 }
      }
   }
   */

   inline void printGraph(const Graph & g) const {
      print_graph(g);
   }

   inline void printDotFile(const string & fileName,
                            const int    * vertexWt,
                            const Graph  & g) const {

      ofstream os;
      UtilOpenFile(os, fileName);

      printf("PrintDotFile fileName = %s\n", fileName.c_str());

      map<string, string> graph_attr, vertex_attr, edge_attr;
      graph_attr ["fontsize"] = "20";

      vertex_attr["fontsize"] = "20";
      vertex_attr["shape"]    = "circle";

      edge_attr  ["style"]    = "dotted";
      edge_attr  ["len"]      = "1.7";
      edge_attr  ["fontsize"] = "20";
      
      //property_map<Graph, edge_weight_t>::type e_weight
      //   = get(edge_weight, g);
      
      write_graphviz(os, g, 
                     sg_make_label_writer_vertex(vertexWt),
                     sg_make_label_writer(get(edge_weight, g)),
                     make_graph_attributes_writer(graph_attr,
                                                  vertex_attr, 
                                                  edge_attr));
      os.close();
   }

   inline void copyGraph(const Graph & gFrom,
			 Graph       & gTo){
      gTo = gFrom;
   }

   inline void clearSubGraph(){
      m_sg.clear();
   }
   
   inline void clearSubGraph(Graph & g){
      g.clear();
   }

   inline void clearVertex(const int vertex,
			   Graph &   g){
      clear_vertex(0, g);
   }

   inline int findConnectedComponents(vector<int> & component){
      return connected_components(m_sg, &component[0]);
   }

   inline int findConnectedComponents(Graph       & g,
				      vector<int> & component){
      return connected_components(g, &component[0]);
   }

   inline int getDegree(const int nodeIndex){
      return degree(nodeIndex, m_sg);
   }

   inline double getDepotDegree(){
      return m_depotDegree;
   }
   
};

#endif
