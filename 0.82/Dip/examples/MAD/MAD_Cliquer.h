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

#ifndef MAD_CLIQUER_INCLUDED
#define MAD_CLIQUER_INCLUDED

//---
//--- The MAD Cliquer -- it sounds like a college football mascot.
//---
//---    Ironic, since my college baseball team's mascot was the 
//---    "MAD Hatter" (Stetson University '94-'98). They use to pay
//---    a guy to get drunk and run around the games dresesd as a 
//---    a big cowboy hat. We were the coolest team in the league.
//---

// --------------------------------------------------------------------- //
extern "C"{
#include "MAD_CliquerI.h"
}

// --------------------------------------------------------------------- //
static int record_clique_func(set_t            s,
                              graph_t        * g,
                              clique_options * opts) {
   
   //---
   //--- make a copy of clique and append to clique_list
   //---
   int * clique_countPtr = (int*) opts->user_data;
   if ((*clique_countPtr) >= opts->clique_list_length) {
      opts->clique_list
         = (set_t*) realloc(opts->clique_list,
                            (opts->clique_list_length+512)*sizeof(set_t));
      opts->clique_list_length += 512;
   }
   opts->clique_list[*clique_countPtr] = set_duplicate(s);
   (*clique_countPtr)++;
   return 1;
}

// --------------------------------------------------------------------- //
static int set_clique_func(set_t            s,
                           graph_t        * g,
                           clique_options * opts) {
   
   //---
   //--- append clique to clique_list
   //---

   //---
   //--- store clique_list_size in user_data
   //---
   int * clique_countPtr = (int*) opts->user_data;
   if ((*clique_countPtr) >= opts->clique_list_length) {
      opts->clique_list
         = (set_t*) realloc(opts->clique_list,
                            (opts->clique_list_length+512)*sizeof(set_t));
      opts->clique_list_length += 512;
   }
   opts->clique_list[*clique_countPtr] = s;
   (*clique_countPtr)++;
   return 1;
}


// --------------------------------------------------------------------- //
class MAD_Cliquer{
 public:
   int              m_clique_count;
   graph_t        * m_g;
   clique_options * m_copts;

   int qualex_count;
   
 public:
   inline graph_t * graphNew(const int n_verts){
      return graph_new(n_verts);
   }

   inline void graphFree(graph_t * g){
      graph_free(g);
   }
   
   inline void copyGraphNonPos(graph_t       * gDest,
                               const graph_t * gSource,
                               const double  * weights){
      int i,j;
      CoinAssert(gSource->n == gDest->n);
      for(i = 0; i < gSource->n; i++) {
         if(weights[i] < -DecompEpsilon){
            for (j = 0; j < gSource->n; j++) {
               if(weights[j] < -DecompEpsilon){
                  if (SET_CONTAINS_FAST(gSource->edges[i], j)) {
                     addEdge(gDest, i, j);
                  }                  
               }
            }
         }   
      }
   }
   
   inline void addEdge(graph_t   * g,
                       const int   i, 
                       const int   j){
      GRAPH_ADD_EDGE(g, i, j);
   }   
   inline void delEdge(graph_t   * g,
                       const int   i, 
                       const int   j){
      GRAPH_DEL_EDGE(g, i, j);
   }   
   inline int getNumVertices(graph_t * g){
      return g->n;
   }
   inline int * getVertWeight(graph_t * g){
      return g->weights;
   }
   inline void setVertWeight(graph_t   * g,
                             const int * vertWeight){
      memcpy(g->weights, vertWeight, g->n * sizeof(int));
   }

   inline void cliqueFindAll(graph_t * g,
                             int       minWeight,
                             int       maxWeight,
                             int       maximal){
      clique_find_all(g, minWeight, maxWeight, maximal, m_copts);
      m_clique_count = *(int*)m_copts->user_data;
   }

   inline void cliqueUnweightedFindAll(graph_t * g,
                                       int       minWeight,
                                       int       maxWeight,
                                       int       maximal){
      clique_unweighted_find_all(g, minWeight, maxWeight, maximal, m_copts);
      m_clique_count = *(int*)m_copts->user_data;
   }

   inline void cliqueFindOne(graph_t * g,
                             int       minWeight,
                             int       maxWeight,
                             int       maximal){
      set_t s = clique_find_single(g, minWeight, maxWeight, maximal, m_copts);
      set_clique_func(s, g, m_copts);
      m_clique_count = *(int*)m_copts->user_data;
   }
   
   inline void cliquePrint(graph_t * g,
                           set_t     s) {
      unsigned int i;      
      printf("size=%d, weight=%d:  ",
             set_size(s), graph_subgraph_weight(g, s));
      for (i = 0; i < SET_MAX_SIZE(s); i++) {
         if (SET_CONTAINS(s,i)) {
            printf(" %d",i);
         }
      }
      printf("\n");
      return;
   }

   inline void cliquePopulateVector(const int     index,
                                    vector<int> & clique) {
      unsigned int i;
      set_t        s = m_copts->clique_list[index];
      for (i = 0; i < SET_MAX_SIZE(s); i++) {
         if (SET_CONTAINS(s,i)) {
            clique.push_back(i);
         }
      }
   }
   
   inline void cliquePrintAll(graph_t * g){
      int i;
      for(i = 0; i < *(int*)m_copts->user_data; i++) {
         cliquePrint(g, m_copts->clique_list[i]);
      }
   }
   
   inline void cliqueFreeMemory(){
      int i;
      for(i = 0; i < m_clique_count; i++) {         
         set_free(m_copts->clique_list[i]);
      }
      m_clique_count = 0;
      m_copts->user_data          = (void*) &m_clique_count;      
      m_copts->clique_list_length = 0;         
   }
   
   inline void printGraph(graph_t * g){
      return graph_print(g);
   }

   inline void cliqueFindOneQualex(vector<int> & ind){
     string gfilename = "dimacs_g.";
     string wfilename = "dimacs_w.";
     
     gfilename += UtilIntToStr(qualex_count);
     wfilename += UtilIntToStr(qualex_count);
     
#ifdef _MSC_VER
     string cmd = "C:\\cygwin\\home\\magala\\COIN\\coin-Decomp\\TempFix\\qualex-ms\\qualex-ms ";
#else
     string cmd = "~/COIN/coin-Decomp/TempFix/qualex-ms/qualex-ms ";
#endif
     cmd += gfilename + " -w" + wfilename;
     printf("\ncmd = %s", cmd.c_str());
     fflush(stdout);
     system(cmd.c_str());
     
     readDimacsSolution(ind);
     
     qualex_count++;     
   }

   inline void printGraphDimacs(graph_t * g){
  
     string filename = "dimacs_g.";
     filename += UtilIntToStr(qualex_count);
     ofstream os(filename.c_str());
     int i,j;
     os << "p edge " << g->n << " " << graph_edge_count(g) << endl;
     for(i = 0; i < g->n; i++) {       
       for (j = 0; j <= i; j++) {          
         if (SET_CONTAINS_FAST(g->edges[i], j)) {
           os << "e " << i+1 << " " << j+1 << endl;
         }                  
       }
     }     
     os.close();
   }

   inline void readDimacsSolution(vector<int> & ind){
     string filename = "dimacs_g.";
     filename += UtilIntToStr(qualex_count);
     filename += ".sol";
     ifstream is(filename.c_str());
     CoinAssert(is);
     char dummy[1000];
     is.getline(dummy,1000);
     is.getline(dummy,1000);
     int i;
     string dummyStr;
     while(!is.eof()){
       is >> dummyStr >> i;
       //printf("\ni = %d", i);
       //fflush(stdout);
       if(is.eof())
         break;
       ind.push_back(i);
     }
     is.close();
   }
   

   inline void printWeightDimacs(const int      n_verts,
                                 const double * weight){    
     string filename = "dimacs_w.";
     filename += UtilIntToStr(qualex_count);
     ofstream os(filename.c_str());
     int i;
     for(i = 0; i < n_verts; i++) {
       os << weight[i] << endl;
     }
     os.close();
   }
   
 public:
   MAD_Cliquer(const int n_verts) :
     m_clique_count(0),
     m_g           (0),      
     m_copts       (0),     
     qualex_count(0)
     {
         m_g     = graphNew(n_verts);         
         CoinAssertHint(m_g, "Error: Out of Memory");

         m_copts = clique_default_options;
         m_copts->user_data          = (void*) &m_clique_count;
         m_copts->user_function      = record_clique_func;
         m_copts->clique_list_length = 0;         
      }
   ~MAD_Cliquer() {
      cliqueFreeMemory();
      graphFree(m_g);    
   }
};


#endif
