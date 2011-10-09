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

// ------------------------------------------------------------------------- //
#include "TSP_Concorde.h"
#include "TSP_DecompApp.h"
#include "TSP_SubtourCut.h"

/*--------------------------------------------------------------------------*/
int TSP_DecompApp::generateCutsSubtour(DecompCutList & newCuts){

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "generateCutsSubtour()", m_param.LogDebugLevel, 2);
   
   //TODO: use access methods
   TSP_Concorde     & tspConcorde      = m_tsp.m_concorde;

   vector<ConcordeSubtourCut> subtourCuts;

   int c;
   int n_subtour  = tspConcorde.generateCutsSubtour(subtourCuts);
   int n_prevcuts = static_cast<int>(newCuts.size());
   
   for(c = 0; c < n_subtour; c++){
      vector<int>  & S   = subtourCuts[c].S;
      vector<bool> & inS = subtourCuts[c].inS;

      if(S.size() >= 2){
	 TSP_SubtourCut * sec_cut = new TSP_SubtourCut(inS, S);
	 
	 UTIL_DEBUG(m_param.LogDebugLevel, 3,
		    sec_cut->print();
		    );
	 newCuts.push_back(sec_cut);	 
      }
      else{
	 cout << "ERROR size of S < 2 (not adding)" << endl;
      }
   }

   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "generateCutsSubtour()", m_param.LogDebugLevel, 2);

   return static_cast<int>(newCuts.size()) - n_prevcuts;
}


#if 0
//TODO
/*--------------------------------------------------------------------------*/
int TSP_DecompApp::generateCutsBlossoms(vector<int>    & edge_list,
					vector<double> & edge_x,
					DecompCutList  & newCuts){

   //Timer timer_comb;
   CCtsp_lpcut_in * tsp_cuts = NULL;
   CCtsp_lpcut_in * tsp_cut;

   int n_comb = 0;
   int seed = (int) CCutil_real_zeit ();
   CCrandstate rstate;
   CCutil_sprand(seed, &rstate);
   int rval = CCtsp_exactblossom(&tsp_cuts, &n_comb, m_tsp.n_vertices, 
				 edge_x.size(), &edge_list[0], &edge_x[0], 
				 &rstate);

   if(m_param.LogLevel > 5)
      cout << "CONCORDE found blossoms : " << n_comb << endl;
   if (rval) {
      cerr << "CCtsp_exactblossom failed" << endl;
      abort();
   }

   tsp_cut = tsp_cuts;
   for(int c = 0; c < n_comb; c++){
      int cliquecount = tsp_cut->cliquecount;
    
      //i am assuming the first clique is the handle and the rest are the teeth
      int j, tmp;
      vector<int> H;
      dynamic_bitset<> inH(m_tsp.n_vertices);
      CC_FOREACH_NODE_IN_CLIQUE(j, tsp_cut->cliques[0], tmp){
	 H.push_back(j);
	 inH.set(j);
      }

      vector< vector<int> > teeth;
      vector< dynamic_bitset<> > inTeeth;
      for (int i = 1; i < cliquecount; i++){
	 vector<int> T;
	 dynamic_bitset<> inT(m_tsp.n_vertices);
	 CC_FOREACH_NODE_IN_CLIQUE(j, tsp_cut->cliques[i], tmp){
	    T.push_back(j);
	    inT.set(j);
	 }
	 teeth.push_back(T);
	 inTeeth.push_back(inT);
      }

      TSP_CombCut * comb_cut = new TSP_CombCut(inH, H, inTeeth, teeth);
      newCuts.push_back(comb_cut);
      tsp_cut = tsp_cut->next;
   }
   //double t_comb = timer_comb.Interval(Timer::Milliseconds)/1000.0;
   //cut_stats.insert(make_pair("comb   :",make_pair(n_comb,t_comb)));
   CCtsp_free_lpcut_in(tsp_cuts);
   tsp_cuts = NULL;
   return n_comb;
}

/*--------------------------------------------------------------------------*/
int TSP_DecompApp::generateCutsComb(vector<int>    & edge_list,
				    vector<double> & edge_x,
				    DecompCutList  & newCuts){

   //Timer timer_comb;
   CCtsp_lpcut_in * tsp_cuts = NULL;
   CCtsp_lpcut_in * tsp_cut;

   int n_comb = 0;
   int rval = CCtsp_block_combs(&tsp_cuts, &n_comb, m_tsp.n_vertices, 
				edge_x.size(), &edge_list[0], &edge_x[0], 1);

   if(m_param.LogLevel > 5)
      cout << "CONCORDE found n_combs : " << n_comb << endl;
   if (rval) {
      cerr << "CCtsp_block_combs failed" << endl;
      abort();
   }

   tsp_cut = tsp_cuts;
   for(int c = 0; c < n_comb; c++){
      int cliquecount = tsp_cut->cliquecount;
    
      //i am assuming the first clique is the handle and the rest are the teeth
      int j, tmp;
      vector<int> H;
      dynamic_bitset<> inH(m_tsp.n_vertices);
      CC_FOREACH_NODE_IN_CLIQUE(j, tsp_cut->cliques[0], tmp){
	 H.push_back(j);
	 inH.set(j);
      }

      vector< vector<int> > teeth;
      vector< dynamic_bitset<> > inTeeth;
      for (int i = 1; i < cliquecount; i++){
	 vector<int> T;
	 dynamic_bitset<> inT(m_tsp.n_vertices);
	 CC_FOREACH_NODE_IN_CLIQUE(j, tsp_cut->cliques[i], tmp){
	    T.push_back(j);
	    inT.set(j);
	 }
	 teeth.push_back(T);
	 inTeeth.push_back(inT);
      }

      TSP_CombCut * comb_cut = new TSP_CombCut(inH, H, inTeeth, teeth);
      newCuts.push_back(comb_cut);
      tsp_cut = tsp_cut->next;
   }
   //double t_comb = timer_comb.Interval(Timer::Milliseconds)/1000.0;
   //cut_stats.insert(make_pair("comb   :",make_pair(n_comb,t_comb)));
   CCtsp_free_lpcut_in(tsp_cuts);
   tsp_cuts = NULL;
   return n_comb;
}
#endif
