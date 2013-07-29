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

//===========================================================================//
#include "TSP_DecompApp.h"
#include "TSP_Concorde.h"
#include "TSP_Boost.h"
//===========================================================================//
#include "DecompVar.h"

//===========================================================================//
void TSP_DecompApp::initializeApp(UtilParameters & utilParam)  {

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "initializeApp()", m_appParam.LogLevel, 2);
   
   //---
   //--- get application parameters
   //---   
   m_appParam.getSettings(utilParam);   
   if(m_appParam.LogLevel >= 1)
      m_appParam.dumpSettings(m_osLog);

   //---
   //--- read TSPLIB instance
   //---
   string fileName = m_appParam.DataDir 
      + UtilDirSlash() + m_appParam.Instance + ".tsp";
   UtilGraphLib & graphLib = m_tsp.m_graphLib;
   graphLib.read_data(fileName.c_str());

   //---
   //--- read best known lb/ub 
   //---
   string bestKnownFile  = m_appParam.DataDir + UtilDirSlash();
   bestKnownFile        += ".." + UtilDirSlash();
   bestKnownFile        += "TSPLIB_opt";
   
   ifstream is;
   string   instanceName;
   double   bestBound, bestKnownLB, bestKnownUB;
   UtilOpenFile(is, bestKnownFile);
   while(!is.eof()){
      is >> instanceName >> bestBound;
      instanceName = UtilStrTrim(instanceName);
      if(instanceName == m_appParam.Instance){
         bestKnownLB = bestBound;
         bestKnownUB = bestBound;
         break;
      }
   }
   setBestKnownLB(bestKnownLB);
   setBestKnownUB(bestKnownUB);

   //---
   //--- build complete graph over V \ {m_vert} = m_boost.m_cgV
   //---   for use with one-tree solver
   //---
   m_tsp.m_boost.buildCompleteGraphMinusVert(m_tsp.m_vert, 
					     graphLib.n_vertices);

   //---
   //--- create models
   //---
   createModels();

   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "initializeApp()", m_appParam.LogLevel, 2);   
}

//===========================================================================//
void TSP_DecompApp::createModels(){

   //---
   //--- This function does the work to create the different models 
   //---  that will be used. This memory is owned by the user. It will
   //---  be passed to the application interface and used by the algorithms.
   //---
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createModels()", m_appParam.LogLevel, 2);
     
   //--- 
   //--- Symmetric Travelling Salesman Problem (TSP)  
   //--- min  sum{e in E}     c_e x_e
   //--- s.t. x(delta(v))   = 2, for v in V       (1) 
   //---      x(delta(S))  >= 2, for S subseteq V (2) (gen-dynamic)
   //---      x in {0,1},        for all e in E   (3) 
   //---    
   UtilGraphLib & graphLib = m_tsp.m_graphLib;
   int            numCols  = graphLib.n_edges; 
   
   //---
   //--- create the cost vector c
   //---
   m_objective = new double[numCols];
   if(!m_objective)
      throw UtilExceptionMemory("createModels", "TSP_DecompApp");
   copy(graphLib.edge_wt, graphLib.edge_wt + numCols, m_objective);
   setModelObjective(m_objective);

   //---
   //--- Two matching relaxation.
   //---  core  = { some subset of (2) }
   //---  relax = { (1) and (3) }         => 2-matching relaxation
   //---   NOTE: for validity, must generate rest of (2)
   //---

   //---
   //--- One tree relaxation.
   //---  core  = { 1 }
   //---  relax = { one-tree TODO write algebra } => 1-tree relaxation
   //---   NOTE: for validity, must generate rest of (2)
   //---
   if(m_appParam.ModelNameCore == "SUBTOUR"){
      DecompConstraintSet * model = new DecompConstraintSet();  
      createModelTrivialSEC(model);
      m_models.push_back(model);
      setModelCore(model, m_appParam.ModelNameCore);      
   }
   if(m_appParam.ModelNameRelax == "2MATCH" ||
      m_appParam.ModelNameCore  == "2MATCH"){
      DecompConstraintSet * model = new DecompConstraintSet();  
      createModel2Match(model);
      m_models.push_back(model);
      if(m_appParam.ModelNameRelax == "2MATCH"){
         assert(m_appParam.ModelNameCore == "SUBTOUR");
	 setModelRelax(model, m_appParam.ModelNameRelax);	 
      }
      else{
         assert(m_appParam.ModelNameRelax == "SUBTOUR");
	 setModelRelax(NULL);
         setModelCore (model, m_appParam.ModelNameCore);
      }
   }  
   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "createModels()", m_appParam.LogLevel, 2);
}

//===========================================================================//
void TSP_DecompApp::createModel2Match(DecompConstraintSet * modelCS){
  
   //---
   //--- s.t. x(delta(v))   = 2, for v in V       (1) [A', b']
   //---      x in {0,1},        for all e = E    (3) [A', b']
   //---
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createModel2Match()", m_appParam.LogLevel, 2);

   UtilGraphLib & graphLib = m_tsp.m_graphLib;
   int      n_cols         = graphLib.n_edges; 
   int    * rowIndices     = new int[2 * n_cols];
   int    * colIndices     = new int[2 * n_cols];
   double * elements       = new double[2 * n_cols];
   int edge_index    = 0;
   int triplet_index = 0;
   for(int u = 1; u < graphLib.n_vertices; u++){
      for(int v = 0; v < u; v++){
	 rowIndices[triplet_index]     = u;
	 rowIndices[triplet_index + 1] = v;
	 colIndices[triplet_index]     = edge_index;
	 colIndices[triplet_index + 1] = edge_index;
	 triplet_index += 2;
	 edge_index++;
      }
   }  
   UtilFillN(elements,             2 * n_cols,           1.0);
   UtilFillN(modelCS->colLB,       n_cols,               0.0);
   UtilFillN(modelCS->colUB,       n_cols,               1.0);
   UtilFillN(modelCS->rowLB,       graphLib.n_vertices,  2.0);
   UtilFillN(modelCS->rowUB,       graphLib.n_vertices,  2.0);
   UtilIotaN(modelCS->integerVars, n_cols,               0);
  
   modelCS->M = new CoinPackedMatrix(true, rowIndices, colIndices,
				     elements, 2 * n_cols);
  
   UTIL_DELARR(rowIndices);
   UTIL_DELARR(colIndices);
   UTIL_DELARR(elements); 

   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "createModel2Match()", m_appParam.LogLevel, 2);
}

//===========================================================================//
void TSP_DecompApp::createModelTrivialSEC(DecompConstraintSet * modelCS){

   //--- 
   //--- Since we generate SECs (2) dynamically, we will simply start
   //---  with some trivial SECs (|S| = 2).
   //---    ?? todo - make it |S|=2, actually using redudant |S|=1
   //---    ?? what you have below is actually a relaxation of (1)!! >= 
   int u, v;
   UtilGraphLib & graphLib = m_tsp.m_graphLib;
   int            n_cols   = graphLib.n_edges; 

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createModelTrivialSEC()", m_appParam.LogLevel, 2);

   modelCS->M = new CoinPackedMatrix(false, 0, 0);
   for(u = 0; u < graphLib.n_vertices; u++){
      CoinPackedVector row;
      for(v = 0; v < graphLib.n_vertices; v++){
	 if(u == v)
	    continue;
	 row.insert(UtilIndexU(u,v), 1.0);
      }
      modelCS->rowLB.push_back(2.0);
      modelCS->rowUB.push_back(DecompInf);
      modelCS->M->appendRow(row);
   }
   UtilFillN(modelCS->colLB,       n_cols,  0.0);
   UtilFillN(modelCS->colUB,       n_cols,  1.0);
   UtilIotaN(modelCS->integerVars, n_cols,    0);

   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "createModelTrivialSEC()", m_appParam.LogLevel, 2);
}

//===========================================================================//
int TSP_DecompApp::generateCuts(const double              * x, 
                                DecompCutList             & newCuts){
  

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "generateCuts()", m_appParam.LogLevel, 2);
   
   int                n_cuts           = 0;
   UtilGraphLib     & graphLib         = m_tsp.m_graphLib;
   TSP_Concorde     & tspConcorde      = m_tsp.m_concorde;   
   tspConcorde.buildSubGraph(graphLib.n_vertices,
			     graphLib.n_edges, x);

   if(m_appParam.CutSubtoursX)
      n_cuts += generateCutsSubtour(newCuts);

   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "generateCuts()", m_appParam.LogLevel, 2);   
   return n_cuts;   
}

//===========================================================================//
DecompSolverStatus 
TSP_DecompApp::solveRelaxed(const int          whichBlock,
                            const double     * redCostX,
                            const double       convexDual,
                            DecompVarList    & varList){
   
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "solveRelaxed()", m_appParam.LogLevel, 2);

   DecompSolverStatus   solverStatus = DecompSolStatNoSolution;
   UtilGraphLib       & graphLib     = m_tsp.m_graphLib;
   Graph              & cgV          = m_tsp.m_boost.m_cgV;   
   int                  n_vertices   = graphLib.n_vertices;
   int                  u, index;

   //TODO: BranchEnforceInSubProb option?

   if(m_appParam.ModelNameRelax == "SUBTOUR"){
      vector< pair<int, double> > edge_cost;
      edge_cost.reserve(n_vertices);	 
      for(u = 0; u < n_vertices; u++){
         if(u != m_tsp.m_vert){
            index = UtilIndexU(m_tsp.m_vert, u);
            edge_cost.push_back(make_pair(index, redCostX[index]));
         }
      }   
      solveOneTree(redCostX, convexDual, edge_cost, varList, cgV);      
      solverStatus = DecompSolStatOptimal;
   }

   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "solveRelaxed()", m_appParam.LogLevel, 2);
   return solverStatus;
}

//===========================================================================//
void TSP_DecompApp::solveOneTree(const double               * cost, 
                                 const double                 alpha,
                                 vector< pair<int,double> > & edge_cost,
                                 DecompVarList              & vars,
                                 Graph                      & g) {

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "solveOneTree()", m_appParam.LogLevel, 2);

   property_map<Graph, edge_index_t>::type e_index_g   = get(edge_index,  g);
   property_map<Graph, edge_weight_t>::type e_weight_g = get(edge_weight, g);

   UtilGraphLib & graphLib = m_tsp.m_graphLib;
   const int max_exchanges = 4; 

   //---
   //--- (partial) sort in increasing order
   //---
   int sort_size = std::min<int>(static_cast<int>(edge_cost.size()), 
				 2 + max_exchanges);
   partial_sort(edge_cost.begin(), 
		edge_cost.begin() + sort_size,
		edge_cost.end(), UtilIsLessThan<int,double>());

   if(m_appParam.LogLevel > 4){
      (*m_osLog) << "Partial Sorted List [size = " << sort_size << "]" << endl;
      for(vector< pair<int,double> >::iterator tmp = edge_cost.begin();
	  tmp != edge_cost.end(); tmp++){
	 (*m_osLog) << "\nsorted edge_cost: " << (*tmp).first; 
	 UtilPrintEdge((*tmp).first);
	 (*m_osLog) << " cost: " << (*tmp).second;
      }
      (*m_osLog) << "\n";
   }  
   
   //---
   //--- update the edge weights in boost object
   //---
   int index;
   graph_traits<Graph>::edge_iterator ei, ei_end;
   for(tie(ei,ei_end) = edges(g); ei != ei_end; ei++){
      index = e_index_g(*ei);
      e_weight_g[*ei] = cost[index];
      if(m_appParam.LogLevel > 4)
         (*m_osLog) << "cost edge( " 
                    << source(*ei,g) << "," << target(*ei,g) << "): " 
                    << cost[index];
      assert(index == UtilIndexU(static_cast<int>(source(*ei,g)),
				 static_cast<int>(target(*ei,g))));
   }
   

   
   int n_edges    = graphLib.n_edges;
   int n_vertices = graphLib.n_vertices;
   vector<bool> inMST(n_edges, false);
   
   vector<int>    inds; 
   vector<double> els(n_vertices, 1.0);
   double rc  = 0.0;
   double obj = 0.0;
   inds.reserve(graphLib.n_vertices);

   //TODO: this should all be in boostGraphI
   //---
   //--- find the minimum spanning tree
   //---

   //boost::print_graph(g);
   
   vector<graph_traits < Graph >::edge_descriptor> spanning_tree;
   vector<graph_traits < Graph >::edge_descriptor>::iterator vei;
   kruskal_minimum_spanning_tree(g, back_inserter(spanning_tree));

   if(m_appParam.LogLevel > 4)
      (*m_osLog) << "Spanning Tree:" << endl;

   int edge_index;
   for (vei  = spanning_tree.begin(); vei != spanning_tree.end(); ++vei) {
      edge_index         = e_index_g[*vei];      
      rc                += cost[edge_index];
      obj               += graphLib.edge_wt[edge_index];      
      inMST[edge_index]  = true;
      inds.push_back(edge_index);

      if(m_appParam.LogLevel > 4){
	 UtilPrintEdge(edge_index);
	 (*m_osLog) << " -> " << cost[edge_index] << " rc : " << rc << endl; 
      }

   }
   
   const double bigM     = DecompInf;   
   int          exchange = 0;
   vector< pair<int,double> >::iterator vpi = edge_cost.begin();
   
   //---
   //--- add the cheapest edge to vert
   //---
   inds.push_back((*vpi).first);
   rc  += (*vpi).second;
   obj += graphLib.edge_wt[(*vpi).first];
   
   if(m_appParam.LogLevel > 4){
      (*m_osLog) << "Adding edge:" << endl;
      UtilPrintEdge((*vpi).first);
      (*m_osLog) << " -> " << cost[(*vpi).first] << " rc : " << rc << endl;
   }

   vpi++;
   for(; vpi != edge_cost.end(); vpi++){
      if(exchange >= max_exchanges)
	 break;
      if(cost[(*vpi).first] >= bigM/2.0)
	 break;

      //---
      //--- add the 2nd cheapest and keep exchanging this one
      //---
      inds.push_back((*vpi).first);    
      rc  += (*vpi).second;
      obj += graphLib.edge_wt[(*vpi).first];

      if(m_appParam.LogLevel > 4){
	 (*m_osLog) << "Adding edges:" << endl;
	 UtilPrintEdge((*vpi).first);
	 (*m_osLog) << " -> " << cost[(*vpi).first] << " rc : " << rc << endl;
	 (*m_osLog) << "Creating var with reduced = " << rc - alpha 
                    << " obj = " << obj << endl;
      }

      DecompVar * oneTree = new DecompVar(inds, els, rc - alpha, obj);
      //oneTree->setBlockId(0);//this will happen by default
      vars.push_back(oneTree);

      exchange++;
      inds.pop_back();
      rc  -= (*vpi).second;
      obj -= graphLib.edge_wt[(*vpi).first];
   }

   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "solveOneTree()", m_appParam.LogLevel, 2);
}

//===========================================================================//
bool TSP_DecompApp::APPisUserFeasible(const double * x, 
				      const int      n_cols,
				      const double   tolZero){

   //---
   //--- A feasible TSP route:
   //---   a.) is binary (assume it already is when it gets here)
   //---   b.) is connected
   //---   c.) all nodes have degree 2
   //---

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "APPisUserFeasible()", m_appParam.LogLevel, 2);
   
   const int n_vertices = m_tsp.m_graphLib.n_vertices;
   TSP_Boost & tspBoost = m_tsp.m_boost;
   
   //wipe out the old support graph
   tspBoost.clearSubGraph();
   
   //construct the current support graph
   tspBoost.buildSubGraph(n_cols, x);
   
   //find the connected components of support_graph
   vector<int> component(n_vertices);
   int n_comp = tspBoost.findConnectedComponents(component);
   
   
   //STOP? is binary already?
   //if (b) n_comp == 1, and (a) binary, then construct a feasible solution
   if(n_comp == 1){

      //(c) all nodes have degree 2
      //  NOTE: for CPM with 2-matching constraints this should already 
      //        be true already
      for(int i = 0; i < n_vertices; i++){
	 if(tspBoost.getDegree(i) != 2){
	    return false;
	 }
      }
      return true;
      
      /*
        if(isBinary(x, n_cols, tolZero)){
        constructRoute(m_sg); 
        return true; 
        }
        else{
        return false; //it is connected, but fractional
        }
      */
   }

   //TODO: feasibility cuts?
   //(*m_osLog) << "Not Feasible: disconnected, n_comp : " << n_comp << endl;

   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "APPisUserFeasible()", m_appParam.LogLevel, 2);

   return false;
}

//===========================================================================//
void TSP_DecompApp::printOriginalColumn(const int   index, 
                                        ostream   * os) const {
   UtilPrintEdge(index, os);
}
