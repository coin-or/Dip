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
#include "UtilMacros.h"
//===========================================================================//
#include "VRP_CVRPsep.h"
#include "VRP_DecompApp.h"
//===========================================================================//
#include "DecompAlgo.h"

//===========================================================================//
void VRP_DecompApp::initializeApp(UtilParameters & utilParam)  {

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "initializeApp()", m_appParam.LogLevel, 2);
   
   //---
   //--- get application parameters
   //---   
   m_appParam.getSettings(utilParam);   
   if(m_appParam.LogLevel >= 1)
      m_appParam.dumpSettings(m_osLog);

   //---
   //--- read VRPLIB instance
   //---
   string fileName = m_appParam.DataDir 
      + UtilDirSlash() + m_appParam.Instance + ".vrp";
   UtilGraphLib & graphLib = m_vrp.m_graphLib;
   graphLib.read_data(fileName.c_str());
   m_vrp.m_numRoutes = m_appParam.NumRoutes;

   //---
   //--- set pointer for CVRPsep object
   //---
   m_cvrpSep.init(&m_vrp);

   //---
   //--- init Concorde object
   //---
#ifdef VRP_DECOMPAPP_USECONCORDE
   m_concorde.init(&m_vrp);
#endif

   //---
   //--- read best known lb/ub 
   //---
   //---   Columns in vrplib.opt:
   //---    1. Problem Instance
   //---    2. # of Customers
   //---    3. # of Vehicles
   //---    4. Vehicle Capacity
   //---    5. Tightness (Demand/Capacity)
   //---    6. Gap %
   //---    7. Upper Bound
   //---
   string   bestKnownFile 
      = m_appParam.DataDir + UtilDirSlash() + "vrplib.opt";
   ifstream is;
   string   instanceName;
   int      numCustomers;
   int      numVehicles;
   int      capacity;
   double   tightness;
   double   gap;
   double   upperBound;
   UtilOpenFile(is, bestKnownFile);
   while(!is.eof()){
      is >> instanceName 
         >> numCustomers
         >> numVehicles
         >> capacity
         >> tightness
         >> gap
         >> upperBound;
      instanceName = UtilStrTrim(instanceName);
      if(instanceName == m_appParam.Instance){
         //gap = 100*(ub-lb)/ub
         //lb  = ub - gap*ub/100 
         m_bestKnownLB = upperBound - (gap * upperBound / 100.0);
         m_bestKnownUB = upperBound;
         printf("Instance = %s, BestLB = %g, BestUB = %g\n",
                instanceName.c_str(), m_bestKnownLB, m_bestKnownUB);
         break;
      }
   }
   is.close();

   //---
   //--- create models
   //---
   createModels();

   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "initializeApp()", m_appParam.LogLevel, 2);   
}

//===========================================================================//
void VRP_DecompApp::createModels(){
   //---
   //--- This function does the work to create the different models 
   //---  that will be used. This memory is owned by the user. It will
   //---  be passed to the application interface and used by the algorithms.
   //---

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createModels()", m_appParam.LogLevel, 2);
   
   //--- 
   //--- Vehicle Routing Problem (VRP)
   //---
   //---  min  sum{e in E}    c_e x_e
   //---  s.t. sum{e in delta(0)} x_e  = 2k                             (1)
   //---       sum{e in delta(i)} x_e  = 2,     forall i in  V          (2)
   //---       sum{e in delta(S)} x_e >= 2b(S), forall S sub N, |S| > 1 (3)
   //---       x_e = {0,1}, for all e = {i,j}, i,j!=0                   (4)
   //---       x_e = {0,2}, for all e = {0,j}                           (5)
   //---
   //--- Some notes on conventions: 
   //---   a.) j = 0 is depot 
   //---   b.) |E| = V(V-1)/2
   //---   c.) order of edges: i = 1 -> V-1, j = 0 -> i-1
   //---
   UtilGraphLib & graphLib = m_vrp.m_graphLib;
   int            n_cols   = graphLib.n_edges; 
   
   //---
   //--- create the cost vector c
   //---
   m_objective = new double[n_cols];
   if(!m_objective)
      throw UtilExceptionMemory("createModels", "VRP_DecompApp");
   copy(graphLib.edge_wt, graphLib.edge_wt + n_cols, m_objective);
   setModelObjective(m_objective);

   //--- 
   //--- Possible decompositions:
   //---
   //--- CPM:
   //---  A'' = 2-degree (on customers) and 2k-degree on depot 
   //---        with GSECs generated dynamically.
   //---
   //--- PC/RC:
   //---  A'  = Multiple Traveling Salesman Problem
   //---  A'' = 2-degree (on customers) and 2k-degree on depot 
   //---        with GSECs generated dynamically.
   //--- NOTE: In this case, 2-degree and 2k-degree is actually
   //---   redundant with MTSP polytope. All we really have to do 
   //---   is generated violated GSECs. But, starting with nothing
   //---   would cause a Catch-22 in framework. So, start with at 
   //---   least something. Maybe start with trivial GSECs? or just
   //---   the 2k-degree constraint?
   //---
   //TODO: k-tree, b-matching?
   // mTSP is actually a nesting of k-tree? could be interesting... 
   // mTSP is also a nesting of b-matching... 
   // but what makes vrp hard is the GSECs... ?
   // some relaxation we can try that includes GSECs? 
   //   solve GSECs with relaxed =2 with Decomp called recursively??
   if(m_appParam.ModelNameCore == "2DEGREE"){
      DecompConstraintSet * model = new DecompConstraintSet();
      createModelTwoDegree(model);   
      m_models.push_back(model);
      setModelCore(model, m_appParam.ModelNameCore);
   }
   if(m_appParam.ModelNameRelax == "2DEGREE"){
      DecompConstraintSet * model = new DecompConstraintSet();
      createModelTwoDegree(model);   
      m_models.push_back(model);
      setModelRelax(model, m_appParam.ModelNameRelax);
   }
   
   //current design - must set empty constraint set
   if(m_appParam.ModelNameRelax == "MTSP"){
      DecompConstraintSet * model = new DecompConstraintSet();
      m_models.push_back(model);
      setModelRelax(model, m_appParam.ModelNameRelax);
   }

   if(m_appParam.ModelNameRelax == "ESPPRCC"){
      DecompConstraintSet * model = new DecompConstraintSet();
      createModelESPPCC(model);   
      m_models.push_back(model);
      m_modelESPPRC = model;
      setModelRelax(model, m_appParam.ModelNameRelax);
   }
   
   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "createModels()", m_appParam.LogLevel, 2);
}

//===========================================================================//
void VRP_DecompApp::createModelTwoDegree(DecompConstraintSet * model){
  
   //---
   //--- s.t. sum{e in delta(0)} x_e  = 2k                             (1)
   //---      sum{e in delta(i)} x_e  = 2,     forall i in  V          (2)
   //---

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createModelTwoDegree()", m_appParam.LogLevel, 2);

   UtilGraphLib & graphLib   = m_vrp.m_graphLib;   
   const int      n_cols     = graphLib.n_edges;   
   const int      n_rows     = graphLib.n_vertices; //includes depot
   const int      n_vertices = graphLib.n_vertices;
   const int      n_edges    = graphLib.n_edges;

   //---
   //--- reserve some space for efficient fill-in of col and row bounds
   //---
   model->colLB.reserve(n_cols);
   model->colUB.reserve(n_cols);
   model->rowLB.reserve(n_rows);
   model->rowUB.reserve(n_rows);

   //---
   //--- set the column lower and upper bounds
   //---
   UtilFillN(model->colLB, n_cols, 0.0);
   UtilFillN(model->colUB, n_cols, 1.0); 

   //---
   //--- edges coming from depot can have value 2 (adjust upper bounds)
   //---
   int i, j;
   int depot_index;
   for(i = 1; i < n_vertices; i++){
      depot_index = UtilIndexU(0, i);
      model->colUB[depot_index] = 2.0;
   }

   //---
   //--- set the row lower and upper bounds (lb=ub=2)
   //---   the first row (vertex = 0) represents the depot (lb=ub=2k)
   //---
   fill_n(back_inserter(model->rowLB), n_rows, 2.0);
   fill_n(back_inserter(model->rowUB), n_rows, 2.0);   
   model->rowLB[0] = 2.0 * m_vrp.m_numRoutes;
   model->rowUB[0] = 2.0 * m_vrp.m_numRoutes;

   //---
   //--- two row entries (coeff = 1.0) for each column (edge = (u,v))
   //---     at row u and row v
   //---
   int      nEdgesX2    = 2 * n_edges;
   int    * rowIndices  = new int[nEdgesX2];
   int    * colIndices  = new int[nEdgesX2];
   double * elements    = new double[nEdgesX2];
   CoinAssertHint(rowIndices && colIndices && elements,
		  "Error: Out of Memory");
   
   int edge_index       = 0;
   int triplet_index    = 0;
   for(i = 1; i < n_vertices; i++){
      for(j = 0; j < i; j++){
	 rowIndices[triplet_index]      = i;
	 rowIndices[triplet_index + 1]  = j;
	 colIndices[triplet_index]      = edge_index;
	 colIndices[triplet_index + 1]  = edge_index;
	 triplet_index                 += 2;
	 edge_index++;
      }
   }  
   UtilFillN(elements, nEdgesX2, 1.0);

   //---
   //--- create a column-ordered CoinPackedMatrix
   //---
   model->M = new CoinPackedMatrix(true,
                                   rowIndices,
                                   colIndices,
                                   elements,
                                   nEdgesX2);

   //---
   //--- mark the integers
   //---
   UtilIotaN(model->integerVars, n_cols, 0);
   
   //---
   //--- free local memory
   //---
   UTIL_DELARR(rowIndices);
   UTIL_DELARR(colIndices);
   UTIL_DELARR(elements);
   
   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "createModelTwoDegree()", m_appParam.LogLevel, 2);
}

//===========================================================================//
int VRP_DecompApp::generateCuts(const double              * x, 
                                DecompCutList             & newCuts){

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "generateCuts()", m_appParam.LogLevel, 2);

   UtilGraphLib  & graphLib  = m_vrp.m_graphLib;
   const int       nEdges    = graphLib.n_edges;
   
   //---
   //--- calculate the number of nonzeros in vector x
   //---
   //CAREFUL: tolerance used here must match
   //  tolerance use in buildLpSol!
   int nNzs = UtilNumNonzeros(x, nEdges, DecompEpsilon);
   
   //---
   //--- build CVRPsep LP solution object
   //---
   m_cvrpSep.buildLpSol(x, nNzs,DecompEpsilon);

   //---
   //--- separate capacity cuts
   //---
   m_cvrpSep.sepCapacityCuts();
   
   //---
   //---  create DecompCuts from CVRPsep cuts
   //--- 
   m_cvrpSep.createVrpCuts(newCuts);

   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "generateCuts()", m_appParam.LogLevel, 2);      

   return newCuts.size();
}

//===========================================================================//
DecompSolverStatus VRP_DecompApp::solveRelaxed(const int          whichBlock,
                                               const double     * redCostX,
                                               const double       convexDual,
                                               DecompVarList    & varList){

   
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "solveRelaxed()", m_appParam.LogLevel, 2);
   
   int            i, status   = 0;
   double         coeff;
   int            index;
   vector<int>    vrpRouteInd;
   vector<double> vrpRouteEls;
   
   DecompSolverStatus solverStatus = DecompSolStatNoSolution;
   if(m_appParam.ModelNameRelax == "MTSP"){
#ifdef VRP_DECOMPAPP_USECONCORDE
      //---
      //--- translate the reduced cost to the expanded graph for MTSP
      //---	 

      m_concorde.setExpandedCost(redCostX);
      
      //---
      //--- solve the TSP
      //---	 
      status = m_concorde.solveTSP(vrpRouteInd, vrpRouteEls); 
      assert(!status);

      double varRedCost  = 0.0;
      double varOrigCost = 0.0;
      for(i = 0; i < static_cast<int>(vrpRouteInd.size()); i++){
         coeff        = vrpRouteEls[i];
         index        = vrpRouteInd[i];
         varRedCost  += coeff * redCostX[index];
         varOrigCost += coeff * m_objective[index];
      }
      UTIL_DEBUG(m_appParam.LogLevel, 5,
		 (*m_osLog) << "VAR varRedCost=" << varRedCost-convexDual;
		 (*m_osLog) << "varOrigCost=" << varOrigCost << endl;
		 );

      DecompVar * var = new DecompVar(vrpRouteInd,
                                      vrpRouteEls,
                                      varRedCost - convexDual,
                                      varOrigCost);
      var->setBlockId(0);
      varList.push_back(var);      
      solverStatus = DecompSolStatOptimal;
#endif
   }


   if(m_appParam.ModelNameRelax == "ESPPRCC"){
      const int          blockId         = 0;
      const DecompAlgo * decompAlgo      = getDecompAlgo();
      //DecompAlgoModel  & decompAlgoModel 
        // = decompAlgo->getModelRelax(blockId);




      //solve then map back and set to exact - cannot let 
      // built in solve it
      exit(1);
   }


   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "solveRelaxed()", m_appParam.LogLevel, 2);
   return solverStatus;
}

//===========================================================================//
bool VRP_DecompApp::APPisUserFeasible(const double * x, 
				      const int      n_cols,
				      const double   tolZero){

   //---
   //--- A feasible VRP solution:
   //---  a.) integral (assume it already is when it gets here)
   //---  b.) 2 degree at each node i in V
   //---  c.) 2k degree at depot
   //---  d.) remove the depot, then each component must satisfy capacity
   //---
   //TODO: back to retuning cuts from here? disconnected gsecs... 

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "APPisUserFeasible()", m_appParam.LogLevel, 2);

   UtilGraphLib & graphLib = m_vrp.m_graphLib;
   bool           feasible = true;
   
   //---
   //--- wipe out the old support graphs
   //---
   m_boost.clearSubGraph(m_boost.m_sg);
   m_boost.clearSubGraph(m_boost.m_sg0);

   //---
   //--- construct the current support graph to m_sg
   //---
   m_boost.buildSubGraph(n_cols, x);

   //---
   //---  b.) 2 degree at each node i in V
   //---  c.) 2k degree at depot (since we can have single customer
   //---      routes where x[i,0] = 2.0, we cannot just check degree
   //---      of depot in support graph)
   //---
   int       v;
   const int n_vertices = graphLib.n_vertices;
#if 0
   //this is always true because always part of LP....  think... 
   if(!UtilIsZero(m_boost.getDepotDegree() - (2 * m_vrp.m_numRoutes))){
      printf("depot degree is wrong: %g\n", m_boost.getDepotDegree());      
      UtilPrintFuncEnd(m_osLog, m_classTag,
		       "APPisUserFeasible()", m_appParam.LogLevel, 2);
      return false;
   }
   for(v = 1; v < n_vertices; v++){
      if(m_boost.getDegree(v) != 2){
	 //this is still not a good check in the case of 1-customer routes
	 printf("degree of %d is wrong: %d\n", v, m_boost.getDegree(v));
	 UtilPrintFuncEnd(m_osLog, m_classTag,
			  "APPisUserFeasible()", m_appParam.LogLevel, 2);
	 return false;
      }
   }
#endif

   //---
   //--- copy the support graph to m_sg0 and remove vertex 0
   //---
   m_boost.copyGraph(m_boost.m_sg, m_boost.m_sg0);
   m_boost.clearVertex(0, m_boost.m_sg0);

   UTIL_DEBUG(m_appParam.LogLevel, 5,
              (*m_osLog) << "m_sg:\n";
              m_boost.printGraph(m_boost.m_sg);
              (*m_osLog) << "m_sg:\n";
              m_boost.printGraph(m_boost.m_sg0);
              
              string             baseName   = "isFeas";
              const DecompAlgo * decompAlgo = getDecompAlgo();
              string fileNameDot = baseName
              + ".n" + UtilIntToStr(decompAlgo->getNodeIndex())
              + ".c" + UtilIntToStr(decompAlgo->getCutCallsTotal())
              + ".p" + UtilIntToStr(decompAlgo->getPriceCallsTotal()) + ".dot";
              m_boost.printDotFile(fileNameDot,
                                   graphLib.vertex_wt,
                                   m_boost.m_sg);
              );


   //---
   //---  d.) remove the depot, then each component must satisfy capacity
   //---
   vector<int> component(n_vertices);
   int c;
   int n_comp
      = m_boost.findConnectedComponents(m_boost.m_sg0, component);

   //---
   //--- NOTE:
   //---  we might have the case that even though solution is integral
   //---  and degree of all verts = 2 and depot degree = 2k but not all
   //---  vertices have a path to depot (i.e., we have subtours)
   //---
   //---  in this case, we will have more components than routes
   //---
   //---  since we find components of sg0, 0 will be a component itself
   //---  so we are looking for numRoutes+1 components
   //---
   if(n_comp != (m_vrp.m_numRoutes+1)){
      UTIL_DEBUG(m_appParam.LogLevel, 2,
                 (*m_osLog) 
                 << "not feasible, n_comp=" << n_comp 
                 << " n_routes="             << m_vrp.m_numRoutes << endl; 
                 );
      return false;
   }

   vector< vector<bool > > comp;
   for(c = 0; c < n_comp; c++){
      vector<bool> compV(n_vertices, false);
      comp.push_back(compV);
   }

   vector<int>    comp_count(n_comp, 0);
   vector<int>    comp_demand(n_comp, 0);
   //vector<double> comp_cut(n_comp, 0.0);
   for (v = 1; v < n_vertices; ++v){
      UTIL_DEBUG(m_appParam.LogLevel, 2,
                 (*m_osLog) << "component["      << v << "] = " << component[v]
                 << " comp_demand = " << comp_demand[component[v]]
                 << " vrp.vertex_wt[" << v << "] = " << graphLib.vertex_wt[v]
                 << endl;
                 );
      c               = component[v];
      comp[c][v]      = true;
      comp_demand[c] += graphLib.vertex_wt[v];
      comp_count[c]  ++;
      //comp_cut[c]    += solution[INDEX_U(0,v)];
   }

   
  for(c = 0; c < n_comp; c++){
     //each component's cut value must be either 0 or 2 since we are integral
     //if the comp_cut    = 0,        then we have a violated SEC
     //if the comp_demand > capacity, then we have a violated GSEC 

     //---
     //--- this is ok, it just means a vehicle only serves one client
     //---   TODO: sanity check, in this case should have x[0,j] = 2
     //---
     if(comp_count[c] <= 1) 
	continue;

     //this cannot happen given current core...?
     //if(comp_cut[c] < m_param.tolerance){
     //new_cuts.push_back(new VRP_GSECCut(comp[c], m_vrp.vertex_wt, 
     //m_vrp.capacity));
     //	feasible = false;
     //}

     if(comp_demand[c] > (graphLib.capacity + tolZero)){
	//new_cuts.push_back(new VRP_GSECCut(comp[c], m_vrp.vertex_wt, 
	//			   m_vrp.capacity));
        UTIL_DEBUG(m_appParam.LogLevel, 2,
                   (*m_osLog) 
                   << "demand of comp " << c << " is " << comp_demand[c]
                   << " and exceeds cap " << graphLib.capacity;
                   );
	UtilPrintFuncEnd(m_osLog, m_classTag,
			 "APPisUserFeasible()", m_appParam.LogLevel, 2);
	return false;
     }
  }

  return feasible;
}

//===========================================================================//
void VRP_DecompApp::printOriginalColumn(const int   index, 
                                        ostream   * os) const {
   UtilPrintEdge(index, os);
}
