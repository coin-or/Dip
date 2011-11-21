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

//===========================================================================//
//  SWITCHED DISPATCH PROBLEM WITH UNIT COMMITMENT
//===========================================================================//
#include "DecompVar.h"
#include "SDPUC_DecompApp.h"

//===========================================================================//
void SDPUC_DecompApp::initializeApp(UtilParameters & utilParam) {
      
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "initializeApp()", m_appParam.LogLevel, 2);

   //---
   //--- get application parameters
   //---
   m_appParam.getSettings(utilParam);
   if(m_appParam.LogLevel >= 1)
      m_appParam.dumpSettings();

   //---
   //--- read problem instance
   //---   
   string instanceFile   = m_appParam.DataDir 
      + UtilDirSlash() + m_appParam.Instance;
   int rc = m_instance.readInstance(instanceFile, false);
   if(rc)
      throw UtilException("Error in readInstance",
                          "initializeApp", "MCF_DecompApp");
   //---
   //--- create models
   //---
   createModels();
   
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "initializeApp()", m_appParam.LogLevel, 2);
}

//===========================================================================//
void SDPUC_DecompApp::createModels(){

   //---
   //--- This function does the work to create the different models 
   //---  that will be used. This memory is owned by the user. It will
   //---  be passed to the application interface and used by the algorithms.
   //---
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createModels()", m_appParam.LogLevel, 2);

   //---
   //---  Switched Dispatch Problem with Unit Commitment (SDPUC).
   //---
   //--- We are given: 
   //---    (1) a directed graph G=(N,A),
   //---    (2) a set of time periods T, 
   //---
   //--- min  sum{(i,j) in A} f1[i,j] y1[i,j] 
   //---		  + sum{t in T} sum{(i,j) in A} f2[i,j] y2[i,j,t] + c[i,j,t] x[i,j,t]
   //--- s.t. sum{(j,i) in A} x[i,j,t] - 
   //---        sum{(i,j) in A} x[i,j,t] = d[i,t],  for all i in N, t in T
   //---      x[i,j,t] >= l[i,j,t] z[i,j,t],       for all (i,j) in A, t in T
   //---      x[i,j,t] <= u[i,j,t] z[i,j,t],       for all (i,j) in A, t in T
   //---	  r[i,j] x[i,j,t] - theta[j] + theta[i] <=  M (1 - z[i,j,t]) for all i,j,t
   //---	  r[i,j] x[i,j,t] - theta[j] + theta[i] >= -M (1 - z[i,j,t]) for all i,j,t
   //---	  z[i,j,t] <= y1[i,j]  for all (i,j) in A, t in T				  //arc-investment
   //---	  z[i,j,t] - z[i,j,t-1] <= y2[i,j,t]  for all (i,j) in A, t in T  //arc(unit) commitment
   //---	  y[i,j] binary for all (i,j) in A
   //---
   //--- NOTE: to make sure the problem is always feasible, 
   //----      demand may have to be modelled as arcs with large negative costs
   //---
   //---
   //--- The decomposition is formed as:
   //---
   //--- MASTER (A''):
   //---	  z[i,j,t] <= y1[i,j]  for all (i,j) in A, t in T				  //arc-investment
   //---	  z[i,j,t] - z[i,j,t-1] <= y2[i,j,t]  for all (i,j) in A, t in T  //arc(unit) commitment
   //---	  y[i,j] binary for all (i,j) in A
   //---
   //--- SUBPROBLEM (A'): (one block for each t in T)
   //---     sum{(j,i) in A} x[i,j,t] - 
   //---        sum{(i,j) in A} x[i,j,t] = d[i,t],  for all i in N
   //---      x[i,j,t] >= l[i,j,t] z[i,j,t],       for all (i,j) in A
   //---      x[i,j,t] <= u[i,j,t] z[i,j,t],       for all (i,j) in A
   //---	  r[i,j] x[i,j,t] - theta[j] + theta[i] <=  M (1 - z[i,j,t]) for all i,j
   //---	  r[i,j] x[i,j,t] - theta[j] + theta[i] >= -M (1 - z[i,j,t]) for all i,j

   //---

   //---
   //--- Get information about this problem instance.
   //---
   int   i, t, a, colIndex;
   int   numTimeperiods = m_instance.m_numTimeperiods;
   int   numArcs        = m_instance.m_numArcs;
   int	 numNodes		= m_instance.m_numNodes;
   int   numCols        = numArcs							//y1-vars
      + 3 * numTimeperiods * numArcs	//y2-, z-, and x-vars
      + numTimeperiods * numNodes;		//theta-vars
   SDPUC_Instance::arc * arcs = m_instance.m_arcs;
   SDPUC_Instance::timeseries * ts = m_instance.m_timeseries;
   cout << "\nnumCols=" << numCols << " numTimePeriods=" << numTimeperiods;
   cout << "numNodes=" << numNodes << " numArcs=" << numArcs << endl;

   //---
   //--- Construct the objective function and set it
   //---    y1-var columns indexed as [a]   = a  in [0 ; numArcs-1]
   //---	y2-var columns indexed as [a,t]  = a + numArcs in [numArcs ; numArcs * (1 + numTimeperiods) - 1]
   //---	z-var columns indexed as [a,t] = t*numArcs + a +  numArcs * (1 + numTimeperiods) 
   //---										in [numArcs*(1+numTimeperiods); numArcs*(1 + 2*numTimeperiods) - 1]
   //---	x-var columns indexed as [a,t] = t*numArcs + a + numArcs*(1 + 2*numTimeperiods) , 
   //---										in [numArcs*(1 + 2*numTimeperiods) ; numArcs*(1 + 3*numTimeperiods) - 1]
   //---	theta-var columns indexed as [i,t] = t*numNodes + i + numArcs*(1 + 3*numTimeperiods) , 
   //---										   in [numArcs*(1 + 3*numTimeperiods) ; numArcs*(1 + 3*numTimeperiods) + numNodes*numTimeperiods - 1]
   //--
   m_objective = new double[numCols];
   //initialise to 0
   for(i = 0; i < numCols; i++){
      m_objective[i] = 0;
   }
   if(!m_objective)
      throw UtilExceptionMemory("createModels", "MCF_DecompApp");
   colIndex = 0;
   for(a = 0; a < numArcs; a++) {
      m_objective[colIndex++] = arcs[a].fcost1;  //fixed arc investment cost 
   }
   for(t = 0; t < numTimeperiods; t++) {
      for(a = 0; a < numArcs; a++) {
	 m_objective[colIndex++] = arcs[a].fcost2;  //fixed arc "start-up" cost 
      }
   }
   colIndex = numArcs*(1 + 2*numTimeperiods);  //start-index for x-vars
   for(t = 0; t < numTimeperiods; t++) {
      for(a = 0; a < numArcs; a++) {
         m_objective[colIndex++] = arcs[a].mcost * ts[0].values[t] ;  //arc cost * probability (assume ts[0] indicate timeperiod probabilities)
      }
   }	 
   //---
   //--- set the objective 
   //---        
   setModelObjective(m_objective);
   /*cout << "obj = " ;
     for(i = 0; i < numCols; i++){
     cout << m_objective[i] << " ";
     }
     cout << endl;*/
   
   //---
   //--- create the core/master model and set it
   //---
   DecompConstraintSet * modelCore = new DecompConstraintSet();      
   createModelCore(modelCore);
   setModelCore(modelCore, "core");

   //---
   //--- create the relaxed/subproblem models and set them
   //---
   for(t = 0; t < numTimeperiods; t++){
      DecompConstraintSet * modelRelax = new DecompConstraintSet();
      string                modelName  = "relax" + UtilIntToStr(t);
      if(m_appParam.UseSparse)
         createModelRelaxSparse(modelRelax, t);
      else
         createModelRelax(modelRelax, t);
      
      setModelRelax(modelRelax, modelName, t);
   }

   //---
   //--- create an extra "empty" block for the master-only vars
   //---   since I don't know what OSI will do with empty problem
   //---   we will make column bounds explicity rows
   //---
   int nMasterOnlyCols = static_cast<int>(modelCore->masterOnlyCols.size());
   if(nMasterOnlyCols){
      if(m_appParam.LogLevel >= 1)
         (*m_osLog) << "Create model part Master-Only." << endl;

      createModelMasterOnlys(modelCore->masterOnlyCols);
   }
   
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "createModels()", m_appParam.LogLevel, 2);
}

//===========================================================================//
void SDPUC_DecompApp::createModelCore(DecompConstraintSet * model){

   //---
   //--- MASTER (A''):
   //---	  z[i,j,t] <= y1[i,j]  for all (i,j) in A, t in T		  //arc-investment
   //---	  z[i,j,t] - z[i,j,t-1] <= y2[i,j,t]  for all (i,j) in A, t in T  //arc(unit) commitment
   //---	  y[i,j] binary for all (i,j) in A
   //---                                    
   int   i, t, a, colIndex;
   int   numTimeperiods =		 m_instance.m_numTimeperiods;
   int   numArcs        =		 m_instance.m_numArcs;
   int   numNodes		=		 m_instance.m_numNodes;
   int   numCols        =		 numArcs   //y1-vars
      + 3 * numTimeperiods * numArcs	           //y2-, z-, and x-vars
      + numTimeperiods * numNodes;		   //theta-vars
   int   numRows		  =		 numArcs * numTimeperiods;
   int	 col_yStartIndex  =		 0;
   int	 col_zStartIndex  =		 numArcs*(1+numTimeperiods);
   int	 col_xStartIndex  =		 numArcs * (1 + 2*numTimeperiods);
   int	 col_thetaStartIndex  =	 numArcs*(1 + 3*numTimeperiods) ;

   SDPUC_Instance::arc * arcs = m_instance.m_arcs;

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createModelCore()", m_appParam.LogLevel, 2);

   //---
   //--- create space for the model matrix (row-majored)
   //---
   model->M = new CoinPackedMatrix(false, 0.0, 0.0);
   if(!model->M)
      throw UtilExceptionMemory("createModelCore", "SDPUC_DecompApp");
   model->M->setDimensions(0, numCols);
   model->reserve(numRows, numCols);
   
   //---
   //--- create the rows and set the col/row bounds
   //---
   UtilFillN(model->colLB, numCols,  -DecompInf);
   UtilFillN(model->colUB, numCols,  DecompInf);
   
   for(a = 0; a < numArcs; a++){
      colIndex = a;
      //add y1-vars   
      model->colLB[colIndex] = 0;
      model->colUB[colIndex] = 1;
	 
      for(t = 0; t < numTimeperiods; t++){
	 int t_prev = 0;
	 if(t == 0) {
            t_prev = numTimeperiods - 1;
         }
	 else {
            t_prev = t - 1;
         }

	 colIndex = a + t * numArcs + numArcs;
	 //add y2-vars   
	 model->colLB[colIndex] = 0;
	 model->colUB[colIndex] = 1;

	 CoinPackedVector row1;        // 0 <= y1[i,j] - z[i,j,t]  for all (i,j) in A, t in T
	 CoinPackedVector row2;        // 0 <= y2[i,j,t] - z[i,j,t] + z[i,j,t-1]  for all (i,j) in A, t in T 
	 CoinPackedVector rowFix_y1;   //fix y1=1
	 CoinPackedVector y2upper1;   // y2(t) <= z(t)  : if off in t then we dont start-up
	 CoinPackedVector y2upper2;   // y2(t) <= 1-z(t-1)  : if on in t-1 then dont start up in t

	 //insert y1-var coefficient
	 row1.insert(a, 1.0);
	 rowFix_y1.insert(a, 1.0);
	 //insert y2-var coefficient
	 colIndex = t * numArcs + a + numArcs;
	 row2.insert(colIndex, 1.0);
	 y2upper1.insert(colIndex, -1.0);
	 y2upper2.insert(colIndex, -1.0);
	 //add z-vars   
	 colIndex = t * numArcs + a + col_zStartIndex;
	 model->colLB[colIndex] = 0;
	 model->colUB[colIndex] = 1;
	 //insert z-var coefficient
	 row1.insert(colIndex, -1.0);
	 row2.insert(colIndex, -1.0);
	 y2upper1.insert(colIndex, 1.0);
	 colIndex = t_prev * numArcs + a + col_zStartIndex;
	 row2.insert(colIndex, 1.0);
	 y2upper2.insert(colIndex, -1.0);

         std::string rowName1_1 = "MP1_1_" + UtilIntToStr(a) + "_" + UtilIntToStr(t);
         std::string rowName1_2 = "MP1_2_" + UtilIntToStr(a) + "_" + UtilIntToStr(t);
         std::string rowName2_1 = "MP2_1_" + UtilIntToStr(a) + "_" + UtilIntToStr(t);
         std::string rowName2_2 = "MP2_2_" + UtilIntToStr(a) + "_" + UtilIntToStr(t);
         std::string rowNameFix = "fix_y1_" + UtilIntToStr(a) + "_" + UtilIntToStr(t);
 
	 //TODO: any issue with range constraint
	 model->appendRow(row1, 0.0, DecompInf, rowName1_1);  //add MP1_constraints (arc investments)
         
         model->appendRow(row1, -DecompInf, 1.0, rowName1_2);  //add MP1_constraints (arc investments)
	 if(arcs[a].tail == 0) {   //ONLY for supply arcs (!!)
	    model->appendRow(row2, 0.0, DecompInf, rowName2_1);  //add MP2_constraints (arc commitment) 
	    model->appendRow(row2, -DecompInf, 1.0, rowName2_2);  //add MP2_constraints (arc commitment) 
	 }
	 model->appendRow(rowFix_y1, 1.0, DecompInf, rowNameFix);  //add fix y1 vars
	 //model->appendRow(y2upper1, 0.0, DecompInf, std::string("y2-upperbound-1"));  //add upperbounds on y2
	 //model->appendRow(y2upper2, -1.0, DecompInf, std::string("y2-upperbound-2"));  //..to strengthen formulation
      }
   }

   //---
   //--- create column names (helps with debugging)
   //---
   //y-vars
   for(a = 0; a < numArcs; a++){
      std::string colName = "y1(a" + UtilIntToStr(a) + "(" +
	 UtilIntToStr(arcs[a].tail) + "," +
	 UtilIntToStr(arcs[a].head) + "))";
      model->colNames.push_back(colName);
   }
   for(t = 0; t < numTimeperiods; t++){
      for(a = 0; a < numArcs; a++){
	 std::string colName = "y2(t" + UtilIntToStr(t) + ",a" + UtilIntToStr(a) + "(" +
            UtilIntToStr(arcs[a].tail) + "," +
            UtilIntToStr(arcs[a].head) + "))";
	 model->colNames.push_back(colName);
      }
   }
   //z-vars
   for(t = 0; t < numTimeperiods; t++){
      for(a = 0; a < numArcs; a++){
         std::string colName = "z(t" + UtilIntToStr(t) + ",a" + UtilIntToStr(a) + "(" +
            UtilIntToStr(arcs[a].tail) + "," +
            UtilIntToStr(arcs[a].head) + "))";
         model->colNames.push_back(colName);
      }
   }
   //x-vars
   for(t = 0; t < numTimeperiods; t++){
      for(a = 0; a < numArcs; a++){
         std::string colName = "x(t" + UtilIntToStr(t) + ",a" + UtilIntToStr(a) + "(" +
            UtilIntToStr(arcs[a].tail) + "," +
            UtilIntToStr(arcs[a].head) + "))";
         model->colNames.push_back(colName);
      }
   }
   //theta-vars
   for(t = 0; t < numTimeperiods; t++){
      for(i = 0; i < numNodes; i++){
         std::string colName = "theta(t" + UtilIntToStr(t) + ",n" + UtilIntToStr(i) + ")";
         model->colNames.push_back(colName);
      }
   }

   //---
   //--- create a list of the "active" columns (those related 
   //---   to this commmodity) all other columns are fixed to 0
   //---
   UtilFillN(model->colLB, numCols,  0.0);
   UtilFillN(model->colUB, numCols,  0.0);
   colIndex = 0;
   for(a = 0; a < numArcs; a++){
      //set y-columns active
	  
      //model->colLB[colIndex] = 0;
      //model->colUB[colIndex] = 1;
      //model->activeColumns.push_back(colIndex);
      //set y-columns as master-only columns
      colIndex = col_yStartIndex + a;
      model->masterOnlyCols.push_back(colIndex);  //y1-vars
      for(t = 0; t < numTimeperiods; t++){
	 colIndex = col_yStartIndex + t * numArcs + a + numArcs;
	 model->masterOnlyCols.push_back(colIndex);  //y2-vars
      }
   }

   if(m_appParam.LogLevel >= 3){
      (*m_osLog) << "Master only columns:" << endl;
      UtilPrintVector(model->masterOnlyCols, m_osLog);
      if(model->getColNames().size() > 0)
	 UtilPrintVector(model->masterOnlyCols, 
			 model->getColNames(), m_osLog);
   }
      
   //---
   //--- set the indices of the integer variables of model
   //---
   UtilIotaN(model->integerVars, col_xStartIndex, col_yStartIndex);
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "createModelCore()", m_appParam.LogLevel, 2);
   //---
   //--- display problem
   //---
   //int j = 0;
   //cout << "find: \nActive cols: ";
   //std::vector<int>::const_iterator it;
   //for(it = model->getActiveColumns().begin(); it != model->getActiveColumns().end(); it++){
   //cout << *it << " ";
   //}
   //
   //cout << "\ns.t.\n";
   //for(j = 0; j < model->getNumRows(); j++){
   //cout << model->rowNames[j] << " : \t\t";
   //cout << model->rowLB[j] << " \t <= \t" ;
   //for (i = 0; i < model->getNumCols(); i++){
   ////cout << "numCols=" << model->getNumCols() << endl;
   //if(model->getMatrix()->getCoefficient(j,i) != 0) {
   ////cout << "i" << i << " ";
   //cout << " " << model->M->getCoefficient(j,i) << " " ;
   //cout << model->getColNames()[i] ;
   //}
   //else {cout << "" ;}
   //}
   //cout << " \t <= \t " << model->rowUB[j] << endl ;
   //
   //}

}


//===========================================================================//
void SDPUC_DecompApp::createModelRelax(DecompConstraintSet * model,
				       int                   tpId){

   //---
   //--- SUBPROBLEM (A'): (one block for each k in K)
   //---      sum{(j,i) in A} x[k,i,j] - 
   //---         sum{(i,j) in A} x[k,i,j] = d[i,k], for all i in N
   //---      x[k,i,j] integer >= l[i,j] <= u[i,j], for all (i,j) in A
   //--- For k=(s,t) in K,
   //---    d[i,k] = -d[k] if i=s
   //---           =  d[k] if i=t
   //---           =  0, otherwise

   //--- SUBPROBLEM (A'): (one block for each t in T)
   //---     sum{(j,i) in A} x[i,j,t] - 
   //---        sum{(i,j) in A} x[i,j,t] = d[i,t],  for all i in N\{s}, where s is the super source
   //---      x[i,j,t] >= l[i,j,t] z[i,j,t],       for all (i,j) in A
   //---      x[i,j,t] <= u[i,j,t] z[i,j,t],       for all (i,j) in A
   //---	  r[i,j] x[i,j,t] - theta[j] + theta[i] <=  M (1 - z[i,j,t]) for all i,j in A
   //---	  r[i,j] x[i,j,t] - theta[j] + theta[i] >= -M (1 - z[i,j,t]) for all i,j in A
   //---
   int         a, i, j, head, tail, colIndex, source;
   int         numTimeperiods = m_instance.m_numTimeperiods;
   int         numArcs        = m_instance.m_numArcs;
   int		   numSwitchings  = m_instance.m_numSwitchings;
   int         numNodes       = m_instance.m_numNodes;   
   int         numCols        = numArcs						    //y1-vars
      + 3 * numTimeperiods * numArcs	//y2-, z-, and x-vars
      + numTimeperiods * numNodes;		//theta-vars

   SDPUC_Instance::arc        * arcs     = m_instance.m_arcs;
   SDPUC_Instance::node       * nodes    = m_instance.m_nodes;
   SDPUC_Instance::timeseries * ts		  = m_instance.m_timeseries;
  
   int		   numACArcs = 0;
   for(a = 0; a < numArcs; a++){
      if(arcs[a].acline == 1) { numACArcs++; }
   }

   int         numRows        = numNodes - 1    // balance
      + 2 * numArcs     // capacity 
      + 2 * numACArcs   // and kirchoffs constraints
      + 1;				// max. allowed no. of switches employed  sum{z} <= k

   int	 col_yStartIndex  =		 0;
   int	 col_zStartIndex  =		 numArcs*(1+numTimeperiods);
   int	 col_xStartIndex  =		 numArcs * (1 + 2*numTimeperiods);
   int	 col_thetaStartIndex  =	 numArcs*(1 + 3*numTimeperiods) ;

   double   bigM = 100;

  


   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createModelRelax()", m_appParam.LogLevel, 2);

  


   //---
   //--- create space for the model matrix (row-majored)
   //---
   model->M = new CoinPackedMatrix(false, 0.0, 0.0);
   if(!model->M)
      throw UtilExceptionMemory("createModelCore", "MCF_DecompApp");
   model->M->setDimensions(0, numCols);
   model->reserve(numRows, numCols);

   //---
   //--- set super source node
   //---
   source = 0;
   
   //---
   //--- create the rows 
   //---   NOTE: this is somewhat inefficient (but simple)
   //---
   int hej = 0;
   cout << "Generating sub-problem " << tpId << endl; hej++;
   //cout << "y_index = " << col_yStartIndex << endl;
   //cout << "z_index = " << col_zStartIndex << endl;
   //cout << "x_index = " << col_xStartIndex << endl;
   //cout << "theta_index = " << col_thetaStartIndex << endl;
   //--- create balance constraints
   for(i = 0; i < numNodes; i++){
      CoinPackedVector row;
      for(a = 0; a < numArcs; a++){
         tail = arcs[a].tail;
         head = arcs[a].head;
         if(head == i){
	    colIndex = col_xStartIndex + tpId * numArcs + a;
            row.insert(colIndex, 1.0);
         }
         else if(tail == i){
            colIndex = col_xStartIndex + tpId * numArcs + a;
            row.insert(colIndex, -1.0);
         }
      }

      //set demand d
      double d = nodes[i].demand * ts[nodes[i].tsdemand].values[tpId];
      std::string rowName = "balance_" + UtilIntToStr(i) + "_" + UtilIntToStr(tpId);
      if(i == source) {
         model->appendRow(row, -DecompInf, 0.0, rowName);
      }
      else {
         model->appendRow(row, d, d, rowName);
      }
   }

   //--- create capacity constraints and kirchoffs constraints
   for(a = 0; a < numArcs; a++){
      CoinPackedVector rowLB, rowUB; //lower-, and upperbound
      //for(a = 0; a < numArcs; a++){
      tail = arcs[a].tail;
      head = arcs[a].head;
      double lb = arcs[a].lb * ts[arcs[a].tscap].values[tpId];
      double ub = arcs[a].ub * ts[arcs[a].tscap].values[tpId];
      double reactance = arcs[a].weight;
      colIndex = col_zStartIndex + tpId * numArcs + a;		 
      rowLB.insert(colIndex, -lb);
      rowUB.insert(colIndex, -ub);
 		
      colIndex = col_xStartIndex + tpId * numArcs + a;
      rowLB.insert(colIndex, 1.0);
      rowUB.insert(colIndex, 1.0);
		 
		 
      //set flow lower and upperbound
      std::string rowNameLB = "lb_" + UtilIntToStr(a) + "_" + UtilIntToStr(tpId);
      std::string rowNameUB = "ub_" + UtilIntToStr(a) + "_" + UtilIntToStr(tpId);
      model->appendRow(rowLB, 0.0, DecompInf, rowNameLB);
      model->appendRow(rowUB, -DecompInf, 0.0, rowNameUB);

      //set kirchoffs voltage constraints for ac-arcs
      if(arcs[a].acline == 1) {
	 CoinPackedVector rowK1, rowK2; //Kirchoff1, and Kirchoff2
	 colIndex = col_zStartIndex + tpId * numArcs + a;
	 rowK1.insert(colIndex, bigM);
	 rowK2.insert(colIndex, -bigM);
	 colIndex = col_xStartIndex + tpId * numArcs + a;
	 rowK1.insert(colIndex, reactance);
	 rowK2.insert(colIndex, reactance);
	 for(i = 0; i < numNodes; i++){
	    if(head == i){
	       colIndex = col_thetaStartIndex + tpId * numNodes + i;
	       rowK1.insert(colIndex, -1.0);
	       rowK2.insert(colIndex, -1.0);
	    }
	    else if(tail == i){
	       colIndex = col_thetaStartIndex + tpId * numNodes + i;
	       rowK1.insert(colIndex, 1.0);
	       rowK2.insert(colIndex, 1.0);
	    }
	 }
         std::string rowNameK1 = "k1_" + UtilIntToStr(a) + "_" + UtilIntToStr(tpId);
         std::string rowNameK2 = "k2_" + UtilIntToStr(a) + "_" + UtilIntToStr(tpId);
	 model->appendRow(rowK1, -DecompInf, bigM, rowNameK1);
	 model->appendRow(rowK2, -bigM, DecompInf, rowNameK2);
      }

   }

   //--- create max. no. of switchings-constraint
   CoinPackedVector row;
   for(a = 0; a < numArcs; a++){
      colIndex = col_zStartIndex + tpId * numArcs + a;
      if(arcs[a].acline == 1) {  //only for arcs with voltage constraints
	 row.insert(colIndex, 1.0);
      }
   }
   //model->appendRow(row, numACArcs-numSwitchings, numACArcs, std::string("sum{1-z}<=k"));
   std::string rowNameSumK = "sumk_" + UtilIntToStr(tpId);
   model->appendRow(row, numACArcs-numSwitchings, numACArcs, rowNameSumK);


   //---
   //--- create a list of the "active" columns (those related 
   //---   to this commmodity) all other columns are fixed to 0
   //---
   UtilFillN(model->colLB, numCols,  0.0);
   UtilFillN(model->colUB, numCols,  0.0);
   colIndex = 0;
   for(a = 0; a < numArcs; a++){
      //set y-columns active
      //if(tpId == 0) {
      //colIndex = col_yStartIndex + a;
      //   model->colLB[colIndex] = 0;
      //   model->colUB[colIndex] = 1;
      //   model->activeColumns.push_back(colIndex);
      //cout << "" << colIndex << ", " ;
      //}
      //set z-columns active
      colIndex = col_zStartIndex + tpId * numArcs + a;
      model->colLB[colIndex] = 0; //1 - arcs[a].switchable;  //if arc not switchable then fix z=1
      model->colUB[colIndex] = 1;
      model->activeColumns.push_back(colIndex);
      //set x-columns active
      colIndex = col_xStartIndex + tpId * numArcs + a;
      double           arcLB = min(0.0, arcs[a].lb * ts[arcs[a].tscap].values[tpId]); //!!****
      double           arcUB = arcs[a].ub * ts[arcs[a].tscap].values[tpId];
      model->colLB[colIndex] = arcLB;
      model->colUB[colIndex] = arcUB;
      model->activeColumns.push_back(colIndex);
   }
   //set theta-columns active
   for(i = 0; i < numNodes; i++){
      colIndex = col_thetaStartIndex + tpId * numNodes + i;
      model->colLB[colIndex] = -DecompInf;
      model->colUB[colIndex] =  DecompInf;
      model->activeColumns.push_back(colIndex);
   }  

   //---
   //--- set the indices of the integer variables of model
   //---
   UtilIotaN(model->integerVars, col_xStartIndex, col_yStartIndex);
  
   ////Display problem
   //cout << "find: \nActive cols: ";
   //std::vector<int>::const_iterator it;
   //for(it = model->getActiveColumns().begin(); it < model->getActiveColumns().end(); it++){
   //cout << *it << " ";
   //}
   //cout << "\ns.t.\n";
   //for(j = 0; j < model->getNumRows(); j++){
   //cout << model->rowNames[j] << " : \t\t";
   //cout << model->rowLB[j] << " \t <= \t" ;
   //for (i = 0; i < model->getNumCols(); i++){
   ////cout << "numCols=" << model->getNumCols() << endl;
   //if(model->getMatrix()->getCoefficient(j,i) != 0) {
   ////cout << "i" << i << " ";
   //cout << " " << model->M->getCoefficient(j,i) << " (" << i << ") + ";
   ////cout << model->getColNames()[i] ;
   //}
   //else {cout << "" ;}
   //}
   //cout << " \t <= \t " << model->rowUB[j] << endl ;
   //
   //}

   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "createModelRelax()", m_appParam.LogLevel, 2);
}

//===========================================================================//
void SDPUC_DecompApp::createModelRelaxSparse(DecompConstraintSet * model,
					     int                   commId){

   ////---
   ////--- SUBPROBLEM (A'): (one block for each k in K)
   ////---      sum{(j,i) in A} x[k,i,j] - 
   ////---         sum{(i,j) in A} x[k,i,j] = d[i,k], for all i in N
   ////---      x[k,i,j] integer >= l[i,j] <= u[i,j], for all (i,j) in A
   ////--- For k=(s,t) in K,
   ////---    d[i,k] = -d[k] if i=s
   ////---           =  d[k] if i=t
   ////---           =  0, otherwise
   ////---
   //int         a, i, head, tail, origColIndex, source, sink;
   //int         numArcs        = m_instance.m_numArcs;
   //int         numNodes       = m_instance.m_numNodes;   
   //int         numCommodities = m_instance.m_numCommodities;
   //int         numCols        = numArcs;
   //int         numRows        = numNodes;
   //int         numColsOrig    = numArcs * numCommodities;
   //SILCEP_Instance::arc       * arcs        = m_instance.m_arcs;
   //SILCEP_Instance::commodity * commodities = m_instance.m_commodities;

   //UtilPrintFuncBegin(m_osLog, m_classTag,
   //     "createModelRelaxSparse()", m_appParam.LogLevel, 2);

   ////---
   ////--- create space for the model matrix (row-majored)
   ////---
   //model->M = new CoinPackedMatrix(false, 0.0, 0.0);
   //if(!model->M)
   //   throw UtilExceptionMemory("createModelCore", "MCF_DecompApp");
   //model->M->setDimensions(0, numCols);
   //model->reserve(numRows, numCols);
   //model->setSparse(numColsOrig);

   ////---
   ////--- get this commodity's source and sink node
   ////---
   //source = commodities[commId].source;
   //sink   = commodities[commId].sink;
   //
   ////---
   ////--- create the rows 
   ////---   NOTE: this is somewhat inefficient (but simple)
   ////---
   //for(i = 0; i < numNodes; i++){
   //   CoinPackedVector row;
   //   for(a = 0; a < numArcs; a++){
   //      tail = arcs[a].tail;
   //      head = arcs[a].head;
   //      if(head == i)
   //         row.insert(a, 1.0);
   //      else if(tail == i)
   //         row.insert(a, -1.0);
   //   }
   //   if(i == source)
   //      model->appendRow(row, 
   //                       -commodities[commId].demand, 
   //                       -commodities[commId].demand);
   //   else if(i == sink)
   //      model->appendRow(row, 
   //                       commodities[commId].demand, 
   //                       commodities[commId].demand);
   //   else 
   //      model->appendRow(row, 0.0, 0.0);
   //}

   ////---
   ////--- set the colLB, colUB, integerVars and sparse mapping
   ////---
   //origColIndex = commId * numArcs;
   //for(a = 0; a < numArcs; a++){
   //   double           arcLB = arcs[a].lb;
   //   double           arcUB = arcs[a].ub;
   //   model->pushCol(arcLB, arcUB, true, origColIndex);
   //   origColIndex++;      
   //}
   //   
   //UtilPrintFuncEnd(m_osLog, m_classTag,
   //                 "createModelRelaxSparse()", m_appParam.LogLevel, 2);
}

//===========================================================================//
void
SDPUC_DecompApp::createModelMasterOnlys(vector<int> & masterOnlyCols){

   int   numTimeperiods = m_instance.m_numTimeperiods;
   int   numArcs        = m_instance.m_numArcs;
   int	 numNodes		= m_instance.m_numNodes;
   int   numCols        = numArcs						//y1-vars
      + 3 * numTimeperiods * numArcs	//y2-, z-, and x-vars
      + numTimeperiods * numNodes;		//theta-vars
   SDPUC_Instance::arc * arcs = m_instance.m_arcs;
   
   int            nBlocks     = m_instance.m_numTimeperiods; //static_cast<int>(m_blocks.size());
   const int      nCols       = numCols;
   //const double * colLB[nCols]	;//	  = ;
   //const double * colUB[nCols]; //       = m_mpsIO.getColUpper();
   //const char   * integerVars = m_mpsIO.integerColumns();
   int            nMasterOnlyCols =     numArcs * (1 + numTimeperiods); //static_cast<int>(m_instance.m_numArcs);
   

   if(m_appParam.LogLevel >= 1){
      (*m_osLog) << "nCols           = " << nCols << endl;
      (*m_osLog) << "nMasterOnlyCols = " << nMasterOnlyCols << endl;
   }

   if(nMasterOnlyCols == 0)
      return;


   int i;
   vector<int>::iterator vit;
   for(vit = masterOnlyCols.begin(); vit != masterOnlyCols.end(); vit++){
      i = *vit;

      //set upper and lower bounds 
      //	  colUB[i] = 1;
      //	  colLB[i] = 0;

      //THINK:
      //  what-if master-only var is integer and bound is not at integer
            
      DecompConstraintSet * model = new DecompConstraintSet();
      model->m_masterOnly      = true;
      model->m_masterOnlyIndex = i;
      model->m_masterOnlyLB    = 0; //colLB[i];
      model->m_masterOnlyUB    = 1; //colUB[i];
      //0=cont, 1=integer
      model->m_masterOnlyIsInt = true; //integerVars[i] ? true : false;

      /*   if(m_appParam.ColumnUB <  1.0e15)
	   if(colUB[i] >  1.0e15)
	   model->m_masterOnlyUB = m_appParam.ColumnUB;
	   if(m_appParam.ColumnLB > -1.0e15)
	   if(colLB[i] < -1.0e15)
	   model->m_masterOnlyLB = m_appParam.ColumnLB;
      */
      //m_modelR.insert(make_pair(nBlocks, model));
      setModelRelax(model, 
                    "master_only" + UtilIntToStr(i), nBlocks);
      nBlocks++;
   }

   return;   
}


//===========================================================================//
DecompSolution
SDPUC_DecompApp::createInitialSolution(){
   
  

   int         a, i, t,  colIndex;
   int         numTimeperiods = m_instance.m_numTimeperiods;
   int         numArcs        = m_instance.m_numArcs;
   int		   numSwitchings  = m_instance.m_numSwitchings;
   int         numNodes       = m_instance.m_numNodes;   
   int         numCols        = numArcs						   //y1-vars
      + 3 * numTimeperiods * numArcs	//y2-, z-, and x-vars
      + numTimeperiods * numNodes;		//theta-vars

   SDPUC_Instance::arc        * arcs     = m_instance.m_arcs;
   SDPUC_Instance::node       * nodes    = m_instance.m_nodes;
   SDPUC_Instance::timeseries * ts		  = m_instance.m_timeseries;
  
   int		   numACArcs = 0;
   for(a = 0; a < numArcs; a++){
      if(arcs[a].acline == 1) { numACArcs++; }
   }

   int         numRows        = numNodes - 1    //balance
      + 2 * numArcs     // capacity 
      + 2 * numACArcs   // and kirchoffs constraints
      + 1;				// max. allowed no. of switches employed  sum{z} <= k

   int	 col_yStartIndex  =		 0;
   int	 col_zStartIndex  =		 numArcs*(1+numTimeperiods);
   int	 col_xStartIndex  =		 numArcs * (1 + 2*numTimeperiods);
   int	 col_thetaStartIndex  =	 numArcs*(1 + 3*numTimeperiods) ;
   double   bigM = 100;

   const int   size = numCols;
   double values_[100000];
   const double quality = 1e75;
   //const double * cost[size];

   //double my_non_const_array[100000];
   //int dummy = read_array_from_file(my_non_const_array);
   //double const (&array)[100000] = my_non_const_array;
   
   
   //initialise all values to 0
   for(i = 0; i < numCols; i++){
      values_[i] = 0;
   }
   // set z = 1
   for(a = 0; a < numArcs; a++){
      for(t = 0; t < numTimeperiods; t++){
	 colIndex = col_zStartIndex + t * numArcs + a;
	 values_[colIndex] = 1.0;
      }  
   }

   double const (&values)[100000] = values_;  

   DecompSolution sol(size, values, quality);

   return sol;
}
