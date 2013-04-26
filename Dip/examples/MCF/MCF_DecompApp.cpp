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
#include "DecompVar.h"
#include "MCF_DecompApp.h"

//===========================================================================//
void MCF_DecompApp::initializeApp(UtilParameters & utilParam) {
      
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
   int rc = m_instance.readInstance(instanceFile);
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
void MCF_DecompApp::createModels(){

   //---
   //--- This function does the work to create the different models 
   //---  that will be used. This memory is owned by the user. It will
   //---  be passed to the application interface and used by the algorithms.
   //---
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createModels()", m_appParam.LogLevel, 2);

   //---
   //--- (Integer) Multi-Commodity Flow Problem (MCF).
   //---
   //--- We are given: 
   //---    (1) a directed graph G=(N,A),
   //---    (2) a set of commodities K, where each commodity is 
   //---         a source-sink pair.
   //---
   //--- min  sum{k in K} sum{(i,j) in A} w[i,j] x[k,i,j]
   //--- s.t. sum{(j,i) in A} x[k,i,j] - 
   //---        sum{(i,j) in A} x[k,i,j] = d[i,k],  for all i in N, k in K
   //---      sum{k in K} x[k,i,j] >= l[i,j],       for all (i,j) in A 
   //---      sum{k in K} x[k,i,j] <= u[i,j],       for all (i,j) in A 
   //---      x[k,i,j] integer >= l[i,j] <= u[i,j], for all (i,j) in A
   //--- For k=(s,t) in K,
   //---    d[i,k] = -d[k] if i=s
   //---           =  d[k] if i=t
   //---           =  0, otherwise
   //---
   //--- NOTE: to make sure the problem is always feasible, dummy arcs
   //---   have been added between all source-sink commodity pairs and have
   //---   been given a 'big' weight.
   //---
   //---
   //--- The decomposition is formed as:
   //---
   //--- MASTER (A''):
   //---      sum{k in K} x[k,i,j] >= l[i,j],       for all (i,j) in A 
   //---      sum{k in K} x[k,i,j] <= u[i,j],       for all (i,j) in A 
   //---      x[k,i,j] integer >= l[i,j] <= u[i,j], for all (i,j) in A
   //---
   //--- SUBPROBLEM (A'): (one block for each k in K)
   //---      sum{(j,i) in A} x[k,i,j] - 
   //---         sum{(i,j) in A} x[k,i,j] = d[i,k], for all i in N
   //---      x[k,i,j] integer >= l[i,j] <= u[i,j], for all (i,j) in A
   //---

   //---
   //--- Get information about this problem instance.
   //---
   int   k, a, colIndex;
   int   numCommodities = m_instance.m_numCommodities;
   int   numArcs        = m_instance.m_numArcs;
   int   numCols        = numCommodities * numArcs;
   MCF_Instance::arc * arcs = m_instance.m_arcs;

   //---
   //--- Construct the objective function and set it
   //---    columns indexed as [k,a]= k*numArcs + a
   //---
   m_objective = new double[numCols];
   if(!m_objective)
      throw UtilExceptionMemory("createModels", "MCF_DecompApp");
   colIndex = 0;
   for(k = 0; k < numCommodities; k++)
      for(a = 0; a < numArcs; a++)
         m_objective[colIndex++] = arcs[a].weight;

   //---
   //--- set the objective 
   //---        
   setModelObjective(m_objective);

   //---
   //--- create the core/master model and set it
   //---
   DecompConstraintSet * modelCore = new DecompConstraintSet();      
   createModelCore(modelCore);
   setModelCore(modelCore, "core");

   //---
   //--- create the relaxed/subproblem models and set them
   //---
   for(k = 0; k < numCommodities; k++){
      DecompConstraintSet * modelRelax = new DecompConstraintSet();
      string                modelName  = "relax" + UtilIntToStr(k);
      if(m_appParam.UseSparse)
         createModelRelaxSparse(modelRelax, k);
      else
         createModelRelax(modelRelax, k);
      
      setModelRelax(modelRelax, modelName, k);
   }
   
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "createModels()", m_appParam.LogLevel, 2);
}

//===========================================================================//
void MCF_DecompApp::createModelCore(DecompConstraintSet * model){

   //---
   //--- MASTER (A''):
   //---      sum{k in K} x[k,i,j] >= l[i,j],       for all (i,j) in A 
   //---      sum{k in K} x[k,i,j] <= u[i,j],       for all (i,j) in A 
   //---      x[k,i,j] integer >= l[i,j] <= u[i,j], for all (i,j) in A
   //---                                    
   int   k, a, colIndex;
   int   numCommodities = m_instance.m_numCommodities;
   int   numArcs        = m_instance.m_numArcs;
   int   numCols        = numCommodities * numArcs;
   int   numRows        = 2              * numArcs;
   MCF_Instance::arc * arcs = m_instance.m_arcs;

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createModelCore()", m_appParam.LogLevel, 2);

   //---
   //--- create space for the model matrix (row-majored)
   //---
   model->M = new CoinPackedMatrix(false, 0.0, 0.0);
   if(!model->M)
      throw UtilExceptionMemory("createModelCore", "MCF_DecompApp");
   model->M->setDimensions(0, numCols);
   model->reserve(numRows, numCols);
   
   //---
   //--- create the rows and set the col/row bounds
   //---
   UtilFillN(model->colLB, numCols,  0.0);
   UtilFillN(model->colUB, numCols,  DecompInf);
   for(a = 0; a < numArcs; a++){
      CoinPackedVector row;
      double           arcLB = arcs[a].lb;
      double           arcUB = arcs[a].ub;
      for(k = 0; k < numCommodities; k++){
         colIndex = k * numArcs + a;
         model->colLB[colIndex] = arcLB;
         model->colUB[colIndex] = arcUB;
         row.insert(colIndex, 1.0);
      }
      //TODO: any issue with range constraints?
      model->appendRow(row, -DecompInf, arcUB);
      string rowNameUB = "capUB(" +
	 UtilIntToStr(a)            + "_" + 
	 UtilIntToStr(arcs[a].tail) + "," +
	 UtilIntToStr(arcs[a].head) + ")";
	 model->rowNames.push_back(rowNameUB);

      model->appendRow(row, arcLB,  DecompInf);
      string rowNameLB = "capLB(" +
	 UtilIntToStr(a)            + "_" + 
	 UtilIntToStr(arcs[a].tail) + "," +
	 UtilIntToStr(arcs[a].head) + ")";
	 model->rowNames.push_back(rowNameLB);
   }

   //---
   //--- create column names (helps with debugging)
   //---
   for(k = 0; k < numCommodities; k++){
      for(a = 0; a < numArcs; a++){
         string colName = "x(comm_" + UtilIntToStr(k) + "," +
            UtilIntToStr(a)            + "_" +
            UtilIntToStr(arcs[a].tail) + "," +
            UtilIntToStr(arcs[a].head) + ")";
         model->colNames.push_back(colName);
      }
   }
      
   //---
   //--- set the indices of the integer variables of model
   //---
   UtilIotaN(model->integerVars, numCols, 0);
   
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "createModelCore()", m_appParam.LogLevel, 2);
}


//===========================================================================//
void MCF_DecompApp::createModelRelax(DecompConstraintSet * model,
                                     int                   commId){

   //---
   //--- SUBPROBLEM (A'): (one block for each k in K)
   //---      sum{(j,i) in A} x[k,i,j] - 
   //---         sum{(i,j) in A} x[k,i,j] = d[i,k], for all i in N
   //---      x[k,i,j] integer >= l[i,j] <= u[i,j], for all (i,j) in A
   //--- For k=(s,t) in K,
   //---    d[i,k] = -d[k] if i=s
   //---           =  d[k] if i=t
   //---           =  0, otherwise
   //---
   int         a, i, head, tail, colIndex, source, sink;
   int         numCommodities = m_instance.m_numCommodities;
   int         numArcs        = m_instance.m_numArcs;
   int         numNodes       = m_instance.m_numNodes;   
   int         numCols        = numCommodities * numArcs;
   int         numRows        = numNodes;
   MCF_Instance::arc       * arcs        = m_instance.m_arcs;
   MCF_Instance::commodity * commodities = m_instance.m_commodities;

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
   //--- get this commodity's source and sink node
   //---
   source = commodities[commId].source;
   sink   = commodities[commId].sink;
   
   //---
   //--- create the rows 
   //---   NOTE: this is somewhat inefficient (but simple)
   //---
   for(i = 0; i < numNodes; i++){
      CoinPackedVector row;
      for(a = 0; a < numArcs; a++){
         tail = arcs[a].tail;
         head = arcs[a].head;
         if(head == i){
            colIndex = commId * numArcs + a;
            row.insert(colIndex, 1.0);
         }
         else if(tail == i){
            colIndex = commId * numArcs + a;
            row.insert(colIndex, -1.0);
         }
      }
      if(i == source)
         model->appendRow(row, 
                          -commodities[commId].demand, 
                          -commodities[commId].demand);
      else if(i == sink)
         model->appendRow(row, 
                          commodities[commId].demand, 
                          commodities[commId].demand);
      else 
         model->appendRow(row, 0.0, 0.0);

      string rowName = "flow(" +	 
	 UtilIntToStr(commId) + "_" +
	 UtilIntToStr(i)      + "_" +
	 UtilIntToStr(source) + "," + 
	 UtilIntToStr(sink)   + ")";
      model->rowNames.push_back(rowName);      
   }

   //---
   //--- create a list of the "active" columns (those related 
   //---   to this commmodity) all other columns are fixed to 0
   //---
   UtilFillN(model->colLB, numCols,  0.0);
   UtilFillN(model->colUB, numCols,  0.0);
   colIndex = commId * numArcs;
   for(a = 0; a < numArcs; a++){
      double           arcLB = arcs[a].lb;
      double           arcUB = arcs[a].ub;
      model->colLB[colIndex] = arcLB;
      model->colUB[colIndex] = arcUB;
      model->activeColumns.push_back(colIndex);
      colIndex++;      
   }
      
   //---
   //--- set the indices of the integer variables of model
   //---
   UtilIotaN(model->integerVars, numCols, 0);
   
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "createModelRelax()", m_appParam.LogLevel, 2);
}

//===========================================================================//
void MCF_DecompApp::createModelRelaxSparse(DecompConstraintSet * model,
                                           int                   commId){

   //---
   //--- SUBPROBLEM (A'): (one block for each k in K)
   //---      sum{(j,i) in A} x[k,i,j] - 
   //---         sum{(i,j) in A} x[k,i,j] = d[i,k], for all i in N
   //---      x[k,i,j] integer >= l[i,j] <= u[i,j], for all (i,j) in A
   //--- For k=(s,t) in K,
   //---    d[i,k] = -d[k] if i=s
   //---           =  d[k] if i=t
   //---           =  0, otherwise
   //---
   int         a, i, head, tail, origColIndex, source, sink;
   int         numArcs        = m_instance.m_numArcs;
   int         numNodes       = m_instance.m_numNodes;   
   int         numCommodities = m_instance.m_numCommodities;
   int         numCols        = numArcs;
   int         numRows        = numNodes;
   int         numColsOrig    = numArcs * numCommodities;
   MCF_Instance::arc       * arcs        = m_instance.m_arcs;
   MCF_Instance::commodity * commodities = m_instance.m_commodities;

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createModelRelaxSparse()", m_appParam.LogLevel, 2);

   //---
   //--- create space for the model matrix (row-majored)
   //---
   model->M = new CoinPackedMatrix(false, 0.0, 0.0);
   if(!model->M)
      throw UtilExceptionMemory("createModelCore", "MCF_DecompApp");
   model->M->setDimensions(0, numCols);
   model->reserve(numRows, numCols);
   model->setSparse(numColsOrig);

   //---
   //--- get this commodity's source and sink node
   //---
   source = commodities[commId].source;
   sink   = commodities[commId].sink;
   
   //---
   //--- create the rows 
   //---   NOTE: this is somewhat inefficient (but simple)
   //---
   for(i = 0; i < numNodes; i++){
      CoinPackedVector row;
      for(a = 0; a < numArcs; a++){
         tail = arcs[a].tail;
         head = arcs[a].head;
         if(head == i)
            row.insert(a, 1.0);
         else if(tail == i)
            row.insert(a, -1.0);
      }
      if(i == source)
         model->appendRow(row, 
                          -commodities[commId].demand, 
                          -commodities[commId].demand);
      else if(i == sink)
         model->appendRow(row, 
                          commodities[commId].demand, 
                          commodities[commId].demand);
      else 
         model->appendRow(row, 0.0, 0.0);
   }

   //---
   //--- set the colLB, colUB, integerVars and sparse mapping
   //---
   origColIndex = commId * numArcs;
   for(a = 0; a < numArcs; a++){
      double           arcLB = arcs[a].lb;
      double           arcUB = arcs[a].ub;
      model->pushCol(arcLB, arcUB, true, origColIndex);
      origColIndex++;      
   }
      
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "createModelRelaxSparse()", m_appParam.LogLevel, 2);
}
