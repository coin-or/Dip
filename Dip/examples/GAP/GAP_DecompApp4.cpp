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
#include "GAP_Status.h"
#include "GAP_DecompApp3.h"
//===========================================================================//
#include "DecompVar.h"

//===========================================================================//
void GAP_DecompApp::initializeApp(UtilParameters & utilParam) {
   

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "initializeApp()", m_appParam.LogLevel, 2);
   
   //---
   //--- get application parameters
   //---
   m_appParam.getSettings(utilParam);
   if(m_appParam.LogLevel >= 1)
      m_appParam.dumpSettings(m_osLog);
   
   //---
   //--- read instance
   //
   string instanceFile = m_appParam.DataDir
      + UtilDirSlash() + m_appParam.Instance;
   m_instance.readInstance(instanceFile);

   //---
   //--- read best known lb/ub (for debugging)
   //---
   string bestKnownFile = m_appParam.DataDir + UtilDirSlash() + "gap.opt";
   m_instance.readBestKnown(bestKnownFile, m_appParam.Instance);
   setBestKnownLB(m_instance.getBestKnownLB());
   setBestKnownUB(m_instance.getBestKnownUB());

   //---
   //--- create models
   //---
   createModels();
   
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "initializeApp()", m_appParam.LogLevel, 2);
}

// --------------------------------------------------------------------- //
int GAP_DecompApp::createModelPartAP(DecompConstraintSet * model){
   
   int             i, j, colIndex;
   int             status     = GAPStatusOk;
   int             nTasks     = m_instance.getNTasks();    //n
   int             nMachines  = m_instance.getNMachines(); //m
   int             nCols      = nTasks * nMachines;
   int             nRows      = nTasks;
   
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createModelPartAP()", m_appParam.LogLevel, 2);

   //---
   //--- Build the core model constraints (AP = assignment problem).
   //---
   //--- m is number of machines (index i)
   //--- n is number of tasks    (index j)
   //---
   //--- sum{i in 1..m} x[i,j] = 1, j in 1..n
   //---
   //--- Example structure: m=3, n=4
   //---  x   x   x     = 1 [j=1]
   //---   x   x   x    = 1 [j=2]
   //---    x   x   x   = 1 [j=3]
   //---     x   x   x  = 1 [j=4]
   //---

   //---
   //--- Allocate an empty row-ordered CoinPackedMatrix. Since we plan
   //---   to add rows, set the column dimension and let the row dimension
   //---   be set dynamically.
   //---   
   model->M = new CoinPackedMatrix(false, 0.0, 0.0);
   CoinAssertHint(model->M, "Error: Out of Memory");
   model->M->setDimensions(0, nCols);

   //---
   //--- we know the sizes needed, so reserve space for them (for efficiency)
   //---
   model->reserve(nRows, nCols);

   //---
   //--- create one row per task
   //---   rowNames are not needed, they are used for debugging
   //---
   for(j = 0; j < nTasks; j++){
      CoinPackedVector row;
      string           rowName = "a(j_" + UtilIntToStr(j) + ")";
      for(i = 0; i < nMachines; i++){
         colIndex = getIndexIJ(i,j);
         row.insert(colIndex, 1.0);
      }
      model->appendRow(row, 1.0, 1.0, rowName);
   }
   
   //---   
   //--- set the col upper and lower bounds (all in [0,1])
   //---
   UtilFillN(model->colLB, nCols,  0.0);
   UtilFillN(model->colUB, nCols,  1.0);

   //---
   //--- set the indices of the integer variables of model 
   //---   (all vars are binary)
   //---
   UtilIotaN(model->integerVars, nCols, 0);
   
   //---
   //--- set column names for debugging
   //---
   colIndex = 0;
   for(i = 0; i < nMachines; i++){
      for(j = 0; j < nTasks; j++){
         string colName = "x("
            + UtilIntToStr(colIndex) + "_"
            + UtilIntToStr(i) + "," + UtilIntToStr(j) + ")";
         model->colNames.push_back(colName);
         colIndex++;
      }      
   }
   
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "createModelPartAP()", m_appParam.LogLevel, 2);

   return status;
}

//===========================================================================//
int GAP_DecompApp::createModelPartKP(DecompConstraintSet * model){
   //---
   //--- helper method - create model with nMachines KPs
   //---
   vector<int> whichKnaps;
   int         nMachines = m_instance.getNMachines();
   UtilIotaN(whichKnaps, nMachines, 0);
   return createModelPartKP(model, whichKnaps);
}

//===========================================================================//
int GAP_DecompApp::createModelPartKP(DecompConstraintSet * model, 
                                     int                   whichKnap){
   //---
   //--- helper method - create model with one (whichKnap) KP
   //---
   vector<int> whichKnaps;
   whichKnaps.push_back(whichKnap);
   return createModelPartKP(model, whichKnaps);
}

//===========================================================================//
int GAP_DecompApp::createModelPartKP(DecompConstraintSet * model, 
                                     vector<int>         & whichKnaps){
   
   int          i, j, b, colIndex;
   int          status     = GAPStatusOk;
   int          nTasks     = m_instance.getNTasks();    //n
   int          nMachines  = m_instance.getNMachines(); //m
   int          nKnaps     = static_cast<int>(whichKnaps.size());
   const int *  weight     = m_instance.getWeight();   
   const int *  capacity   = m_instance.getCapacity();
   int          nCols      = nTasks * nKnaps;
   int          nRows      = nKnaps;

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createModelPartKP()", m_appParam.LogLevel, 2);

   //---
   //--- Build the relax model constraints (KP = assignment problem).
   //---
   //--- m is number of machines (index i)
   //--- n is number of tasks    (index j)
   //---
   //--- sum{j in 1..n} w[i,j] x[i,j] <= b[i], i in 1..m
   //--- x[i,j] in {0,1}, i in 1..m, j in 1..n
   //---
   //--- Example structure: m=3, n=4
   //---  xxxx         <= b[i=1]
   //---      xxxx     <= b[i=2]
   //---          xxxx <= b[i=3]
   //---

   //---
   //--- Allocate an empty row-ordered CoinPackedMatrix. Since we plan
   //---   to add rows, set the column dimension and let the row dimension
   //---   be set dynamically.
   //--- NOTE: this matrix is sparse version. So, nCols is just the number
   //---   of active columns. That is, if we are generating one KP per block,
   //---   then there are only nTasks columns in this model.
   //---  
   model->M = new CoinPackedMatrix(false, 0.0, 0.0);
   CoinAssertHint(model->M, "Error: Out of Memory");
   model->M->setDimensions(0, nCols);

   //---
   //--- we know the sizes needed, so reserve space for them (for efficiency)
   //---
   model->reserve(nRows, nCols);

   //---
   //--- tell the solver that this block is represented sparsely
   //---   as an argument, tell the method how many columns are 
   //---   in the original compact formulation
   //---
   model->setSparse(nTasks * nMachines);

   //---
   //--- create the columns using the pushCol interface
   //---    and setup the mapping between sparse and dense models
   //---
   //---   pushCol(
   //---           column lower bound
   //---           column upper bound
   //---           is it integer?
   //---           index for the original compact formulation )
   //---
   vector<int>::iterator it;
   for(it = whichKnaps.begin(); it != whichKnaps.end(); it++){
      b = *it;
      for(i = 0; i < nMachines; i++)
         for(j = 0; j < nTasks; j++)
            if(i == b){
	       //---
	       //--- set column names for debugging
	       //---
	       string colName = "x("
		  + UtilIntToStr(getIndexIJ(i,j)) + "_"
		  + UtilIntToStr(i) + "," + UtilIntToStr(j) + ")";
	       model->colNames.push_back(colName);	       
	       model->pushCol(0.0, 1.0, true, getIndexIJ(i,j));
	    }
   }
   
   //---
   //--- create one row per knapsack
   //----  this is a sparse matrix, so use mapping between original
   //---   rowNames are not needed, they are used for debugging
   //---    
   const map<int,int> &  origToSparse = model->getMapOrigToSparse();
   map<int,int>::const_iterator mit;
   for(it = whichKnaps.begin(); it != whichKnaps.end(); it++){
      i = *it;
      CoinPackedVector row;
      string           rowName = "k(i_" + UtilIntToStr(i) + ")";
      for(j = 0; j < nTasks; j++){
	 colIndex = getIndexIJ(i,j);//dense
	 mit      = origToSparse.find(colIndex);
	 assert(mit != origToSparse.end());         
         row.insert(mit->second, weight[colIndex]);
      }      
      model->appendRow(row, -DecompInf, capacity[i], rowName);
   }
      
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "createModelPartKP()", m_appParam.LogLevel, 2);

   return status;
}

// --------------------------------------------------------------------- //
int GAP_DecompApp::createModels(){

   //---
   //--- This function does the work to create the different models 
   //---  that will be used. This memory is owned by the user. It will
   //---  be passed to the application interface and used by the algorithms.
   //---
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createModels()", m_appParam.LogLevel, 2);
   

   //---
   //--- Generalized Assignment Problem (GAP)
   //---   m is number of machines (index i)
   //---   n is number of tasks    (index j)
   //---
   //--- min  sum{i in 1..m, j in 1..n} p[i,j] x[i,j]
   //--- s.t. sum{           j in 1..n} w[i,j] x[i,j] <= b[i], i in 1..m
   //---      sum{i in 1..m           }        x[i,j]  = 1   , j in 1..n
   //---      x[i,j] in {0,1}, i in 1..m, j in 1..n
   //---
   //--- Example structure: m=3, n=4
   //---  xxxx         <= b[i=1]
   //---      xxxx     <= b[i=2]
   //---          xxxx <= b[i=3]
   //---  x   x   x     = 1 [j=1]
   //---   x   x   x    = 1 [j=2]
   //---    x   x   x   = 1 [j=3]
   //---     x   x   x  = 1 [j=4]
   //---

   //---
   //--- Get information about this problem instance.
   //--
   int          i;
   string       modelName;
   int          status     = GAPStatusOk;
   int          nTasks     = m_instance.getNTasks();    //n
   int          nMachines  = m_instance.getNMachines(); //m
   const int *  profit     = m_instance.getProfit();
   int          nCols      = nTasks * nMachines;

   //---
   //--- Construct the objective function (the original problem is 
   //---  a maximization, so we flip the sign to make it minimization).
   //---
   m_objective = new double[nCols];
   assert(m_objective);
   if(!m_objective)
      return GAPStatusOutOfMemory;
   for(i = 0; i < nCols; i++)
      m_objective[i] = profit[i];

   //---
   //--- A'[i] for i=1..m: m independent knapsacks
   //---  sum{j in 1..n}  w[i,j] x[i,j] <= b[i]
   //---  x[i,j] in {0,1}, i in 1..m, j in 1..n
   //---
   //--- A'':
   //---  sum{i in 1..m} x[i,j] = 1, j in 1..n
   //---
   //--- Example structure: m=3, n=4
   //--- A'[i=1]:
   //---  xxxx         <= b[i=1]
   //--- A'[i=2]:
   //---      xxxx     <= b[i=2]
   //--- A'[i=3]:
   //---          xxxx <= b[i=3]
   //---
   //--- A'':
   //---  x   x   x     = 1 [j=1]
   //---   x   x   x    = 1 [j=2]
   //---    x   x   x   = 1 [j=3]
   //---     x   x   x  = 1 [j=4]
   //---
   setModelObjective(m_objective);
   
   DecompConstraintSet * modelCore = new DecompConstraintSet();
   status = createModelPartAP(modelCore);
   if(status) return status;
   
   setModelCore(modelCore, "AP");
   m_models.push_back(modelCore);
   
   for(i = 0; i < nMachines; i++){
      DecompConstraintSet * modelRelax = new DecompConstraintSet();
      status = createModelPartKP(modelRelax, i);
      
      modelName = "KP" + UtilIntToStr(i);
      setModelRelax(modelRelax, modelName, i);
      m_models.push_back(modelRelax);
   }
   
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "createModels()", m_appParam.LogLevel, 2);
   return status;
}

//--------------------------------------------------------------------- //
void GAP_DecompApp::printOriginalColumn(const int   index, 
					ostream   * os) const {
   pair<int,int> p = m_instance.getIndexInv(index);
   (*os) << "x[ " << index << " : " << p.first << " , " << p.second << " ]";
}
