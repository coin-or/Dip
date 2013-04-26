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
#include "GAP_DecompApp.h"
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
   //--- read best known lb/ub
   //---
   string bestKnownFile = m_appParam.DataDir + UtilDirSlash() + "gap.opt";
   m_instance.readBestKnown(bestKnownFile, m_appParam.Instance);
   setBestKnownLB(m_instance.getBestKnownLB());
   setBestKnownUB(m_instance.getBestKnownUB());

   //---
   //--- open space for GAP_Knapsack objects
   //---
   int          k;
   const int    nTasks    = m_instance.getNTasks();
   const int    nMachines = m_instance.getNMachines();
   const int  * capacity  = m_instance.getCapacity();
   const int  * weight    = m_instance.getWeight();
   const int  * profit    = m_instance.getProfit();
   GAP_KnapPisinger * knapK  = 0;
   
   m_knap.reserve(nMachines);
   for(k = 0; k < nMachines; k++){
      knapK = new GAP_KnapPisinger(nTasks, 
				   capacity[k], 
				   weight + (k * nTasks),
				   profit + (k * nTasks));      
      m_knap.push_back(knapK);
   }

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
   
   model->M = new CoinPackedMatrix(false, 0.0, 0.0);
   CoinAssertHint(model->M, "Error: Out of Memory");
   model->M->setDimensions(0, nCols);
   model->reserve(nRows, nCols);

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
   //--- set the col upper and lower bounds
   //---
   UtilFillN(model->colLB, nCols,  0.0);
   UtilFillN(model->colUB, nCols,  1.0);

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
   
   //---
   //--- set the indices of the integer variables of model
   //---
   UtilIotaN(model->integerVars, nCols, 0);
   
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "createModelPartAP()", m_appParam.LogLevel, 2);

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
   m_models.insert(make_pair("AP", modelCore));
   
   
   for(i = 0; i < nMachines; i++){
      modelName = "KP" + UtilIntToStr(i);
      setModelRelax(NULL, modelName, i);
   }

   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "createModels()", m_appParam.LogLevel, 2);
   return status;
}

//--------------------------------------------------------------------- //
DecompSolverStatus 
GAP_DecompApp::solveRelaxed(const int             whichBlock,
                            const double        * redCostX,
                            const double          convexDual,
                            list<DecompVar*>    & vars){

   if(!m_appParam.UsePisinger)
      return DecompSolStatNoSolution;
   
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "solveRelaxed()", m_appParam.LogLevel, 2);
   
   vector<int>      solInd;
   vector<double>   solEls;
   double           varRedCost  = 0.0;
   double           varOrigCost = 0.0;         
   double         * origCost    = m_objective;
   const double   * redCostXB   = redCostX + getOffsetI(whichBlock);
   const double   * origCostB   = origCost + getOffsetI(whichBlock);
   
  
      
   //---
   //--- print out red cost 
   //---
   /*{
     int j;	    
     const int    nTasks    = m_instance.getNTasks();
     const int *  weight    = m_instance.getWeight() + getOffsetI(b);
     for(j = 0; j < nTasks; j++){	    
     printf("RedCost[j=%d, wt=%d]: %g\n", j, weight[j], redCostXB[j]);
     }      
     }*/

   m_knap[whichBlock]->solve(whichBlock, 
                             redCostXB,
                             origCostB,
                             solInd,
                             solEls,
                             varRedCost,
                             varOrigCost);
   
   //printf("b=%d alpha              = %g\n", b, alpha);
   //printf("b=%d varRedCost - alpha = %g\n", b, varRedCost - alpha);
   //printf("b=%d varOrigCost        = %g\n", b, varOrigCost);

   
   UTIL_DEBUG(m_appParam.LogLevel, 4,
	      printf("PUSH var with RC = %g\n", varRedCost - convexDual);
	      );
   
   DecompVar * var = new DecompVar(solInd, solEls, 
				   varRedCost - convexDual, varOrigCost);
   var->setBlockId(whichBlock);
   vars.push_back(var);
   
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "APPsolveRelaxed()", m_appParam.LogLevel, 2);   
   return DecompSolStatOptimal;
}


//--------------------------------------------------------------------- //
void GAP_DecompApp::printOriginalColumn(const int   index, 
					ostream   * os) const {
   pair<int,int> p = m_instance.getIndexInv(index);
   (*os) << "x[ " << index << " : " << p.first << " , " << p.second << " ]";
}
