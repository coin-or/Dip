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
#include "MMKP_DecompApp.h"

//===========================================================================//
void MMKP_DecompApp::initializeApp(UtilParameters & utilParam) {
      
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
   string dataFormat     = m_appParam.DataFormat;
   if(dataFormat == "khan" || dataFormat == "hifi")
      m_instance.readInstance(instanceFile, dataFormat);      
   else if(dataFormat == "gsimon")
      m_instance.readInstanceSimon(instanceFile);
   else
      throw UtilException("Unknown data format", 
                          "initializeApp", "MMKP_DecompApp");
   
   //---
   //--- read best known lb/ub
   //---
   if(dataFormat == "khan"){
      string bestKnownFile = m_appParam.DataDir + UtilDirSlash() + "mmkp.opt";
      m_instance.readBestKnown(bestKnownFile, m_appParam.Instance);
      setBestKnownLB(m_instance.getBestKnownLB());
      setBestKnownUB(m_instance.getBestKnownUB());
   }

   //---
   //--- open space for MMKP_MCKnap objects
   //---
   int                    nKnapRows  = m_instance.getNKnapRows();
   int                    nGroupRows = m_instance.getNGroupRows();
   int                    nGroupCols = m_instance.getNGroupCols();   
   const double *         capacity   = m_instance.getCapacity();
   const double * const * weight     = m_instance.getWeight();
   MMKP_MCKnap  *         mcknapK    = NULL;
   
   int k;
   m_mcknap.reserve(nKnapRows);
   for(k = 0; k < nKnapRows; k++){
      mcknapK = new MMKP_MCKnap(nGroupRows, nGroupCols);      
      mcknapK->setMCKnapData(capacity[k], weight[k]);
      m_mcknap.push_back(mcknapK);
   }
   
   //---
   //--- open memory for auxiliary memory pool - need?
   //---
   m_auxMemPool.allocateMemory(nGroupRows * nGroupCols);

   //---
   //--- create models
   //---
   createModels();
   
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "initializeApp()", m_appParam.LogLevel, 2);
}

//===========================================================================//
void MMKP_DecompApp::createModels(){

   //---
   //--- This function does the work to create the different models 
   //---  that will be used. This memory is owned by the user. It will
   //---  be passed to the application interface and used by the algorithms.
   //---

   //---
   //--- Get information about this problem instance.
   //---
   int             i;
   const int       nGroupRows = m_instance.getNGroupRows();
   const int       nGroupCols = m_instance.getNGroupCols();
   const int       nKnapRows  = m_instance.getNKnapRows();
   const double *  value      = m_instance.getValue();
   const int       numCols    = nGroupRows * nGroupCols;  

   //---
   //--- Multi-Dimensional Multi-Choice Knapsack Problem (MMKP).
   //---
   //---  max  sum{i in 1..n, j in 1..l[i]}  v[i,j]   x[i,j] <==>
   //---  min  sum{i in 1..n, j in 1..l[i]} -v[i,j]   x[i,j]
   //---  s.t. sum{i in 1..n, j in 1..l[i]}  r[k,i,j] x[i,j] <= b[k], k in 1..m
   //---       sum{j in 1..l[i]}                      x[i,j]  = 1   , i in 1..n
   //---       x[i,j] in {0,1}, i in 1..n, j in 1..l[i]
   //---
   //--- Multi-Choice Knapsack Polytope (MCKP) [for a fixed k]
   //---  s.t. sum{i in 1..n, j in 1..l[i]}  r[k,i,j] x[i,j] <= b[k]
   //---       sum{j in 1..l[i]}                      x[i,j]  = 1   , i in 1..n
   //---       x[i,j] in {0,1}, i in 1..n, j in 1..l[i]
   //---
   //--- Multi-Dimensional Knapsack Polytope (MDKP) 
   //---  s.t. sum{i in 1..n, j in 1..l[i]}  r[k,i,j] x[i,j] <= b[k], k in 1..m
   //---       x[i,j] in {0,1}, i in 1..n, j in 1..l[i]
   //---

   //--- 
   //--- MMKP example structure (n=4, m=3)
   //---   xxxx <= b[1]
   //---   xxxx <= b[2]
   //---   xxxx <= b[3]
   //---   x     = 1 (i=1)
   //---    x    = 1 (i=2)
   //---     x   = 1 (i=3)
   //---      x  = 1 (i=4)
   //---
   //--- MDKP example structure
   //---   xxxx <= b[1]
   //---   xxxx <= b[2]
   //---   xxxx <= b[3]
   //---
   //--- MCKP example structure
   //---   xxxx <= b[1]
   //---   x     = 1 (i=1)
   //---    x    = 1 (i=2)
   //---     x   = 1 (i=3)
   //---      x  = 1 (i=4)
   //---
   //--- MC2KP example structure
   //---   xxxx <= b[1]
   //---   xxxx <= b[2]
   //---   x     = 1 (i=1)
   //---    x    = 1 (i=2)
   //---     x   = 1 (i=3)
   //---      x  = 1 (i=4)
   //---

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createModels()", m_appParam.LogLevel, 2);

   //---
   //--- Construct the objective function (the original problem is 
   //---  a maximization, so we flip the sign to make it minimization).
   //---
   m_objective = new double[numCols];
   if(!m_objective)
      throw UtilExceptionMemory("createModels", "MMKP_DecompApp");
   for(i = 0; i < numCols; i++)
      m_objective[i] = -value[i];
   setModelObjective(m_objective);
   
   //---
   //--- Model
   //---    relax = MCKP (with knapsack constraint  1   )
   //---    core  = MDKP (with knapsack constraints 2..m)
   //---
   //--- A' (relax) = MCKP
   //---   sum{i in 1..n, j in 1..l[i]}  r[1,i,j] x[i,j] <= b[1]
   //---   sum{j in 1..l[i]}                      x[i,j]  = 1   , i in 1..n
   //---   x[i,j] in {0,1}, i in 1..n, j in 1..l[i]
   //--- A'' (core)
   //---   sum{i in 1..n, j in 1..l[i]}  r[k,i,j] x[i,j] <= b[k], k in 2..m
   //---   x[i,j] in [0,1], i in 1..n, j in 1..l[i]
   //---
   if(m_appParam.ModelNameCore == "MDKP0"){
      DecompConstraintSet * modelCore  = new DecompConstraintSet();      
      createModelPartMDKPCompl(modelCore);
      m_models.push_back(modelCore);
      setModelCore(modelCore, m_appParam.ModelNameCore);
   }
   if(m_appParam.ModelNameRelax == "MCKP0"){
      DecompConstraintSet * modelRelax = new DecompConstraintSet();      
      createModelPartMCKP(modelRelax);
      m_models.push_back(modelRelax);
      setModelRelax(modelRelax, m_appParam.ModelNameRelax);
   }
   
   //---
   //--- Model (nested oracles)
   //---    relax     = MCKP (with knapsack constraint  1   )
   //---    core      = MDKP (with knapsack constraints 2..m)
   //---    relax*[b] = MCKP (with knapsack constraints 1,b ), b=2..m
   //---
   //--- A'[b] (relax*) = MC2KP0
   //---   sum{i in 1..n, j in 1..l[i]}  r[1,i,j] x[i,j] <= b[k], k in {1,b}
   //---   sum{j in 1..l[i]}                      x[i,j]  = 1   , i in 1..n
   //---   x[i,j] in {0,1}, i in 1..n, j in 1..l[i]
   //---
   if(m_appParam.ModelNameRelaxNest == "MC2KP0"){
      for(i = 1; i < nKnapRows; i++){
         DecompConstraintSet * modelRelax = new DecompConstraintSet();      
         createModelPartMC2KP(modelRelax, 0, i);         
         m_models.push_back(modelRelax);
         setModelRelaxNest(modelRelax, 
			   m_appParam.ModelNameRelaxNest + UtilIntToStr(i));
      }
   }
   
   //---
   //--- Model
   //---    relax = MDKP (with knapsack constraints 1..m)
   //---    core  = MCP  (with choice constraints)
   //---
   //--- A' (relax) = MDKP
   //---   sum{i in 1..n, j in 1..l[i]}  r[k,i,j] x[i,j] <= b[k], k in 1..m
   //---   x[i,j] in {0,1}, i in 1..n, j in 1..l[i]
   //--- A'' (master)
   //---   sum{j in 1..l[i]}                      x[i,j]  = 1   , i in 1..n
   //---   x[i,j] in [0,1], i in 1..n, j in 1..l[i]
   //---
   if(m_appParam.ModelNameCore == "MCP"){
      DecompConstraintSet * modelCore  = new DecompConstraintSet();      
      createModelPartMCP(modelCore);
      m_models.push_back(modelCore);
      setModelCore(modelCore, m_appParam.ModelNameCore);
   }
   if(m_appParam.ModelNameRelax == "MDKP"){      
      DecompConstraintSet * modelRelax = new DecompConstraintSet();
      createModelPartMDKP(modelRelax);
      m_models.push_back(modelRelax);
      setModelRelax(modelRelax, m_appParam.ModelNameRelax);
   }

   //---
   //--- Model
   //---    relax = MDKP (with knapsack constraints 1..floor(m/2))
   //---    relax = MDKP (with knapsack constraints 1..m) -> nested
   //---    core  = MCP  (with choice constraints)
   //---
   //--- A' (relax) = MDKP
   //---   sum{i in 1..n, j in 1..l[i]}  r[k,i,j] x[i,j] <= b[k], k in 1..m
   //---   x[i,j] in {0,1}, i in 1..n, j in 1..l[i]
   //--- A'' (master) + missing from MDKP
   //---   sum{j in 1..l[i]}                      x[i,j]  = 1   , i in 1..n
   //---   x[i,j] in [0,1], i in 1..n, j in 1..l[i]
   //---
   //--- which half? those with tight knapsack
   //---
   if(m_appParam.ModelNameCore == "MMKPHalf"){
      DecompConstraintSet * modelCore  = new DecompConstraintSet();      
      createModelPartMMKPHalf(modelCore);
      m_models.push_back(modelCore);
      setModelCore(modelCore, m_appParam.ModelNameCore);
   }
   if(m_appParam.ModelNameRelax == "MDKPHalf"){      
      DecompConstraintSet * modelRelax = new DecompConstraintSet();
      createModelPartMDKPHalf(modelRelax);
      m_models.push_back(modelRelax);
      setModelRelax(modelRelax, m_appParam.ModelNameRelax);
   }
   if(m_appParam.ModelNameRelaxNest == "MDKP"){      
      DecompConstraintSet * modelRelax = new DecompConstraintSet();
      createModelPartMDKP(modelRelax);
      m_models.push_back(modelRelax);
      setModelRelaxNest(modelRelax, m_appParam.ModelNameRelax);
   }
   if(m_appParam.ModelNameRelaxNest == "MMKP"){      
      DecompConstraintSet * modelRelax = new DecompConstraintSet();
      createModelPartMMKP(modelRelax);
      m_models.push_back(modelRelax);
      setModelRelaxNest(modelRelax, m_appParam.ModelNameRelax);
   }
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "createModels()", m_appParam.LogLevel, 2);
}

//===========================================================================//
void MMKP_DecompApp::createModelPartMDKP(DecompConstraintSet * model){
   vector<int> whichKnaps;
   int         nKnapRows = m_instance.getNKnapRows();
   UtilIotaN(whichKnaps, nKnapRows, 0);
   createModelPartMDKP(model, whichKnaps);
}

//===========================================================================//
void MMKP_DecompApp::createModelPartMDKPCompl(DecompConstraintSet * model,
                                              int                   whichKnap){
   int         i;
   vector<int> whichKnaps;
   const int   nKnapRows  = m_instance.getNKnapRows();
   for(i = 0; i < nKnapRows; i++){
      if(i == whichKnap)
         continue;
      whichKnaps.push_back(i);
   }
   createModelPartMDKP(model, whichKnaps);
}

//===========================================================================//
void MMKP_DecompApp::createModelPartMDKPHalf(DecompConstraintSet * model){
   int         i;
   vector<int> whichKnaps;
   const int   nKnapRows  = m_instance.getNKnapRows();
   const int   nHalfRows  = static_cast<int>(std::floor(nKnapRows/2.0));
   for(i = 0; i < nHalfRows; i++){
      whichKnaps.push_back(i);
   }
   createModelPartMDKP(model, whichKnaps);
}


//===========================================================================//
void MMKP_DecompApp::createModelPartMDKP(DecompConstraintSet * model,
                                         vector<int>         & whichKnaps){

   //---
   //--- Multi-Dimensional Knapsack Polytope
   //---  sum{i in 1..n, j in 1..l[i]}  r[k,i,j] x[i,j] <= b[k], k in K
   //---   
   int             i, j, colIndex;
   int             nGroupRows    = m_instance.getNGroupRows();
   int             nGroupCols    = m_instance.getNGroupCols();
   const double *  capacity      = m_instance.getCapacity();
   const double * const * weight = m_instance.getWeight();
   int             nKnapRows     = static_cast<int>(whichKnaps.size());
   int             numCols       = nGroupRows * nGroupCols;
   int             numRows       = nKnapRows;

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createModelPartMDKP()", m_appParam.LogLevel, 2);

   model->M = new CoinPackedMatrix(false, 0.0, 0.0);
   if(!model->M)
      throw UtilExceptionMemory("createModelPartMDKP", "MMKP_DecompApp");
   model->M->setDimensions(0, numCols);
   model->reserve(numRows, numCols);

   //---
   //---  sum{i in 1..n, j in 1..l[i]}  r[k,i,j] x[i,j] <= b[k], k in K
   //---      
   vector<int>::const_iterator vi;
   for(vi = whichKnaps.begin(); vi != whichKnaps.end(); vi++){
      CoinPackedVector   rowK;
      const double     * weightK = weight[*vi];
      colIndex = 0;
      for(i = 0; i < nGroupRows; i++){
         for(j = 0; j < nGroupCols; j++){
            rowK.insert(colIndex, weightK[colIndex]);
            colIndex++;
         }
      }
      model->appendRow(rowK, -DecompInf, capacity[*vi]);
   }
      
   //---   
   //--- set the col upper and lower bounds
   //---
   UtilFillN(model->colLB, numCols,  0.0);
   UtilFillN(model->colUB, numCols,  1.0);

   //---
   //--- set the indices of the integer variables of model
   //---
   UtilIotaN(model->integerVars, numCols, 0);
   
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "createModelPartMDKP()", m_appParam.LogLevel, 2);
}

//===========================================================================//
void MMKP_DecompApp::createModelPartMCP(DecompConstraintSet * model){
   vector<int> whichKnaps;   
   createModelPartMCKP(model, whichKnaps);
}

//===========================================================================//
void MMKP_DecompApp::createModelPartMCKP(DecompConstraintSet * model,
                                         int                   whichKnap){
   vector<int> whichKnaps;
   whichKnaps.push_back(whichKnap);
   createModelPartMCKP(model, whichKnaps);
}

//===========================================================================//
void MMKP_DecompApp::createModelPartMMKPHalf(DecompConstraintSet * model){
   vector<int> whichKnaps;
   int         i;
   const int   nKnapRows  = m_instance.getNKnapRows();
   const int   nHalfRows  = static_cast<int>(std::floor(nKnapRows/2.0));
   for(i = nHalfRows; i < nKnapRows; i++)
      whichKnaps.push_back(i);
   createModelPartMCKP(model, whichKnaps);
}

//===========================================================================//
void MMKP_DecompApp::createModelPartMMKP(DecompConstraintSet * model){
   vector<int> whichKnaps;
   int         i;
   const int   nKnapRows  = m_instance.getNKnapRows();
   for(i = 0; i < nKnapRows; i++)
      whichKnaps.push_back(i);
   createModelPartMCKP(model, whichKnaps);
}

//===========================================================================//
void MMKP_DecompApp::createModelPartMC2KP(DecompConstraintSet * model,
                                          int                   whichKnap1,
                                          int                   whichKnap2){
   vector<int> whichKnaps;
   whichKnaps.push_back(whichKnap1);
   whichKnaps.push_back(whichKnap2);
   createModelPartMCKP(model, whichKnaps);
}

//===========================================================================//
void MMKP_DecompApp::createModelPartMCKP(DecompConstraintSet * model,
                                         vector<int>         & whichKnaps){
   
   
   //---
   //--- Multi-Choice Knapsack Polytope [for a fixed k] (Subproblem[k]):
   //---  sum{i in 1..n, j in 1..l[i]}  r[k,i,j] x[i,j] <= b[k]
   //---  sum{j in 1..l[i]}                      x[i,j]  = 1   , i in 1..n
   //---
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createModelPartMCKP()", m_appParam.LogLevel, 2);
   
   
   int             i, j, colIndex;
   int             nGroupRows    = m_instance.getNGroupRows();
   int             nGroupCols    = m_instance.getNGroupCols();
   const double *  capacity      = m_instance.getCapacity();
   const double * const * weight = m_instance.getWeight();
   int             nKnapRows     = static_cast<int>(whichKnaps.size());
   int             numCols       = nGroupRows * nGroupCols;
   int             numRows       = nGroupRows + nKnapRows;

   //TODO: should this all be more opaque?
   model->M = new CoinPackedMatrix(false, 0.0, 0.0);
   if(!model->M)
      throw UtilExceptionMemory("createModelPartMDKP", "MMKP_DecompApp");
   model->M->setDimensions(0, numCols);   
   model->reserve(numRows, numCols);

   //---
   //---  s.t. sum{i in 1..n, j in 1..l[i]}  r[k,i,j] x[i,j] <= b[k]
   //---   
   vector<int>::const_iterator vi;
   for(vi = whichKnaps.begin(); vi != whichKnaps.end(); vi++){
      CoinPackedVector   rowK;
      const double     * weightK = weight[*vi];
      colIndex  = 0;
      for(i = 0; i < nGroupRows; i++){
         for(j = 0; j < nGroupCols; j++){
            rowK.insert(colIndex, weightK[colIndex]);
            colIndex++;
         }
      }
      model->appendRow(rowK, -DecompInf, capacity[*vi]);
   }

   //---
   //---       sum{j in 1..l[i]} x[i,j]  = 1, i in 1..n
   //---
   colIndex = 0;
   for(i = 0; i < nGroupRows; i++){
      CoinPackedVector row;
      for(j = 0; j < nGroupCols; j++){
         row.insert(colIndex, 1.0);
         colIndex++;
      } 
      model->appendRow(row, 1.0, 1.0);
   }
   
   //---   
   //--- set the col upper and lower bounds
   //---
   UtilFillN(model->colLB, numCols,  0.0);
   UtilFillN(model->colUB, numCols,  1.0);

   //---
   //--- set the indices of the integer variables of model
   //---
   UtilIotaN(model->integerVars, numCols, 0);
   
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "createModelPartMCKP()", m_appParam.LogLevel, 2);
}

//===========================================================================//
DecompSolverStatus 
MMKP_DecompApp::solveRelaxed(const int          whichBlock,
                             const double     * redCostX,
                             const double       convexDual,
                             DecompVarList    & varList){
                                            
   if(!m_appParam.UsePisinger)
      return DecompSolStatNoSolution;

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "solveRelaxed()", m_appParam.LogLevel, 2);
   
   //---
   //--- this allows user direct access to access methods in 
   //---   algorithm interface (in case they want to use any 
   //---   of its data)
   //---
   //--- for example, if the user wants to enforce the branching
   //---   decisions in the oracle
   //--- TODO: can this be done using mcknap solver?
   //---
   //const DecompAlgo * decompAlgo = getDecompAlgo();
   //const double     * colLBNode  = decompAlgo->getColLBNode();
   //const double     * colUBNode  = decompAlgo->getColUBNode();

   //---
   //--- in the case where oracle=MCKP0, we have a specialized solver
   //--- in the case where oracle=MDKP,  we do not have a specialized solver
   //---   so, we just return with no solution and the framework will 
   //---   attempt to use the built-in MILP solver
   //---
   DecompSolverStatus solverStatus = DecompSolStatNoSolution;
   if(m_appParam.ModelNameRelax == "MCKP0"){
      vector<int>           solInd;
      vector<double>        solEls;
      double                varRedCost    = 0.0;
      double                varOrigCost   = 0.0;
      MMKP_MCKnap         * mcknapK       = m_mcknap[whichBlock];
      //TODO: check status return codes here
      mcknapK->solveMCKnap(redCostX, m_objective,
                           solInd, solEls, varRedCost, varOrigCost);
      assert(static_cast<int>(solInd.size()) == m_instance.getNGroupRows());
      
      UTIL_DEBUG(m_param.LogDebugLevel, 4,
                 printf("PUSH var with k = %d RC = %g origCost = %g\n", 
                        whichBlock, varRedCost - convexDual, varOrigCost);
                 );
      
      //the user should not have to calculate orig cost too
      //  the framework can do this... in fact the framework shoudl
      //  calculate red-cost too... but the user might want to check this stuff
      
      //this is way too confusing for user to remember they need -alpha!
      //  let framework do that - also setting the block id - framework!
      DecompVar * var = new DecompVar(solInd, solEls, 
                                      varRedCost - convexDual, varOrigCost);
      var->setBlockId(whichBlock);
      varList.push_back(var);      
      solverStatus = DecompSolStatOptimal;
   }
           
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "solveRelaxed()", m_appParam.LogLevel, 2);
   return solverStatus;
}


//===========================================================================//
void MMKP_DecompApp::printOriginalColumn(const int   index, 
					 ostream   * os) const {
   pair<int,int> p = m_instance.getIndexInv(index);
   (*os) << "x[ " << index << " : " << p.first << " , " << p.second << " ]";
}

