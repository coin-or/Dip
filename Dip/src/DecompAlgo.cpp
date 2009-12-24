//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)       //
//                                                                           //
// Copyright (C) 2002-2009, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

//===========================================================================//
#include "DecompApp.h"
#include "DecompAlgo.h"
#include "DecompAlgoD.h"
#include "DecompAlgoC.h"
#include "DecompCutOsi.h"
#include "DecompAlgoCGL.h"
#include "DecompSolverResult.h"
//===========================================================================//
//#define DEBUG_SOLVE_RELAXED
//#define   RELAXED_THREADED
#define   NTHREADS   2
//#define   DO_INTERIOR //also in DecompAlgoPC

//===========================================================================//
//#define STAB_DUMERLE

//===========================================================================//
#include "OsiClpSolverInterface.hpp"
#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglMixedIntegerRounding2.hpp"

#ifdef __DECOMP_LP_CPX__
#define USE_MULTI_RAY
#endif


//===========================================================================//
#ifdef RELAXED_THREADED
static void * solveRelaxedThread(void * args);
#endif
typedef struct SolveRelaxedThreadArgs SolveRelaxedThreadArgs;
struct SolveRelaxedThreadArgs {
   DecompAlgo         * algo;
   int                   threadId;
   int                   batchSz;//blocks
   int                 * batch;  //blocks
   vector<DecompAlgoModel*>   * algoModel;
   
   int                   nBaseCoreRows;
   double              * u;
   double              * redCostX;
   const double              * origCost;
   int                   n_origCols;
   bool                  checkDup;
   bool                  doExact;
   bool                  doCutoff;
   list<DecompVar*>    * vars;
};
   
//===========================================================================//
void DecompAlgo::checkBlocksColumns(){

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "checkBlocksColumns()", m_param.LogDebugLevel, 2);

   if(m_modelRelax.size() == 0)
      return; 

   //---
   //--- sanity check that the blocks are column disjoint
   //---
   map<int, DecompAlgoModel>::iterator mid1;
   map<int, DecompAlgoModel>::iterator mid2;
   for(mid1 = m_modelRelax.begin(); mid1 != m_modelRelax.end(); mid1++){
      DecompAlgoModel     & modelRelax1 = (*mid1).second;
      DecompConstraintSet * model       = modelRelax1.getModel();
      if(!model->getMatrix())
         return;
      set<int>            & activeCols1 
         = modelRelax1.getModel()->activeColumnsS;
      for(mid2 = m_modelRelax.begin(); mid2 != m_modelRelax.end(); mid2++){
         if(mid1 == mid2)
            continue;
         DecompAlgoModel & modelRelax2 = (*mid2).second;
         set<int>        & activeCols2
            = modelRelax2.getModel()->activeColumnsS;
         set<int>          activeCols1inter2;
         set_intersection(activeCols1.begin(), activeCols1.end(),
                          activeCols2.begin(), activeCols2.end(), 
                          inserter(activeCols1inter2,
                                   activeCols1inter2.begin()));
         if(activeCols1inter2.size() > 0){
            cerr << "NOTE: the columns in block " << modelRelax1.getBlockId()
                 << " -> " << modelRelax1.getModelName() << " and block "
                 << modelRelax2.getBlockId() 
                 << " -> " << modelRelax2.getModelName() << " overlap." 
                 << endl;
            throw UtilException("Columns in some blocks overlap.",
				"checkBlocksColumns", "DecompAlgo");
         }
      }
   }
   
   //---
   //--- sanity check that the union of active columns in blocks
   //---   should cover all columns in core - if not, these are 'master-only'
   //---   columns which can be dealt with using either LD or the using the 
   //---   ideas of Rob Pratt discussion (9/27/09), or defined explicitly 
   //---   by user
   //---
   set<int> activeColsUnion;
   set<int>::iterator sit;
   for(mid1 = m_modelRelax.begin(); mid1 != m_modelRelax.end(); mid1++){
      DecompAlgoModel & modelRelax = (*mid1).second;
      set<int>        & activeCols 
         = modelRelax.getModel()->activeColumnsS;
      set_union(activeCols.begin(),       activeCols.end(),
                activeColsUnion.begin(),  activeColsUnion.end(),
                inserter(activeColsUnion, activeColsUnion.begin()));
   }

   int                     i;
   const DecompAlgoModel & modelCore      = getModelCore();
   bool                    allColsCovered = true;
   for(i = 0; i < modelCore.getModel()->getNumCols(); i++){
      sit = activeColsUnion.find(i);
      if(sit == activeColsUnion.end()){
         (*m_osLog) << "Column " << setw(5) << i << " -> " 
                    << setw(25) << modelCore.getModel()->colNames[i] 
                    << " is missing from union of blocks." << endl;
         allColsCovered = false;
      }
   }
   if(!allColsCovered)
      throw UtilException("Some columns not covered in blocks",
                          "checkBlocksColumns", "DecompAlgo");

   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "checkBlocksColumns()", m_param.LogDebugLevel, 2);
}

//===========================================================================//
void DecompAlgo::initSetup(UtilParameters * utilParam,
                           string         & sectionParam){

   //---
   //--- read in algorithm parameters
   //---
   m_param.getSettings(*utilParam, sectionParam);
   UTIL_MSG(m_param.LogLevel, 2,
	    (*m_osLog)
	    << "Initial Algo Setup"
	    << " (algo = " << DecompAlgoStr[m_algo] << ")\n";
	    );
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "initSetup()", m_param.LogDebugLevel, 2);
   if(m_param.LogLevel > 1){
      m_param.dumpSettings(sectionParam, m_osLog);        
   }

   //---
   //--- create DecompAlgoModel objects from DecompAppModel objects
   //---   these just store pointers to the models provided by user
   //---   and will store pointers to the approriate OSI objects
   //---
   getModelsFromApp();

   m_numConvexCon = static_cast<int>(m_modelRelax.size());
   assert(m_numConvexCon >= 1);

   DecompConstraintSet * modelCore = m_modelCore.getModel();   
   assert(modelCore);
   UTIL_DEBUG(m_param.LogDebugLevel, 1,	      
	      if(modelCore){
                 (*m_osLog)<< "ModelCore   cols: " << modelCore->getNumCols()
                           << " rows: "            << modelCore->getNumRows()
                           << "\n";
	      }
	      else{
                 (*m_osLog) << "ModelCore is Empty.\n";
	      }
	      );

   //---
   //--- sanity checks on user input
   //---
   checkBlocksColumns();

   //---
   //--- if we have a core, allocate a pool of memory for re-use
   //---
   if(modelCore){
      m_memPool.allocateMemory(modelCore->getNumCols(),
			       modelCore->getNumRows());
   }

   //---
   //--- By default the relaxation can be solved using a generic IP solver.
   //---
   //--- Here, for each relaxation, we initialize an OSI interface and load
   //--- the problem data.
   //---
   map<int, DecompAlgoModel>         ::iterator mit;
   map<int, vector<DecompAlgoModel> >::iterator mivt;
   vector<DecompAlgoModel>           ::iterator vit;
   for(mit = m_modelRelax.begin(); mit != m_modelRelax.end(); mit++)
      createOsiSubProblem((*mit).second);                          
   for(mivt  = m_modelRelaxNest.begin(); 
       mivt != m_modelRelaxNest.end(); mivt++){
      for(vit  = (*mivt).second.begin(); 
          vit != (*mivt).second.end(); vit++)
         createOsiSubProblem((*vit));
   }
   UTIL_DEBUG(m_param.LogDebugLevel, 1,
	      (*m_osLog) << "Number of Convexity Constraints: " 
	      << m_numConvexCon << endl;              
	      for(mit  = m_modelRelax.begin(); 
                  mit != m_modelRelax.end(); mit++){
		 DecompConstraintSet * model = (*mit).second.getModel();
                 if(model->M){
		    (*m_osLog)
		       << "ModelRelax  cols: " << model->getNumCols()
		       << " rows: "            << model->getNumRows()
		       << endl;
		 }
	      }
	      );
   
     
   //---
   //--- open memory to store the current solution (in terms of x)
   //---
   const int             nCols     = modelCore->getNumCols();
   const double        * colLB     = modelCore->getColLB();
   const double        * colUB     = modelCore->getColUB();
   assert(nCols > 0);     
   
   m_xhat      = new double[nCols];
   m_colLBNode = new double[nCols];
   m_colUBNode = new double[nCols];
   assert(m_xhat && m_colLBNode && m_colUBNode);          
   memcpy(m_colLBNode, colLB, nCols * sizeof(double));
   memcpy(m_colUBNode, colUB, nCols * sizeof(double));



   //---
   //--- PC: create an initial set of points F'[0] subseteq F' (c    + eps)
   //--- DC: create an initial set of points F'[0] subseteq F' (xhat + eps)
   //--- RC: do nothing - DecompAlgo base?? WHY - need an shat to get going
   //---  C: do nothing - DecompAlgo base
   //---
   DecompVarList initVars;
   m_nodeStats.varsThisCall += generateInitVars(initVars);
   
   //---
   //--- create the master OSI interface
   //---
   m_masterSI = new OsiLpSolverInterface();
   CoinAssertHint(m_masterSI, "Error: Out of Memory");
   m_masterSI->messageHandler()->setLogLevel(m_param.LogLpLevel);
#ifdef __DECOMP_LP_CLP__
   OsiClpSolverInterface * masterSIClp = dynamic_cast<OsiClpSolverInterface*>(m_masterSI);
   masterSIClp->getModelPtr()->setLogLevel(m_param.LogLpLevel);
#endif

   //---
   //--- init CGL object
   //---  NOTE: do not allow PC gomory cuts for now
   //---
   m_cgl = new DecompAlgoCGL(m_param.LogDebugLevel,
			     m_algo);
   m_cgl->setLogStream(m_osLog);
   m_cgl->setLogLevel (m_param.LogDebugLevel);
   m_cgl->initGenerators(m_param.CutCglClique,
			 m_param.CutCglOddHole,
			 m_param.CutCglFlowC,
			 m_param.CutCglKnapC,
			 m_param.CutCglMir,
			 m_param.CutCglGomory);
   //(m_algo == CUT) ? 
   //                    m_param.CutCglGomory : 0);
      
   //---
   //--- create master problem
   //---
   createMasterProblem(initVars);
   UTIL_MSG(m_param.LogLevel, 2,
            (*m_osLog) 
            << "Model core nCols= " << modelCore->getNumCols()
            << " nRows = "          << modelCore->getNumRows() << "\n";
            ); 

   //---
   //--- construct cutgen solver interface
   //---
   if(m_param.CutCGL){
      m_cutgenSI = new OsiClpSolverInterface();
      assert(m_cutgenSI);
      loadSIFromModel(m_cutgenSI, true);

      //---
      //--- add an objective cut to the cut generator LP
      //---      obj >= globalLB
      //--- 
      if(m_algo == PRICE_AND_CUT){
         //---
         //--- THINK:
         //--- this is causing an issue later - because packs 0's in matrix
         //---   once gets to cut generator
         //---         
         //CoinPackedVector objCut(nCols, getOrigObjective());
         CoinPackedVector objCut;
         const double * objCoeff = getOrigObjective();
         int i;
		 for(i = 0; i < m_cutgenSI->getNumCols(); i++){
            if(!UtilIsZero(objCoeff[i]))
               objCut.insert(i, objCoeff[i]);
         }         
         m_cutgenObjCutInd = m_cutgenSI->getNumRows();
         m_cutgenSI->addRow(objCut, -DecompInf, DecompInf);      
      }
   }

   //---
   //--- construct auxillary compact lp interface
   //---
   if(m_param.InitCompactSolve){
      //TODO: would be nice if we could utilize IP presolve here?
      m_auxSI = new OsiLpSolverInterface();
      assert(m_auxSI);
      loadSIFromModel(m_auxSI);
   }

   /*#ifdef STAB_DUMERLE
   //---
   //--- using the cut gen OSI, solve the initial LP
   //---   for the compact formulation to get starting duals
   //--- TODO: what if no CutCGL - need its own object
   //---
   //--- we are only going to use the duals from core as estimates
   //---   of the duals for master
   //---
   assert(m_param.CutCGL);
   m_cutgenSI->initialSolve();
   assert(m_cutgenSI->isProvenOptimal());

   const double         * dualSol  = m_cutgenSI->getRowPrice();
   const vector<string> & rowNames = m_cutgenSI->getRowNames();
   
   int r;
   for(r = 0; r < modelCore->getNumRows(); r++){
      if(fabs(dualSol[r]) > DecompEpsilon){
         if(r < static_cast<int>(rowNames.size())){
            printf("INIT DUAL FOR CORE ROW[%6d -> %25s] = %12.10f\n", 
                   r, rowNames[r].c_str(), dualSol[r]);
         }
         else
            printf("INIT DUAL[%6d] = %12.10f\n", r, dualSol[r]);
      }
   }
   #endif*/
   
   
   
   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "initSetup()", m_param.LogDebugLevel, 2);
}

//===========================================================================//
void DecompAlgo::createOsiSubProblem(DecompAlgoModel & algoModel){

   //TODO: design question, we are assuming that master solver is
   //  an LP solver and relaxed solver is an IP - it really should
   //  be a generic object and an LP or IP solver is just one option
   //  for a solver

   OsiSolverInterface  * subprobSI = NULL;
   DecompConstraintSet * model     = algoModel.getModel();
   if(!model || !model->M)
      return;

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createOsiSubProblem()", m_param.LogDebugLevel, 2);

   int nInts = model->getNumInts();
   int nCols = model->getNumCols();
   int nRows = model->getNumRows();
   
   subprobSI = new OsiLpSolverInterface();//THINK: not IpOsi?
   assert(subprobSI);
   subprobSI->messageHandler()->setLogLevel(m_param.LogLpLevel);
   
   //TODO: use assign vs load? just pass pointers?
   subprobSI->loadProblem(*model->getMatrix(), 
                          model->getColLB(), 
                          model->getColUB(),
                          NULL, //null objective
                          model->getRowLB(),
                          model->getRowUB());
   if(nInts > 0)
      subprobSI->setInteger(model->getIntegerVars(), nInts);
         
   //---
   //--- set column and row names (if they exist)
   //---  
   vector<string> & colNames = model->colNames;
   vector<string> & rowNames = model->rowNames;
   subprobSI->setIntParam(OsiNameDiscipline, 1);//1=Lazy
   if(colNames.size())
      subprobSI->setColNames(colNames, 0, nCols, 0);
   if(rowNames.size())
      subprobSI->setRowNames(rowNames, 0, nRows, 0);   

   //---
   //--- set subproblem pointer
   //---
   algoModel.setOsi(subprobSI);   

   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "createOsiSubProblem()", m_param.LogDebugLevel, 2);

}

//===========================================================================//
void DecompAlgo::getModelsFromApp(){
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "getModelsFromApp()", m_param.LogDebugLevel, 2);

   //---
   //--- check to make sure certain things are set
   //---
   if(!m_app->m_objective)
      throw UtilException("Application objective has not been set", 
                          "getModelsFromApp", "DecompAlgo");
   if(!m_app->m_modelCore.getModel())
      throw UtilException("Application core constraint set has not been set", 
                          "getModelsFromApp", "DecompAlgo");

   m_objective = m_app->m_objective;
   m_modelCore = m_app->m_modelCore;

   map<int, DecompAppModel>         ::iterator mit;
   map<int, vector<DecompAppModel> >::iterator mivt;
   vector<DecompAppModel>           ::iterator vit;
   for(mit  = m_app->m_modelRelax.begin(); 
       mit != m_app->m_modelRelax.end(); mit++){
      DecompAlgoModel algoModel = (*mit).second;
      m_modelRelax.insert(make_pair((*mit).first, algoModel));
   }                   
   for(mivt  = m_app->m_modelRelaxNest.begin(); 
       mivt != m_app->m_modelRelaxNest.end(); mivt++){
      vector<DecompAlgoModel> v;
      for(vit  = (*mivt).second.begin(); 
          vit != (*mivt).second.end(); vit++)
         v.push_back(*vit);
      m_modelRelaxNest.insert(make_pair((*mivt).first, v));
   }
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "getModelsFromApp()", m_param.LogDebugLevel, 2);
}

//===========================================================================//
void DecompAlgo::loadSIFromModel(OsiSolverInterface           * si,
                                 bool                           doInt){
   
   //---
   //--- Initialize the solver interface.
   //---    min c(x)
   //---     A' x   >= b'  [optional]
   //---     A''x   >= b'' 
   //---     l <= x <= u
   //---
   //--- relaxV contains [A',  b' ] in terms of x (if explicit)
   //--- core   contains [A'', b''] in terms of x
   //---

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "loadSIFromModel()", m_param.LogDebugLevel, 2);
  
   DecompConstraintSet * core  = m_modelCore.getModel();
   DecompConstraintSet * relax = NULL;
   int nCols  = core->getNumCols();
   int nRowsC = core->getNumRows();
   int nRowsR = 0;
   int nRows  = nRowsC;

   //---
   //--- create matrix from core matrix
   //---
   CoinPackedMatrix * M = new CoinPackedMatrix(*core->M);    
   assert(M);

   //---
   //--- append to bottom the relax matrix/matrices
   //---  create block file (for use in MILPBlock app)
   //---
#ifdef CREATE_BLOCKFILE
   ofstream os("blockFile.txt");
#endif
   map<int, DecompAlgoModel>::iterator mit;
   for(mit  = m_modelRelax.begin(); mit != m_modelRelax.end(); mit++){
      relax = (*mit).second.getModel();
      if(relax->M){
	 nRowsR  = relax->getNumRows();

#ifdef CREATE_BLOCKFILE         
         {
            int r;
            os << (*mit).second.getBlockId();
            os << " " << nRowsR << endl;
            for(r = 0; r < nRowsR; r++){
               os << nRows + r << " ";
            }
            os << endl;
         }
#endif

	 nRows  += nRowsR;
         M->bottomAppendPackedMatrix(*relax->M);
      }
   }
#ifdef CREATE_BLOCKFILE         
   os.close();
#endif

   //---
   //--- set column and row bounds and objective coeffs
   //---   
   double * colLB    = new double[nCols];
   double * colUB    = new double[nCols];
   double * objCoeff = new double[nCols];
   double * rowLB    = new double[nRows];
   double * rowUB    = new double[nRows];
   assert(colLB && colUB && objCoeff && rowLB && rowUB);
   
   memcpy(colLB,    core->getColLB(),   nCols  * sizeof(double));
   memcpy(colUB,    core->getColUB(),   nCols  * sizeof(double));
   memcpy(objCoeff, getOrigObjective(), nCols  * sizeof(double));  
   memcpy(rowLB,    core->getRowLB(),   nRowsC * sizeof(double));
   memcpy(rowUB,    core->getRowUB(),   nRowsC * sizeof(double));
   int rowIndex = nRowsC;
   for(mit  = m_modelRelax.begin(); mit != m_modelRelax.end(); mit++){
      relax = (*mit).second.getModel();
      if(relax->M){
	 nRowsR  = relax->getNumRows();
	 memcpy(rowLB + rowIndex, relax->getRowLB(),
		nRowsR * sizeof(double));
	 memcpy(rowUB + rowIndex, relax->getRowUB(), 
		nRowsR * sizeof(double));
	 rowIndex += nRowsR;
      }
   }

   //---
   //--- assign problem pointers to OSI object (OSI will delete this memory)
   //---
   si->assignProblem(M, colLB, colUB, objCoeff, rowLB, rowUB);

   //---
   //--- set integer variables
   //---
   //---  Due to a design issue with OsiCpx/Xxx, if we declare 
   //---  integers in the master formulation, we have issues retrieving
   //---  any information. So, the master formulation must be the 
   //---  continuous relaxation (even in the standard CPM case).
   //---  
   //---  However, the CGL framework expects integer arguments. So,
   //---  even in the CPM case (like the PC case) we are going to need
   //---  to carry around a second copy of the core in an OSI object
   //---  we call the cutgenSI. This will have to be kept up to date in
   //---  both cases CPM and PC. This also is a problem for gomory cuts
   //---  or any cut generator that depends on the LP solver specific 
   //---  information (like basis for gomory).
   //---
   //---  TODO: we might be able to get around gomory in future by 
   //---  setting the basis in cutgenSI (from masterSI) during CPM and
   //---  taking current recomposed point and doing a crossover to basis
   //---  in PC case.
   //---
   if(doInt){
      int nInts = core->getNumInts();
      if(nInts > 0)
	 si->setInteger(&core->getIntegerVars()[0], nInts);
   }

   //---
   //--- set column and row names
   //---
   si->setIntParam(OsiNameDiscipline, 1);//1=Lazy
   vector<string> & colNames  = core->colNames;
   vector<string> & rowNamesC = core->rowNames;
   if(colNames.size())
      si->setColNames(colNames,  0, nCols, 0);
   if(rowNamesC.size())
      si->setRowNames(rowNamesC, 0, nRowsC, 0);
   rowIndex = nRowsC;
   for(mit  = m_modelRelax.begin(); mit != m_modelRelax.end(); mit++){
      relax = (*mit).second.getModel();
      if(relax->M){
         vector<string> & rowNamesR = relax->rowNames;
	 nRowsR = relax->getNumRows();
	 if(rowNamesR.size())
	    si->setRowNames(rowNamesR, 0, nRowsR, rowIndex);
         rowIndex += nRowsR;    
      }
   }
   
   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "loadSIFromModel()", m_param.LogDebugLevel, 2);
}

//===========================================================================//
//think - if need for RC?
//especially with Phase I this is very PC specific
void DecompAlgo::createMasterProblem(DecompVarList & initVars){

   //---
   //--- Initialize the solver interface for the master problem.
   //---
   //--- For the master constraint system:
   //---  m_modelCore  contains [A'', b''], in terms of x.
   //---  m_modelRelax contains [A', b'],   and might contiain multiple blocks.
   //---
   //--- For each block we must add a convexity constraint. Let K be the set
   //--- of blocks.
   //---
   //--- Notation:
   //---  n       = orig number of vars
   //---  m''     = orig number of rows in A'', b''
   //---  |K|     = the number of blocks that defines [A' b']
   //---  s       = a solution (e.p.) to the relaxed problem, size (1xn)
   //---  c       = original cost, size (1xn)
   //---  F'[k]   = the current set of relaxed e.p.'s for block k in K
   //--- a''[i,j] = entry at row i, column j for A'' matrix
   //---  C       = original set of columns (n = |C|)
   //---  R''     = original set of rows in A'' (m''=|R''|)
   //---  
   //--- The Dantzig-Wolfe LP:
   //---
   //--- min  sum{k in K, s in F'[k]}
   //---        sum{j in C}(c[j]   * s[j]) * lambda[k][s]
   //--- s.t. sum{k in K, s in F'[k]}
   //---        sum{j in C}(a''[i,j] * s[j])* lambda[k][s] ~ b''[i],  i in R''
   //---      sum{s in F'[k]} lambda[k][s] = 1, k in K
   //---      lambda[k][s]                >= 0, k in K, s in F'[k]
   //---
   //--- NOTE: if 0 is feasible to subproblem, we can relax convexity to <= 1
   //---   

   //---
   //--- Change for Phase I model.
   //---   Add a slack and/or surplus variable to each master constraint
   //---   including the bounds for branching?? THINK... 
   //---
   //--- THINK:
   //--- Do we bother removing these vars once feasible? What about the 
   //--- fact that adding cuts could once again cause infeasible.... 
   //---
   //--- What do we do after a branching? jump back to Phase I?
   //---
   //---
   //--- Phase I:
   //--- min  sum{i in R''} (splus[i] + sminus[i])
   //---
   //--- Phase II:
   //--- min  sum{k in K, s in F'[k]}
   //---        sum{j in C}(c[j]   * s[j]) * lambda[k][s]
   //---
   //--- s.t. sum{k in K, s in F'[k]}
   //---        sum{j in C}(a''[i,j] * s[j])* lambda[k][s] +
   //---          splus[i] - sminus[i] ~ b''[i],  i in R''
   //---      sum{s in F'[k]} lambda[k][s] = 1, k in K
   //---      lambda[k][s]                >= 0, k in K, s in F'[k]
   //---      splus[i]                    >= 0, i in R''
   //---      sminus[i]                   >= 0, i in R''
   //---

  
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createMasterProblem()", m_param.LogDebugLevel, 2);

   DecompConstraintSet * modelCore = m_modelCore.getModel();   
   assert(modelCore); //TODO - can core be empty?   


   int r, c, startRow, endRow;
   int nColsCore = modelCore->getNumCols();
   int nRowsCore = modelCore->getNumRows();
   int nIntVars  = modelCore->getNumInts();
   int nInitVars = static_cast<int>(initVars.size());
   assert(initVars.size() > 0);
   
   double * dblArrNCoreCols = new double[nColsCore];
   assert(dblArrNCoreCols);

   //---
   //--- set the row counts
   //---
   m_nRowsOrig   = nRowsCore;
   m_nRowsBranch = 2 * nIntVars;
   m_nRowsConvex = m_numConvexCon;
   m_nRowsCuts   = 0;
   
   //---
   //--- set the row types of the rows
   //---   original rows, branch rows, convexity rows
   //---
   UtilFillN(m_masterRowType, m_nRowsOrig,   DecompRow_Original);
   UtilFillN(m_masterRowType, m_nRowsBranch, DecompRow_Branch);
   UtilFillN(m_masterRowType, m_nRowsConvex, DecompRow_Convex);

   //---
   //--- In order to implement simple branching, we are going to
   //--- treat all column bounds as explicit constraints. Then branching
   //--- for DW can be done in the same way it is done for regular CPM.
   //---   We want to add these directly to the core so as to facilitate
   //---   operations to expand rows. Basically, we treat these just like cuts.
   //---
   coreMatrixAppendColBounds();
   nRowsCore = modelCore->getNumRows();
   
   //THINK - what need this for?
   //number of original core rows
   modelCore->nBaseRowsOrig = modelCore->nBaseRows;
   //number of original core rows plus branching rows
   modelCore->nBaseRows     = modelCore->getNumRows();

   //---
   //--- create a matrix for the master LP
   //---  make room for original rows, branching rows and convexity rows
   //---
   int                nRows    = m_nRowsOrig + m_nRowsBranch + m_nRowsConvex;
   int                nColsMax = nInitVars  
      + 2 * (m_nRowsOrig + m_nRowsBranch); 
   double           * colLB    = new double[nColsMax];
   double           * colUB    = new double[nColsMax];
   double           * objCoeff = new double[nColsMax];   
   double           * denseCol = new double[nRows];
   CoinPackedMatrix * masterM  = new CoinPackedMatrix(true, 0, 0);
   vector<string>     colNames;
   assert(colLB && colUB && objCoeff && denseCol && masterM);
   
   //---
   //--- set the number of rows, we will add columns
   //---
   masterM->setDimensions(nRows, 0);
   
   //---
   //--- create artifical columns in master LP for:
   //---  original rows
   //---  branching rows
   //---
   startRow = 0;
   endRow   = m_nRowsOrig;  
   masterMatrixAddArtCols(masterM, 
			  colLB, 
			  colUB, 
			  objCoeff, 
			  colNames,
			  startRow, endRow, 'O');
   startRow = m_nRowsOrig;
   endRow   = m_nRowsOrig + m_nRowsBranch;
   masterMatrixAddArtCols(masterM, 
			  colLB, 
			  colUB, 
			  objCoeff, 
			  colNames,
			  startRow, endRow, 'B');
   
      
   int colIndex     = 0;
   int blockIndex   = 0;
   DecompVarList::iterator li;
   //TODO:
   //  this should be calling a function to add var to lp so don't dup code
   //TODO:
   //  check for duplicates in initVars
   for(li = initVars.begin(); li != initVars.end(); li++){
      
      //---
      //--- appending these variables (lambda) to end of matrix
      //---   after the artificials
      //---
      colIndex         = masterM->getNumCols();
      m_colIndexUnique = colIndex;

      //---
      //--- store the col index for this var in the master LP
      //---   NOTE: if we remove columns, this will be wrong
      //---
      (*li)->setColMasterIndex(colIndex);
      
      //---
      //--- we expect the user to define the block id in the var object
      //---
      blockIndex = (*li)->getBlockId();

      //---
      //--- give the column a name
      //---
      string colName = "lam(c_" + UtilIntToStr(m_colIndexUnique)
	 + ",b_" + UtilIntToStr(blockIndex) + ")";
      colNames.push_back(colName);
      
      UTIL_DEBUG(m_param.LogDebugLevel, 5,
		 (*li)->print(m_osLog, m_app);
		 );
	      
      //--- 
      //--- get dense column = A''s, append convexity constraint on end 
      //---   this memory is re-used, so be sure to clear out
      //---
      //STOP: see addVarsToPool - here init, so no cuts to deal with
      // but we don't use dense array here - make a difference?
      //modelCore->M->times((*li)->m_s, denseCol);
      (*li)->fillDenseArr(modelCore->getNumCols(),
			  dblArrNCoreCols);
      modelCore->M->times(dblArrNCoreCols, denseCol);


      UtilFillN(denseCol + nRowsCore, m_numConvexCon, 0.0);
      assert(blockIndex >= 0);
      assert(blockIndex < m_numConvexCon);
      denseCol[nRowsCore + blockIndex] = 1.0;

      //---
      //--- create a sparse column from the dense column
      //---
      // THINK: do i need a DecompCol?
      // THINK: does this allocate memory for coinpackedvec twice?
      CoinPackedVector * sparseCol 
	 = UtilPackedVectorFromDense(nRowsCore + m_numConvexCon,
				     denseCol, m_param.TolZero);
      UTIL_DEBUG(m_param.LogDebugLevel, 5,
		 (*m_osLog) << "\nSparse Col: \n";
		 UtilPrintPackedVector(*sparseCol, m_osLog);
		 );
      
      //TODO: check for duplicates (against m_vars)
      //      or force initVars to be sent in with no dups?	 	 
      //--- 
      //--- append the sparse column to the matrix
      //---  
      masterM->appendCol(*sparseCol);
      colLB[colIndex]    = 0.0;      
      colUB[colIndex]    = DecompInf;
      objCoeff[colIndex] = 0.0;       //PHASE I

      //---
      //--- set master column type
      //---
      m_masterColType.push_back(DecompCol_Structural);
	 
      //---
      //--- clean-up 
      //---
      UTIL_DELPTR(sparseCol); 
   } //END: for(li = initVars.begin(); li != initVars.end(); li++)
   
   //---
   //--- insert the initial set of variables into the master variable list
   //---
   //THINK: now doing in loop, so can check for dups
   appendVars(initVars);
      
   //---
   //--- THINK: do we want to adjust modelCore directly here?
   //---   adjust row bounds for convexity constraint 
   //---
   //TODO: in memory
   double * zeroSol    = new double[nColsCore];
   assert(zeroSol);

   UtilFillN(zeroSol, nColsCore, 0.0);
   //TODO - REVISIT - that's not the right check
   //  needs to be feasible to subproblem?
   bool     isZeroFeas = isIPFeasible(zeroSol);     
   UTIL_DEBUG(m_param.LogDebugLevel, 5,
	      if(isZeroFeas)
		 (*m_osLog) << "Zero Sol is Feasible - relax convexity con.\n";
	      );
   
   //---
   //--- row bounds from core inclding
   //---   original rows
   //---   branching rows
   //---
   vector<double> masterRowLB(modelCore->rowLB);
   vector<double> masterRowUB(modelCore->rowUB);
   //---
   //--- row bounds for convexity constraints
   //---
   if(isZeroFeas){
      for(r = 0; r < m_numConvexCon; r++){
	 masterRowLB.push_back(-DecompInf);
	 masterRowUB.push_back(1.0);     
      }
   }
   else{
      for(r = 0; r < m_numConvexCon; r++){
	 masterRowLB.push_back(1.0);	   
	 masterRowUB.push_back(1.0);
      }
   }

   //---
   //--- load the problem into master's solver interface
   //---
   assert(masterM->getNumRows() == static_cast<int>(masterRowLB.size()));
   assert(masterM->getNumRows() == static_cast<int>(masterRowUB.size()));
   assert(masterM->getNumRows() == static_cast<int>(m_masterRowType.size()));
   assert(masterM->getNumCols() == static_cast<int>(m_masterColType.size()));
   m_masterSI->loadProblem(*masterM,  
			   colLB, colUB, objCoeff, 
			   &masterRowLB[0], 
			   &masterRowUB[0]);
          
   //---
   //--- load column and row names to OSI
   //---
   int nRowNames = static_cast<int>(modelCore->rowNames.size());
   int nColNames = static_cast<int>(colNames.size());
   if(nRowNames || nColNames)
      m_masterSI->setIntParam(OsiNameDiscipline, 2); //Full-Names
   if(nRowNames > 0){      
      assert(nRowNames == modelCore->getNumRows());
      m_masterSI->setRowNames(modelCore->rowNames, 0, nRowNames, 0);
      vector<string> conRowNames;
      for(r = 0; r < m_numConvexCon; r++){
	 string rowName = "conv(b_" + UtilIntToStr(r) + ")";         
	 conRowNames.push_back(rowName);
      }
      m_masterSI->setRowNames(conRowNames, 0, m_numConvexCon, nRowNames);
   }
   if(nColNames > 0)
      m_masterSI->setColNames(colNames, 0, nColNames, 0);

   //TODO: make a function
   UTIL_DEBUG(m_param.LogDebugLevel, 4,                 
              for(r = 0; r < m_masterSI->getNumRows(); r++){	
                 const string rowN = m_masterSI->getRowName(r);
                 (*m_osLog) << "Row[" << setw(4) << r << "] Name: "
                            << setw(30) << rowN << " Type: "
                            << setw(20) << DecompRowTypeStr[m_masterRowType[r]]
                            << endl;
              }
              for(c = 0; c < m_masterSI->getNumCols(); c++){	
                 const string colN = m_masterSI->getColName(c);
                 (*m_osLog) << "Col[" << setw(4) << c << "] Name: "
                            << setw(30) << colN << " Type: "
                            << setw(20) << DecompColTypeStr[m_masterColType[c]]
                            << endl;
              }
              );
              
   //---
   //--- reset unique col index id
   //---
   m_colIndexUnique = masterM->getNumCols();
      
   //---
   //--- free local memory
   //---  
   UTIL_DELPTR(masterM);
   UTIL_DELARR(denseCol);
   UTIL_DELARR(colLB);
   UTIL_DELARR(colUB);
   UTIL_DELARR(objCoeff);
   UTIL_DELARR(zeroSol);
   UTIL_DELARR(dblArrNCoreCols);
   
   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "createMasterProblem()", m_param.LogDebugLevel, 2);
}

//===========================================================================//
void DecompAlgo::masterMatrixAddArtCol(CoinPackedMatrix * masterM,
                                       char              LorG,
                                       int               rowIndex,
                                       int               colIndex,
                                       DecompColType     colType,
                                       double          & colLB,
                                       double          & colUB,
                                       double          & objCoeff){
   CoinPackedVector artCol;
   if(LorG == 'L')
      artCol.insert(rowIndex, -1.0);
   else
      artCol.insert(rowIndex,  1.0);
   masterM->appendCol(artCol);
   colLB    = 0.0;
   colUB    = DecompInf;
   objCoeff = 1.0;       
   m_masterColType.push_back(colType);
   m_masterArtCols.push_back(colIndex);
}

//===========================================================================//
void DecompAlgo::masterMatrixAddArtCols(CoinPackedMatrix * masterM,
                                        double          * colLB,
                                        double          * colUB,
                                        double          * objCoeff,
                                        vector<string>  & colNames,
                                        int               startRow,
                                        int               endRow,
                                        char              origOrBranch){
   
   //---
   //--- min sp + sm
   //---
   //--- ax  = b --> ax + sp - sm  = b, sp >= 0, sm >= 0 
   //--- ax <= b --> ax      - sm <= b,          sm >= 0 
   //--- ax >= b --> ax + sp      >= b, sp >= 0
   //---
   int              r, colIndex;
   DecompConstraintSet * modelCore = m_modelCore.getModel();
   vector<char>   & rowSense  = modelCore->rowSense;
   vector<string> & rowNames  = modelCore->rowNames;
   bool             hasNames  = rowNames.empty() ? false : true;
   string           colName;
   string           strIndex;
   string           colNameL  = origOrBranch == 'O' ? "sOL(c_" : "sBL(c_";
   string           colNameG  = origOrBranch == 'O' ? "sOG(c_" : "sBG(c_";
   DecompColType    colTypeL  = origOrBranch == 'O' ? 
      DecompCol_ArtForRowL : DecompCol_ArtForBranchL;
   DecompColType    colTypeG  = origOrBranch == 'O' ? 
      DecompCol_ArtForRowG : DecompCol_ArtForBranchG;
   
   colIndex = masterM->getNumCols();
   for(r = startRow; r < endRow; r++){
      if(hasNames)
	 strIndex = UtilIntToStr(colIndex);
      switch(rowSense[r]){
      case 'L':
	 masterMatrixAddArtCol(masterM, 'L', r, colIndex, colTypeL,
			       colLB[colIndex], colUB[colIndex], 
			       objCoeff[colIndex]);
	 if(hasNames){
	    colName = colNameL + strIndex + "_" + rowNames[r] + ")";
	    colNames.push_back(colName);
	 }
         m_artColIndToRowInd.insert(make_pair(colIndex, r));
	 colIndex++;
	 break;
      case 'G':
	 masterMatrixAddArtCol(masterM, 'G', r, colIndex, colTypeG,
			       colLB[colIndex], colUB[colIndex],
			       objCoeff[colIndex]);
	 if(hasNames){
	    colName = colNameG + strIndex + "_" + rowNames[r] + ")";
	    colNames.push_back(colName);
	 }
         m_artColIndToRowInd.insert(make_pair(colIndex, r));
	 colIndex++;
	 break;
      case 'E':
	 masterMatrixAddArtCol(masterM, 'L', r, colIndex, colTypeL,
			       colLB[colIndex], colUB[colIndex], 
			       objCoeff[colIndex]);
	 if(hasNames){
	    colName = colNameL + strIndex + "_" + rowNames[r] + ")";
	    colNames.push_back(colName);
	 }
         m_artColIndToRowInd.insert(make_pair(colIndex, r));
	 colIndex++;
	 masterMatrixAddArtCol(masterM, 'G', r, colIndex, colTypeG,
			       colLB[colIndex], colUB[colIndex],
			       objCoeff[colIndex]);
	 if(hasNames){
	    colName = colNameG + strIndex + "_" + rowNames[r] + ")";
	    colNames.push_back(colName);
	 }
         m_artColIndToRowInd.insert(make_pair(colIndex, r));
	 colIndex++;
	 break;
      default:
	 assert(0);
      }
   }   
}

//===========================================================================//
void DecompAlgo::coreMatrixAppendColBounds(){
   //---
   //--- In order to implement simple branching, we are going to
   //--- treat all column bounds as explicit constraints. Then branching
   //--- for DW can be done in the same way it is done for regular CPM.
   //---
   //--- THINK: this needs some investigation. In some cases, this is not a 
   //--- great idea for performance. But, the advantage is in ease of 
   //--- implementation. The user does not need to do any sort of specialzed 
   //--- branching for DW.
   //---
   //--- NOTE: this idea won't work for identical subproblem case
   //---
   int       i, j;
   char      sense;
   double    rhs;
   bool      doNames = true; //TODO: make an option
   
   DecompConstraintSet * modelCore   = m_modelCore.getModel();
   const int             nIntVars    = modelCore->getNumInts();
   const double        * colLBCore   = modelCore->getColLB();
   const double        * colUBCore   = modelCore->getColUB(); 
   const int           * integerVars = modelCore->getIntegerVars();
   vector<string>      & colNames    = modelCore->getColNamesMutable();
   vector<string>      & rowNames    = modelCore->getRowNamesMutable();


   //TODO: use mem pool? or just create block (identity) if doing PC?   
   const int   numRows   = 2 * nIntVars;
   int       * rowStarts = new int[numRows + 1];   
   int       * rowInd    = new int[numRows];
   double    * rowEls    = new double[numRows];
   assert(rowStarts && rowInd && rowEls);

   //---
   //--- first  nColsCore rows are x <= u
   //--- second nColsCore rows are x >= l
   //---
   rowStarts[0] = 0;
   for(i = 0; i < numRows; i++){
      if(i < nIntVars){
	 j = integerVars[i];
	 //x <= u
	 rowStarts[i+1] = rowStarts[i] + 1;
	 rowInd[i]      = j;
	 rowEls[i]      = 1.0;
      }
      else{
	 //x >= l
	 j = integerVars[i - nIntVars];
	 rowStarts[i+1] = rowStarts[i] + 1;
	 rowInd[i]      = j;
	 rowEls[i]      = 1.0;
      }
   }

   //---
   //--- append as actual rows to A'' (duals used in pricing)
   //---
   modelCore->M->appendRows(numRows, rowStarts, rowInd, rowEls);

   //---
   //--- now convert to sense for hashing
   //---
   for(i = 0; i < numRows; i++){
      if(i < nIntVars){
	 //x <= u
	 j = modelCore->integerVars[i];
	 modelCore->rowLB.push_back(-DecompInf);
	 modelCore->rowUB.push_back(colUBCore[j]);
	 sense = 'L';
	 rhs   = colUBCore[j];
	 if(doNames){
	    string rowName = "ub(" + colNames[j] + ")";
	    rowNames.push_back(rowName);
	 }
      }
      else{
	 //x >= l
	 j = modelCore->integerVars[i - nIntVars];
	 modelCore->rowLB.push_back(colLBCore[j]);
	 modelCore->rowUB.push_back(DecompInf);
	 sense = 'G';
	 rhs   = colLBCore[j];
	 if(doNames){
	    string rowName = "lb(" + colNames[j] + ")";
	    rowNames.push_back(rowName);
	 }
      }
      modelCore->rowRhs.push_back(rhs);
      modelCore->rowSense.push_back(sense);     
      assert(sense != 'R');
      assert(sense != 'N');
     
      string rowHash = UtilCreateStringHash(1,
					    rowInd + i,
					    rowEls + i,
					    sense, rhs);      
      modelCore->rowHash.push_back(rowHash);      
   }
   
   UTIL_DELARR(rowStarts);
   UTIL_DELARR(rowInd);
   UTIL_DELARR(rowEls);   
}

//===========================================================================//
DecompStatus DecompAlgo::processNode(const int    nodeIndex,
                                     const double globalLB,
                                     const double globalUB){  
   
   
   double                mostNegRC = 0.0;
   DecompConstraintSet * modelCore = m_modelCore.getModel();
   m_stabEpsilon = 0.0;

   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "processNode()", m_param.LogDebugLevel, 1);
   
   UTIL_MSG(m_param.LogLevel, 2,
	    double gap = DecompInf;
	    if(globalLB > -DecompInf && globalUB < DecompInf){
	       if(globalUB != 0.0)
		  gap = fabs(globalUB-globalLB)/fabs(globalUB);
	       else
		  gap = fabs(globalUB-globalLB);
	    }
            (*m_osLog)
            << "Process Node " << nodeIndex
            << " (algo = "     << DecompAlgoStr[m_algo] 
            << ", phaseLast = "<< DecompPhaseStr[m_phaseLast]
	    << ") gLB = "      << UtilDblToStr(globalLB) 
	    << " gUB = "       << UtilDblToStr(globalUB) 
	    << " gap = "       << UtilDblToStr(gap,5) 
            << " time = "      << UtilDblToStr(globalTimer.getCpuTime(), 3)
            << endl;
	    );
   
   //---
   //--- init status
   //---
   m_useInitLpDuals = true;
   m_status         = STAT_UNKNOWN;
   m_globalLB       = globalLB;
   m_globalUB       = globalUB;

   //---
   //--- check SolveMasterAsIp setting
   //---   on by default, but if only one block, turn off
   //---
   if(m_numConvexCon == 1)
      m_param.SolveMasterAsIp = 0;
   
   //---
   //--- init stats and timer
   //--- 
   m_stats.timerDecomp.reset();   
   m_nodeStats.init();
   m_nodeStats.nodeIndex      = nodeIndex;
   m_nodeStats.objBest.first  = globalLB;
   m_nodeStats.objBest.second = globalUB;
   m_compressColsLastPrice    = 0;
   m_compressColsLastNumCols  = m_masterSI->getNumCols();
   m_phaseIObj.clear();
   
   //---
   //--- get initial phase
   //--- 
   //--- CPM            <-- CUT
   //--- PC  (node = 0) <-- PRICEI 
   //---     (node > 0) <-- <same as last node>
   //---
   phaseInit(m_phaseLast);
   m_phase = m_phaseLast;

   //---
   //--- it is possible that phaseInit can find
   //---  the node infeasible
   //---
   if(m_phase == PHASE_DONE)
      m_status = STAT_INFEASIBLE;
   else{
      if(m_param.LogDumpModel > 1)
         printCurrentProblem(m_masterSI,
                             "masterProb",
                             m_nodeStats.nodeIndex,
                             m_nodeStats.cutCallsTotal,
                             m_nodeStats.priceCallsTotal);
      
      //---
      //--- find the initial solution (dual and/or primal)
      //---
      m_status = solutionUpdate(m_phase);
   }

   if(m_status != STAT_INFEASIBLE){
      //for CPM, can't this just access from m_masterSI?
      recomposeSolution(getMasterColSolution(), m_xhat);
      UTIL_DEBUG(m_param.LogDebugLevel, 4,
                 m_app->printOriginalSolution(modelCore->getNumCols(),
                                              modelCore->getColNames(),
                                              m_xhat);
                 );

      //TODO: solution pool?
      //TODO: check if this is IP feasible
      //   make that a function
      if(isIPFeasible(m_xhat)){
	 if(m_app->APPisUserFeasible(m_xhat,
				     modelCore->getNumCols(),
				     m_param.TolZero)){
	    //printf("m_xhat is APP FEASIBLE, m_xhatIPFeas size = %d\n",
            //   (int)m_xhatIPFeas.size());
	    
	    //check for dup sol                  
	    bool isDup = m_xhatIPFeas.size() > 0 ? true : false;
	    vector<DecompSolution*>::iterator vit;
	    for(vit  = m_xhatIPFeas.begin();
		vit != m_xhatIPFeas.end(); vit++){
	       const DecompSolution * xhatIPFeas = *vit;
	       const double         * values
		  = xhatIPFeas->getValues();
	       for(int c = 0; c < modelCore->getNumCols(); c++){
		  if(!UtilIsZero(values[c] - m_xhat[c])){
		     isDup = false;
		     break;
		  }                           
	       }
	    }
	    
	    if(isDup){
	       printf("IS DUP, not pushing\n");
	    }
	    else{
	       DecompSolution * decompSol
		  = new DecompSolution(modelCore->getNumCols(),
				       m_xhat,
				       getOrigObjective());
	       //getMasterObjValue());
	       //solution pool?
	       m_xhatIPFeas.push_back(decompSol);
	       //printf("m_xhatIPFeas size = %d\n", 
               //     (int)m_xhatIPFeas.size());
	    }
	 }
	 vector<DecompSolution*>::iterator vi;
	 DecompSolution * viBest = NULL;
	 double bestBoundUB = m_nodeStats.objBest.second;
	 for(vi = m_xhatIPFeas.begin(); vi != m_xhatIPFeas.end(); vi++){
	    const DecompSolution * xhatIPFeas = *vi; 
	    if(xhatIPFeas->getQuality() <= bestBoundUB){
	       bestBoundUB = xhatIPFeas->getQuality();
	       viBest = *vi;
	    }
	 }
	 if(viBest){
	    //save the best
	    setObjBoundUB(bestBoundUB);
	    m_xhatIPBest = viBest;
	 }
      }
   
      //for CPM, dont' we need to update obj lb here? in case that no cuts
      //  are found, then node is done and we need to update boundx
      if(m_algo == CUT)
	 updateObjBoundLB();
   }

   
   //---
   //--- main processing loop
   //---
   while(m_phase != PHASE_DONE){

      UTIL_MSG(m_param.LogLevel, 2,
	       (*m_osLog)
	       << "Processing Node " << nodeIndex
	       << " (algo = "     << DecompAlgoStr[m_algo] 
	       << ", phase = "    << DecompPhaseStr[m_phase]
	       << ") thisLB = "      
	       << UtilDblToStr(m_nodeStats.objBest.first)
	       << " globalUB = " 
	       << UtilDblToStr(m_nodeStats.objBest.second) 
	       << endl;
	       );

      
      //---
      //--- update the phase based on parms, status and current phase
      //---
      phaseUpdate(m_phase, m_status);

      //---
      //--- check if we have exceeded time
      //---   THINK: should this check be in phaseUpdate?
      //---
      if(m_stats.timerOverall.isPast(m_param.LimitTime)){
	 UTIL_MSG(m_param.LogLevel, 2,
		  (*m_osLog)
		  << "Node " << nodeIndex << " process stopping on time."
		  << endl;);
	 m_stopCriteria = DecompStopTime;
	 m_phase        = PHASE_DONE;
      }
      
      //---
      //--- if the lower bound meets the global ub, we are done
      //---   
      if(m_nodeStats.objBest.first >= 
	 (m_nodeStats.objBest.second - DecompEpsilon)){
	 UTIL_MSG(m_param.LogLevel, 2,
		  (*m_osLog)
		  << "Node " << nodeIndex << " process stopping on bound."
		  << " This LB= " 
		  << UtilDblToStr(m_nodeStats.objBest.first)
		  << " Global UB= " 
		  << UtilDblToStr(m_nodeStats.objBest.second) << "." << endl;);
	 m_stopCriteria = DecompStopBound;
	 m_phase        = PHASE_DONE;	 	 
      }
      if(m_phase == PHASE_DONE)
	 break;

      bool          isGapTight = false;      
      DecompVarList newVars;
      DecompCutList newCuts;
      switch(m_phase){
      case PHASE_PRICE1:
      case PHASE_PRICE2:
	 m_nodeStats.priceCallsRound++;
	 m_nodeStats.priceCallsTotal++;
	 
	 //---
	 //--- after adding some rows, the columns in the var pool
	 //--- might no longer be valid, so we need to re-expand everything
	 //---
	 if(m_varpool.size() > 0){
	    if(!m_varpool.colsAreValid()){
	       UTIL_MSG(m_param.LogDebugLevel, 3,
			(*m_osLog) << "EXPANDING varpool.\n";);
	       m_varpool.reExpand(*modelCore, m_param.TolZero);
	    }
	 
	    //---
	    //--- THINK....
	    //---
	    if(m_status == STAT_FEASIBLE)
	       m_varpool.setReducedCosts(getRowPrice(), m_status);
	    else{
	       //if doing RC, never called??
	       const double * u = getDualRays(1)[0];
	       m_varpool.setReducedCosts(u, m_status);
	       UTIL_DELARR(u);
	    }
	 }
	 
	 //---
	 //--- attempt to generate some new variables with rc < 0
	 //---
	 mostNegRC                  = 0.0;
	 m_nodeStats.varsThisCall   = generateVars(m_status,
						   newVars, mostNegRC);
	 m_nodeStats.varsThisRound += m_nodeStats.varsThisCall;
	 m_nodeStats.cutsThisCall   = 0;
	 
	 if(m_nodeStats.varsThisCall > 0){
	    //--- 
	    //--- add the newly generated variables to the var pool
	    //---
	    addVarsToPool(newVars);	    
	    //---
	    //--- add variables from the variable pool to the master problem
	    //---
	    addVarsFromPool();            
	 }
	 if(m_isColGenExact           && 
	    m_rrIterSinceAll == 0     &&
	    m_status == STAT_FEASIBLE &&
	    m_phase  == PHASE_PRICE2)
	    isGapTight = updateObjBoundLB(mostNegRC);


	 //---
	 //--- update stab parameters delta=duals, epsilon reduced
	 //---
         /*#ifdef STAB_DUMERLE
           m_stabEpsilon *= 0.95;
	 dualSol = m_masterSI->getRowPrice();	 
	 int i, r;
	 for(i = 0; i < m_masterSI->getNumCols(); i++){
         DecompColType type = m_masterColType[i];
         if(isMasterColArtificial(i)){
         if(type == DecompCol_ArtForRowL    ||
         type == DecompCol_ArtForBranchL ||
         type == DecompCol_ArtForCutL){
         r = m_artColIndToRowInd[i];
         printf("Master Col i=%d type=%s r=%d dual=%g\n",
         i, DecompColTypeStr[type].c_str(), r, dualSol[r]);
         m_masterSI->setObjCoeff(i, -dualSol[r]);
         }
         else if(type == DecompCol_ArtForRowG    ||
         type == DecompCol_ArtForBranchG ||
         type == DecompCol_ArtForCutG){
         r = m_artColIndToRowInd[i];
         printf("Master Col i=%d type=%s r=%d dual=%g\n",
         i, DecompColTypeStr[type].c_str(), r, dualSol[r]);
         m_masterSI->setObjCoeff(i, dualSol[r]);         
         }
         //CAN'T DO THIS IF IN PHASEI!
         //m_masterSI->setColBounds(i, 0.0, 0.0);//TODO
         m_masterSI->setColBounds(i, 0.0, m_stabEpsilon);//TODO
         }
	 }
         #endif*/


	 break;
      case PHASE_CUT:
	 m_nodeStats.cutCallsRound++;
	 m_nodeStats.cutCallsTotal++;

	 //---
	 //--- after adding some cols, the rows in the cut pool
	 //--- might no longer be valid, so we need to re-expand everything
	 //---
	 if(!m_cutpool.rowsAreValid() && (m_cutpool.size()>0)){
	    UTIL_MSG(m_param.LogDebugLevel, 3,
		     (*m_osLog) << "EXPANDING cutpool.\n";);
	    m_cutpool.reExpand(m_vars, 
			       modelCore->getNumCols(),
			       m_nArtCols);
	 }
	 
	 //THINK: here is where you will do sep of m_xhat vs shat
	 m_cutpool.calcViolations(m_xhat);
         
	 //---
	 //--- attempt to generate some new cuts with vio > 0
	 //---
	 m_nodeStats.cutsThisCall   = generateCuts(m_xhat, newCuts);
	 m_nodeStats.cutsThisRound += m_nodeStats.cutsThisCall;
	 m_nodeStats.varsThisCall   = 0;
	 
	 if(m_nodeStats.cutsThisCall > 0){
	    //this updates the lb based on last solve, not this solve!
	    // gen cut doesn't change bound until we resolve	    
	    //if(m_algo == CUT)
	    //   updateObjBoundLB();
	    
	    //--- 
	    //--- add the newly generated cuts to the cut pool
	    //---            
	    addCutsToPool(m_xhat, newCuts, m_nodeStats.cutsThisCall);
	    
	    //---
	    //--- add cuts from the cut pool to the master problem
	    //---
	    addCutsFromPool();        
	 }
	 break;
      case PHASE_DONE:
	 break;
      default:
	 assert(0);
      }

      //---
      //--- Careful here -- for many apps the user will use heuristics
      //--- during column generation phase. If we update LB after each 
      //--- column added we might stop too early if this LB exceeds the
      //--- tree's global upper bound.
      //---
      //--- We need the user to tell us if they solved it exactly or not.      
      //---
      //this should be in phaseUpdate?

      
      //TODO: moved this into phaseUpdate for PC - need to revisit CPM!
      //TODO: now moved to phaseUpdate - what about case of no branch object!?
      if(m_phase != PHASE_DONE){
	 //---
	 //--- perform a solution update
	 //---  PC: take PARM steps of simplex
	 //---  ?? DC: take PARM steps of simplex (INF case?)
	 //---  RC: take PARM steps of subgradient
	 //---  VC: take PARM steps of volume
	 //---
	 if(m_param.LogDumpModel > 1)
	    printCurrentProblem(m_masterSI,
				"masterProb",
				m_nodeStats.nodeIndex,
				m_nodeStats.cutCallsTotal,
				m_nodeStats.priceCallsTotal);


	 //---
	 //--- check to see if something got added
	 //---
	 int nChanges = m_nodeStats.cutsThisCall + m_nodeStats.varsThisCall;
	 UTIL_DEBUG(m_param.LogDebugLevel, 3,
		    (*m_osLog) << "nNewVars = " << m_nodeStats.varsThisCall
		    << endl;
		    (*m_osLog) << "nNewCuts = " << m_nodeStats.cutsThisCall
		    << endl;
		    );
		    
	 if(!isDone()){
	    
	    if(nChanges){
	       //why is this not in the switch statement?
	       m_status = solutionUpdate(m_phase);
	       //make this some update that can be override for CPM vs PC
	       // or move this update to phaseUpdate??? 
	       if(m_algo == CUT && m_nodeStats.cutsThisCall > 0){
		  updateObjBoundLB();
	       }
	    }
	 }


	 ////////////////THINK
	 //???? shouldn't we have recompose and look for IP feasible
	 //  inside of solutionUpdate - as it should be checked in every
	 //  case


	 //what happens often in pricing - you find integer point,
	 //but lb can be improved so you price, but still get integer point
	 //as min, so you keep adding the same ip point - need cmp to
	 //check for dups - esp against previous one

	 //---
	 //--- check if IP feasible (are we done?) 
	 //--- TODO: for nonexplicity, also check user app isfeasible
	 //---
	 //TODO: should this whole section be phaseDone?
	 if(nChanges && m_status != STAT_INFEASIBLE){
	    recomposeSolution(getMasterColSolution(), m_xhat);
	    UTIL_DEBUG(m_param.LogDebugLevel, 4,
		       m_app->printOriginalSolution(modelCore->getNumCols(),
                                                    modelCore->getColNames(),
                                                    m_xhat);
		       );
	    //TODO: solution pool?
	    //this is checked again in phase update...
	    //first, check to see if LP solution is already ip and user feas
	    if(isIPFeasible(m_xhat)){
	       //printf("m_xhat is IP FEASIBLE\n");
	       if(m_app->APPisUserFeasible(m_xhat,
					   modelCore->getNumCols(),
					   m_param.TolZero)){
		  //printf("m_xhat is APP FEASIBLE, m_xhatIPFeas size = %d\n",
                  // (int)m_xhatIPFeas.size());
		  
		  //check for dup sol                  
		  bool isDup = m_xhatIPFeas.size() > 0 ? true : false;
		  vector<DecompSolution*>::iterator vit;
		  for(vit  = m_xhatIPFeas.begin();
		      vit != m_xhatIPFeas.end(); vit++){
		     const DecompSolution * xhatIPFeas = *vit;
		     const double         * values
			= xhatIPFeas->getValues();
		     for(int c = 0; c < modelCore->getNumCols(); c++){
			if(!UtilIsZero(values[c] - m_xhat[c])){
			   isDup = false;
			   break;
			}                           
		     }
		  }

		  if(isDup){
		     printf("IS DUP, not pushing\n");
		  }
		  else{
		     DecompSolution * decompSol
			= new DecompSolution(modelCore->getNumCols(),
					     m_xhat,
					     getOrigObjective());
		     //getMasterObjValue());
		     //solution pool?
		     m_xhatIPFeas.push_back(decompSol);
		     //printf("m_xhatIPFeas size = %d\n", 
                     //(int)m_xhatIPFeas.size());
		  }
	       }
	    }
        
	    //heuristics(m_xhat, m_xhatIPFeas);
	    m_app->APPheuristics(m_xhat, getOrigObjective(), m_xhatIPFeas);


	    //TODO: make this a function!
	    
	    vector<DecompSolution*>::iterator vi;
	    DecompSolution * viBest = NULL;
	    double bestBoundUB = m_nodeStats.objBest.second;
	    for(vi = m_xhatIPFeas.begin(); vi != m_xhatIPFeas.end(); vi++){
	       const DecompSolution * xhatIPFeas = *vi; 
	       //if i am in this list, i must already be feasible... 
	       //if(isIPFeasible(xhatIPFeas->getValues())){
	       //printf("bestBound = %g, xhatIPFeas qual = %g\n",
	       //     bestBoundUB, xhatIPFeas->getQuality());
	       if(xhatIPFeas->getQuality() <= bestBoundUB){
		  bestBoundUB = xhatIPFeas->getQuality();
		  viBest = *vi;
	       }
	       //}         
	    }
	    if(viBest){
	       //save the best
	       setObjBoundUB(bestBoundUB);
	       m_xhatIPBest = viBest;
	    }
	 }
     	 
	 if(nChanges && m_phase != PHASE_PRICE1){
	    //---
	    //--- check on tailoff
	    //---	 
	    if(isTailoffLB(m_param.TailoffLength,
			   m_param.TailoffPercent)){
	       printf("TAILING OFF -- STOP processNode\n");
	       m_stopCriteria = DecompStopTailOff;
	       m_phaseLast    = m_phase;
	       m_phase        = PHASE_DONE;
	    }
	    
            
	    //---
	    //--- did the master objective change?
	    //---   if not, make sure the columns just added cannot be 
	    //---   deleted
	    //---
	    int i;
            UTIL_DEBUG(m_param.LogDebugLevel, 3,
                       (*m_osLog) 
                       << "m_masterObjLast = " << setw(10) 
                       << UtilDblToStr(m_masterObjLast)
                       << " thisMaster = " << setw(10)
                       << UtilDblToStr(getMasterObjValue())
                       << " varsThisCall = " << setw(5) 
                       << m_nodeStats.varsThisCall;);
	    //what if it just changed due to cuts?
	    if(UtilIsZero(m_masterObjLast - getMasterObjValue(), 1.0e-4)){
	       m_objNoChange = true;
               
               UTIL_DEBUG(m_param.LogDebugLevel, 3,
                          (*m_osLog) << "No objective change" << endl;);           
	       //0 1 2 3 4
	       // w new vars
	       //4x 3x 2 1 0
	       // "DecompCol_Structural_NoDelete",
	       if(m_nodeStats.varsThisCall > 0){
		  int sz = static_cast<int>(m_masterColType.size());
                  for(i  = sz - 1;
		      i >= sz - m_nodeStats.varsThisCall; i--){
                     UTIL_DEBUG(m_param.LogDebugLevel, 3,
                                (*m_osLog) << "Col " << i << " has type "
                                << DecompColTypeStr[m_masterColType[i]] 
                                << endl;);
		     assert(m_masterColType[i] == DecompCol_Structural);
		     m_masterColType[i] = DecompCol_Structural_NoDelete;
		  }
	       }
	    }
	    else{
	       m_objNoChange = false;
               UTIL_DEBUG(m_param.LogDebugLevel, 3, (*m_osLog) << endl;);
	    }   
	    m_masterObjLast = getMasterObjValue();
            
            
	    if(m_phase != PHASE_DONE && m_param.CompressColumns){
	       //---
	       //--- adjust columns effectiveness count
	       //---	    
	       adjustColumnsEffCnt();
	       
	       //---
	       //--- periodically, get rid of ineffective columns
	       //--- periodic:
	       //---    every K iterations OR
	       //---    numCols has doubled since last compression
	       //---
	       compressColumns();
	    }
	 }
      }
   } //while(phase != PHASE_DONE)

   phaseDone();
   
   //need to check again, if we get ip feasible in first LP
   //but this will cause dups... if we also find above?
   if(m_xhatIPFeas.size() == 0 && m_status != STAT_INFEASIBLE){
      //this is checked again in phase update...
      //first, check to see if LP solution is already ip and user feas
      if(isIPFeasible(m_xhat)){               
	 if(m_app->APPisUserFeasible(m_xhat,
				     modelCore->getNumCols(),
				     m_param.TolZero)){
	    DecompSolution * decompSol
	       = new DecompSolution(modelCore->getNumCols(),
				    m_xhat,
				    getOrigObjective());
	    m_xhatIPFeas.push_back(decompSol);
	    m_xhatIPBest = decompSol;
	 }
      }
   }

   if(m_xhatIPBest){
      UTIL_DEBUG(m_param.LogLevel, 3,
		 (*m_osLog) << "Best Feasible Solution with Quality = " 
		 << UtilDblToStr(m_xhatIPBest->getQuality(), 6) << "\n";
		 m_app->printOriginalSolution(modelCore->getNumCols(),
					      modelCore->getColNames(),
					      m_xhatIPBest->getValues());
		 );
   }


   UTIL_DEBUG(m_param.LogDebugLevel, 3,
	      (*m_osLog) << "StatOut     : " 
	      << DecompStatusStr[m_status] << "\n";
	      (*m_osLog) << "StopCriteria: " 
	      << DecompAlgoStopStr[m_stopCriteria] << "\n";
	      (*m_osLog) << "RelGap      : " 
	      << UtilDblToStr(m_relGap, 6) << "\n";
	      );

   m_stats.thisDecomp.push_back(m_stats.timerDecomp.getRealTime());


   //if i am root and doing price and cut, solve this IP to get ub...
   //  e.g., cutting stock works well -> better to do at AlpsDecompTreeNode

   UTIL_MSG(m_param.LogDebugLevel, 3,
	    m_stats.printOverallStats(m_osLog);
	    m_nodeStats.printObjHistory(m_osLog);
	    m_nodeStats.printObjHistoryLB(m_osLog);
	    m_nodeStats.printObjHistoryUB(m_osLog);
	    );
   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "processNode()", m_param.LogDebugLevel, 1);

   return m_status;
}

//--------------------------------------------------------------------- //
void DecompAlgo::setMasterBounds(const double * lbs,
                                 const double * ubs){
   
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "setMasterBounds()", m_param.LogDebugLevel, 2);
   
   int       c, coreColIndex;
   DecompConstraintSet * modelCore = m_modelCore.getModel();   
   const int nIntVars = modelCore->getNumInts();
   const int nCols    = modelCore->getNumCols();
   const int beg      = modelCore->nBaseRowsOrig;

   //---
   //--- make copy so we can enforce in subproblems
   //---
   memcpy(m_colLBNode, lbs, nCols * sizeof(double));
   memcpy(m_colUBNode, ubs, nCols * sizeof(double));
   
   //TODO: can reuse this memory
   int      nRows   = 2 * nIntVars;
   int    * index   = new int[nRows];
   char   * sense   = new char[nRows];
   double * rhs     = new double[nRows];
   double * range   = new double[nRows];
   const int * integerVars = modelCore->getIntegerVars();


   //lbs,ubs is indexed on core column index
   // but c is being looped over integers here... 

   //---
   //--- the row index for column c's UB (x <= u) is: beg            + c
   //--- the row index for column c's LB (x >= l) is: beg + nIntVars + c
   //---
   for(c = 0; c < nIntVars; c++){
      //x <= u
      coreColIndex = integerVars[c];
      index[c]     = beg + c; //row index into master
      sense[c]     = 'L';
      rhs  [c]     = ubs[coreColIndex];
      range[c]     = 0.0;
      if(m_masterRowType[beg+c] != DecompRow_Branch){
	 printf("ERROR: row %d type: %s\n",
		beg+c,
		DecompRowTypeStr[m_masterRowType[beg+c]].c_str());
      }
      assert(m_masterRowType[beg+c] == DecompRow_Branch);
   }
   for(c = nIntVars; c < (2*nIntVars); c++){
      //x >= l
      coreColIndex = integerVars[c - nIntVars];
      index[c]     = beg + c; 
      sense[c]     = 'G';
      rhs  [c]     = lbs[coreColIndex];
      range[c]     = 0.0;
      if(m_masterRowType[beg+c] != DecompRow_Branch){
	 printf("ERROR: row %d type: %s\n",
		beg+c,
		DecompRowTypeStr[m_masterRowType[beg+c]].c_str());
      }
      assert(m_masterRowType[beg+c] == DecompRow_Branch);
   }
  
   m_masterSI->setRowSetTypes(index, index+(2*nIntVars), sense, rhs, range);
   UTIL_DELARR(index);
   UTIL_DELARR(sense);
   UTIL_DELARR(rhs);
   UTIL_DELARR(range);

   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "setMasterBounds()", m_param.LogDebugLevel, 2);   
}

//===========================================================================//
DecompStatus DecompAlgo::solutionUpdate(const DecompPhase phase,
                                        //TODO: not current luse?
                                        const int         maxInnerIter,
                                        const int         maxOuterIter){
   
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "solutionUpdate()", m_param.LogDebugLevel, 2);
   
   m_stats.timerOther1.reset();
   
   int i, cpxStat=0, cpxMethod=0;
   DecompStatus status = STAT_UNKNOWN;

   //---
   //--- was missing all along? 9/28/09
   //---
   //#ifdef __DECOMP_LP_CLP__
   //m_masterSI->setHintParam(OsiDoPresolveInResolve, false, OsiHintDo);
   //#else
   m_masterSI->setHintParam(OsiDoPresolveInResolve, true, OsiHintDo);
   //#endif
   
   //m_masterSI->setIntParam(OsiMaxNumIteration, maxInnerIter);
  
   
   //THINK:
   //if we allow for interior, need crossover too?

#ifdef __DECOMP_LP_CPX__
   OsiCpxSolverInterface * masterCpxSI 
      = dynamic_cast<OsiCpxSolverInterface*>(m_masterSI);
   CPXENVptr env = masterCpxSI->getEnvironmentPtr();
#ifdef DO_INTERIOR
   CPXLPptr  lp  = masterCpxSI->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL);
#endif
   CPXsetintparam( env, CPX_PARAM_PREIND, CPX_ON );
   CPXsetintparam( env, CPX_PARAM_SCRIND, CPX_ON );
   CPXsetintparam( env, CPX_PARAM_SIMDISPLAY, 2 );


   int preInd = 0;
   CPXgetintparam(env, CPX_PARAM_PREIND, &preInd);
   printf("preind=%d\n",preInd);
#endif

  
   switch(phase){
   case PHASE_PRICE1:
   case PHASE_PRICE2:
      m_masterSI->setDblParam(OsiDualObjectiveLimit, DecompInf);

      //m_masterSI->setHintParam(OsiDoDualInResolve, false, OsiHintDo);
      m_masterSI->setHintParam(OsiDoDualInResolve, true, OsiHintDo);


      //if(m_algo == DECOMP)//THINK!
      // m_masterSI->setHintParam(OsiDoPresolveInResolve, false, OsiHintDo);
#ifdef DO_INTERIOR //TODO: option, not compile time
      //CPXhybbaropt(env, lp, 0);//if crossover, defeat purpose
      CPXbaropt(env, lp);
      cpxMethod = CPXgetmethod(env, lp);
      cpxStat = CPXgetstat(env, lp);
      if(cpxStat)
	 printf("cpxMethod=%d, cpxStat = %d\n", cpxMethod, cpxStat);
#else
      m_masterSI->resolve();
#endif
      break;
   case PHASE_CUT:
      m_masterSI->setHintParam(OsiDoDualInResolve, true, OsiHintDo);
      m_masterSI->resolve();
      break;
   default:
      assert(0);
   }
   
   UTIL_MSG(m_param.LogDebugLevel, 3,             
	    (*m_osLog) << "Solution update n_cols:"
	    << setw(10) << m_masterSI->getNumCols() << " n_rows: " 
	    << setw(10) << m_masterSI->getNumRows() << " n_iter: "
	    << setw(10) << m_masterSI->getIterationCount() << " time: "
	    << setw(10) << m_stats.timerOther1.getRealTime()
	    << endl;
	    );

#ifdef __DECOMP_LP_CLP__
   UTIL_DEBUG(m_param.LogDebugLevel, 4,
              {
                 OsiClpSolverInterface * osiClp  = dynamic_cast<OsiClpSolverInterface*>(m_masterSI);
                 printf("clp status        = %d\n", osiClp->getModelPtr()->status());
                 printf("clp prob status   = %d\n", osiClp->getModelPtr()->problemStatus());
                 printf("clp second status = %d\n", osiClp->getModelPtr()->secondaryStatus());
              }
              );
#endif

   
   UTIL_DEBUG(m_param.LogDebugLevel, 4,
	      (*m_osLog)
	      << "Iteration Count               : "
	      << m_masterSI->getIterationCount() << "\n"
	      << "isAbandoned()                 : " 
	      << m_masterSI->isAbandoned() << "\n"
	      << "isProvenOptimal()             : " 
	      << m_masterSI->isProvenOptimal() << "\n"
	      << "isProvenPrimalInfeasible()    : " 
	      << m_masterSI->isProvenPrimalInfeasible() << "\n"
	      << "isProvenDualInfeasible()      : " 
	      << m_masterSI->isProvenDualInfeasible() << "\n"
	      << "isPrimalObjectiveLimitReached : " 
	      << m_masterSI->isDualObjectiveLimitReached() << "\n"
	      << "isDualObjectiveLimitReached   : " 
	      << m_masterSI->isDualObjectiveLimitReached() << "\n"
	      << "isIterationLimitReached       : " 
	      << m_masterSI->isIterationLimitReached() << "\n";
	      );
                
   if(m_masterSI->isProvenOptimal()){
      status = STAT_FEASIBLE;  

      //if we are using cpx, we need to save the 
      //solution and we cannot use getColSolution() later on
      //for example, after addCols is called, cache is lost
      const int      nCols  = m_masterSI->getNumCols();
      const double * colSol = m_masterSI->getColSolution();
      colSolution.clear();
      colSolution.reserve(nCols);
      for(i = 0; i < nCols; i++)
	 colSolution.push_back(colSol[i]);
      assert((int)colSolution.size() == nCols);
      
      UTIL_DEBUG(m_param.LogDebugLevel, 4,
		 (*m_osLog)
		 << "MasterObj                     : "
		 << UtilDblToStr(getMasterObjValue()) << "\n";
		 );
      
      //sanity check
      if(m_algo != CUT)
	 checkMasterDualObj(); 


      //---
      //--- HACK: there is some bug in CLP where infeasible is declared optimal
      //---   but then we get back solution at state when it internally gave up
      //---   
      //--- Check to see if some lambda < 0 - i.e., junk. If so, assume that
      //---  it meant to return infeasible.
      //---
      for(i = 0; i < nCols; i++){
         if(colSol[i] < -1){
            (*m_osLog) << "ERROR: NEGATIVE LAMBDA, but Osi returns as optimal"
                       << " assume it was meant to be infeasible." << endl;
            status = STAT_INFEASIBLE;
         }
      }
   }
   else if(m_masterSI->isProvenPrimalInfeasible() ||
	   m_masterSI->isProvenDualInfeasible()){
      //for interior, if infeasible, the status is not
      //  getting picked up properly by OSI


      status = STAT_INFEASIBLE;
      //---
      //--- it is possible that presolver determined infeasibility
      //--- but, we will need a dual ray, so we should resolve with
      //--- presolve off
      //---
      m_masterSI->setDblParam(OsiDualObjectiveLimit, DecompInf);
      m_masterSI->setHintParam(OsiDoPresolveInResolve, false, OsiHintDo);
      m_masterSI->resolve();
      m_masterSI->setHintParam(OsiDoPresolveInResolve, true,  OsiHintDo);
   }
   else
      {
#ifdef DO_INTERIOR
	 if(m_masterSI->isDualObjectiveLimitReached()){
	    status = STAT_INFEASIBLE;
	 }
	 else
#endif
	    {
	       assert(0);
	    }
      }


   //---
   //--- HACK: there is some bug in CLP where infeasible is declared optimal
   //---   but then we get back solution at state when it internally gave up
   //---   
   //--- Check to see if some lambda < 0 - i.e., junk. If so, assume that
   //---  it meant to return infeasible.
   //---
   

   m_stats.thisSolUpdate.push_back(m_stats.timerOther1.getRealTime()); 

   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "solutionUpdate()", m_param.LogDebugLevel, 2);   
   return status;
}

//===========================================================================//
int DecompAlgo::generateInitVars(DecompVarList & initVars){

   int          c, attempts;
   double       aveC;

   DecompConstraintSet * modelCore = m_modelCore.getModel();   
   const int      limit      = m_param.LimitInitVars;
   const int      limit2     = 2 * limit;
   const int      nCoreCols  = modelCore->getNumCols();
   const double * objCoeff   = getOrigObjective();

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "generateInitVars()", m_param.LogDebugLevel, 2);

   //---
   //--- APP: create an initial set of points F'[0] subseteq F'
   //---  The base implementation of this function does nothing.
   //---  This is the user's chance to implement something application 
   //---  specific.
   //---
   m_app->generateInitVars(initVars);
   
   int nInitVars = static_cast<int>(initVars.size());
   UTIL_DEBUG(m_param.LogDebugLevel, 4,
	      (*m_osLog) 
	      << "nInitVars from app = " << nInitVars 
	      << " userLimit = " << limit << endl;
	      );
   
   if(nInitVars < limit){      
      //---
      //--- create an initial set of points F'[0] subseteq F'
      //---   randomly by solving zSP(c + eps), eps = U[0,ave(c)]
      //---

      //---
      //--- NOTE: in GAP case, subproblem is knapsack, if use orig cost
      //---   all cost > 0, so will get NULL column, later on reduced costs
      //---   will give negative values, so this is not a problem
      //---

      double * costeps = new double[nCoreCols];
      assert(objCoeff);
      aveC = UtilAve(objCoeff, nCoreCols);
      
      attempts = 0;
      DecompSolverResult subprobResult(nCoreCols);
      while((nInitVars < limit) && (attempts < limit2)){         
	 //---
	 //--- perturb the cost vector
	 //---
	 srand(attempts);
	 for(c = 0; c < nCoreCols; c++){
	    double r = 0.0;
	    if(attempts != 0)
	       r = UtilURand(-aveC, aveC);
	    costeps[c] = objCoeff[c] + r;
	 }
	 	 
	 //---
	 //--- APP: solve zSP(c + eps) 
	 //---
         map<int, DecompAlgoModel>::iterator mit; 
         for(mit = m_modelRelax.begin(); mit != m_modelRelax.end(); mit++){
            DecompAlgoModel & algoModel = (*mit).second;
            solveRelaxed(costeps,               //reduced cost (fake here)
                         objCoeff,              //original cost vector 
                         9e15,                  //alpha        (fake here)
                         nCoreCols,             //num core columns
                         false,                 //isNested
                         algoModel,
			 &subprobResult,        //results
			 initVars);             //var list to populate
	 }

         
         map<int, vector<DecompAlgoModel> >::iterator mivt;
         vector<DecompAlgoModel>           ::iterator vit; 
         for(mivt  = m_modelRelaxNest.begin(); 
             mivt != m_modelRelaxNest.end(); mivt++){
            for(vit  = (*mivt).second.begin(); 
                vit != (*mivt).second.end(); vit++){
               solveRelaxed(
                            costeps,               //reduced cost (fake here)
                            objCoeff,              //original cost vector 
                            9e15,                  //alpha        (fake here)
                            nCoreCols,             //num core columns
                            true,                  //isNested
                            (*vit),
                            &subprobResult,        //results
                            initVars);             //var list to populate
            }
         }
         
         //---
	 //--- THINK: check for duplicate variables - done in solveRelaxed
	 //---   don't assume the user does the duplicate check - should be 
         //---   done by col pool also
         //---
         nInitVars = static_cast<int>(initVars.size());
	 attempts++;
      }

      UTIL_DEBUG(m_param.LogDebugLevel, 4,
		 (*m_osLog) << "\nm_varsThisCall = "
		 << initVars.size() << "\n";
		 );
      
      
      //--- 
      //--- TODO: solve a few iterations of subgradient to get init vars?
      //---    
      //--- TODO: put them in the var pool??
      //---
      UTIL_DELARR(costeps);//TODO: use mem-pool
   }



   //---
   //--- generate init vars by solving root LP and 
   //---   running DC at each iteration
   //---
   if(m_param.InitVarsWithCutDC){
      printf("======= BEGIN Gen Init Vars - call CPM process root node\n");
      DecompAlgoC cpm(m_app, m_utilParam);
      cpm.m_param.CutDC = 2;
      cpm.processNode(0, -DecompInf, DecompInf);

      //---
      //--- copy the vars generated in passes of DC into initVars
      //---   to warm-start DW master formulation
      //---
      //m_vars.insert(m_vars.end(), cpm.m_vars.begin(), cpm.m_vars.end());
      initVars.splice(initVars.end(), cpm.m_vars);
      printf("VARS moved into PC object initVars.size=%d\n",
	     static_cast<int>(initVars.size()));
      //printVars(m_osLog);//use this to warm start DW



      //a hidden advantage of decomp in DC?
      DecompSolution * bestSol = NULL;
      vector<DecompSolution*>::iterator it;
      
      //there will be just one, i think, just need to copy it over here
      double thisBound;
      double bestBoundUB = m_nodeStats.objBest.second;
      for(it  = cpm.m_xhatIPFeas.begin();
	  it != cpm.m_xhatIPFeas.end(); it++){
         thisBound = (*it)->getQuality();
         printf("From init vars, IP Feasible with Quality = %g\n", thisBound);
	 if((*it)->getQuality() <= bestBoundUB){
            bestBoundUB = (*it)->getQuality();
            bestSol     = (*it);
	 }
      }
      //need to make copy of solution, since D.m_xhatIpFeas goes out of scope
      if(bestSol){
	 DecompSolution * bestSolCp = new DecompSolution(*bestSol);
	 m_xhatIPFeas.push_back(bestSolCp);	 
	 setObjBoundUB(bestSolCp->getQuality());
	 m_xhatIPBest = bestSolCp;
	 m_xhatIPBest->print();
      }
      


      printf("======= END   Gen Init Vars - call CPM process root node\n");
   }


   if(m_param.InitVarsWithIP){
      printf("======= BEGIN Gen Init Vars - call Direct IP solver\n");
      DecompAlgoC          direct(m_app, m_utilParam);
      DecompSolverResult * result = new DecompSolverResult(nCoreCols);
      assert(result);

      direct.solveDirect(m_param.InitVarsWithIPLimitTime, result);
      if(result->m_solution){
	 //---
	 //--- if an incumbent was found, create a var(s) from it
	 //---
	 if(m_numConvexCon == 1){
	    DecompVar * directVar = new DecompVar(nCoreCols,
						  result->m_solution,
						  0.0, result->m_objUB);
	    initVars.push_back(directVar);
	 }
	 else{
	    map<int, DecompAlgoModel>::iterator mid;
	    for(mid = m_modelRelax.begin(); mid != m_modelRelax.end(); mid++){
	       int               blockId       = (*mid).first;
	       DecompAlgoModel & modelRelax    = (*mid).second;
	       vector<int>     & activeColumns 
		  = modelRelax.getModel()->activeColumns;
	       vector<int>       ind;
	       vector<double>    els;
	       double            origCost = 0.0;    
	       vector<int>::iterator it;
	       for(it  = activeColumns.begin(); 
		   it != activeColumns.end(); it++){
		  if(!UtilIsZero(result->m_solution[*it])){
		     ind.push_back(*it);
		     els.push_back(result->m_solution[*it]);
		     origCost += objCoeff[*it] * result->m_solution[*it];
		  }
	       }
	       DecompVar * directVar 
		  = new DecompVar(ind, els, 0.0, origCost);
	       directVar->setBlockId(blockId);
	       initVars.push_back(directVar);
	    }
	 }
	 
	 //---
	 //--- update the upper bound
	 double bestBoundUB = m_nodeStats.objBest.second;
	 if(result->m_objUB < bestBoundUB){
	    DecompSolution * directSol = new DecompSolution(nCoreCols,
							    result->m_solution,
							    result->m_objUB);
	    m_xhatIPFeas.push_back(directSol);
	    m_xhatIPBest = directSol;
	    setObjBoundUB(result->m_objUB);	    
	 }
      }
      
      printf("======= END   Gen Init Vars - call Direct IP solver\n");      
   }


   //---
   //--- check init vars for incumbent
   //---
   DecompVarList::iterator vli;
   for(vli = initVars.begin(); vli != initVars.end(); vli++){
      //---
      //--- unlikey to happen - but we should check ALL columns
      //---  to see if they are IP feasible 
      //---
      (*vli)->fillDenseArr(modelCore->getNumCols(),
                           m_memPool.dblArrNCoreCols);            
      
      if(isIPFeasible(m_memPool.dblArrNCoreCols)){               
	 if(m_app->APPisUserFeasible(m_memPool.dblArrNCoreCols,
				     modelCore->getNumCols(),
				     m_param.TolZero)){
	    DecompSolution * decompSol
	       = new DecompSolution(modelCore->getNumCols(),
				    m_memPool.dblArrNCoreCols,
				    (*vli)->getOriginalCost());
            m_xhatIPBest = decompSol;
	    m_xhatIPFeas.push_back(decompSol);
            //printf("var is ip feas with obj = %g\n",
            //     (*vli)->getOriginalCost());
	    setObjBoundUB((*vli)->getOriginalCost());
	 }
      }
   }

   //---
   //--- this will update the global UB before we start processing
   //---   if we were lucky enough to find an incumbent in init vars
   //---
   //setCutoffUB(getCutoffUB());

   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "generateInitVars()", m_param.LogDebugLevel, 2);
   nInitVars = static_cast<int>(initVars.size());
   return nInitVars;
}

//===========================================================================//
//once we do RC, this probably won't be in base anyway
bool DecompAlgo::updateObjBoundLB(const double mostNegRC){
   //---
   //--- C    : LB = masterLP obj
   //--- PC   : LB = zDW_RMP + RC* <= zDW <= zDW_RMP
   //---    where RC* is the most negative reduced cost
   //---    assuming the relaxation subproblem was solved exactly
   //---
   //--- Careful here -- for many apps the user will use heuristics
   //--- during column generation phase. If we update LB after each 
   //--- column added we might stop too early if this LB exceeds the
   //--- tree's global upper bound.
   //---
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "updateObjBoundLB()", m_param.LogDebugLevel, 2);

#ifdef STAB_DUMERLE
   const double * dualSol = m_masterSI->getRowPrice();
   const double * rowRhs  = m_masterSI->getRightHandSide();
   double         zDW_UB  = 0.0;
   int r;
   for(r = 0; r < m_masterSI->getNumRows(); r++){
      zDW_UB += dualSol[r] * rowRhs[r];
   }
   printf("MasterObj = %g SumPiB = %g\n",
	  getMasterObjValue(), zDW_UB);
#else
   double zDW_UB = getMasterObjValue();
#endif


   double zDW_LB = zDW_UB + mostNegRC;
   setObjBoundLB(zDW_LB);

   //---
   //--- THINK: if do stabilization
   //---  zDW_RMP = b*u ?
   //--- NO - because zDW_RMP = cx + delta*slacks

   //TODO: stats - we want to play zDW_LB vs UB... 

   UTIL_MSG(m_param.LogDebugLevel, 3,
	    (*m_osLog)
	    << "MasterObj = "   << UtilDblToStr(zDW_UB)
	    << "\tmostNegRC = " << UtilDblToStr(mostNegRC) << "\n"
	    << "ThisLB = "      << UtilDblToStr(zDW_LB) << "\t"
	    << "BestLB = "      << UtilDblToStr(m_nodeStats.objBest.first)
	    << "\n";
	    );

   UTIL_DEBUG(m_param.LogDebugLevel, 1,
	      (*m_osLog)
	      << "PriceCallsRound= " 
	      << setw(3)  << m_nodeStats.priceCallsRound
	      << setw(13) << "\tmostNegRC=" 
	      << setw(13) << UtilDblToStr(mostNegRC, 4) 
	      << setw(13) << "\tthisLB="    
	      << setw(13) << UtilDblToStr(zDW_LB, 4) 
	      << endl;
	      ); 

   if((getNodeIndex() == 0) &&
      (zDW_LB > (m_app->getBestKnownUB() + DecompEpsilon))){
      (*m_osLog) << "ERROR: in root node, bestKnownUB = " 
		 << m_app->getBestKnownUB() << " thisBoundLB = " 
		 << zDW_LB << endl;
      assert(0);
   }
	      
   //---
   //--- if RC* is relatively small (i.e., zDW_LB-zDW_UB is small)
   //---  then consider it converged - stop and move on to branching
   //--- TODO: make this an option
   //---
   bool   isGapTight = false;
   double tightGap   = m_param.MasterGapLimit;
   double relGap     = 0.00;
   if(UtilIsZero(zDW_UB)){
      relGap = fabs(mostNegRC);
   }
   else{
      relGap = fabs(mostNegRC) / fabs(zDW_UB);
   }
   if(relGap <= tightGap)
      isGapTight = true;
   if(m_param.LogDebugLevel >= 2){
      (*m_osLog) << "DW relGap = " << UtilDblToStr(relGap) 
		 << " isTight = " << isGapTight << "\n";
   }
   m_relGap = relGap;

   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "updateObjBoundLB()", m_param.LogDebugLevel, 2);
   return isGapTight;
}

//===========================================================================//
void DecompAlgo::masterPhaseItoII(){
   
   //---
   //--- switch from Phase I to Phase II
   //---
   UTIL_MSG(m_app->m_param.LogDebugLevel, 4,
	    (*m_osLog) << "Switching from PhaseI to PhaseII\n";
	    );      

   int i;
   int nMasterCols = m_masterSI->getNumCols();
   //---
   //--- set obj for all columns to original cost
   //--- set obj for artificial columns to 0
   //--- fix column bounds for artificial columns to 0
   //---

#ifdef STAB_DUMERLE
   //---
   //--- set cost on slacks to delta, where delta is init'd to 
   //---  initial duals (this is to be done in Phase 2 only)
   //---    
   //--- min deltap sp - deltam sm
   //---
   //--- ax  = b --> ax + sp - sm  = b, sp >= 0 <= epsp, sm >= 0 <= epsm
   //--- ax <= b --> ax      - sm <= b,                  sm >= 0 <= epsm
   //--- ax >= b --> ax + sp      >= b, sp >= 0 <= epsp
   //---
   int            r;
   const double * dualSol = NULL;
   if(m_useInitLpDuals){
      dualSol = m_cutgenSI->getRowPrice();
      m_useInitLpDuals = false;
      m_stabEpsilon    = 0.0;
   }
   else{
      dualSol = m_masterSI->getRowPrice();
   }
   assert(nMasterCols == static_cast<int>(m_masterColType.size()));
   for(i = 0; i < nMasterCols; i++){
      DecompColType type = m_masterColType[i];
      if(type == DecompCol_ArtForRowL    ||
         type == DecompCol_ArtForBranchL ||
         type == DecompCol_ArtForCutL){
         r = m_artColIndToRowInd[i];
	 printf("Master Col i=%d type=%s r=%d dual=%g\n",
		i, DecompColTypeStr[type].c_str(), r, dualSol[r]);
         m_masterSI->setObjCoeff(i, -dualSol[r]);
      }
      else if(type == DecompCol_ArtForRowG    ||
              type == DecompCol_ArtForBranchG ||
              type == DecompCol_ArtForCutG){
         r = m_artColIndToRowInd[i];
	 printf("Master Col i=%d type=%s r=%d dual=%g\n",
		i, DecompColTypeStr[type].c_str(), r, dualSol[r]);
         m_masterSI->setObjCoeff(i, dualSol[r]);         
      }
      else{
         m_masterSI->setObjCoeff(i, 0.0);
      }
      if(isMasterColArtificial(i)){
         //m_masterSI->setColBounds(i, 0.0, 0.0);//TODO
         m_masterSI->setColBounds(i, 0.0, m_stabEpsilon);//TODO
      }
   }   
   DecompVarList::iterator li;            
   for(li = m_vars.begin(); li != m_vars.end(); li++){
      assert(isMasterColStructural((*li)->getColMasterIndex()));
      m_masterSI->setObjCoeff((*li)->getColMasterIndex(),
			      (*li)->getOriginalCost());
   }   
   if(m_param.LogDumpModel > 1)
      printCurrentProblem(m_masterSI,
			  "masterProb_switchItoII",
			  m_nodeStats.nodeIndex,
			  m_nodeStats.cutCallsTotal,
			  m_nodeStats.priceCallsTotal);
#else
   assert(nMasterCols == static_cast<int>(m_masterColType.size()));
   for(i = 0; i < nMasterCols; i++){
      m_masterSI->setObjCoeff(i, 0.0);
      if(isMasterColArtificial(i))
	 m_masterSI->setColBounds(i, 0.0, 0.0);
   }   
   DecompVarList::iterator li;            
   for(li = m_vars.begin(); li != m_vars.end(); li++){
      assert(isMasterColStructural((*li)->getColMasterIndex()));
      m_masterSI->setObjCoeff((*li)->getColMasterIndex(),
			      (*li)->getOriginalCost());
   }   
   if(m_param.LogDumpModel > 1)
      printCurrentProblem(m_masterSI,
			  "masterProb_switchItoII",
			  m_nodeStats.nodeIndex,
			  m_nodeStats.cutCallsTotal,
			  m_nodeStats.priceCallsTotal);
#endif
}

//===========================================================================//
void DecompAlgo::masterPhaseIItoI(){
   
   //---
   //--- switch from Phase II to Phase I
   //---
   UTIL_MSG(m_app->m_param.LogDebugLevel, 4,
	    (*m_osLog) << "Switching from PhaseII to PhaseI\n";
	    );      

   int i;
   int nMasterCols = m_masterSI->getNumCols();
   //---
   //--- set obj for all columns to 0
   //--- set obj for artificial columns to 1
   //--- unfix column bounds for artificial columns
   //---
   assert(nMasterCols == static_cast<int>(m_masterColType.size()));
   for(i = 0; i < nMasterCols; i++){
      if(isMasterColStructural(i))
	 m_masterSI->setObjCoeff(i, 0.0);
      else{
	 m_masterSI->setObjCoeff(i, 1.0);
	 m_masterSI->setColBounds(i, 0.0, DecompInf);
      }
   }   
   if(m_param.LogDumpModel > 1)
      printCurrentProblem(m_masterSI,
			  "masterProb_switchIItoI",
			  m_nodeStats.nodeIndex,
			  m_nodeStats.cutCallsTotal,
			  m_nodeStats.priceCallsTotal);
}

//how does this logic work when dealing with TSP where
// cuts define validity? must continue if user supplies
// cut generator - think?
//===========================================================================//
void DecompAlgo::phaseUpdate(DecompPhase  & phase,
			     DecompStatus & status){
   
   bool              mustSwitch, considerSwitch;
   bool              isCutPossible, isPricePossible, gapTight;

   DecompPhase           nextPhase  = PHASE_UNKNOWN;
   DecompStatus          nextStatus = status;
   pair<double,double> & objBest    = m_nodeStats.objBest;
   int & priceCallsTotal            = m_nodeStats.priceCallsTotal;
   int & cutCallsTotal              = m_nodeStats.cutCallsTotal;
   int & priceCallsRound            = m_nodeStats.priceCallsRound;
   int & cutCallsRound              = m_nodeStats.cutCallsRound;
   int & varsThisCall               = m_nodeStats.varsThisCall;
   int & cutsThisCall               = m_nodeStats.cutsThisCall;
   int & varsThisRound              = m_nodeStats.varsThisRound;
   int & cutsThisRound              = m_nodeStats.cutsThisRound;

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "phaseUpdate()", m_param.LogDebugLevel, 2);
   m_phaseLast = phase;

   UTIL_MSG(m_param.LogDebugLevel, 3,  
	    (*m_osLog) << "cutsThisRound  : " << cutsThisRound   << "\n";
	    (*m_osLog) << "varsThisRound  : " << varsThisRound   << "\n";
	    (*m_osLog) << "cutsThisCall   : " << cutsThisCall    << "\n";
	    (*m_osLog) << "varsThisCall   : " << varsThisCall    << "\n";
	    (*m_osLog) << "cutCallsTotal  : " << cutCallsTotal   << "\n";
	    (*m_osLog) << "priceCallsTotal: " << priceCallsTotal << "\n";
	    (*m_osLog) << "cutCallsRound  : " << cutCallsRound   << "\n";
	    (*m_osLog) << "priceCallsRound: " << priceCallsRound << "\n";
	    (*m_osLog) << "PHASEIN        : "
	    << DecompPhaseStr[phase] << "\n";
	    (*m_osLog) << "STATIN         : "
	    << DecompStatusStr[status] << "\n";
	    (*m_osLog) << "BestLB         : "
	    << UtilDblToStr(objBest.first) << "\n";
	    (*m_osLog) << "BestUB         : "
	    << UtilDblToStr(objBest.second) << "\n";
	    );


   //---
   //--- was the current model found to be infeasible?
   //---
   if(status == STAT_INFEASIBLE){
      //---
      //--- otherwise, switch to PHASEI
      //---   NOTE: this can happen when a new cut (or branch cut) is added
      //---      
      masterPhaseIItoI();
      nextPhase  = PHASE_PRICE1;
      nextStatus = solutionUpdate(nextPhase);
      goto PHASE_UPDATE_FINISH;
   }

   //---
   //--- check to see if we have exceeded the total iter limits
   //---
   isCutPossible   =
      (m_param.LimitRoundCutIters   > 0) &&
      (cutCallsTotal                < m_param.LimitTotalCutIters);
   isPricePossible =
      (m_param.LimitRoundPriceIters > 0) &&
      (priceCallsTotal              < m_param.LimitTotalPriceIters);

   switch(phase){
   case PHASE_PRICE1:
      {
	 //---
	 //--- we are in PHASEI, check to see if solution is feasible
	 //---   to original by checking to see if all artificials are 0
	 //---
	 const double phaseIObj = m_masterSI->getObjValue();
	 UTIL_MSG(m_app->m_param.LogDebugLevel, 3,  
		  (*m_osLog) 
		  << "PhaseIObj= " << UtilDblToStr(phaseIObj) << endl;);
         m_phaseIObj.push_back(phaseIObj);
	 //if(phaseIObj <= DecompZero){
	 if(phaseIObj <= DecompEpsilon){ //11/09/09
	    //---
	    //--- switch to PHASE II (art=0)
	    //---
	    masterPhaseItoII();
	    m_nodeStats.resetCutRound();
	    m_nodeStats.resetPriceRound();
	    if(m_algo == DECOMP){
               nextPhase  = PHASE_DONE;
               nextStatus = STAT_FEASIBLE;
	    }
	    else{
               nextPhase  = PHASE_PRICE2;
               nextStatus = solutionUpdate(nextPhase);
	    }
	    goto PHASE_UPDATE_FINISH;
	 }
	 else{
	    //---
	    //--- we still have active artificials
	    //---    if this is the first call, just continue
	    //---    otherwise, check to see if any new (redCost<0)
	    //---    columns were found to help break infeasibility
	    //---
	    //--- if no vars were found, we are really infeasible
	    //---  else, repeat PhaseI
	    //---
	    if(priceCallsTotal == 0 || varsThisCall > 0){
	       nextPhase = PHASE_PRICE1;
	       goto PHASE_UPDATE_FINISH;	       
	    }
	    else{
	       UTIL_MSG(m_app->m_param.LogDebugLevel, 3,  
			(*m_osLog) 
			<< "Node " << getNodeIndex() 
			<< " is Infeasible." << endl;);
	       m_stopCriteria = DecompStopInfeasible;
	       nextPhase      = PHASE_DONE;
	       nextStatus     = STAT_INFEASIBLE;
	       goto PHASE_UPDATE_FINISH;
	    }
	 }
      }//END: case PHASE_PRICE1
      break;
   case PHASE_PRICE2:
      {
	 assert(status == STAT_FEASIBLE || status == STAT_UNKNOWN);
	 //---
	 //--- if this is the first call, just continue
	 //---
	 if(priceCallsTotal == 0 && cutCallsTotal == 0){
	    nextPhase = PHASE_PRICE2;
	    goto PHASE_UPDATE_FINISH;
	 }

	 //---
	 //--- Princess Bride (1987):
	 //---   "truly, you have a dizzying intellect"
	 //---
	 mustSwitch     = false;
	 considerSwitch = false;    

	 //---
	 //--- if we hit the total limit or
	 //---    we found no new vars this call
	 //--- then we must switch (or we are done)
	 //---
	 if(!isPricePossible || (varsThisCall == 0) || (varsThisRound == 0))
	    mustSwitch = true;

	 //---
	 //--- if we hit the round limit, we must consider switching
	 //---
	 if(priceCallsRound >= m_param.LimitRoundPriceIters)
	    considerSwitch = true;
	 
	 if(mustSwitch){
	    //---
	    //--- we must switch from pricing
	    //---
	    if(!isCutPossible){
	       //---
	       //--- if we exceed the cut iter limit, we are done
	       //---
	       nextPhase = PHASE_DONE;
	    }
	    else{
	       if((cutCallsTotal > 0)  &&
		  (cutsThisRound == 0) &&
		  (varsThisRound == 0)){
		  //---
		  //--- nothing new happened, so we are done
		  //---
		  nextPhase = PHASE_DONE;
	       }
	       else{
		  //---
		  //--- something new happened, so try cuts again
		  //---
		  nextPhase = PHASE_CUT;
		  m_nodeStats.resetCutRound();
	       }
	    }
	 }//END: if(mustSwitch)
	 else if(considerSwitch){
	    //---
	    //--- we consider switching from pricing
	    //---
	    if(!isCutPossible){
	       if(!isPricePossible){
		  //---
		  //--- if we exceed both iter limits, we are done
		  //---
		  nextPhase = PHASE_DONE;
	       }
	       else{
		  //---
		  //--- if we exceed the cut iter limit, but not the price lim
		  //--- since we are not in mustSwitch, m_varsThisRound > 0,
		  //--- so we can go back to pricing, even though it violates
		  //--- the round counter, because we have no other choice
		  //---
		  nextPhase = PHASE_PRICE2;
	       }
	    }
	    else{
	       if((cutsThisRound == 0) && (varsThisRound == 0)){
		  //---
		  //--- nothing new happened, so we are done
		  //---
		  nextPhase = PHASE_DONE;
	       }
	       else{
		  //---
		  //--- something new happened, so try cuts again
		  //---
		  nextPhase = PHASE_CUT;
		  m_nodeStats.resetCutRound();
	       }
	    }
	 } //END: else if(considerSwitch)
	 else
	    nextPhase = PHASE_PRICE2;
      } //END: case PHASE_PRICE2

      //---
      //--- are we suggesting another price phase but gap is tight?
      //---
      //--- ppp, gaptight, ccppp, gaptight c
      gapTight = isGapTight();
      if(gapTight && isCutPossible)
	 if(cutCallsTotal == 0 ||  //haven't even tried cuts yet
	    varsThisRound > 0){    //some new vars, try cut again
	    //---
	    //--- we haven't even tried cuts yet, give it a try
	    //---
	    nextPhase = PHASE_CUT;
	 }
      if(nextPhase == PHASE_PRICE2 && gapTight){
         m_stopCriteria = DecompStopGap;
         
         //---
         //--- if branching candidate does not exist
         //---  do nothing, proceed however you would normally
         //--- else
         //---  force a switch off PRICE2 - go to cuts if possible
         //---
         
         //---
         //--- Even if we are stop on gap, we need to be careful of 
         //--- the following: If the last solution was integral (no 
         //--- branching candidates) but we are not done pricing out 
         //--- (i.e., a column with negative RC still exist) and we 
         //--- declare that we are tailing off then the node will get 
         //--- put back in the node work queue. This can lead to that 
         //--- node being repeatedly stopped and reseted. It is better 
         //--- to just price it out since we cannot branch on it in 
         //--- this state.
         //---
         UTIL_DEBUG(m_param.LogDebugLevel, 3,
                    (*m_osLog) << "Gap is tight" << endl;);
         int    branchedOnIndex = -1;
         double branchedOnValue =  0;
         chooseBranchVar(branchedOnIndex, branchedOnValue);
         if(branchedOnIndex != -1){
            UTIL_DEBUG(m_param.LogDebugLevel, 3,
                       (*m_osLog) << "Gap is tight and we have a " 
                       << "branch candidate" << endl;);
            nextPhase = PHASE_DONE;
            break;
         }
         else{
            UTIL_DEBUG(m_param.LogDebugLevel, 3,
                       (*m_osLog) << "Gap is tight and we have NO " 
                       << "branch candidate" << endl;);
         }
      }
      
      break;
   case PHASE_CUT:
      {
	 //---
	 //--- if this is the first call, just continue
	 //---
	 if(priceCallsTotal == 0 && cutCallsTotal == 0){
	    nextPhase = PHASE_CUT;
	    goto PHASE_UPDATE_FINISH;
	 }

	 //---
	 //--- if tight was gap, the we went to cuts and found none,
	 //---   then stop on gap
	 //---
	 gapTight = isGapTight();
	 if(priceCallsTotal > 0 && cutsThisCall == 0 && gapTight){
	    m_stopCriteria = DecompStopGap;
	    
	    //---
	    //--- Even if we are stop on gap, we need to be careful of 
	    //--- the following: If the last solution was integral (no 
	    //--- branching candidates) but we are not done pricing out 
	    //--- (i.e., a column with negative RC still exist) and we 
	    //--- declare that we are tailing off then the node will get 
	    //--- put back in the node work queue. This can lead to that 
	    //--- node being repeatedly stopped and reseted. It is better 
	    //--- to just price it out since we cannot branch on it in 
	    //--- this state.
	    //---
	    UTIL_DEBUG(m_param.LogDebugLevel, 3,
		       (*m_osLog) << "Gap is tight" << endl;);
	    int    branchedOnIndex = -1;
	    double branchedOnValue =  0;
	    chooseBranchVar(branchedOnIndex, branchedOnValue);
	    if(branchedOnIndex != -1){
	       UTIL_DEBUG(m_param.LogDebugLevel, 3,
			  (*m_osLog) << "Gap is tight and we have a " 
			  << "branch candidate" << endl;);
	       nextPhase = PHASE_DONE;
	       goto PHASE_UPDATE_FINISH;
	    }
	    else{
	       UTIL_DEBUG(m_param.LogDebugLevel, 3,
			  (*m_osLog) << "Gap is tight and we have NO " 
			  << "branch candidate" << endl;);
	    }
	 }
	 
      
	 //---
	 //--- Princess Bride (1987):
	 //---   "truly, you have a dizzying intellect"
	 //---
	 mustSwitch     = false;
	 considerSwitch = false;    	 
	 	
	 //---
	 //--- if we hit the total limit or
	 //---    we found no new cuts this call
	 //--- then we must switch (or we are done)
	 //---
	 if(!isCutPossible || (cutsThisCall == 0) || (cutsThisRound == 0))
	    mustSwitch = true;
	 
	 //---
	 //--- if we hit the round limit, we must consider switching
	 //---
	 if(cutCallsRound >= m_param.LimitRoundCutIters)
	    considerSwitch = true;	 
	 
	 if(mustSwitch){
	    //---
	    //--- we must switch from cutting
	    //---
	    if(!isPricePossible){
	       //---
	       //--- if we exceed the price iter limit, we are done
	       //---
	       nextPhase = PHASE_DONE;
	    }
	    else{
	       if((priceCallsTotal >  0) &&
		  (cutsThisRound   == 0) && 
		  (varsThisRound   == 0)){
		  //---
		  //--- nothing new happened, so we are done
		  //---
		  nextPhase = PHASE_DONE;
	       }
	       else{
		  //---
		  //--- something new happened, so try price again
		  //---
		  nextPhase = PHASE_PRICE2;
		  m_nodeStats.resetPriceRound();
	       }
	    }
	 }//END: if(mustSwitch)
	 else if(considerSwitch){
	    //---
	    //--- we consider switching from cutting
	    //---
	    if(!isPricePossible){
	       if(!isCutPossible){
		  //---
		  //--- if we exceed both iter limits, we are done
		  //---
		  nextPhase = PHASE_DONE;
	       }
	       else{
		  //---
		  //--- if we exceed the price iter limit, but not the cut lim
		  //--- since we are not in mustSwitch, m_cutsThisRound > 0,
		  //--- so we can go back to cutting, even though it violates
		  //--- the round counter, because we have no other choice
		  //---
		  nextPhase = PHASE_CUT;
	       }
	    }
	    else{
	       if((cutsThisRound == 0) && (varsThisRound == 0)){
		  //---
		  //--- nothing new happened, so we are done
		  //---
		  nextPhase = PHASE_DONE;
	       }
	       else{
		  //---
		  //--- something new happened, so try price again
		  //---
		  nextPhase = PHASE_PRICE2;
		  m_nodeStats.resetPriceRound();
	       }
	    }
	 } //END: else if(considerSwitch)
	 else
	    nextPhase = PHASE_CUT;
      } //END: case PHASE_CUT
      break;
   default:
      assert(0);
      //UtilAssert(0, "Bad Phase in phaseUpdate!", m_osLog);      
   } //END: switch(phase)
   
      
 PHASE_UPDATE_FINISH:
   UTIL_MSG(m_param.LogDebugLevel, 3,  
	    (*m_osLog) << "PhaseOut: " << DecompPhaseStr[nextPhase]; 
	    (*m_osLog) << " StatusOut: " << DecompStatusStr[nextStatus]; 
	    (*m_osLog) << endl;
	    );
   phase  = nextPhase;
   status = nextStatus;
   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "phaseUpdate()", m_param.LogDebugLevel, 2);
}

//NOTE: not ok for CPX... do self?

#ifdef __DECOMP_LP_CPX__

#ifdef USE_MULTI_RAY
//------------------------------------------------------------------------- //
vector<double*> DecompAlgo::getDualRays(int maxNumRays){
   OsiCpxSolverInterface * siCpx
      = dynamic_cast<OsiCpxSolverInterface*>(m_masterSI);
   const int m = m_masterSI->getNumRows();
   const int n = m_masterSI->getNumCols();
   const double * rowRhs    = m_masterSI->getRightHandSide();
   const char   * rowSense  = m_masterSI->getRowSense();
   int r, b, c;
   vector<double*> rays;


   //Ax + Is = b
   // ax     <= b
   // ax + s  = b, s >= 0
   // ax     >= b
   // ax + s  = b, s <= 0


   UTIL_DEBUG(m_param.LogDebugLevel, 5,  
	      for(r = 0; r < m; r++){
		 (*m_osLog) << "Row r: " << r << " sense: " << rowSense[r]
			    << " rhs: " << rowRhs[r] << endl;
	      }
	      );

   m_masterSI->enableSimplexInterface(false);

   double * tabRhs   = new double[m];
   int    * basics   = new int[m];
   double * yb       = new double[m];
   double * bInvRow  = new double[m];
   double * bInvARow = new double[n];

   //STOP ============================================
   //tabRhs and yb do NOT match up.... is this an issue?
   //have to hand adjust or use tabRhs since proof is based on B-1
   //which matches up with bhead - what to do in the case of CLP?

   //but, we are multiplying this by A'' later on which is based on 
   //original variable space, not the one adjusted by simplex - so if
   //we return the dual ray directly from B-1 then do B-1A by hand -
   //do we have a problem?

   //need to add a check that B-1A matches my dualray.A calculation
   //in generate vars... it might be ok and yb not ok, because the
   //adjustments in simplex might only be related to rhs... 

   
   //i don't think Osi returns tabRhs... that should be changed
   CPXgetbhead(siCpx->getEnvironmentPtr(), 
	       siCpx->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL), 
	       basics, tabRhs);

   //as a sanity check print out the basis status next to the yb vs tabRhs
   //calculation.... let's see why and where things don't match up... 
   
   //yb, where y is a row of B-1 (note, can get from bhead?)
   UTIL_DEBUG(m_param.LogDebugLevel, 6,  
	      (*m_osLog) << "\nB-1:";
	      for(r = 0; r < m; r++){
		 yb[r] = 0.0;
		 m_masterSI->getBInvRow(r, bInvRow);
		 (*m_osLog) << "\nB-1Row r: " << r << ": " << endl;
		 for(b = 0; b < m; b++){
		    yb[r] += bInvRow[b] * rowRhs[b];
		    (*m_osLog) << setw(6) << "bind: "
			       << setw(4) << basics[b]
			       << setw(12) << bInvRow[b] 
			       << " [" 
			       << setw(12) << rowRhs[b] 
			       << "] "
			       << setw(8) << " +=: " 
			       << setw(12) << bInvRow[b] * rowRhs[b] 
			       << setw(8) << " yb: "
			       << setw(12) << yb[r]
			       << setw(8) << " tabRhs: " 
			       << setw(12) << tabRhs[r] << endl;
		 }
		 if(!UtilIsZero(yb[r] - tabRhs[r]))
		    (*m_osLog) << " DIFF is " << yb[r] - tabRhs[r] << endl;
		 assert(UtilIsZero(yb[r] - tabRhs[r], 1.0e-4));
	      }
	      );


   for(r = 0; r < m; r++){
      yb[r] = 0.0;
      m_masterSI->getBInvRow(r, bInvRow);
      for(b = 0; b < m; b++){
	 yb[r] += bInvRow[b] * rowRhs[b];//(B-1)_r.b
      }
      if(!UtilIsZero(yb[r] - tabRhs[r])){
	 (*m_osLog) << " DIFF is " << yb[r] - tabRhs[r] << endl;
	 (*m_osLog) << "\nB-1Row r: " << r << ": basics[r]=" << basics[r] 
		    << endl;
	 yb[r] = 0.0;
	 for(b = 0; b < m; b++){
	    if(UtilIsZero(bInvRow[b]))
	       continue;
	    yb[r] += bInvRow[b] * rowRhs[b];
	    (*m_osLog) << setw(6) << "bind: "
		       << setw(4) << basics[b]
		       << setw(12) << bInvRow[b] 
		       << " [" 
		       << setw(12) << rowRhs[b]; 
	    if(basics[b] < 0){//== -rowIndex-1
	       (*m_osLog) << " sense = " << rowSense[-(basics[b]+1)];
	    }
	    (*m_osLog) << "] "
		       << setw(8) << " +=: " 
		       << setw(12) << bInvRow[b] * rowRhs[b] 
		       << setw(8) << " yb: "
		       << setw(12) << yb[r]
		       << setw(8) << " tabRhs: " 
		       << setw(12) << tabRhs[r] << endl;
	 }
      }
      
      //assert(UtilIsZero(yb[r] - tabRhs[r], 1.0e-4));
   }
		 

   for(r = 0; r < m; r++){
      if(UtilIsZero(tabRhs[r]))
	 continue;
   
      //all pos case? if yb < 0 (then we want to minimize B-1Ax, x in P')
      //all neg case? if yb > 0 (then we want to maximize B-1Ax, x in P')
      UTIL_DEBUG(m_param.LogDebugLevel, 6,  
		 (*m_osLog) << "\nB-1A:";
		 );
      if(tabRhs[r] > 0){    //instead of yb
	 //Ted also checks that it is a slack var here - why?  
	 bool allneg = true;      
	 m_masterSI->getBInvARow(r, bInvARow);
	 UTIL_DEBUG(m_param.LogDebugLevel, 6,  
		    (*m_osLog) << "\nB-1ARow r: " << r << ": ";
		    );
	 allneg = true;
	 for(c = 0; c < n; c++){
	    UTIL_DEBUG(m_param.LogDebugLevel, 6,  
		       (*m_osLog) << bInvARow[c] << " ";	
		       );
	    if(bInvARow[c] >= DecompEpsilon){
	       allneg = false;
	       break;
	    }
	 }
	 if(allneg){
	    UTIL_DEBUG(m_param.LogDebugLevel, 6,  
		       (*m_osLog) << " ---> allneg";
		       );
	    double * dualRay  = new double[m];	    
	    m_masterSI->getBInvRow(r, dualRay);
	    transform(dualRay, dualRay + m, dualRay, negate<double>());   
	    rays.push_back(dualRay);	 
	 }
      }
      else{
	 bool allpos = true;
	 m_masterSI->getBInvARow(r, bInvARow);
	 UTIL_DEBUG(m_param.LogDebugLevel, 6,  
		    (*m_osLog) << "\nB-1ARow r: " << r << ": ";
		    );
	 allpos = true;
	 for(c = 0; c < n; c++){
	    UTIL_DEBUG(m_param.LogDebugLevel, 6,  
		       (*m_osLog) << bInvARow[c] << " ";
		       );
	    if(bInvARow[c] <= -DecompEpsilon){
	       allpos = false;
	       break;
	    }
	 }
	 if(allpos){
	    UTIL_DEBUG(m_param.LogDebugLevel, 6,  
		       (*m_osLog) << " ---> allpos";
		       );
	    double * dualRay  = new double[m];
	    m_masterSI->getBInvRow(r, dualRay);
	    rays.push_back(dualRay);	 
	 }
      }
   }

   UTIL_DELARR(tabRhs);
   UTIL_DELARR(basics);
   UTIL_DELARR(yb);
   UTIL_DELARR(bInvRow);
   UTIL_DELARR(bInvARow);
   m_masterSI->disableSimplexInterface();   

   printf("rays.size = %d\n", static_cast<int>(rays.size()));   
   if(rays.size() <= 0){
      printf("NO RAYS using standard lookup - try dualfarkas\n");
      double   proof_p;
      double * dualRay = new double[m];
      CPXdualfarkas(siCpx->getEnvironmentPtr(),
		    siCpx->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL), 
		    dualRay, &proof_p);
      (*m_osLog) << "After dual farkas proof_p = " << proof_p << "\n";
      transform(dualRay, dualRay + m, dualRay, negate<double>());   
      for(int i = 0; i < m; i++){
	 printf("dualRay[%d]: %g\n", i, dualRay[i]);
      }
      rays.push_back(dualRay);
   }
   

   //NOTE: you will have dup rays here - need to filter out... 
   printf("rays.size = %d", static_cast<int>(rays.size()));
   for(size_t i = 0; i < rays.size(); i++){
      bool isProof = isDualRayInfProof(rays[i],
				       m_masterSI->getMatrixByRow(),
				       m_masterSI->getColLower(),
				       m_masterSI->getColUpper(),
				       m_masterSI->getRightHandSide(),
				       NULL);      
      if(!isProof){
	 isDualRayInfProof(rays[i],
			   m_masterSI->getMatrixByRow(),
			   m_masterSI->getColLower(),
			   m_masterSI->getColUpper(),
			   m_masterSI->getRightHandSide(),
			   m_osLog);
      }
      assert(isProof);
   }
   assert(rays.size() > 0);
   return rays;
}
#else
//------------------------------------------------------------------------- //
//TEST THIS
vector<double*> DecompAlgo::getDualRays(int maxNumRays){
   OsiCpxSolverInterface * siCpx
      = dynamic_cast<OsiCpxSolverInterface*>(m_masterSI);
   const int m = m_masterSI->getNumRows();
   const int n = m_masterSI->getNumCols();
   double proof_p;
   bool   isProof;
  
   vector<double*> rays;
   double * ray = new double[m];
   int err 
      = CPXdualfarkas(siCpx->getEnvironmentPtr(),
		      siCpx->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL), 
		      ray, &proof_p);//proof_p
   if(err){
      cerr << "CPXdualfarkas returns err " << err << endl;
      abort();
   }
   
   cout << "After dual farkas proof_p = " << proof_p << "\n";



   //We have to flip because in this context we want to max B-1Ax, x in P'
   double * pneg = new double[m];
   transform(ray, ray + m, pneg, negate<double>());   
   rays.push_back(pneg);


#if 1
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 5,
	      bool isProof = isDualRayInfProof(rays[0],
					       m_masterSI->getMatrixByRow(),
					       m_masterSI->getColLower(),
					       m_masterSI->getColUpper(),
					       m_masterSI->getRightHandSide(),
					       NULL);
	      printf("isProof = %d\n", isProof);
	      printBasisInfo(m_masterSI, m_osLog);
	      fflush(stdout);
	      if(!isProof){
		 isDualRayInfProof(ray,
				   m_masterSI->getMatrixByRow(),
				   m_masterSI->getColLower(),
				   m_masterSI->getColUpper(),
				   m_masterSI->getRightHandSide(),
				   m_osLog);
		 printBasisInfo(m_masterSI, m_osLog);
		 fflush(stdout);
	      }
	      );
   assert(isDualRayInfProof(ray,
			    m_masterSI->getMatrixByRow(),
			    m_masterSI->getColLower(),
			    m_masterSI->getColUpper(),
			    m_masterSI->getRightHandSide(),
			    NULL));
#endif

  

   return rays;
}
#endif
#endif



#ifdef __DECOMP_LP_CLP__


#ifdef USE_MULTI_RAY
//------------------------------------------------------------------------- //
//STOP - try this... 
vector<double*> DecompAlgo::getDualRays(int maxNumRays){
   const int m = m_masterSI->getNumRows();
   const int n = m_masterSI->getNumCols();
   const double * rowRhs    = m_masterSI->getRightHandSide();
   const char   * rowSense  = m_masterSI->getRowSense();
   int i, r, b, c;
   vector<double*> rays;

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "getDualRays()", m_param.LogDebugLevel, 2);


   UTIL_DEBUG(m_param.LogDebugLevel, 5,  
	      for(r = 0; r < m; r++){
		 (*m_osLog) << "Row r: " << r << " sense: " << rowSense[r]
			    << " rhs: " << rowRhs[r] << endl;
	      }
	      );

   m_masterSI->enableSimplexInterface(false);

   //with simplex interface, this is slightly different... 
   const double * colSolution = m_masterSI->getColSolution();
   const double * rowAct      = m_masterSI->getRowActivity();//==slacks?

   double * tabRhs   = new double[m];//osi_clp does not give this?
   //B-1b just equals x, but what if art column then is slack var
   
   int    * basics   = new int[m];
   double * yb       = new double[m];
   double * bInvRow  = new double[m];
   double * bInvARow = new double[n];

   m_masterSI->getBasics(basics);

   for(r = 0; r < m; r++){
      i = basics[r];
      if(i < n){
	 tabRhs[r] = colSolution[i]; //should == B-1b
	 //printf("tabRhs[c:%d]: %g\n", i, tabRhs[r]);
      }
      else{
	 //this really should be slack vars...
	 //assuming clp does Ax-Is = b, s = ax-b ??? nope... 
	 //tabRhs[r] = rowAct[i - n] - rowRhs[i - n];
	 tabRhs[r] = rowRhs[i - n] - rowAct[i - n];
	 //printf("tabRhs[r:%d]: %g [act: %g rhs: %g sense: %c]\n",
	 //	i-n, tabRhs[r], rowAct[i-n], rowRhs[i-n], rowSense[i-n]);
      }
   }
   

   //as a sanity check print out the basis status next to the yb vs tabRhs
   //calculation.... let's see why and where things don't match up... 
   
   //yb, where y is a row of B-1 (note, can get from bhead?)

   //B-1b is tab rhs, is this equivalent to x for struct columns?


   UTIL_DEBUG(m_param.LogDebugLevel, 6,  
	      (*m_osLog) << "\nB-1:";
	      for(r = 0; r < m; r++){		 
		 if(UtilIsZero(tabRhs[r]))
		    continue;
		 yb[r] = 0.0;
		 m_masterSI->getBInvRow(r, bInvRow);
		 (*m_osLog) << "\nB-1Row r: " << r << ": " << endl;
		 for(b = 0; b < m; b++){
		    yb[r] += bInvRow[b] * rowRhs[b];		    
		    (*m_osLog) << setw(6) << "bind: "
			       << setw(4) << basics[b]
			       << setw(12) << bInvRow[b] 
			       << " [" 
			       << setw(12) << rowRhs[b] 
			       << "] "
			       << setw(8) << " +=: " 
			       << setw(12) << bInvRow[b] * rowRhs[b] 
			       << setw(8) << " yb: "
			       << setw(12) << yb[r]
			       << setw(8) << " tabRhs: " 
			       << setw(12) << tabRhs[r]
			       << endl;
		 }
		 if(!UtilIsZero(yb[r] - tabRhs[r]))
		    (*m_osLog) << " DIFF is " << yb[r] - tabRhs[r] << endl;
		 assert(UtilIsZero(yb[r] - tabRhs[r], 1.0e-4));
	      }
	      );


   for(r = 0; r < m; r++){
      if(UtilIsZero(tabRhs[r]))
	 continue;
   
      //all pos case? if yb < 0 (then we want to minimize B-1Ax, x in P')
      //all neg case? if yb > 0 (then we want to maximize B-1Ax, x in P')
      if(tabRhs[r] > 0){ //instead of yb
	 //Ted also checks that it is a slack var here - why?  
	 bool allneg = true;
	 //not getting back slacks part here... need?
	 m_masterSI->getBInvARow(r, bInvARow);
	 UTIL_DEBUG(m_param.LogDebugLevel, 6,  
		    (*m_osLog) << "B-1ARow r: " << r << ": ";
		    );
	 allneg = true;
	 for(c = 0; c < n; c++){
	    UTIL_DEBUG(m_param.LogDebugLevel, 6,  
		       (*m_osLog) << bInvARow[c] << " ";	
		       );
	    if(bInvARow[c] >= DecompEpsilon){
	       allneg = false;
	       break;
	    }
	 }
	 if(allneg){
	    UTIL_DEBUG(m_param.LogDebugLevel, 6,  
		       (*m_osLog) << " ---> allneg";
		       );
	    double * dualRay  = new double[m];	    
	    m_masterSI->getBInvRow(r, dualRay);
	    transform(dualRay, dualRay + m, dualRay, negate<double>());   
	    rays.push_back(dualRay);	 
	 }
      }
      else{
	 bool allpos = true;
	 m_masterSI->getBInvARow(r, bInvARow);
	 UTIL_DEBUG(m_param.LogDebugLevel, 6,  
		    (*m_osLog) << "B-1ARow r: " << r << ": ";
		    );
	 allpos = true;
	 for(c = 0; c < n; c++){
	    UTIL_DEBUG(m_param.LogDebugLevel, 6,  
		       (*m_osLog) << bInvARow[c] << " ";
		       );
	    if(bInvARow[c] <= -DecompEpsilon){
	       allpos = false;
	       break;
	    }
	 }
	 if(allpos){
	    UTIL_DEBUG(m_param.LogDebugLevel, 6,  
		       (*m_osLog) << " ---> allpos";
		       );
	    double * dualRay  = new double[m];
	    m_masterSI->getBInvRow(r, dualRay);
	    rays.push_back(dualRay);	 
	 }
      }
      UTIL_DEBUG(m_param.LogDebugLevel, 6,  
		 (*m_osLog) << endl;
		 );
   }


   UTIL_DELARR(basics);
   UTIL_DELARR(yb);
   UTIL_DELARR(bInvRow);
   UTIL_DELARR(bInvARow);
   m_masterSI->disableSimplexInterface();   

   /*
     if(rays.size() <= 0){
     double   proof_p;
     double * dualRay = new double[m];
     CPXdualfarkas(siCpx->getEnvironmentPtr(),
     siCpx->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL), 
     dualRay, &proof_p);
     (*m_osLog) << "After dual farkas proof_p = " << proof_p << "\n";
     transform(dualRay, dualRay + m, dualRay, negate<double>());   
     for(int i = 0; i < m; i++){
     printf("dualRay[%d]: %g\n", i, dualRay[i]);
     }
     rays.push_back(dualRay);
     }
   */
   


   //NOTE: you will have dup rays here - need to filter out... 
   UTIL_DEBUG(m_param.LogDebugLevel, 5,  
	      (*m_osLog) << "Number of Rays = " << rays.size() << endl;
	      );

   for(int i = 0; i < (int)rays.size(); i++){
      bool isProof = isDualRayInfProof(rays[i],
				       m_masterSI->getMatrixByRow(),
				       m_masterSI->getColLower(),
				       m_masterSI->getColUpper(),
				       m_masterSI->getRightHandSide(),
				       NULL);
      if(!isProof){
	 isDualRayInfProof(rays[i],
			   m_masterSI->getMatrixByRow(),
			   m_masterSI->getColLower(),
			   m_masterSI->getColUpper(),
			   m_masterSI->getRightHandSide(),
			   m_osLog);
      }
      assert(isProof);
   }

   assert(rays.size() > 0);

   UTIL_DELARR(tabRhs);

   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "getDualRays()", m_param.LogDebugLevel, 2);
   
   return rays;
}
#else

//===========================================================================//
vector<double*> DecompAlgo::getDualRays(int maxNumRays){

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "getDualRays()", m_param.LogDebugLevel, 2);
   vector<double*> raysT = m_masterSI->getDualRays(maxNumRays);   
   const double * rayT = raysT[0];
   assert(rayT);


   //stop
   //what is yb, that will tell me if i want to opt over uA or -uA
   //y^T b
   int   i;
   const CoinPackedMatrix * rowMatrix = m_masterSI->getMatrixByRow();
   const double           * rowRhs    = m_masterSI->getRightHandSide();
   const int                m         = rowMatrix->getNumRows();
   double yb = 0.0;
   for(i = 0; i < m; i++){
      yb += rayT[i] * rowRhs[i]; //safe to use rowRhs? or flips in tab going on
   }
   (*m_osLog) << " yb = " << yb << endl;
   //need tabRhs if doing this way? 
   



   //see Clp/examples/decompose.cpp
   //   he flips the infeasibility ray (always...)   
   //---    yA >= 0, yb < 0, or  --> find a yAs <= 0 (min)
   //---    yA <= 0, yb > 0 ??   --> find a yAs >= 0 (max <--> -min)

   vector<double*> rays;
   if(yb > 0){
      double * pneg = new double[m];
      transform(rayT, rayT + m, pneg, negate<double>());
      rays.push_back(pneg);
   }
   else{
      rays.push_back(raysT[0]);
   }
   const double * ray = rays[0];
   assert(ray);

      

#if 1
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 5,
	      bool isProof = isDualRayInfProof(ray,
					       m_masterSI->getMatrixByRow(),
					       m_masterSI->getColLower(),
					       m_masterSI->getColUpper(),
					       m_masterSI->getRightHandSide(),
					       NULL);
	      printf("isProof = %d\n", isProof);
	      fflush(stdout);
	      if(!isProof){
		 isDualRayInfProof(ray,
				   m_masterSI->getMatrixByRow(),
				   m_masterSI->getColLower(),
				   m_masterSI->getColUpper(),
				   m_masterSI->getRightHandSide(),
				   m_osLog);
		 printBasisInfo(m_masterSI, m_osLog);
		 fflush(stdout);
	      }
	      );;
   assert(isDualRayInfProof(ray,
			    m_masterSI->getMatrixByRow(),
			    m_masterSI->getColLower(),
			    m_masterSI->getColUpper(),
			    m_masterSI->getRightHandSide(),
			    NULL));
#endif
   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "getDualRays()", m_param.LogDebugLevel, 2);
   return rays;
}
#endif
#endif

//------------------------------------------------------------------------ //
int DecompAlgo::generateVarsFea(DecompVarList    & newVars, 
                                double           & mostNegReducedCost){
   //---
   //--- solve min{s in F' | RC[s]} to generate new variables
   //--- 
   //--- if LP was feasible,   
   //---    then RC[s] = c.s - (uhat.A''.s + alpha)
   //---               = c.s -  uhat.A''.s - alpha
   //---               = (c - uhat.A'')s - alpha
   //---
   //--- The master LP was formed in the following order
   //---   (A''s) lam[s]  >= b'' - from the original core [A'', b'']
   //---   sum{s} lam[s]   = 1   - convexity constraint
   //--- But, we may have added cuts - which conceptually are added
   //--- into [A'', b'']. But, in reality they are simply appended to
   //--- the end of the LP matrix. So, when we get back the dual vector, 
   //--- we have to be aware of where alpha actually is.
   //---

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "generateVarsFea()", m_param.LogDebugLevel, 2);
   m_stats.timerOther1.reset();

   //---
   //--- TODO:
   //--- Blocks...
   //---  How do we adjust this for different blocks?
   //---  Each block has exactly one convexity constraint. Otherwise,
   //---  the coefficient in the matrix is 0, so there is no effect.
   //---  So, all we need to do is use the approriate block's dual variable
   //---  call it alpha.
   //---
   //--- At this point, we have many strategies -- we could do round robin
   //--- and only generate one var at a time. Or, we can generate vars for
   //--- every block every time (default for now).
   //---

   int                   i, b, r;
   DecompVarList         potentialVars;
   DecompConstraintSet * modelCore   = m_modelCore.getModel();
   const int      m             = m_masterSI->getNumRows();
   int            nBaseCoreRows = modelCore->nBaseRows;
   const int      nCoreCols     = modelCore->getNumCols();
   const double * u             = NULL;
   const double   epsilonRedCost= 1.0e-4;//make option
   const double * origObjective = getOrigObjective();
   double       * redCostX      = NULL;   
   double         alpha         = 0.0;
   int    whichBlock;
   double varRedCost;
   DecompVarList::iterator it;

   if(m_algo == DECOMP)
      nBaseCoreRows = nCoreCols;
   
   assert(!m_masterSI->isProvenPrimalInfeasible());
   u = getRowPrice();


   
   //THINK: includes artificials now
   redCostX = new double[nCoreCols]; // (c - uhat.A") in x-space
   assert(redCostX);

   //---
   //--- PC: get dual vector
   //---   u     --> (A''s) lam[s]  >= b'' --> from core [A'', b'']
   //---   alpha --> sum{s} lam[s]   = 1   --> convexity constraint
   //---
   //--- NOTE: Sign of alpha is negative which is different than infeas case
   //---       the solve relax function adds alpha so no sign switch needed.
   //---
   //---       We flip the sign here.
   //---
      
   //---
   //--- remove the convexity constraint(s) from the dual vector
   //---   TODO/THINK: this is fugly
   //--- this is done because the convexity constraint acts like an
   //--- objective offset in subproblem - it might be cleaner to manage
   //--- it as an offset and just check for < 0 directly, rather than
   //--- less than alpha -- sign switches are a little messy
   //---

   //THINK: we should be checked col-type to make sure we get the
   //  right rows for the convexity constraints

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 4,
	      (*m_osLog) << "m            =" << m << endl;
	      (*m_osLog) << "numConvexCon =" << m_numConvexCon << endl;
	      (*m_osLog) << "nBaseCoreRows=" << nBaseCoreRows  << endl;
	      );

   double * u_adjusted = new double[m - m_numConvexCon];
   assert(u_adjusted);

   //copy the dual vector for original core rows
   CoinDisjointCopyN(u, nBaseCoreRows, u_adjusted);

   //---
   //--- sanity check - make sure we are only skipping convexity rows
   //---    everything else has a dual we want to use
   //---
   if(m_param.DebugLevel >= 1){
      for(r = 0; r < nBaseCoreRows; r++){
	 assert(m_masterRowType[r] == DecompRow_Original ||
		m_masterRowType[r] == DecompRow_Branch);
      }
      for(r = nBaseCoreRows; r < nBaseCoreRows + m_numConvexCon; r++){
	 assert(m_masterRowType[r] == DecompRow_Convex);
      }
      for(r = nBaseCoreRows + m_numConvexCon; r < m; r++){
	 assert(m_masterRowType[r] == DecompRow_Cut);
      }
   }

   // if no cuts, don't need to do any of this
   // if DECOMP, don't need to do any of this
   //append dual vector for any added cuts (skip over convexity constraints)
   assert((m - nBaseCoreRows - m_numConvexCon) == 
	  getNumRowType(DecompRow_Cut));
   CoinDisjointCopyN(u + nBaseCoreRows + m_numConvexCon,  //from
		     m - nBaseCoreRows - m_numConvexCon,  //size
		     u_adjusted + nBaseCoreRows);         //to
   
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 5,
	      for(i = 0; i < m; i++){
		 if(!UtilIsZero(u[i], DecompEpsilon)){
		    (*m_osLog) << "u[" << setw(5) << i << " ]: " 
			       << setw(12) << UtilDblToStr(u[i],3) << " --> " 
			       << DecompRowTypeStr[m_masterRowType[i]] << "\n";
		 }
	      }
	      );

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 5,
	      for(i = 0; i < (m-m_numConvexCon); i++){
		 if(!UtilIsZero(u_adjusted[i], DecompEpsilon)){
		    (*m_osLog) << "uAdj[" << setw(5) << i << " ]: " 
			       << setw(12) << UtilDblToStr(u_adjusted[i],3) 
			       << endl;	
		 }
	      }
	      );

   //---
   //--- calculate reduced costs
   
   if(m_algo == DECOMP){
      //---
      //--- in DW, we use (c-uA'')x, where A'' is the core matrix
      //---    u, in this case has dimension = #core rows
      //--- in D , we use (c-u   )x, we don't use the core matrix
      //---    u, in this case has dimension = #core cols
      //---
      assert((m - m_numConvexCon) == modelCore->M->getNumCols());
      for(i = 0; i < nCoreCols; i++)
	 redCostX[i] = u_adjusted[i];
   }
   else{
      assert((m - m_numConvexCon) == modelCore->M->getNumRows());
      modelCore->M->transposeTimes(u_adjusted, redCostX);
   }
   
   //--- 
   //--- if in Phase I, c=0
   //---   
   if(m_phase == PHASE_PRICE1){
      for(i = 0; i < nCoreCols; i++)
	 redCostX[i] = -redCostX[i];
   }
   else{
      for(i = 0; i < nCoreCols; i++)
	 redCostX[i] = origObjective[i] - redCostX[i];
   }
   

   int var_index = 0;   
   const double * objC = m_masterSI->getObjCoefficients();
   const double * rcLP = m_masterSI->getReducedCost();

   if(m_param.DebugLevel >= 1){   
      //---
      //--- sanity check on duals returned
      //---   complementary slackness (c-uA)x = 0
      //---   also check that the given reduced cost matches the 
      //---   hand calculation
      //---
      const double           * x     = m_masterSI->getColSolution();
      const double           * pi    = m_masterSI->getRowPrice();
      const int                nCols = m_masterSI->getNumCols();
      //const int                nRows = m_masterSI->getNumRows();
      const CoinPackedMatrix * M     = m_masterSI->getMatrixByRow();
      double                 * uA    = new double[nCols];
      M->transposeTimes(pi, uA);
#ifndef DO_INTERIOR 
      for(i = 0; i < nCols; i++){
	 if(!UtilIsZero( x[i], 1.0e-5 ) &&
	    !UtilIsZero( (objC[i] - uA[i]) * x[i], 1.0e-4 ) ){
	    printf("ERR in COMPL-SLACK i:%d objC:%15.10f uA:%15.10f x:%15.10f\n",
		   i, objC[i], uA[i], x[i]); fflush(stdout);
	    assert(0);
	 }
	 if(!UtilIsZero( (objC[i] - uA[i]) - rcLP[i], 1.0e-4 ) ){
	    printf("ERR in RC i:%d objC:%15.10f uA:%15.10f RCLP:%15.10f\n",
		   i, objC[i], uA[i], rcLP[i]); fflush(stdout);
	    assert(0);
	 }      
      }
#endif
      UTIL_DELARR(uA);
   }


   //---
   //--- sanity check - none of the columns currently in master
   //--- should have negative reduced cost
   //---   m_vars contains the variables (in x-space) that have
   //---   been pushed into the master LP (assumes no compression)
   //---
   if(m_param.DebugLevel >= 1){         
      for(it = m_vars.begin(); it != m_vars.end(); it++){
	 double redCost = 0.0; 
	 
	 //m_s      is a sparse vector in x-space (the column)
	 //redCostX is a dense  vector in x-space (the cost in subproblem)
	 b       = (*it)->getBlockId();
	 redCost = (*it)->m_s.dotProduct(redCostX);//??
	 alpha   = u[nBaseCoreRows + b];
	 assert(m_masterRowType[nBaseCoreRows+b] == DecompRow_Convex);	 
	 assert(isMasterColStructural((*it)->getColMasterIndex()));
	 
	 UTIL_DEBUG(m_app->m_param.LogDebugLevel, 5,
		    (*m_osLog) 
		    << "MasterColIndex = " 
		    << setw(6) << (*it)->getColMasterIndex()
		    << "Block = " 
		    << setw(3) << b
		    << "LPRedCost = " << setw(10) 
		    << UtilDblToStr(rcLP[(*it)->getColMasterIndex()],5)
		    << "CalcRedCost = " << setw(10)
		    << UtilDblToStr(redCost-alpha, 5)
		    << "ObjCost = " << setw(10)
		    << UtilDblToStr(objC[(*it)->getColMasterIndex()], 5)
		    << "Alpha = " << setw(10)
		    << UtilDblToStr(alpha, 5)
		    << endl; );
	 
	 //---
	 //--- sanity check - none of the columns currently in master
	 //--- should have negative reduced cost
	 //---
	 double rcLPi = rcLP[(*it)->getColMasterIndex()];
	 if(rcLPi < -epsilonRedCost ){  
	    (*m_osLog) << "VAR v-index:" << var_index++ 
		       << " m-index: " << (*it)->getColMasterIndex()
		       << " b-index: " << b
		       << " rcLP: "   << rcLPi
		       << endl; 
	    (*it)->print(m_osLog, modelCore->colNames, 
			 const_cast<double*>(redCostX));
	    (*m_osLog) << "******** ERROR ********" << endl;	    
	    assert(0);
	 }
      
	 //---
	 //--- check that objective in LP and calculated objective match
	 //---
	 if(m_phase == PHASE_PRICE2){
	    double objCalc = (*it)->m_s.dotProduct(origObjective);
	    if(!UtilIsZero(objCalc - objC[(*it)->getColMasterIndex()], 
			   1.0e-3)){
	       (*m_osLog) << "VAR v-index:" << var_index++ 
			  << " m-index: " << (*it)->getColMasterIndex()
			  << " b-index: " << b
			  << " objLP: "   << objC[(*it)->getColMasterIndex()]
			  << " objCalc: " << objCalc
			  << endl; 
	       (*it)->print(m_osLog, modelCore->colNames, 
			    const_cast<double*>(origObjective));
	       (*m_osLog) << "******** ERROR ********" << endl;	       
	       assert(0);
	    }
	 }
	 
	 //---
	 //--- check that LP reduced cost and calculated reduced cost
	 //---  match up
	 //--- in the case of using dual smoothing, we cannot do this check
	 //---  since the calculated reduced cost is based on the smoothed
	 //---  duals
	 //---  
#ifndef DUAL_SMOOTHING
	 if(!UtilIsZero(rcLP[(*it)->getColMasterIndex()] 
			- (redCost-alpha), 1.0e-3)){
	    //---
	    //--- this whole next section is an expansion of log
	    //---   when there is an issue found that the solver 
	    //---   returns RC that doesn't match the one calculated
	    //---   based on core matrix and duals
	    //---
	    (*m_osLog) << "VAR v-index:" << var_index++ 
		       << " m-index: " << (*it)->getColMasterIndex()
		       << " b-index: " << b
		       << " rc: "      << redCost 
		       << " alpha: "   << alpha
		       << " rc-a: "    << redCost-alpha
		       << " RCLP: "    << rcLP[(*it)->getColMasterIndex()]
		       << endl; 
	    
	    //---
	    //--- this, plus alpha shows the calculation of the red-cost in 
	    //---   x-space, next, look at the same calculation in lambda-space
	    //---
	    (*it)->print(m_osLog, modelCore->colNames, 
			 const_cast<double*>(redCostX));
	    (*m_osLog) << "******** ERROR ********" << endl;
	    
	    //---
	    //--- the rows in lambda that show up should be ones where 
	    //---   these components show up in original A''
	    //---
	    double* uA2    = new double[nCoreCols];
	    modelCore->M->transposeTimes(u_adjusted, uA2);
	    (*it)->print(m_osLog, modelCore->colNames, 
			 const_cast<double*>(uA2));
	    UTIL_DELARR(uA2);
	    
	    
	    //---
	    //--- RC of a col of lambda-space is u.A[i]
	    //---
	    (*m_osLog) << " objLP: "   
		       << UtilDblToStr(objC[(*it)->getColMasterIndex()],4) 
		       << endl;
	    
	    //---
	    //--- recalc column of master A.s for this s
	    //---
	    if(m_algo != DECOMP){
	       double * denseS = new double[modelCore->getNumCols()];
	       (*it)->fillDenseArr(modelCore->getNumCols(), denseS);
	       int r;
	       const CoinPackedMatrix * Mr = modelCore->getMatrix();	       
	       for(r = 0; r < modelCore->getNumRows(); r++){
		  printf("\nROW %d\n", r);
		  CoinShallowPackedVector vec = Mr->getVector(r);
		  UtilPrintPackedVector(vec, m_osLog, 
					modelCore->getColNames(), 
					denseS);
		  
	       }
	       UTIL_DELARR(denseS);
	    }
	    
	    
	    const CoinPackedMatrix * Mc = m_masterSI->getMatrixByCol();
	    CoinShallowPackedVector vec 
	       = Mc->getVector((*it)->getColMasterIndex());
	    UtilPrintPackedVector(vec, m_osLog, m_masterSI->getColNames(), u);
	    double uA = vec.dotProduct(u);
	    (*m_osLog) << " objLP: "   
		       << UtilDblToStr(objC[(*it)->getColMasterIndex()],4) 
		       << endl;
	    (*m_osLog) << " uA   : "   << UtilDblToStr(uA,4) << endl;
	    (*m_osLog) << " RC   : "   
		       << UtilDblToStr(objC[(*it)->getColMasterIndex()] - uA,
				       4) << endl;
	    (*m_osLog) << " RCLP : "   
		       << UtilDblToStr(rcLP[(*it)->getColMasterIndex()],4) 
		       << endl;
	    assert(0);
	    (*m_osLog) << endl;
	 } //END: if(!UtilIsZero(rcLP[(*it)->getColMasterIndex()] ...
#endif      
      } //END: for(it = m_vars.begin(); it != m_vars.end(); it++)
   } //END: if(m_param.DebugLevel >= 1)

  




   //---
   //--- if doing round-robin, solve just one block unless
   //---   in PhaseI (do all blocks)
   //---   in PhaseII every n iterations (so we can get a valid update to LB)
   //---
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
	      (*m_osLog) 
	      << "RoundRobin iterSinceAll= " << m_rrIterSinceAll 
	      << " lastBlock= "              <<  m_rrLastBlock << endl;
	      );
   int doAllBlocks = false;
   if(m_rrIterSinceAll >= (m_param.RoundRobinInterval * m_numConvexCon)){
      doAllBlocks      = true;
      m_rrIterSinceAll = 0;
   }
   
   vector<double>     mostNegRCvec(m_numConvexCon, DecompInf);  
   DecompSolverResult solveResult(nCoreCols);
   OsiSolverInterface          * subprobSI   = NULL;

   //---
   //--- solve min{ (c - u.A'')x - alpha |  x in F'}
   //---         
   if(doAllBlocks){

      //---
      //--- here is where you would thread it
      //---
#ifdef RELAXED_THREADED
      //round-robin assign blocks to threads
      printf("===== START Threaded solve of subproblems. =====\n");
      pthread_t                  threads[NTHREADS];
      pthread_attr_t             pthread_custom_attr;
      SolveRelaxedThreadArgs     arg[NTHREADS];
      DecompSolverResult      ** solverResultT = new DecompSolverResult*[NTHREADS]; 
      DecompVarList              potentialVarsT[NTHREADS];
      int * batch[NTHREADS];
      int   batchSize[NTHREADS];

      int t, b;
      for(t = 0; t < NTHREADS; t++){
	 batch[t]         = new int[m_numConvexCon];
	 batchSize[t]     = 0;
	 solverResultT[t] = new DecompSolverResult(nCoreCols);
      }

      t = 0;
      for(b = 0; b < m_numConvexCon; b++){
	 batch[t][batchSize[t]++] = b;
	 if((++t) % NTHREADS == 0)
	    t = 0;
      }

      pthread_attr_init(&pthread_custom_attr);

      bool useCutoff = false;      
      if(m_phase == PHASE_PRICE2)
	 useCutoff = m_param.SubProbUseCutoff;
      m_isColGenExact = false; //so it goes to IP solver - ugh

      for(t = 0; t < NTHREADS; t++){

	 arg[t].algo          = this;
	 arg[t].threadId      = t;
	 arg[t].batchSz       = batchSize[t];
	 arg[t].batch         = batch[t];	 
	 arg[t].nBaseCoreRows = nBaseCoreRows;
	 arg[t].u             = const_cast<double*>(u); //duals (for alpha)
	 arg[t].redCostX      = redCostX; 
	 arg[t].origCost      = origObjective;
	 arg[t].n_origCols    = nCoreCols;
	 arg[t].vars          = &potentialVarsT[t];

	 pthread_create(&threads[t],
			&pthread_custom_attr,
			solveRelaxedThread, (void*)&arg[t]);
      }
      for(t = 0; t < NTHREADS; t++){
	 pthread_join(threads[t], NULL);
      }
      m_isColGenExact = true;

      //clean-up memory
      for(t = 0; t < NTHREADS; t++){
	 printf("arg[%d].vars size=%d\n", 
		t, static_cast<int>(arg[t].vars->size()));
	 for(it  = arg[t].vars->begin();
	     it != arg[t].vars->end(); it++){
	    varRedCost = (*it)->getReducedCost();
	    whichBlock = (*it)->getBlockId();
	    alpha      = u[nBaseCoreRows + whichBlock];
	    UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
		       (*m_osLog) 
		       << "alpha[block=" << whichBlock << "]:" << alpha 
		       << " varRedCost: " << varRedCost << "\n";
		       );
	    
	    
	 }
      }
      printf("===== END   Threaded solve of subproblems. =====\n");
      for(t = 0; t < NTHREADS; t++){
	 //one function to do this?
	 for(it  = arg[t].vars->begin();
	     it != arg[t].vars->end(); it++){	    
	    potentialVars.push_back(*it);
	 }
      }
      //put the vars from all threads into one vector


#else
      bool useCutoff = false;
      bool useExact  = true;
      if(m_phase == PHASE_PRICE2)
         useCutoff = m_param.SubProbUseCutoff;
      
      map<int, DecompAlgoModel>::iterator mit; 
      for(mit = m_modelRelax.begin(); mit != m_modelRelax.end(); mit++){
         DecompAlgoModel & algoModel = (*mit).second;
         subprobSI = algoModel.getOsi();
         b         = algoModel.getBlockId();
         UTIL_DEBUG(m_app->m_param.LogDebugLevel, 4,
                    (*m_osLog) << "solve relaxed model = "
                    << algoModel.getModelName() << endl;);         
	 //---
	 //--- PC: get dual vector
	 //---   alpha --> sum{s} lam[s]   = 1 - convexity constraint
	 //---
	 alpha = u[nBaseCoreRows + b];
         
	 //TODO: stat return, restrict how many? pass that in to user?	 
	 //---
	 //--- NOTE: the variables coming back include alpha in 
	 //---       calculation of reduced cost
	 //---
	 solveRelaxed(redCostX, 
                      origObjective,
		      alpha, 
		      nCoreCols, 
		      false, //isNested
                      algoModel,
		      &solveResult,
		      potentialVars);
	 //if cutoff delcares infeasible, we know subprob >= 0
	 //  we can use 0 as valid (but possibly weaker bound)
	 if(solveResult.m_isCutoff){
	    mostNegRCvec[b] = min(mostNegRCvec[b],0.0);
	 }
      }

      
      map<int, vector<DecompAlgoModel> >::iterator mivt;
      vector<DecompAlgoModel>           ::iterator vit; 
      //TODO: need time limit and gap limit?
      useExact  = false;
      useCutoff = true;
      for(mivt  = m_modelRelaxNest.begin(); 
          mivt != m_modelRelaxNest.end(); mivt++){
         for(vit  = (*mivt).second.begin(); 
             vit != (*mivt).second.end(); vit++){
            subprobSI = (*vit).getOsi();
            b         = (*vit).getBlockId();
	    alpha     = u[nBaseCoreRows + b];
            UTIL_DEBUG(m_app->m_param.LogDebugLevel, 4,
                       (*m_osLog) << "solve relaxed nested model = "
                       << (*vit).getModelName() << endl;);
            solveRelaxed(redCostX,
                         origObjective,         //original cost vector 
                         alpha,
                         nCoreCols,             //num core columns
			 true,                  //isNested
                         (*vit),
                         &solveResult,          //results
                         potentialVars);        //var list to populate
            if(solveResult.m_isCutoff){
               mostNegRCvec[b] = min(mostNegRCvec[b],0.0);
            }            
         }
      } 





#endif
   }
   else{

      //---
      //--- keep trying until find a block with candidate variable
      //---
      for(i = 0; i < m_numConvexCon; i++){
	 b = (m_rrLastBlock + 1) % m_numConvexCon;

         map<int, DecompAlgoModel>::iterator mit; 
         mit = m_modelRelax.find(b);
         assert(mit != m_modelRelax.end());

         DecompAlgoModel & algoModel = (*mit).second;
         subprobSI = algoModel.getOsi();

	 //if(static_cast<int>(m_subprobSI.size()) > b)
         // subprobSI = m_subprobSI[b];


	 //---
	 //--- PC: get dual vector
	 //---   alpha --> sum{s} lam[s]   = 1 - convexity constraint
	 //---
	 alpha = u[nBaseCoreRows + b];
	 
	 //TODO: stat return, restrict how many? pass that in to user?	 
	 //---
	 //--- NOTE: the variables coming back include alpha in 
	 //---       calculation of reduced cost
	 //---
	 solveRelaxed(redCostX, 
                      origObjective,
		      alpha, 
		      nCoreCols, 
		      false,//isNested
                      algoModel,
		      &solveResult,
		      potentialVars);

	 if(solveResult.m_isCutoff){
	    mostNegRCvec[b] = min(mostNegRCvec[b],0.0);
	 }

	 
	 m_rrIterSinceAll++;
	 m_rrLastBlock = b;

	 bool foundNegRC = false;
	 for(it = potentialVars.begin(); it != potentialVars.end(); it++){
	    varRedCost = (*it)->getReducedCost();
	    if(varRedCost < -epsilonRedCost){ //TODO: strict, -dualTOL?
	       foundNegRC = true;
	       break;
	    }
	 }
	 if(foundNegRC)
	    break;
      }
   }


   for(it = potentialVars.begin(); it != potentialVars.end(); it++){
      varRedCost = (*it)->getReducedCost();
      whichBlock = (*it)->getBlockId();
      alpha      = u[nBaseCoreRows + whichBlock];
      UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
		 (*m_osLog) 
		 << "alpha[block=" << whichBlock << "]:" << alpha 
		 << " varRedCost: " << varRedCost << "\n";
		 );

      //---
      //--- unlikey to happen - but we should check ALL columns
      //---  to see if they are IP feasible - whether or not the
      //---  column has negative red-cost or not
      //---
      (*it)->fillDenseArr(modelCore->getNumCols(),
			  m_memPool.dblArrNCoreCols);            

      //---
      //--- STOP: this isIpFeasible uses isLPFeasible
      //---   isLPFeasible should have 2 settings
      //---    here, it can be true, but might not be
      //---    if checking recompsed point, it must be true, else bug
      //--- 



      if(isIPFeasible(m_memPool.dblArrNCoreCols)){               
	 if(m_app->APPisUserFeasible(m_memPool.dblArrNCoreCols,
				     modelCore->getNumCols(),
				     m_param.TolZero)){
	    DecompSolution * decompSol
	       = new DecompSolution(modelCore->getNumCols(),
				    m_memPool.dblArrNCoreCols,
				    (*it)->getOriginalCost());
	    //TODO: solution pool?
	    m_xhatIPFeas.push_back(decompSol);
	    setObjBoundUB((*it)->getOriginalCost());       
	 }
      }
            
      //---
      //--- for multi-blocks the mostNegReducedCost is the
      //--- sum of the best reduced costs over all blocks
      //---  NOTE: we need all blocks to make it valid
      //---
      if(varRedCost < mostNegRCvec[whichBlock])
	 mostNegRCvec[whichBlock] = varRedCost;
            
      if(varRedCost < -epsilonRedCost){ //TODO: strict, -dualTOL?
	 UTIL_MSG(m_app->m_param.LogDebugLevel, 3,
		  (*m_osLog) << "PUSHING new var with varRedCost= "
		  << UtilDblToStr(varRedCost,5) << endl;);
	 //---
	 //--- the variable has neg reduced cost, push onto list
	 //---
	 newVars.push_back(*it);         
      }
      else{
	 UTIL_DELPTR(*it);	 
      }        
   }

   mostNegReducedCost = 0.0;
   for(b = 0; b < m_numConvexCon; b++){
      mostNegReducedCost += mostNegRCvec[b];
      UTIL_MSG(m_app->m_param.LogDebugLevel, 3,
	       (*m_osLog) 
	       << "mostNegR[block=" << b << "]: " << mostNegRCvec[b]
	       << " mostNegReducedCost: " << mostNegReducedCost << "\n";
	       );      
   }
   
   potentialVars.clear(); //THINK? what does clear do exactly ?   
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 4,   
	      for(it = newVars.begin(); it != newVars.end(); it++){
		 (*it)->print(m_osLog, m_app);
	      }
	      );
   
   
   //---
   //--- free local memory
   //---
   UTIL_DELARR(u_adjusted);
   UTIL_DELARR(redCostX);

   m_stats.thisGenVars.push_back(m_stats.timerOther1.getRealTime()); 
   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "generateVarsFea()", m_param.LogDebugLevel, 2); 
   return static_cast<int>(newVars.size());
}

//------------------------------------------------------------------------ //
int DecompAlgo::generateVars(const DecompStatus   stat,
                             DecompVarList      & newVars, 
                             double             & mostNegReducedCost){
  
   //---
   //--- solve min{s in F' | RC[s]} to generate new variables
   //--- 
   //--- if LP was feasible,   then RC[s] = c.s - (uhat A''s + alpha)
   //--- if LP was infeasible, then RC[s] = (uhat A''s + alpha), dual ray
   //---
   //--- The master LP was formed in the following order
   //---   (A''s) lam[s]  >= b'' - from the original core [A'', b'']
   //---   sum{s} lam[s]   = 1   - convexity constraint
   //--- But, we may have added cuts - which conceptually are added
   //--- into [A'', b'']. But, in reality they are simply appended to
   //--- the end of the LP matrix. So, when we get back the dual vector, 
   //--- we have to be aware of where alpha actually is.
   //---
   
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "generateVars()", m_param.LogDebugLevel, 2);
   return generateVarsFea(newVars, mostNegReducedCost);

   /*#else
     switch(stat){
     case STAT_FEASIBLE:
     return generateVarsFea(newVars, mostNegReducedCost);
     case STAT_INFEASIBLE:
     return generateVarsInf(newVars, mostNegReducedCost);
     default:
     assert(0);
     }*/
   return 0;
}


//TODO - ugh! only PC again?
//------------------------------------------------------------------------ //
//this seems ok for C and PC... but what when we want to do DC within C THINK
int DecompAlgo::generateCuts(double        * xhat,
                             DecompCutList & newCuts){

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "generateCuts()", m_param.LogDebugLevel, 2);
   m_stats.timerOther1.reset();
   DecompConstraintSet          * modelCore   = m_modelCore.getModel();



   m_app->generateCuts(xhat,
		       newCuts);


   //---
   //--- attempt to generate CGL cuts on x??
   //--- the only way this is going to work, is if you carry
   //--- around another OSI problem instance completely in terms of x
   //--- only allow CGL to generate cuts on the original full formulation??
   //--- otherwise, we'd be allowing cuts on cuts - would have to add cuts
   //--- to master and to this version... hmmm... ughh

   //--- for some problems, you can't have the full original formulation
   //--- for PC you probably don't want it anyway... P' comes in from ep's
   //--- you can just generate cuts on Q"? that cut off xhat

   //--- but for C, you need Q' and Q"

   //--- m_masterSI holds the problem in terms of lambda (over Q")
   //--- m_subprobSI holds the problem in terms of x (over P')
   //---
   if(m_param.CutCGL){
      assert(m_cutgenSI);
      if(m_algo == PRICE_AND_CUT){
         //TODO: this could be tighter in the tree
         double gLB = getNodeIndex() ? m_globalLB :
            m_nodeStats.objBest.first;
         m_cutgenSI->setRowLower(m_cutgenObjCutInd, gLB);
      }
      if(m_param.LogDumpModel > 1)
         printCurrentProblem(m_cutgenSI,
                             "CUTGEN",
                             m_nodeStats.nodeIndex,
                             m_nodeStats.cutCallsTotal,
                             m_nodeStats.priceCallsTotal);      
      m_cgl->generateCuts(m_cutgenSI, 
                          m_masterSI,
                          xhat, 
			  modelCore->integerVars,
			  newCuts);
   }



#if 1
   //do DC only if no other cuts were found or if CutDC=2
   //  this is the case of doing for init vars
   if((m_param.CutDC == 1 && newCuts.size() == 0) ||
      (m_param.CutDC == 2)){
     
      printf("\n\n==================================== IN DECOMP\n\n");
      DecompAlgoD D(m_app, m_utilParam,
		    xhat, modelCore->getNumCols());
    
      
      //also might want to use the columns you get here for something... 
      //heur for ubs, etc..

      //either returns a set of cuts or a decomposition, have 
      //that wrap solve()?


      D.solveD(&newCuts);


      //---
      //--- copy the vars generated in passes of DC into initVars
      //---   to warm-start DW master formulation
      //---
      //--- NO: this won't work because it just copies the pointers
      //---  and these will be deleted when D scopes out - see DecompAlgo
      //---  destructor... 
      //m_vars.insert(m_vars.begin(), D.m_vars.begin(), D.m_vars.end());


      //---
      //--- this moves the elements of D.m_vars to m_vars
      //---   this is what we want since D will be deleted after this
      //---
      m_vars.splice(m_vars.end(), D.m_vars);
      //printf("VARS moved into CPM object\n");
      //printVars(m_osLog);//use this to warm start DW



      //a hidden advantage of decomp in BC?
      DecompSolution * bestSol = NULL;
      vector<DecompSolution*>::iterator it;

      double thisBound;
      double bestBoundUB = m_nodeStats.objBest.second;
      for(it  = D.m_xhatIPFeas.begin();
	  it != D.m_xhatIPFeas.end(); it++){
	 thisBound = (*it)->getQuality();
	 printf("From DECOMP, IP Feasible with Quality = %g\n", thisBound);
	 if((*it)->getQuality() <= bestBoundUB){
	    bestBoundUB = (*it)->getQuality();
	    bestSol     = (*it);
	 }
      }
      //need to make copy of solution, since D.m_xhatIpFeas goes out of scope
      if(bestSol){
	 DecompSolution * bestSolCp = new DecompSolution(*bestSol);
	 m_xhatIPFeas.push_back(bestSolCp);	 
	 setObjBoundUB(bestSolCp->getQuality());
	 m_xhatIPBest = bestSolCp;
	 m_xhatIPBest->print();
      }

      //this could also very likely return a new gUB -
      //  for the case when it does find a decomposition 
      //  and luckily it is feasible to original?
     
      //STOP -- 6/6/08
      //if decomp is found, then can't use currently - just looking for
      //farkas -- if decomp is found this means that z_LP = z_DW for that
      //relaxation?? 

      printf("D.m_stopCriteria = %s\n", 
	     DecompAlgoStopStr[D.m_stopCriteria].c_str());
      
      //who deletes this memory? better to pass in newCuts.. 
      printf("\n\n====================================OUT DECOMP\n\n");
      //exit(1);
   }
#endif


  

   m_stats.thisGenCuts.push_back(m_stats.timerOther1.getRealTime()); 

   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "generateCuts()", m_param.LogDebugLevel, 2);
   return static_cast<int>(newCuts.size());
}








//------------------------------------------------------------------------- //
//member of varpool versus algo class? different for DC??
void DecompAlgo::addVarsToPool(DecompVarList & newVars){

   //TODO: how to get rid of duplicate vars - use some kind of hash map!!
   //what we have now is really expensive! will see in profile
   double           * denseCol  = NULL;
   CoinPackedVector * sparseCol = NULL;

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "addVarsToPool()", m_param.LogDebugLevel, 2);

   DecompConstraintSet          * modelCore   = m_modelCore.getModel();
   //---
   //--- sanity check - make sure the number of rows in core is
   //---    num of (orig+branch+cuts) in LP formulation
   //---
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 5,
	      (*m_osLog) 
	      << "num original= " << getNumRowType(DecompRow_Original) << "\n"
	      << "num branch  = " << getNumRowType(DecompRow_Branch)   << "\n"
	      << "num cut     = " << getNumRowType(DecompRow_Cut)      << "\n"
	      << "num convex  = " << getNumRowType(DecompRow_Convex)   << "\n"
	      << "num core    = " << modelCore->getNumRows()         << "\n";
	      );

   if(m_algo != DECOMP){
      assert((getNumRowType(DecompRow_Original) +
	      getNumRowType(DecompRow_Branch)   +
	      getNumRowType(DecompRow_Cut)) == modelCore->getNumRows());
      denseCol = new double[modelCore->getNumRows() + m_numConvexCon];
   }

   int blockIndex;
   DecompVarList::iterator li;//const?
   for(li = newVars.begin(); li != newVars.end(); li++){
      //--- 
      //--- get dense column = A''s, append convexity constraint on end 
      //--- THINK: PC specific
      //---
      //M contains A'', OSI does not - that contains lambda space
      //but in decomp_and_cut, we don't need modelCore->M
      
      //the column that eventually gets added to masterSI, masterSI
      //has A'' (original) + convexity + cuts, in that order
      //so the column must match that order
      //but modelCore->M is only A'' + cuts, no convexity constraint

      //so need something different for RC vs PC vs ... ?

      //TODO - fix this derive method for decomp
      
      if(m_algo == DECOMP){//UGH - derive own!
	 blockIndex = (*li)->getBlockId();
	 sparseCol  = new CoinPackedVector((*li)->m_s);

	 UTIL_DEBUG(m_app->m_param.LogDebugLevel, 5,
		    (*m_osLog) << "\nPRINT m_s\n";
		    UtilPrintPackedVector((*li)->m_s);         
		    );
	 //add in convexity constraint
	 sparseCol->insert(modelCore->getNumCols() + blockIndex, 1.0);

	 UTIL_DEBUG(m_app->m_param.LogDebugLevel, 5,
		    (*m_osLog) << "\nPRINT sparseCol\n";                    
		    UtilPrintPackedVector(*sparseCol);
		    );
      }
      else{
	 //modelCore->M = A'', plus possible cuts, but without convexity
	 //this is WRONG as soon as a cut is added!

	 //this is slow to use sparse m_s because timesMinor needs to do
	 //an expensive lookup x[index[j]]
	 //which is faster, make dense, then this or ...

	 //modelCore holds the original plus cuts - right? yes

	 //---
	 //--- this creates a dense array of the column (in x-space)
	 //---  it is stored as a sparse vector, so we have translate
	 //---  it to do the matrix multiply
	 //---
	 //--- we put the dense array into the mem pool at dblArrNCoreCols
	 //---
	 //printf("addVarsToPool - the var:\n");
	 //(*li)->print();
	 (*li)->fillDenseArr(modelCore->getNumCols(),
			     m_memPool.dblArrNCoreCols);
	 //int i;
	 //printf("the var after making dense:\n");
	 //for(i = 0; i < modelCore->getNumCols();i++){
	 //   printf("denseCol x-space[row i=%5d]: %12.10f\n", 
	 //	   i, m_memPool.dblArrNCoreCols[i]);
	 //}
		

	 //---
	 //--- modelCore->M = A'' + branch-rows + new cuts
	 //---    here, we caculate M.s (where s is a var in x-space)
	 //---
	 //should denseCol be zeor'd out first?
	 modelCore->M->times(m_memPool.dblArrNCoreCols, denseCol);
	 //printf("after coreM.s, dense column of master before convexity =\n");		
	 //for(i = 0; i < modelCore->getNumRows(); i++){
	 //   printf("denseCol lam-space[row i=%5d]: %12.10f\n", 
	 //	   i, denseCol[i]);
	 //}
		

	 //---
	 //--- Let A'' = original rows and branching rows.
	 //---
	 //--- a dense column here gives the coefficients in the
	 //--- reformulation by [A'' + cuts] * column (in x-space)
	 //---
	 //--- dimensions:
	 //---      (m'' + ncuts x n) * (n x 1) -> (m'' + ncuts x 1)	 
	 //---
	 //--- but this is missing the convexity constraints, and these
	 //--- need to be in order as originally setup (after A'', but
	 //--- before cuts)
	 //---
	 //--- shift all the cuts over to make space for convexity cons
	 //---
	 //--- since the convexity constraints are just sum{} lambda = 1
	 //--- we know that there is exactly one entry 1.0 for the approriate
	 //--- convexity row (depends on block id)
	 //---
	 //--- r[0],   r[1],  ..., r[m''-1], 
	 //--- cut[0], cut[1], ... cut[ncuts-1]
	 //---    -->
	 //--- r[0],    r[1],    ..., r[m''-1], 
	 //--- conv[0], conv[1], ..., conv[b-1],
	 //--- cut[0],  cut[1],  ...  cut[ncuts-1]	 
	 //---
	 int r, b;

	 //number of rows in core (A''+cuts)
	 int mpp             = modelCore->getNumRows();
	 //number of rows in original core (before cuts: A'')
	 int convexity_index = modelCore->nBaseRows;
	 //in the master, the convexity constraints are put just 
	 // after A'' (before any cuts were added)
	 assert(m_masterRowType[convexity_index] == DecompRow_Convex);
	 assert(mpp - convexity_index == getNumRowType(DecompRow_Cut));
	 
	 //printf("mpp            =%d\n", mpp);
	 //printf("convexity_index=%d\n", convexity_index);
	 
	 //for(r = 0; r < mpp; r++){
	 // printf("BEFORE dense[%d] = %g\n", r, denseCol[r]);
	 //}

	 //---
	 //--- for each cut row, move it to right/down
	 //---    o=original, b=branch, x=convex, c=cut
	 //---                             0123456789012
	 //---    LP                     : oooobbbbxxccc
	 //---    Core (current denseCol): oooobbbbccc  
	 //--- make room in denseCol, 0 out, then fill in for block
	 //---    nRows=13, nCuts=3, nConv=2 
	 //---       r=12..10 <-- r=10..8
	 //---                             0123456789012
	 //---                           : oooobbbb..ccc
	 //---                           : oooobbbb00ccc
	 //---                           : oooobbbb10ccc
	 //---  
	 int nCuts = mpp - convexity_index;
	 int nRows = m_masterSI->getNumRows();
	 assert(nRows == mpp + m_numConvexCon);
	 
	 //printf("nCuts = %d\n", nCuts);
	 //printf("nRows = %d\n", nRows);
	 //printf("convexity_index = %d\n", convexity_index);
	 for(r = (nRows - 1); r >= (nRows-nCuts); r--){
	    //printf("shift setting denseCol[r=%d] to denseCol[%d] as %12.10f\n",
	    //   r, r-m_numConvexCon, denseCol[r-m_numConvexCon]); 
	    denseCol[r] = denseCol[r - m_numConvexCon];
	 }
	 for(b = 0; b < m_numConvexCon; b++){
	    //printf("setting convex entry %d to 0\n",
	    //   convexity_index+b);
	    denseCol[convexity_index + b] = 0.0;
	 }
	 denseCol[convexity_index + (*li)->getBlockId()] = 1.0;
	 //printf("setting convex entry %d for block %d to 1\n",
	 //convexity_index+(*li)->getBlockId(),
	 //	(*li)->getBlockId());
	 
	 //for(r = 0; r < mpp + m_numConvexCon; r++){
	 // printf("AFTER dense[%d] = %g\n", r, denseCol[r]);
	 //}


	 //---
	 //--- creat a sparse column from the dense column
	 //---
	 sparseCol 
	    = UtilPackedVectorFromDense(modelCore->getNumRows() +
					m_numConvexCon,
					denseCol, m_app->m_param.TolZero);
	 UTIL_DEBUG(m_app->m_param.LogDebugLevel, 5,
		    (*m_osLog) << "\nPRINT sparseCol\n";
		    UtilPrintPackedVector(*sparseCol);
		    );
		    
	 
      }//END: else(m_algo == DECOMP)

      DecompWaitingCol waitingCol(*li, sparseCol);
      //what if we have created a variable that is already in LP!
      //that should not happen the way we currently do things... if it 
      //does, it could flag that we have an error of some sort
    
      //TOOD: since DecompVarList does not have its own class...
      //this is ugly, fix this later... make a helper funciton of DecompVar?

      //TODO: this is very expensive - use hash like in cuts
#define DUAL_SMOOTHING
      if(m_varpool.isDuplicate(m_vars, waitingCol)){
#ifndef DUAL_SMOOTHING
	 UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
		    (*m_osLog) << "Duplicate variable, already in vars!!\n";
		    (*m_osLog) << "\nVAR POOL:\n";
		    m_varpool.print(m_osLog);
		    (*m_osLog) << "\nVARS:\n";
		    printVars(m_osLog);  
		    );
#else
	 UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
		    (*m_osLog) << "Duplicate variable, already in vars!!\n";
		    );
#endif

	 //THIS can happen if coming from dual rays... ????

#ifndef DUAL_SMOOTHING
	 assert(0);
#endif
         
	 waitingCol.deleteVar();
	 waitingCol.deleteCol();
	 //THINK
	 //FUGLY
	 if(m_algo != RELAX_AND_CUT){
	    m_nodeStats.varsThisCall--;
	    m_nodeStats.varsThisRound--;
	 }
	 // assert(0);//this can happen if multi-poly?
	 continue;
      }


      if(m_varpool.isDuplicate(waitingCol)){
	 UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
		    (*m_osLog) << "Duplicate variable, already in var pool\n";
		    );
	 waitingCol.deleteVar();
	 waitingCol.deleteCol();
	 //THINK
	 if(m_algo != RELAX_AND_CUT){
	    m_nodeStats.varsThisCall--;
	    m_nodeStats.varsThisRound--;
	 }
      }
      else{      
	 //didn't you already do this part when you added the var to new_vars?
	 //waitingCol.setReducedCost(pi, phase);
	 m_varpool.push_back(waitingCol);
      }
   } //END: for(li = newVars.begin(); li != newVars.end(); li++)

   UTIL_DELARR(denseCol);
   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "addVarsToPool()", m_param.LogDebugLevel, 2);
}

//------------------------------------------------------------------------- //
//THINK: why are these not member functions of varpool?
//why of algo? because different for different algos...
void DecompAlgo::addVarsFromPool(){

   //TODO: we have checked to make sure there are no dups added to pool
   //do we also need to check that no dup vars are added to LP? for that
   //we'd have to check across m_vars

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "addVarsFromPool()", m_param.LogDebugLevel, 2);
   //const int maxvars_toadd = m_app->m_param.maxvars_periter;  
   //int n_newcols = std::min<int>(m_varpool.size(), maxvars_toadd);

   DecompVarPool::iterator vi;
   DecompVarPool::iterator viLast;

   int n_newcols = static_cast<int>(m_varpool.size());
   if(n_newcols == 0){
      UtilPrintFuncEnd(m_osLog, m_classTag,
		       "addVarsFromPool()", m_param.LogDebugLevel, 2);
      return;
   }

   //---
   //--- sort the pool by increasing reduced cost
   //---
   //THINK: need a copy constructor for object here?
   partial_sort(m_varpool.begin(), 
		m_varpool.begin() + n_newcols,
		m_varpool.end(), 
		is_less_thanD());
  
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 10,
	      (*m_osLog) << "\nVAR POOL BEFORE:\n";
	      m_varpool.print(m_osLog);
	      (*m_osLog) << "\nVARS BEFORE:\n";
	      printVars(m_osLog);  
	      );
  
   //THINK: find_first would be good, but we need the index of the iterator
   //how do i get that?
  
   //---
   //--- never add anything with pos rc
   //---
   int index = 0;
   for(vi = m_varpool.begin(); vi != m_varpool.end(); vi++){
      if(m_algo != RELAX_AND_CUT){//THINK
	 if((*vi).getReducedCost() >  -0.0000001){//TODO
	    //if((*vi).getReducedCost() >  -m_app->m_param.dualTol){
	    break;
	 }
      }
      index++;
   }
   n_newcols = std::min<int>(n_newcols, index);
   
   //TODO
   /*if(n_newcols > 0)
     m_cutpool.setRowsAreValid(false);*/

   //see CoinBuild

   //---
   //--- 1.) build up the block of columns to be added to the master
   //---     create a block for speed, rather than one column at a time
   //--- 2.) copy the var pointers to the DecompModel var list
   //---
   double * clb = new double[n_newcols];
   double * cub = new double[n_newcols];
   double * obj = new double[n_newcols];
   const CoinPackedVectorBase ** colBlock =
      new const CoinPackedVectorBase*[n_newcols];  
   const vector<string> & colNamesM = m_masterSI->getColNames();
   vector<string>   colNames;
   bool             hasNames  = colNamesM.size() > 0 ? true : false;
   const int        colIndex0 = m_masterSI->getNumCols();
   if(hasNames){
      assert(colIndex0 == static_cast<int>(colNamesM.size()));
   }

   index = 0;
   for(vi = m_varpool.begin(); vi != m_varpool.end(); vi++){
      if(index >= n_newcols)
	 break;
      //THINK: bother using CoinPackedVector??
      const CoinPackedVector * col = (*vi).getColPtr();
      DecompVar              * var = (*vi).getVarPtr();
      
      assert(col);
      colBlock[index] = col;
      clb[index]      = (*vi).getLowerBound();
      cub[index]      = (*vi).getUpperBound();
      if(m_phase == PHASE_PRICE1){
	 obj[index] = 0.0;
      }
      else{
	 obj[index] = (*vi).getOrigCost();
      }
      

      
      int blockIndex = var->getBlockId();
      int colIndex   = colIndex0 + index;
      
      var->setColMasterIndex(colIndex);
      m_masterColType.push_back(DecompCol_Structural);
      
      //---
      //--- give the column a name
      //---
      if(hasNames){
	 //since the names use colIndex and we delete, we might
	 // reuse names if we are not careful... we should probably
	 // give each a unique index.... and keep a running count
	 string colName = "lam(c_" + UtilIntToStr(m_colIndexUnique)
	    + ",b_" + UtilIntToStr(blockIndex) + ")";
	 colNames.push_back(colName);
      }
      m_colIndexUnique++;
         
      appendVars(var);
      index++;
   }
   viLast = vi;

   m_masterSI->addCols(n_newcols, colBlock, clb, cub, obj); 
   if(hasNames){      
      m_masterSI->setColNames(colNames, 0, 
			      static_cast<int>(colNames.size()), 
			      colIndex0);
   }
   
   const int n_colsAfter  = m_masterSI->getNumCols();
   assert(colIndex0 + n_newcols == n_colsAfter);
    
   //---
   //--- 3.) delete the col memory and clear the var pointer from varpool
   //---     the column memory is no longer needed, it has been copied into 
   //---     the master object, the variable memory is still needed, its 
   //---     pointer is now in m_vars, and no longer is needed in varpool
   //---
   //THINK is this all neccessary? just to keep memory small? or 
   //doing this for some reason of efficiency?
   for(vi = m_varpool.begin(); vi != viLast; vi++){
      (*vi).deleteCol();    
      (*vi).clearVar(); //needed? dangling pointer if not
   }
   //TODO: is this slow for vector? if so, maybe list is still the way to go
   m_varpool.erase(m_varpool.begin(), viLast);

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 10,
	      (*m_osLog) << "\nVAR POOL AFTER:\n";
	      m_varpool.print(m_osLog);
	      (*m_osLog) << "\nVARS AFTER:\n";
	      printVars(m_osLog);  
	      );

   //---
   //--- free local memory
   //---
   UTIL_DELARR(colBlock);
   UTIL_DELARR(clb);
   UTIL_DELARR(cub);
   UTIL_DELARR(obj);
   
   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "addVarsFromPool()", m_param.LogDebugLevel, 2);
   
}


/*-------------------------------------------------------------------------*/
void DecompAlgo::addCutsToPool(const double  *  x,
                               DecompCutList & newCuts,
                               int           & m_cutsThisCall){

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "addCutsToPool()", m_param.LogDebugLevel, 2);

   //for RC, make sure no cuts we are about to add are already in modelCore
   //also check that we have no duplicate cuts being put in here

   //TODO: do something similiar to check for pos-rc vars

   DecompConstraintSet          * modelCore   = m_modelCore.getModel();
   int  r, cutIndex = 0;
   bool isViolated = false; 
   bool isDupCore;//also check relax?
   bool isDupPool;
   bool addCut;

   DecompCutPool::iterator ci;
   DecompCutList::iterator li = newCuts.begin();
   while(li != newCuts.end()){
      CoinPackedVector * row       = new CoinPackedVector();
        
      //---
      //--- create a row (in terms of original formulation, x), from a cut
      //---
      (*li)->expandCutToRow(row);

      //---
      //--- set the hash string (for quick duplicate checks)
      //---
      (*li)->setStringHash(row);
        
#if 0
      bool isOptViolated = false;
      for(i = 0; i < m_optPoint.size(); i++){
	 isOptViolated = (*li)->calcViolation(row, &m_optPoint[i][0]);
	 if(isOptViolated)
	    (*m_osLog) << "\n\nCUT VIOLATES OPT POINT";
	 (*li)->print();
	 assert(!isOptViolated);
      }
#endif
      //here we will use a hash table - or just STL map
      addCut    = true;   
      isDupCore = false;
      for(r = 0; r < modelCore->getNumRows(); r++){
	 //override isSame( )
	 //in one case you check hash if expanded
	 //in user case you check isSame directly
	 //this will become hash lookup code
	 if(modelCore->rowHash[r] == (*li)->getStrHash()){
	    //---
	    //--- This should not happen, however, it is possible
	    //--- due to roundoff error. Since x = sum{}lambda, 
	    //--- the masterLP might be feasible while an a.x might
	    //--- violate a row bound slightly. This is checked after
	    //--- the recomposition. But, we don't throw an error unless
	    //--- the error is significant. The cut generator might
	    //--- duplicate a cut, because it finds an inequality that
	    //--- does cut off the current point that matches a row/cut
	    //--- already in the LP.
	    //---
	    //--- Like the check in checkPointFeasible, we should check
	    //--- that this duplicated cut violates by only a small 
	    //--- percentage. If not, then it really is an error.
	    //---
            UTIL_MSG(m_app->m_param.LogDebugLevel, 3,
                     (*m_osLog) << "Cut is Duplicate with Core\n";
                     );
            UTIL_MSG(m_app->m_param.LogDebugLevel, 4,
                     (*li)->print();            
                     );
	    isDupCore = true;
	    break;
	 }
      }

      if(isDupCore)
	 addCut = false;
      else{    

	 //---
	 //--- is this cut already in pool 
	 //---
	 isDupPool = false;
	 for(ci = m_cutpool.begin(); ci != m_cutpool.end(); ci++){
	    if((*li)->getStrHash() == (*ci).getCutPtr()->getStrHash()){
               UTIL_MSG(m_app->m_param.LogDebugLevel, 4,
                        (*m_osLog) << "Cut is Duplicate with Pool\n";
                        (*li)->print();
                        );
	       isDupPool = true;
	       break;
	    }           
	 }

	 if(isDupPool)
	    addCut = false;
	 else{
	    isViolated = (*li)->calcViolation(row, x);//also sets it
	    if(!isViolated)
	       addCut = false;
	 }
      }
      
      if(addCut){
	 DecompWaitingRow waitingRow(*li, row);
	 m_cutpool.push_back(waitingRow);
	 li++;
      }
      else{
	 //---
	 //--- cut is not violated, do not put in cut pool, delete memory
	 //---      
	 UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
		    (*m_osLog) << "\nCUT " << cutIndex 
		    << " do not put in pool";
		    );                 
	 UTIL_DELPTR(*li); //need to do?
	 li = newCuts.erase(li); //does this call cut destructor?
	 m_cutsThisCall--;
	 //then don't increment li next iter?? think
      }
      cutIndex++;
   }

   assert(m_cutsThisCall >= 0);
   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "addCutsToPool()", m_param.LogDebugLevel, 2);
}



/*--------------------------------------------------------------------------*/
int DecompAlgo::addCutsFromPool(){

   //this is exactly the same for RC, except that we also have to add 
   //multipliers 

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "addCutsFromPool()", m_param.LogDebugLevel, 2);

   //TODO: do some work here to check for duplicate cuts (actually do that
   //in addCutsToPool) in RC, can add cuts that already did (no "core model"
   //)

   //TODO: do some work here to check for parallel cuts
   DecompConstraintSet          * modelCore   = m_modelCore.getModel();
   //TODO partial sort!
   sort(m_cutpool.begin(),        
	m_cutpool.end(), 
	is_greater_thanD());
  
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 10,
	      (*m_osLog) << "\nCUT POOL BEFORE:\n";
	      m_cutpool.print(m_osLog);
	      (*m_osLog) << "\nCUTS BEFORE:\n";
	      printCuts(m_osLog);  
	      );

   const int maxcuts_toadd = 10000;//m_app->m_param.cut_maxcuts_periter;
  
   int n_newrows = CoinMin(static_cast<int>(m_cutpool.size()), maxcuts_toadd);
   //since we use a list - find_first won't help as it returns an 
   //iterator not an index in the list... UGH
   int index = 0;
   DecompCutPool::iterator li;
   for(li = m_cutpool.begin(); li != m_cutpool.end(); li++){
      if((*li).getViolation() < DecompEpsilon){//PARM
	 break;
      }
      index++;
   }
   //never add anything not violated
   n_newrows = std::min<int>(n_newrows, index);

   //TODO: look into coin build...
   double * rlb = new double[n_newrows];
   double * rub = new double[n_newrows];
   const CoinPackedVectorBase ** rowBlock   =
      new const CoinPackedVectorBase*[n_newrows];
  
   //vector<string> & coreRowNames   = modelCore->getRowNames();
   vector<string>   colNames;
   vector<string>   rowNames;
   string           colName;
   string           rowName;
   int              rowIndex;
   //int              colIndex;
   index = 0;
   for(li = m_cutpool.begin(); li != m_cutpool.end(); li++){
      if(index >= n_newrows)
	 break;
      
      CoinPackedVector * row = (*li).getRowPtr();
      DecompCut        * cut = (*li).getCutPtr();
      rlb[index]      = (*li).getLowerBound();
      rub[index]      = (*li).getUpperBound();
      rowBlock[index] = row;
      rowIndex        = m_masterSI->getNumRows() + index;
      //TODO: allow user to give cut names?
      rowName         = "cut(" + UtilIntToStr(rowIndex) + ")";
      rowNames.push_back(rowName);      

      //---
      //--- add the cut ptr to the list of cuts in masterLP
      //---
      m_cuts.push_back(cut);
            
      //---
      //--- set hash for cut
      //---
      modelCore->rowHash.push_back(cut->getStrHash());
      
      index++;
   }

   if((m_algo == RELAX_AND_CUT)) {
      //this is what we want to update in RC, but not in C
      modelCore->M->appendRows(n_newrows, rowBlock);      
      //create a modelCore->appendRows does all of this direct from 
      //a cut pool    
      //THINK: francois idea, when add cuts to P'?
      char   sense;
      double rhs, range;
      for(index = 0; index < n_newrows; index++){
	 modelCore->rowLB.push_back(rlb[index]);
	 modelCore->rowUB.push_back(rub[index]);
	 UtilBoundToSense(rlb[index], rub[index], DecompInf, 
			  sense, rhs, range);
	 modelCore->rowRhs.push_back(rhs);
	 modelCore->rowSense.push_back(sense);     
      }
   }
   else{
      //---
      //--- add the new rows to master
      //--- add the new rows to core (?) - must in PC, don't really need here... 
      //---
      m_masterSI->addRows(n_newrows, rowBlock, rlb, rub);
      int nRowNames = static_cast<int>(rowNames.size());
      if(nRowNames > 0)
	 m_masterSI->setRowNames(rowNames, 0, nRowNames, 0);
      
      //---
      //--- add to master row types
      //--- add names to modelCore as well (?)
      //---
      int i;
      for(i = 0; i < n_newrows; i++){
	 m_masterRowType.push_back(DecompRow_Cut);
	 //coreRowNames.push_back(rowNames[i]);
      }
   }

   
   
   
   //if(m_isTightenAlgo)
   //   m_masterSI->addRows(n_newrows, rowBlock, rlb, rub);
   //any reason to update this copy?? THINK
   //cuts on cuts... 
   
   
   //clean up
   index = 0;
   for(li = m_cutpool.begin(); li != m_cutpool.end(); li++){
      if(index >= n_newrows)
	 break;
      (*li).deleteRow();
      (*li).clearCut();//need to do this?
      index++;
   }
   m_cutpool.erase(m_cutpool.begin(), li);
  
   //UTIL_DELARR(rowReformBlock);
   UTIL_DELARR(rowBlock);
   UTIL_DELARR(rlb);
   UTIL_DELARR(rub);

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 10,
	      (*m_osLog) << "\nCUT POOL AFTER:\n";
	      m_cutpool.print(m_osLog);
	      (*m_osLog) << "\nCUTS AFTER:\n";
	      printCuts(m_osLog);  //??? add to cuts?
	      (*m_osLog) << "n_newrows = " << n_newrows << "\n";
	      );

   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "addCutsFromPool()", m_param.LogDebugLevel, 2);

   return n_newrows;
}

//------------------------------------------------------------------------- //
bool DecompAlgo::isIPFeasible(const double * x,
                              const double   feasVarTol,
                              const double   feasConTol,
                              const double   intTol){
   
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "isIPFeasible()", m_param.LogDebugLevel, 2);

   DecompConstraintSet * modelCore   = m_modelCore.getModel();
   const int             nInts       = modelCore->getNumInts();
   const int           * integerVars = modelCore->getIntegerVars();   
   const double          intTol10    = 10 * intTol;

   bool ipFeas = true;
   if(!isLPFeasible(x, feasVarTol, feasConTol)){
      ipFeas = false;
      goto FUNC_EXIT;
   }

   int    i, c;
   for(i = 0; i < nInts; i++){
      c = integerVars[i];
      if(!UtilIsIntegral(x[c], intTol)){
	 //Notify, but don't mark in feasible unless 10x worse.
	 UTIL_DEBUG(m_param.LogDebugLevel, 4,
		    (*m_osLog) << "IpFeas Integer Col x[" << c << "]: " 
                    << UtilDblToStr(x[c]) << endl;
		    );
	 
	 if(!UtilIsIntegral(x[c], intTol10)){
	    ipFeas = false;
	    goto FUNC_EXIT;
	 }
      }
   }
   
   
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 4,
	      m_app->printOriginalSolution(modelCore->getNumCols(), 
                                           modelCore->getColNames(),
                                           x);
	      );
   
 FUNC_EXIT:
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 4,
	      (*m_osLog) << "isIPFeasible = " << ipFeas << endl;
	      );
   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "isIPFeasible()", m_param.LogDebugLevel, 2);
   return ipFeas;
}

//------------------------------------------------------------------------- //
bool DecompAlgo::isLPFeasible(const double * x,
                              const double   feasVarTol,
                              const double   feasConTol){


   //---
   //--- The base isFeasible assumes a full explicit description
   //---  in modelCore and m_modelRelax (plus integrality). The user
   //---  app can also define an isFeasible - which will be checked
   //---  first.
   //---

   //TODO: consider different tolerances?
   //---
   //--- There are two cases where this is used which require different
   //---  tolerance levels. If this is to sanity check the recomposed 
   //---  solution x=sum{s}lamabda_s, then there could be alot of round-off
   //---  and we can probably allow a higher level of error. If this is 
   //---  used for checking if we found a new incumbent, there is less additive
   //---  roundoff since usually working with integer values, so tolerances
   //---  can be tighter. For example, in p0033, if we don't enforce a tight
   //---  satisfaction of constraint ax <= -1656, then ax=-1655 might sneak
   //---  by as it is only a 0.06% violation.
   //---

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "isLPFeasible()", m_param.LogDebugLevel, 2);
   
   bool lpFeas = m_modelCore.isPointFeasible(x, 
                                             m_param.LogDebugLevel,
                                             feasVarTol,
                                             feasConTol);
   if(!lpFeas)
      goto FUNC_EXIT;
   

   if(m_modelRelax.size()){
      map<int, DecompAlgoModel>::iterator mit;
      for(mit = m_modelRelax.begin(); mit != m_modelRelax.end(); mit++){
	 lpFeas = (*mit).second.isPointFeasible(x, 
						m_param.LogDebugLevel,
						feasVarTol,
						feasConTol);
         if(!lpFeas)
            goto FUNC_EXIT;
      }
   }
   
 FUNC_EXIT:
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 4,
	      (*m_osLog) << "isLPFeasible = " << lpFeas << endl;
	      );
   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "isLPFeasible()", m_param.LogDebugLevel, 2);
   return lpFeas;
}

//--------------------------------------------------------------------- //
//TODO: also need to accept fixed variables (for consistency
//and for branching), its possible this function can return more
//than one var
#ifdef RELAXED_THREADED
static void * solveRelaxedThread(void * args){
   SolveRelaxedThreadArgs * targs 
      = static_cast<SolveRelaxedThreadArgs*>(args);

   DecompAlgo          * algo         = targs->algo;
   int                   threadId     = targs->threadId;
   int                   batchSz      = targs->batchSz;
   int                 * batch        = targs->batch;   
   int                   nBaseCoreRows= targs->nBaseCoreRows;   
   double              * u            = targs->u;
   double              * redCostX     = targs->redCostX; //read-only
   const double              * origObjective= targs->origCost; //read-only
   int                   nCoreCols    = targs->n_origCols;  
   list<DecompVar*>    * vars         = targs->vars;
   //vector<DecompAlgoModel*> * algoModelV = targs->algoModel;
   
   
   //---
   //--- For pricing,
   //--- redCostX: is the red-cost for each original column  (c - uhat A")_e
   //--- origCost: is the original cost for each original column c_e
   //--- alpha:    is the dual for the convexity constraint 
   //--- 
   //--- The reduced cost of a new variable (column) is the sum of the 
   //--- reduced cost on each of the original columns in the new variable
   //--- minus alpha (this function is responsible for returning the reduced
   //--- cost, which includes alpha).
   //---
   //--- NOTE, redCost does not include alpha as sent in
   //---

   DecompApp * app = algo->getDecompAppMutable();

   int i, b;
   for(i = 0; i < batchSz; i++){
      b = batch[i];
      printf("THREAD %d solving b %d\n", threadId, b); fflush(stdout);

      //is this thread safe??
      DecompAlgoModel & algoModel = algo->getAlgoModelRelax(b);

      int    isExact  = 0;
      int    nVars    = vars->size();
      int    nNewVars = 0;
      double             alpha        = u[nBaseCoreRows + b];
      DecompSolverStatus solverStatus = DecompSolStatNoSolution;
      
      bool useExact  = true;
      bool useCutoff = false;
      DecompSolverResult solveResult(nCoreCols);

      //too many cases where ned locks - make threaded version
      algo->solveRelaxed(redCostX,
			 origObjective,
			 alpha,
			 nCoreCols,
			 true,
			 true,
			 useExact,
			 useCutoff,
			 false,
			 algoModel,
			 &solveResult,
			 *vars);
   }
   return NULL;
}
#endif

//--------------------------------------------------------------------- //
DecompStatus DecompAlgo::solveRelaxed(const double        * redCostX,
                                      const double        * origCost,
                                      const double          alpha,
                                      const int             n_origCols,
				      const bool            isNested,
                                      DecompAlgoModel     & algoModel,
                                      DecompSolverResult  * solveResult,
                                      list<DecompVar*>    & vars){
   
   //---
   //--- For pricing,
   //--- redCostX: is the red-cost for each original column  (c - uhat A")_e
   //--- origCost: is the original cost for each original column c_e
   //--- alpha:    is the dual for the convexity constraint 
   //--- 
   //--- The reduced cost of a new variable (column) is the sum of the 
   //--- reduced cost on each of the original columns in the new variable
   //--- minus alpha (this function is responsible for returning the reduced
   //--- cost, which includes alpha).
   //---
   //--- NOTE, redCost does not include alpha as sent in
   //---

#ifndef RELAXED_THREADED
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "solveRelaxed()", m_param.LogDebugLevel, 2);
#endif

   OsiSolverInterface  * subprobSI  = algoModel.getOsi();
   int                   whichBlock = algoModel.getBlockId();
   bool                  isRoot     = getNodeIndex() ? false : true;
   DecompConstraintSet * modelCore  = getModelCore().getModel();

#ifndef RELAXED_THREADED
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
	      (*m_osLog) << "solve block b = " << whichBlock << endl;
	      (*m_osLog) << "alpha         = " << alpha      << endl;
	      (*m_osLog) << "isNested      = " << isNested   << endl;
	      );
   m_stats.timerOther1.reset();
#endif

   DecompStatus rc     = STAT_UNKNOWN;

#ifndef RELAXED_THREADED
   m_stats.timerOther2.reset();
#endif

   int isExact  = 0;
   int nVars    = static_cast<int>(vars.size());
   int nNewVars = 0;

   //SolveRelaxAsIp
   //0 = If a user function is defined, it will use the user function.
   //    If the user returns an exact solution, it will not run the built-in
   //    IP solve (default).
   //    If a user function is not defined, it will use the built-in IP solve.
   //1 = Use the built-in IP solve, even if there is a user defines a function.
   //2 = Calls the user defined function (if exists) and then calls built-in 
   //    IP solver (use this for debugging).
   

   bool doCutoff = m_param.SubProbUseCutoff;
   bool doExact  = isNested ? false : true;

   //---
   //--- run user method, but disregard so we can compare
   //--- against built-in oracle
   //---
   DecompSolverStatus solverStatus = DecompSolStatNoSolution;
#ifndef RELAXED_THREADED
   if(m_param.SolveRelaxAsIp == 2){
      list<DecompVar*> varsDebug;
      if(isNested)
         solverStatus 
            = m_app->solveRelaxedNest(whichBlock, redCostX, alpha, varsDebug);
      else
         solverStatus 
            = m_app->solveRelaxed(whichBlock, redCostX, alpha, varsDebug);
      solverStatus = DecompSolStatNoSolution;
   }
   
   if(m_param.SolveRelaxAsIp != 1){
      if(isNested)
         solverStatus 
            = m_app->solveRelaxedNest(whichBlock, redCostX, alpha, vars);      
      else
         solverStatus 
            = m_app->solveRelaxed(whichBlock, redCostX, alpha, vars);
      m_stats.thisSolveRelaxApp.push_back(m_stats.timerOther2.getRealTime());
      nNewVars        = static_cast<int>(vars.size()) - nVars;
   }
   m_isColGenExact = (solverStatus == DecompSolStatOptimal);
   UTIL_DEBUG(m_param.LogDebugLevel, 4,
	      (*m_osLog) << "m_isColGenExact = " << m_isColGenExact << endl;
	      );
#endif



   if((!m_isColGenExact && nNewVars <= 0) || (m_param.SolveRelaxAsIp == 2)){ 
      
      //---
      //--- Here, we are going to use the built-in IP solver
      //---  to solve the subproblem. In many cases, the solver
      //---  will get a "good column" quickly (within some gap).
      //---  That is, it is ok to return suboptimal columns that 
      //---  have negative reduced cost.       
      //---
      //--- To possible approaches to speed things up:
      //---   (1) return as soon as an UB < 0 is found
      //---   (2) return when gap is tight
      //---
      //--- However, to prove that the final DW LB is valid, we will need 
      //---  solve the pricing problem to optimaity at some point.
      //---
      
      assert(subprobSI);
      double * milpSolution = NULL;
      double   rcBestCol    = DecompInf;
      double   milpObj      = DecompInf;
      
      //---
      //--- reset the objective to reduced cost
      //---
      algoModel.setOsiObjCoeff(redCostX);
      
      //---
      //--- reset the col lbs/ubs to node bounds
      //---
      //--- for block angular case, the user must tell us 
      //---   the active columns
      //---      
      //--- CAREFUL: this overrides the user subproblem column
      //---  bounds - if for some reason they don't want that 
      //---  to match up with core, this might cause an issue
      //---
      algoModel.setActiveColBounds(m_colLBNode, m_colUBNode);
      
      //---
      //--- dump subproblem model .mps/.lp
      //---
      if(m_param.LogDumpModel > 1){
         if(isNested){
            printCurrentProblem(subprobSI,
                                "subProbN_" + algoModel.getModelName(),
                                m_nodeStats.nodeIndex,
                                m_nodeStats.cutCallsTotal,
                                m_nodeStats.priceCallsTotal,
                                whichBlock);
         }
         else{
            printCurrentProblem(subprobSI,
                                "subProb_" + algoModel.getModelName(),
                                m_nodeStats.nodeIndex,
                                m_nodeStats.cutCallsTotal,
                                m_nodeStats.priceCallsTotal,
                                whichBlock);
         }
      }
      
    
      //---
      //--- solve: min cx, s.t. A'x >= b', x in Z ([A,b] is in modelRelax.M)
      //---    

      //TODO: get best N feasible solutions? is this possible with CBC
      //  this can help if this is the primary form of col-gen     

      algoModel.solveOsiAsIp(solveResult,
                             m_param,
                             doExact,
                             doCutoff,
                             isRoot,
                             alpha-DecompEpsilon);      
      
      rcBestCol = solveResult->m_objLB - alpha; //for sake of bound
      if(solveResult->m_nSolutions)
         milpSolution = solveResult->m_solution;
      
      
#ifndef RELAXED_THREADED      
      //TODO:
      //     z_DW is a LB on z_IP
      //During DW we get... z*_DW + RC* <= z_DW <= z*_DW
      //  but if we have bounds on RC* lbRC <= RC* <= ubRC
      //then z*_DW + lbRC <= z*_DW + RC* <= z_DW
      // we can get a valid LB update just from LB of oracle.. 
      //we can choose to stop and branch any time we want - 
      // we sometimes wait for rc=0 but really we can just look at
      // gap between DW's lb and ub... 
      
      
      m_isColGenExact = solveResult->m_isOptimal;

      UTIL_DEBUG(m_param.LogDebugLevel, 4,
		 (*m_osLog) << "m_isColGenExact = " << m_isColGenExact << endl;
		 );
#endif
      
      // THINK: we really don't want to force the user to create vars
      // and check rc and obj, etc... but they might know how to be smart
      // and produce more than one, etc... THINK

      if(milpSolution){
         //---
         //--- create a DecompVar (shat) from the optimal solution      
         //---
         vector<int>    ind;
         vector<double> els;
         int c;
         double varRedCost  = 0.0; //stupid - == obj ??
         double varOrigCost = 0.0;
         for(c = 0; c < n_origCols; c++){
            if(!UtilIsZero(milpSolution[c], m_app->m_param.TolZero)){
               ind.push_back(c);
               els.push_back(milpSolution[c]);	    
               //the reduced cost of shat: (c-uA").s
               varRedCost  += redCostX[c] * milpSolution[c];	    
               //the original cost of shat: c.s
               varOrigCost += origCost[c] * milpSolution[c];
               UTIL_DEBUG(m_app->m_param.LogDebugLevel, 5,
                          (*m_osLog) << "c: " << c 
                          << " varOrigCost = " << varOrigCost
                          << " origCost = " << origCost[c]
                          << " solution = " << milpSolution[c] << endl;
                          );
            }
         }
         varRedCost -= alpha;//RC = c-uA''s - alpha    
         UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
                    (*m_osLog) << "alpha      = " << alpha << "\n";
                    (*m_osLog) << "varRedCost = " << varRedCost << "\n";
                    (*m_osLog) << "varOrigCost = " << varOrigCost << "\n";
                    );
	    
	    
         /*double bestVarFromApp = DecompInf;
           if(rc != STAT_UNKNOWN){
           if(nNewVars >= 1){
           const DecompVar * appVar = vars.back();
           UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
           (*m_osLog) << "VAR from application varRedCost = "
           << appVar->getReducedCost() << endl;
           );
           bestVarFromApp = appVar->getReducedCost();
           }
           UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
           (*m_osLog) << "VAR from CBC varRedCost         = "
           << varRedCost << endl;
           );
           //assert(varRedCost <= (bestVarFromApp + 1.0e-2));
           if(isExact){
           assert( UtilIsZero(varRedCost-bestVarFromApp, 1.0e-2) );
           }
           }*/
	    
	    
         DecompVar * var = new DecompVar(ind, els, varRedCost, varOrigCost);
         var->setBlockId(whichBlock);
         UTIL_DEBUG(m_app->m_param.LogDebugLevel, 5,
                    var->print(););
         vars.push_back(var);
      }
   } //END:   if((rc == STAT_UNKNOWN)  || (!isExact && nNewVars <= 0)){ 



   //---
   //--- sanity check - if user provides a full description of 
   //---  relaxed problem, make sure the variables from either app or IP
   //---  solver are valid
   //---
   if(vars.size() > 0) {
      //---
      //--- get a pointer to the relaxed model for this block
      //---   even if this check is for a nested model, it should
      //---   be feasible to base relaxed model for this block
      //---	
      double * xTemp = new double[n_origCols];
      assert(xTemp);
      list<DecompVar*>::iterator it;
      for(it = vars.begin(); it != vars.end(); it++){
	 int               whichBlock     = (*it)->getBlockId();
	 if(whichBlock != -1){
	    DecompAlgoModel & algoModelCheck = getModelRelax(whichBlock);
	    UTIL_DEBUG(m_param.LogDebugLevel, 4,
		       (*m_osLog) << "Check that var atisifes relax matrix "
		       << whichBlock << endl;
		       (*it)->print();
		       );
	    (*it)->fillDenseArr(n_origCols, xTemp);
	    //TODO: get rid of this function, use isPointFeasible
	    bool isRelaxFeas 
	       = checkPointFeasible(algoModelCheck.getModel(), xTemp);
	    assert(isRelaxFeas);
	 }
      }
      UTIL_DELARR(xTemp);
   }   
   
#ifndef RELAXED_THREADED   
   m_stats.thisSolveRelax.push_back(m_stats.timerOther1.getRealTime()); 
   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "solveRelaxed()", m_param.LogDebugLevel, 2);
#endif
   return STAT_UNKNOWN;
}


//===========================================================================//
void DecompAlgo::recomposeSolution(const double * solution,
				   double       * rsolution){


   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "recomposeSolution()", m_param.LogDebugLevel, 2);

   //---
   //--- refresh the output array (x-space)
   //---
   DecompConstraintSet          * modelCore   = m_modelCore.getModel();
   UtilFillN(rsolution, modelCore->getNumCols(), 0.0);

   int  j;
   bool isFeas = true;
   for(j = 0; j < m_masterSI->getNumCols(); j++){
      if((fabs(solution[j])   > DecompEpsilon) && isMasterColArtificial(j)){
	 isFeas = false;
	 break;
      }
   }
   if(m_param.LogDebugLevel >= 4){
      int   r;
      const vector<string> & colNames = m_masterSI->getColNames();
      const vector<string> & rowNames = m_masterSI->getRowNames();
      for(j = 0; j < m_masterSI->getNumCols(); j++){
	 if(fabs(solution[j]) > DecompEpsilon){
	    if(j < static_cast<int>(colNames.size()))
	       printf("MASTER %25s PRIM[%6d->%20s] = %12.10f\n",
		      DecompColTypeStr[m_masterColType[j]].c_str(), 
		      j, colNames[j].c_str(), solution[j]);
	    else
	       printf("MASTER %25s PRIM[%6d] = %12.10f\n",
		      DecompColTypeStr[m_masterColType[j]].c_str(), 
		      j, solution[j]);
	    if(isMasterColArtificial(j))
	       isFeas = false;
	 }
      }
      if(m_masterSI->getNumIntegers() == 0){
	 //---
	 //--- in the case of recompose after solutionUpdateAsIp
	 //---  we cannot access valid duals
	 //---
	 const double         * dualSol  = m_masterSI->getRowPrice();
	 for(r = 0; r < m_masterSI->getNumRows(); r++){
	    if(fabs(dualSol[r]) > DecompEpsilon){
	       if(r < static_cast<int>(rowNames.size())){
		  printf("MASTER %25s DUAL[%6d->%20s] = %12.10f\n",
			 DecompRowTypeStr[m_masterRowType[r]].c_str(),
			 r, rowNames[r].c_str(), dualSol[r]);
	       }
	       else
		  printf("MASTER %25s DUAL[%6d] = %12.10f\n",
			 DecompRowTypeStr[m_masterRowType[r]].c_str(),
			 r, dualSol[r]);
	    }
	 }
      }
   }
      
   double lamSol;
   int    i, colIndex = 0;     
   const vector<string> & colNames  = modelCore->getColNames();
   int                    nColNames = static_cast<int>(colNames.size());
   DecompVarList::const_iterator li;
   for(li = m_vars.begin(); li != m_vars.end(); li++){
      colIndex = (*li)->getColMasterIndex();
      lamSol   = solution[colIndex];
      assert(colIndex < m_masterSI->getNumCols());
      assert(isMasterColStructural(colIndex));
      if(lamSol > m_param.TolZero){
	 UTIL_DEBUG(m_param.LogDebugLevel, 4,
		    (*m_osLog) << "LAMBDA[" << colIndex << "]: " << lamSol;
		    if(nColNames)
		       (*li)->print(m_osLog, colNames);
		    else
		       (*li)->print(m_osLog);
		    );
	 CoinPackedVector & v = (*li)->m_s;
	 const int    * inds  = v.getIndices();
	 const double * els   = v.getElements();
	 for(i = 0; i < v.getNumElements(); i++){
	    rsolution[inds[i]] += els[i] * lamSol;
	 }
      }
   }

   UTIL_MSG(m_param.LogDebugLevel, 4,
	    const double * cLB = modelCore->getColLB();
	    const double * cUB = modelCore->getColUB();
	    for(i = 0; i < modelCore->getNumCols(); i++){
	       if(!UtilIsZero(fabs(rsolution[i]))){
		  (*m_osLog) << "x[ " << setw(5) << i << " -> ";
		  if(nColNames)
		     (*m_osLog) << setw(25) << colNames[i];
		  (*m_osLog) << " ] = "  << UtilDblToStr(rsolution[i],6) 
			     << " LB = " << UtilDblToStr(cLB[i],6)
			     << " UB = " << UtilDblToStr(cUB[i],6)
			     << endl;
	       }
	    }
	    );
   
   //---
   //--- if any artificials are positive, then don't check against core
   //---
   if(isFeas){
      //TODO: get rid of this function, use isPointFeasible
      isFeas = checkPointFeasible(modelCore, rsolution);
      assert(isFeas);
   }

   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "recomposeSolution()", m_param.LogDebugLevel, 2);
}

//===========================================================================//
bool DecompAlgo::isTailoffLB(const int    changeLen,
                             const double changePerLimit){

   //---
   //--- check the change percentage in the last changeLen updates
   //---
   assert(changeLen >= 2);
   
   //---
   //--- don't check for tailoff until we have enough iterations
   //---
   if(static_cast<int>(m_nodeStats.objHistoryLB.size()) <= changeLen)
      return false;
   vector< DecompObjBound >::reverse_iterator it 
      = m_nodeStats.objHistoryLB.rbegin();
   int    len       = 0;
   double prevBound = (*it).thisBound;
   double diff      =  DecompInf;
   double sumDiff   = 0.0;
   double aveDiff   = 0.0;
   double perDiff   = 0.0;
   for( ; it != m_nodeStats.objHistoryLB.rend(); it++){
      diff       = fabs(prevBound - (*it).thisBound);
      UTIL_DEBUG(m_param.LogDebugLevel, 3,	      
		 (*m_osLog) 
		 << setw(10) << "prevBound="
		 << setw(10) << UtilDblToStr(prevBound, 2)
		 << setw(10) << ", thisBound="
		 << setw(10) << UtilDblToStr((*it).thisBound) << endl;
		 );
      sumDiff   += diff;
      prevBound  = (*it).thisBound; 
      len++;
      if(len >= changeLen)
	 break;
   }
   aveDiff = sumDiff / len;
   if(UtilIsZero(prevBound)){
      perDiff = aveDiff;
   }
   else{
      perDiff = 100 * aveDiff / fabs(prevBound);
   }
   UTIL_MSG(m_param.LogDebugLevel, 2,	      
	    (*m_osLog) 
	    << setw(10) << "Percentage difference in obj bound="
	    << setw(10) << UtilDblToStr(perDiff, 2) << endl;
	    );

   //---
   //--- if the average percentage difference is more than some threshold
   //---    than we are tailing off
   //---
   if(perDiff > changePerLimit)
      return false;
   else{
      //---
      //--- Even if we are tailing off, we need to be careful of the following:
      //---    If the last solution was integral (no branching candidates)
      //---    but we are not done pricing out (i.e., a column with negative 
      //---    RC still exist) and we declare that we are tailing off then the 
      //---    node will get put back in the node work queue. This can lead
      //---    to that node being repeatedly stopped and reseted. It is
      //---    better to just price it out since we cannot branch on it in 
      //---    this state.
      //---
      int    branchedOnIndex = -1;
      double branchedOnValue =  0;
      chooseBranchVar(branchedOnIndex, branchedOnValue);
      if(branchedOnIndex == -1)
	 return false;
      else return true;   
   }
}
