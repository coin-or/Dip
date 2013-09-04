//===========================================================================//
// This file is part of the DIP Solver Framework.                            //
//                                                                           //
// DIP is distributed under the Eclipse Public License as part of the        //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)       //
//                                                                           //
// Conceptual Design: Matthew Galati, SAS Institute Inc.                     //
//                    Ted Ralphs, Lehigh University                          //
//                                                                           //
// Copyright (C) 2002-2013, Lehigh University, Matthew Galati, Ted Ralphs    //
// All Rights Reserved.                                                      //
//===========================================================================//

//===========================================================================//
#include "DecompApp.h"
#include "DecompVar.h"
#include "DecompAlgoRC.h"

using namespace std;

//===========================================================================//
DecompPhase DecompAlgoRC::phaseInit()
{
   //THINK: should base have some container for master and sub
   //since in RC, don't need OSI??
   //so base should never call m_masterSI, but m_masterContainer
   //or something
   if (m_param.LogDumpModel > 1)
      printCurrentProblem(m_masterSI,
                          "masterProb",
                          m_nodeStats.nodeIndex,
                          m_nodeStats.cutCallsTotal,
                          m_nodeStats.priceCallsTotal);

   //---
   //--- update primal/dual vectors
   //---
   m_status = STAT_FEASIBLE;
   //TODO: what about the INF case?? artificial columns? DC-ABCC version
   // ---
   // --- update the phase
   // ---
   return PHASE_PRICE2;
}

//===========================================================================//
void DecompAlgoRC::phaseDone()
{
   //take the current set of variables and solve DW master to get primal
   //TODO: right now, creating from scratch each time -- really need
   //      to append and warm start - esp if doing alot of branching
   //delete m_masterSI;
   //m_masterSI = NULL;
   //m_masterSI = new OsiLpSolverInterface();
   //CoinAssertHint(m_masterSI, "Error: Out of Memory");
   //m_masterSI->messageHandler()->setLogLevel(m_param.LogLpLevel);
   DecompConstraintSet*           modelCore   = m_modelCore.getModel();
#if 0
   //---
   //--- Initialize the solver interface for the master problem.
   //--- PC: min c(s lam)
   //---     A''s   lam   >= b'',
   //---     sum{s} lam_s  = 1  ,
   //---            lam_s >= 0  , s in F'[0]
   //--- modelCore contains [A'', b''], from the original model in
   //--- terms of x. In this function we create the DW-LP in terms of
   //--- lambda, [[A''s, 1], [b'', 1]] and load that into the OSI
   //--- interface m_masterSI.
   //---
   //--- NOTE: if 0 is feasible to subproblem, we can relax convexity to <= 1
   //---
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "phaseDone()", m_param.LogDebugLevel, 2);
   CoinAssert(m_vars.size() > 0);
   int nColsCore = modelCore->getNumCols();
   int nRowsCore = modelCore->getNumRows();
   modelCore->nBaseRowsOrig = modelCore->nBaseRows;
   modelCore->nBaseRows     = nRowsCore;
   //THINK? should initVars be in pool and do addVarsFromPool here?
   CoinPackedMatrix* M = new CoinPackedMatrix(true, 0, 0);
   M->setDimensions(nRowsCore + 1, 0);
   const int n_cols  = static_cast<int>(m_vars.size());
   double* colLB    = new double[n_cols];
   double* colUB    = new double[n_cols];
   double* obj      = new double[n_cols];
   double* denseCol = new double[nRowsCore + 1];
   CoinAssertHint(colLB && colUB && obj && denseCol, "Error: Out of Memory");
   int col_index     = 0;
   DecompVarList::iterator li;

   for (li = m_vars.begin(); li != m_vars.end(); li++) {
      UTIL_DEBUG(m_param.LogDebugLevel, 5,
                 (*li)->print(m_osLog, m_app);
                );
      //---
      //--- get dense column = A''s, append convexity constraint on end
      //---
      modelCore->M->times((*li)->m_s, denseCol);
      denseCol[nRowsCore] = 1.0;
      //---
      //--- create a sparse column from the dense column
      //---
      // THINK: do i need a DecompCol?
      // THINK: does this allocate memory for coinpackedvec twice?
      CoinPackedVector* sparseCol
         = UtilPackedVectorFromDense(nRowsCore + 1,
                                     denseCol, m_param.TolZero);
      UTIL_DEBUG(m_param.LogDebugLevel, 5,
                 (*m_osLog) << "\nSparse Col: \n";
                 UtilPrintPackedVector(*sparseCol, m_osLog);
                );
      //TODO: do in const blocks
      // ---
      // --- append the sparse column to the matrix
      // ---
      M->appendCol(*sparseCol);
      colLB[col_index] = 0.0; //THINK: (*li)->getLowerBound();
      colUB[col_index] = DecompInf; //THINK: (*li)->getUpperBound(); //FIX!!
      obj[col_index]   = (*li)->getOriginalCost(); //c.s
      col_index++;
      UTIL_DELPTR(sparseCol); //THINK
   }

   //---
   //--- THINK: do we want to adjust modelCore directly here?
   //--- adjust row bounds for convexity constraint
   //---
   //TODO: in memory
   double* zeroSol    = new double[nColsCore];
   CoinAssertHint(zeroSol, "Error: Out of Memory");
   UtilFillN(zeroSol, nColsCore, 0.0);
   //TODO - REVISIT - that's not the right check
   //  needs to be feasible to subproblem?
   bool     isZeroFeas = isIPFeasible(zeroSol);
   UTIL_DEBUG(m_param.LogDebugLevel, 5,

              if (isZeroFeas)
              (*m_osLog) << "Zero Sol is Feasible - relax convexity con.\n";
             );

   vector<double> masterRowLB(modelCore->rowLB);
   vector<double> masterRowUB(modelCore->rowUB);

   if (isZeroFeas) {
      masterRowLB.push_back(-DecompInf);
      masterRowUB.push_back(1.0);
   } else {
      masterRowLB.push_back(1.0);
      masterRowUB.push_back(1.0);
   }

   //---
   //--- load the problem into master's solver interface
   //---
   assert(M->getNumRows() == static_cast<int>(masterRowLB.size()));
   assert(M->getNumRows() == static_cast<int>(masterRowUB.size()));
   m_masterSI->loadProblem(*M,  colLB, colUB, obj,
                           &masterRowLB[0],
                           &masterRowUB[0]);
   // ---
   // --- free local memory
   // ---
   UTIL_DELPTR(M);
   UTIL_DELARR(denseCol);
   UTIL_DELARR(colLB);
   UTIL_DELARR(colUB);
   UTIL_DELARR(obj);
   UTIL_DELARR(zeroSol);
#endif
   m_status = DecompAlgo::solutionUpdate(PHASE_UNKNOWN, 99999, 99999);

   //---
   //--- check if IP feasible (are we done?)
   //--- TODO: for nonexplicity, also check user app isfeasible
   //---
   //TODO: should this whole section be phaseDone?
   if (m_status != STAT_INFEASIBLE) {
      DecompAlgo::recomposeSolution(m_masterSI->getColSolution(), m_xhat);
      UTIL_DEBUG(m_param.LogDebugLevel, 4,
                 m_app->printOriginalSolution(modelCore->getNumCols(),
                                              modelCore->getColNames(),
                                              m_xhat);
                );

      if (isIPFeasible(m_xhat)) {
         if (m_app->APPisUserFeasible(m_xhat,
                                      modelCore->getNumCols(),
                                      m_param.TolZero)) {
            DecompSolution* decompSol
               = new DecompSolution(modelCore->getNumCols(),
                                    m_xhat, m_masterSI->getObjValue());
            m_xhatIPFeas.push_back(decompSol);
         }
      }

      vector<DecompSolution*>::iterator vi;
      DecompSolution* viBest = NULL;
      double bestBoundUB = m_nodeStats.objBest.second;

      for (vi = m_xhatIPFeas.begin(); vi != m_xhatIPFeas.end(); vi++) {
         const DecompSolution* xhatIPFeas = *vi;

         if (isIPFeasible(xhatIPFeas->getValues())) {
            if (xhatIPFeas->getQuality() <= bestBoundUB) {
               bestBoundUB = xhatIPFeas->getQuality();
               viBest = *vi;
            }
         }
      }

      if (viBest) {
         //save the best
         setObjBoundIP(bestBoundUB);
         m_xhatIPBest = viBest;
      }
   }

   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "phaseDone()", m_param.LogDebugLevel, 2);
}

//===========================================================================//
void DecompAlgoRC::recomposeSolution(const double* solution,
                                     double*        rsolution)
{
   printf("RC recomposeSolution does nothing\n");
}

//===========================================================================//
void DecompAlgoRC::createMasterProblem(DecompVarList& initVars)
{
   //---
   //--- there is no master LP in RC, just initialize the dual vector
   //---
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "createMasterProblem()", m_param.LogDebugLevel, 2);
   DecompAlgo::createMasterProblem(initVars);
   CoinAssert(initVars.size() > 0);
   //---
   //--- In order to implement simple branching, we are going to
   //--- treat all column bounds as explicit constraints. Then branching
   //--- for DW can be done in the same way it is done for regular CPM.
   //---
   //TODO: looks like volume doesn't like R rows... change and split up?
   //printf("Volume Algorithm can't work if there is a non ELG row\n");
   //   coreMatrixAppendColBounds();
   //THINK: is this the right place for this
   //TODO: give user option to feed in a good starting dual vector?
   //you might want to do that if switch between DW and RC... etc... THINK
   DecompConstraintSet*           modelCore   = m_modelCore.getModel();
   fill_n(back_inserter(m_u), modelCore->getNumRows(), 0.0);
   //TODO
   m_rc = new double[modelCore->getNumCols()]; //better in constructor?
   CoinAssertHint(m_rc, "Error: Out of Memory");
   //m_vars here will contain all the shat's used - we can only
   //use one of them in subgradient, so just save the last one
   //DecompVarList::iterator it = initVars.begin();
   //assert(*it);
   //think
   //double redCost = (*it)->getOriginalCost() - (*it)->m_s.dotProduct(&m_u[0]);
   //(*it)->setReducedCost(redCost);
   //for rc we want to calc the reduced cost
   //m_vars.push_back(*it);
   //UTIL_DEBUG(m_param.LogDebugLevel, 3,
   //	      (*it)->print(m_osLog);
   //	      );
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "createMasterProblem()", m_param.LogDebugLevel, 2);
}

// ------------------------------------------------------------------------- //
bool DecompAlgoRC::isDone()
{
   //iter count is checked by phaseUpdate
   //need to check step limit in here
   //need to check if ub-lb gap is small? isn't that always checked?

   //printf("\nm_UB: %12.10f, m_LB: %12.10f", m_UB, m_LB);
   if ((m_step < 1.0e-3) ||             //step length too small
         m_zeroSub         ||              //0 subgradient
         UtilIsZero(m_UB - m_LB, 1.0e-3)) { //gap is small
      return true;
   }

   return false;
}

//do we need a var pool at all here?

// ------------------------------------------------------------------------- //
int DecompAlgoRC::addCutsFromPool()
{
   int nNewRows = DecompAlgo::addCutsFromPool();
   m_u.reserve(m_u.size() + nNewRows);
   UtilFillN(m_u, nNewRows, 0.0);

   //is this the right place to do this?
   //we were pricing out, then cutting, but when we price out,
   //step=0, so we need to start over with step size
   //best place for this would be in phaseUpdatE?

   if (nNewRows > 0) {
      m_step = 2.0;
   }

   return nNewRows;
}

// ------------------------------------------------------------------------- //
int DecompAlgoRC::generateVars(const DecompStatus   stat,
                               DecompVarList&     newVars,
                               double&            mostNegReducedCost)
{
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "generateVars()", m_param.LogDebugLevel, 2);
   DecompConstraintSet*           modelCore   = m_modelCore.getModel();
   //really only returning one var here...
   //for RC, doesn't have to be negative??
   mostNegReducedCost = DecompInf;//bad name here
   //TODO: whenever a cut is added, if doing RC, you need to add an
   //element to u.... do we need to override gen cuts just for that?
   //gen cuts (x) is wrong here anyway... need gen cuts (s)
   assert(static_cast<int>(m_u.size()) == modelCore->getNumRows());
   //THINK:
   // if we overload getDual method, then this can be same as PC genvars?
   // in PC, we are letting OSI house the duals, in RC we do it ourself
   // seems silly to use OSI at all for RC?
   //
   // almost want Decomp to house dual even for PC...
   //TODO: m_app->m_model seems dumb
   //reduced cost = c - uA
   const double* origObjective = getOrigObjective();
   modelCore->M->transposeTimes(&m_u[0], m_rc);

   for (int c = 0; c < modelCore->getNumCols(); c++) {
      printf("RC[%d] -> c: %g - uA: %g = m_rc: %g\n",
             c, origObjective[c], m_rc[c],
             origObjective[c] - m_rc[c]);
      m_rc[c] = origObjective[c] - m_rc[c];
   }

   //double alpha = 0.0;
   DecompVarList potentialVars;
   //TODO: stat return, restrict how many? pass that in to user?
   //only take those with negative reduced cost?
   //check for dups here
   //TODO: blocks!
   //WRONG...
#if 0
   solveRelaxed(
      0,
      m_rc, origObjective, alpha,
      modelCore->getNumCols(), false, true, true, false,
      m_subprobSI[0],
      NULL,//WILL FAIL
      potentialVars);//NO CHECK RC??
#endif
   //another way to do this is to just collect all m_vars,
   //not worrying about duplicates -- and when we get to DW
   //strip out the dups before constructing the master formulation
   DecompVarList::iterator it;
   double varRedCost;

   for (it = potentialVars.begin(); it != potentialVars.end(); it++) {
      varRedCost = (*it)->getReducedCost();
      newVars.push_back(*it);

      if (varRedCost < mostNegReducedCost) {
         mostNegReducedCost = varRedCost;
         //TODO: if this winds up being a dup in addVarsToPool
         //this memory will be erased -- which will cause a problem
         //because we need this var
         //TODO: FUGLY - make a copy so if dup doesn't cause prob
         //i don't like this! ugh
         //m_shatVar = *it; //THINK
         m_shatVar = *(*it);
         m_shatVar.fillDenseArr(modelCore->getNumCols(), m_xhat);
      }
   }

   potentialVars.clear(); //THINK? what does clear do exactly ?
#if 0

   for (it = newVars.begin(); it != newVars.end(); it++)
      UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
                 (*it)->print(m_osLog);
                );

#endif
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "generateVars()", m_param.LogDebugLevel, 2);
   return static_cast<int>(newVars.size());
}

// ------------------------------------------------------------------------- //
DecompStatus DecompAlgoRC::solutionUpdate(const DecompPhase phase,
      const int         maxInnerIter,
      const int         maxOuterIter)
{
   //---
   //--- C, PC: This step solves (or takes a few steps to solve) master LP
   //---        which updates both the primal (x,lambda) and dual(u) vectors.
   //--- RC   : This does one step of subgradient, which updates the dual (u)
   //---        vector, given some shat (the last variable in m_var).
   //generalize this for the many variants of subgradient?
   //TODO: might make life easier to flip all inequalities to <= or >= !
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "solutionUpdate()", m_param.LogDebugLevel, 2);
   m_UB = 1.05 * m_nodeStats.objBest.second; //TODO heuristics 1.05??
   //TODO: tols
   //how to allow hooks to other stabilization methods, bundle, etc... user
   //can simply derive from RC or base and recode this method...
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << "\nVARS m_vars:\n";
              printVars(m_osLog);
             );
   DecompVar* shatVar = &m_shatVar; //m_vars.back();
   //DecompVar * shatVar = m_vars.back();
   CoinAssert(shatVar);
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << "\nshat: ";
              shatVar->print(m_osLog);
             );
   //make this part of class, else realloc every iter...
   //use vector, let it grow?
   int      r;
   //char     sense;
   //double   range;
   DecompConstraintSet*           modelCore   = m_modelCore.getModel();
   int      n_coreRows = modelCore->getNumRows();
   //const double * rhs        = getRightHandSide();
   //const char   * sense      = getRowSense();
   const double* rhs        = &modelCore->rowRhs[0];
   const char*    sense      = &modelCore->rowSense[0];
   double* violation  = new double[n_coreRows];
   double* activity   = new double[n_coreRows];
   modelCore->M->times(shatVar->m_s, activity); //As
   assert(static_cast<int>(m_u.size()) == n_coreRows);
   // =, b  - Ax or Ax - b
   //>=, b  - Ax > 0 is a violation
   //<=, b  - Ax < 0 is a violation
   m_zeroSub = true;

   for (r = 0; r < n_coreRows; r++) {
      violation[r] = rhs[r] - activity[r];
      UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
                 (*m_osLog) << setprecision(8);
                 (*m_osLog) << "r: " << r << " vio: " << violation[r]
                 << " rhs: " << rhs[r] << " act: " << activity[r]
                 << " u: " << m_u[r] << " sense: " << sense[r];
                );

      //beasley suggestion... this is just another basic variant of SG?
      // =, b - Ax != 0 is a violation
      //>=, b - Ax > 0  is a violation
      //<=, Ax - b > 0  is a violation
      switch (sense[r]) {
      case 'E':
         //violation[i] = rhs[i] - activity[i];
         break;

      case 'G':

         //violation[i] = rhs[i] - activity[i];
         //TODO: use tol
         if (violation[r] < 0.0 && m_u[r] >= -1.0e-4 && m_u[r] <= 1.0e-4) {
            violation[r] = 0.0;
         }

         break;

      case 'L':

         //violation[i] = rhs[i] - activity[i];
         if (violation[r] > 0.0 && m_u[r] >= -1.0e-4 && m_u[r] <= 1.0e-4) {
            violation[r] = 0.0;
         }

         break;
      }

      UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
                 (*m_osLog) << " -> vio: " << violation[r] << "\n";
                );

      //?? shouldn't it be if all are feasible?
      if (fabs(violation[r]) > 0.0001) {
         m_zeroSub = false;
      }
   }

   //when to half the step size?
   //same stuff as m_tlb?? setTrueLowerBound?
   double bound    = shatVar->getReducedCost();//c - uA (is this set?)
   //double constant = calcConstant(m_u.size(), &m_u[0]);
   double constant = 0.0;

   for (r = 0; r < n_coreRows; r++) {
      constant += m_u[r] * rhs[r];
   }

   //needs rhs from OSI - this is why it was better to fake it
   //and have the information in OSI... but then carrying around
   //modelCore and OSI for no good reason... but... this is messier

   //LR Bound = (c - uA)shat + ub, assumes u >= 0

   //first iter, LR Bound = cshat - is that a valid LB?
   //yes, actually is a LB for any u >= 0
   //this assumes u is optimal? or just dual feasible
   //is u = 0 dual feasible?

   //but only a valid bound if shat has the lowest reduced cost
   //for a given u... so, for first iter of smallip, should have given
   //(2,1) with LB = 2... RC don't do genInitVars?

   //TODO: think initial dual vector - solve an LP to get started?

   if (bound + constant > m_LB + m_app->m_param.TolZero) {
      m_LB  = bound + constant;
      m_cntSameLB  = 0;
      //count_sameLB = 0;
      //make param to do or not
      //reducedCostFixing(reducedCost);
   } else {
      m_cntSameLB++;
   }

   //if(count_sameLB >= m_app->m_param.RC_sameLBLimit){
   if (m_cntSameLB >= 10) {
      m_step /= 2.0;
      cout << "LB has not changed in " << m_cntSameLB
           << " iterations - halve step: " << m_step << endl;
      m_cntSameLB = 0;
   }

   printf("m_UB: %12.10f, m_LB: %12.10f\n", m_UB, m_LB);
   assert((m_UB - m_LB) > -0.0001);
   double theta = 0.0;
   double denom = 0.0;

   for (r = 0; r < n_coreRows; r++) {
      denom += violation[r] * violation[r];
   }

   if (denom > 0.0) {
      theta = m_step * ((1.00 * m_UB) - m_LB) / denom;
   }

   //TODO: debug util to print a list of values in a nice format to log?
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << "m_UB: " << m_UB << " m_LB: " << m_LB
              << " denom: " << denom << " m_step: " << m_step
              << " theta: " << theta << "\n";
             );

   //STOP 10/6/07
   //How do we deal with range constraints? What does volume do, for example?

   for (r = 0; r < n_coreRows; r++) {
      switch (sense[r]) {
      case 'E':
         m_u[r] += theta * violation[r];
         break;

      case 'G':
         //u > 0, g_i > 0 for violations
         m_u[r] = max(0.0, m_u[r] + (theta * violation[r]));
         break;

      case 'L':
         //u < 0, g_i < 0 for violatoins
         m_u[r] = max(0.0, m_u[r] - (theta * violation[r]));
         break;

      default:
         assert(0);
      }

      UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
                 (*m_osLog) << "r: " << r << " m_u: " << m_u[r] << "\n";
                );
   }

   m_iter++;
   //if no var pool buffer is used, make sure you clean up manually
   /*it++;
   UtilDeleteListPtr(potentialVars, it, potentialVars.end());
   potentialVars.clear();*/
   //TODO - back in algo, reduced cost fixing
   //TODO: temp memory don't alloc/free each time!
   //UTIL_DELARR(rhs);
   UTIL_DELARR(violation);
   UTIL_DELARR(activity);
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "solutionUpdate()", m_param.LogDebugLevel, 2);
   return STAT_FEASIBLE;
}

// ------------------------------------------------------------------------- //
bool DecompAlgoRC::updateObjBound(const double mostNegRC)
{
   //---
   //--- C    : LB = masterLP obj
   //--- PC   : LB = zDW_RMP + RC* <= zDW <= zDW_RMP
   //---    where RC* is the most negative reduced cost
   //---    assuming the relaxation subproblem was solved exactly
   //---
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "updateObjBound()", m_param.LogDebugLevel, 2);
   DecompConstraintSet*           modelCore   = m_modelCore.getModel();
   //mostNegRC not used?
   //DecompVar * shatVar = m_shatVar;//m_vars.back();
   //DecompVar * shatVar = m_vars.back();
   //CoinAssert(shatVar);
   int r;
   const int n_coreRows = modelCore->getNumRows();
   //double bound = shatVar->getReducedCost();//c - uA (is this set?)
   double constant = 0.0;
   const double* rhs = &modelCore->rowRhs[0];

   for (r = 0; r < n_coreRows; r++) {
      constant += m_u[r] * rhs[r];
   }

   //double thisBoundLB = shatVar->getReducedCost() + constant;
   double thisBoundLB = mostNegRC + constant;
   setObjBound(thisBoundLB, constant);
   UTIL_DEBUG(m_param.LogDebugLevel, 5,
              (*m_osLog)
              << "ThisLB = " << UtilDblToStr(thisBoundLB) << "\t"
              << "BestLB = " << UtilDblToStr(m_nodeStats.objBest.first)
              << "\n";
             );
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "updateObjBound()", m_param.LogDebugLevel, 2);
   return false;
}
