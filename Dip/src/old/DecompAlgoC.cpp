//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)       //
//                                                                           //
// Copyright (C) 2002-2007, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

#include "DecompApp.h"
#include "DecompVar.h"
#include "DecompAlgoC.h"
#include "DecompPortable.h"
// ------------------------------------------------------------------------- //
const char* DecompAlgoC::m_classTag         = "\nD-ALGOC       : ";

// ------------------------------------------------------------------------- //
void DecompAlgoC::createMasterProblem(DecompVarList& initVars)
{
   //---
   //--- Initialize the solver interface for the master problem.
   //--- C:  min c(x)
   //---     A' x   >= b'  [optional?]
   //---     A''x   >= b''
   //---     l <= x <= u
   //---
   //--- m_modelRelax contains [A',  b' ] in terms of x (if explicit)
   //--- m_modelCore  contains [A'', b''] in terms of x
   //---
   assert(initVars.size() == 0);
   //assume whichModel has been set, or pass it in for most flexible?
   //same question about m_app or pass in app for flexibility - THINK
   //DecompConstraintSet & modelCore  = m_app->modelCore[m_whichModel-1];
   //DecompConstraintSet & modelRelax  = m_app->modelRelax[m_whichModel-1];
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << " <---- createMasterProblem() ---- ";
             );
   //---
   //--- allocate memory here, pass in with assignProblem
   //---  (OSI will delete this memory)
   //---
   int nCols  = m_modelCore->M->getNumCols();
   int nRowsC = m_modelCore->M->getNumRows();
   int nRowsR = 0;
   int nRows  = nRowsC;

   if (m_modelRelax->M) {
      nRowsR  = m_modelRelax->M->getNumRows();
      nRows  += nRowsR;
   }

   double* colLB    = new double[nCols];
   double* colUB    = new double[nCols];
   double* objCoeff = new double[nCols];
   double* rowLB    = new double[nRows];
   double* rowUB    = new double[nRows];
   CoinPackedMatrix* M = new CoinPackedMatrix(*m_modelCore->M);

   if (m_modelRelax->M) {
      M->bottomAppendPackedMatrix(*m_modelRelax->M);
   }

   memcpy(colLB,   &m_modelCore->colLB[0],   nCols  * sizeof(double));
   memcpy(colUB,   &m_modelCore->colUB[0],   nCols  * sizeof(double));
   memcpy(objCoeff, m_app->m_model.objCoeff, nCols  * sizeof(double));
   memcpy(rowLB,   &m_modelCore->rowLB[0],   nRowsC * sizeof(double));
   memcpy(rowUB,   &m_modelCore->rowUB[0],   nRowsC * sizeof(double));

   if (m_modelRelax->M) {
      memcpy(rowLB + nRowsC, &m_modelRelax->rowLB[0], nRowsR * sizeof(double));
      memcpy(rowUB + nRowsC, &m_modelRelax->rowUB[0], nRowsR * sizeof(double));
   }

   //TODO: give user choice to put this in... also, what if not explicit
   //or, could use bottomAppendPackedMatrix
   /*
   if(m_modelRelax->M){
      m_masterSI->addRows(m_modelRelax->getNumRows(),
                          m_modelRelax->M->getVectorStarts(),
                          m_modelRelax->M->getIndices(),
                          m_modelRelax->M->getElements(),
                          &m_modelRelax->rowLB[0],
                          &m_modelRelax->rowUB[0]);
   }
   */
   /*
   m_masterSI->loadProblem(*m_modelCore->M,
                           &m_modelCore->colLB[0],
                           &m_modelCore->colUB[0],
                           m_app->m_model.objCoeff,
                           &m_modelCore->rowLB[0],
                           &m_modelCore->rowUB[0]);
   */
   m_masterSI->assignProblem(M, colLB, colUB, objCoeff, rowLB, rowUB);
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 10,
              m_masterSI->writeLp("master", "lp", 1.0e-6, 10, 2, 1.0, true);
             );
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << "  ---- createMasterProblem() ----> ";
             );
}

/*-------------------------------------------------------------------------*/
void DecompAlgoC::recomposeSolution(const double* solution,
                                    double*        rsolution)
{
   //weird name for this, just setting m_xhat here...
   //util function for memcpy??
   //more C++ version of this?
   //neccessary in Cut?
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << " <---- recomposeSolution()   ---- ";
             );
   //printf("\nnumcols: %d", m_modelCore->getNumCols());
   memcpy(rsolution, solution, m_modelCore->getNumCols() * sizeof(double));
   //for(int i = 0; i < m_modelCore->getNumCols(); i++){
   //  printf("\nrsol[%d]: %g, sol[%d]: %g, m_xhat[%d]: %g",
   //         i, rsolution[i], i, solution[i], i, m_xhat[i]);
   //}
   //m_app->printOriginalSolution(m_modelCore->getNumCols(),
   //                             m_xhat, 0);//m_iteration?
   //what is this for? if CPM?
#if 0
   //keep the vars during tighten
   vector<int>    ind;
   vector<double> els;
   double varRedCost  = 0.0;//need?
   double varOrigCost = 0.0;

   for (int c = 0; c < m_modelCore->getNumCols(); c++) {
      if (fabs(solution[c]) > m_app->m_param.TolZero) {
         ind.push_back(c);
         els.push_back(solution[c]);
         //what is this? m_model...
         varOrigCost += m_app->m_model.objCoeff[c] * solution[c];
      }
   }

   DecompVar* var = new DecompVar(ind, els, varRedCost, varOrigCost);
   m_vars.push_back(var);
#endif
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << "  ---- recomposeSolution()   ---->";
             );
}


// ------------------------------------------------------------------------ //
decompPhase DecompAlgoC::phaseUpdate(const decompPhase phase,
                                     const decompStat  stat)
{
   bool        isCutPossible, mustSwitch, considerSwitch;
   bool        ipFeasible, appFeasible;
   decompPhase nextPhase = PHASE_UNKNOWN;
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << " <---- phaseUpdate() ---- ";
             );
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << "\nm_cutsThisRound  : " << m_cutsThisRound;
              (*m_osLog) << "\nm_varsThisRound  : " << m_varsThisRound;
              (*m_osLog) << "\nm_cutsThisCall  : " << m_cutsThisCall;
              (*m_osLog) << "\nm_varsThisCall  : " << m_varsThisCall;
              (*m_osLog) << "\nm_cutCallsTotal  : " << m_cutCallsTotal;
              (*m_osLog) << "\nm_priceCallsTotal: " << m_priceCallsTotal;
              (*m_osLog) << "\nm_cutCallsRound  : " << m_cutCallsRound;
              (*m_osLog) << "\nm_priceCallsRound: " << m_priceCallsRound;
              (*m_osLog) << "\nPHASEIN          : " << decompPhaseStr[phase];
              (*m_osLog) << "\nSTATIN           : " << decompStatStr[stat];
             );
   //the user must(?) tell us?

   //app user feasible does not need to be virtual
   //the default checks against core and relaxed constraints - see SmallIP
   if (ipFeasible = isIPFeasible(m_xhat)) {
      printf("\nIP FEASIBLE, we are done? there might still be cuts");

      if (appFeasible = m_app->APPisUserFeasible(m_xhat,
                        m_modelCore->getNumCols(),
                        m_app->m_param.TolZero)) {
         printf("\nAPP FEASIBLE");
      }

      if (ipFeasible && !appFeasible) {
         nextPhase = PHASE_CUT;
         goto PHASE_UPDATE_FINISH;
      }
   }

   //what if we are integral but hit cut limit, but are still not userfeas?
   //we can't branch - no fractionals... so is feasible should actually
   //return cuts - it did  in old code... or just keep it simple,
   //if all integral but not feasible - then go to cuts whether the limit
   //has been hit or not... although feas check generates cuts itself
   //so this seems inefficient

   //---
   //--- we have exceeded the cut iter limit and the price iter limit
   //--- we are done
   //---
   if (m_cutCallsTotal   >= m_app->m_param.LimitTotalCutIters) {
      nextPhase = PHASE_DONE;
      goto PHASE_UPDATE_FINISH;
   }

   if (stat == STAT_INFEASIBLE) {
      //what if we are not pricing and it goes infeasible - that should STOP
      //test this with an infeasible milp example
      nextPhase = PHASE_DONE;
      goto PHASE_UPDATE_FINISH;
   }

   //what if user puts in total=1, round=0?
   isCutPossible   = (m_cutCallsTotal   < m_app->m_param.LimitTotalCutIters);
   //TODO: what if goes infeasible in middle?

   switch (phase) {
   case PHASE_INIT: {
      if (isCutPossible) {
         nextPhase = PHASE_CUT;
      } else {
         nextPhase = PHASE_DONE;
      }
   }
   break;
   case PHASE_CUT: {
      mustSwitch     = false;
      considerSwitch = false;

      if (!isCutPossible || (m_cutsThisCall == 0) || (m_cutsThisRound == 0)) {
         mustSwitch = true;
      }

      if (m_cutCallsRound >= m_app->m_param.LimitRoundCutIters) {
         considerSwitch = true;
      }

      if (mustSwitch) {
         //---
         //--- we must switch from cutting
         //---
         nextPhase = PHASE_DONE;
      }//END: if(mustSwitch)
      else if (considerSwitch) {
         //---
         //--- we consider switching from cutting
         //---
         if (!isCutPossible) {
            //---
            //--- if we exceed both iter limits, we are done
            //---
            nextPhase = PHASE_DONE;
         } else {
            //---
            //--- if we exceed the price iter limit, but not the cut limit
            //--- since we are not in mustSwitch, m_cutsThisRound > 0, so
            //--- we can go back to cutting, even though it violates the
            //--- round counter, because we have no other choice
            //---
            nextPhase = PHASE_CUT;
         }
      } //END: else if(considerSwitch)
      else {
         nextPhase = PHASE_CUT;
      }
   }
   break;
   case PHASE_DONE:
      break;
   case PHASE_UNKNOWN:
   default:
      assert(0);//TODO: real error codes?
   }

PHASE_UPDATE_FINISH:
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << "\nPHASEOUT         :"
              << decompPhaseStr[nextPhase] << "\t";
             );
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << "  ---- phaseUpdate() ----> ";
             );
   return nextPhase;
}

//don't need branch - as this will be in alps, but need a way to adjust
//the information needed by each algo given some NodeDesc...


// ------------------------------------------------------------------------ //
int DecompAlgoC::branch(int      branchedOnIndex,
                        double   branchedOnValue)
{
   //int    & statusDown,
   //int    & statusUp){
   //create two childern and solve each - just to get idea of
   //how column bound changes will effect things and be translated
   //need a nodeDesc?
   //need to store global UB, global LB, etc... let this happen at ALPs
   //layer - create an interface routine to handle all this - in between
   //ALPs and DECOMP
   vector<double> & colLB = m_modelCore->colLB;
   vector<double> & colUB = m_modelCore->colUB;
   double oldLB = colLB[branchedOnIndex];
   double oldUB = colUB[branchedOnIndex];
   //branch down
   printf("\n\n===== Branched Down -> resolve");
   m_masterSI->setColBounds(branchedOnIndex,
                            oldLB,
                            floor(branchedOnValue));
   //solve again
   processNode();
   //set back to orig
   //m_masterSI->setColBounds(branchedOnIndex,
   //                         colLB[branchedOnIndex],
   //                         colUB[branchedOnIndex]);
   //branch up
   printf("\n\n===== Branched Up -> resolve");
   m_masterSI->setColBounds(branchedOnIndex,
                            ceil(branchedOnValue),
                            oldUB);
   //solve again
   processNode();
   m_masterSI->setColBounds(branchedOnIndex,
                            oldLB, oldUB);
   return 0;
}
