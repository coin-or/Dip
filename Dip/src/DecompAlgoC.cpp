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
#include "DecompAlgoC.h"
#include "DecompAlgoD.h"
#include "DecompConstraintSet.h"

using namespace std;

//TODO: OsiDualObjLimit = gUB? if LB is higher, then can stop early

//===========================================================================//
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
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "createMasterProblem()", m_param.LogDebugLevel, 2);
   loadSIFromModel(m_masterSI);

   if (m_param.CutCGL) {
      m_cutgenSI = new OsiClpSolverInterface();
      CoinAssertHint(m_cutgenSI, "Error: Out of Memory");
      loadSIFromModel(m_cutgenSI, true);
   }

   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "createMasterProblem()", m_param.LogDebugLevel, 2);
}

//===========================================================================//
void DecompAlgoC::setMasterBounds(const double* lbs,
                                  const double* ubs)
{
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "setMasterBounds()", m_param.LogDebugLevel, 2);
   int c;
   const int n_cols = m_masterSI->getNumCols();
   //TODO: reuse this memory - col size does not change
   int*     index  = new int[n_cols];
   double* bounds = new double[2 * n_cols];

   for (c = 0; c < n_cols; c++) {
      index[c]      = c;
      bounds[2 * c]   = lbs[c];
      bounds[2 * c + 1] = ubs[c];
   }

   m_masterSI->setColSetBounds(index, index + n_cols, bounds);
   UTIL_DELARR(index);
   UTIL_DELARR(bounds);
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "setMasterBounds()", m_param.LogDebugLevel, 2);
}

//===========================================================================//
bool DecompAlgoC::updateObjBound(const double mostNegRC)
{
   //---
   //--- C    : LB = masterLP obj
   //--- PC   : LB = zDW_RMP + RC* <= zDW <= zDW_RMP
   //---    where RC* is the most negative reduced cost
   //---    assuming the relaxation subproblem was solved exactly
   //---
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "updateObjBoundLB()", m_param.LogDebugLevel, 2);
   double thisBoundLB = m_masterSI->getObjValue();
   setObjBound(thisBoundLB, thisBoundLB);
   UTIL_DEBUG(m_param.LogDebugLevel, 5,
              (*m_osLog)
              << "ThisLB = " << UtilDblToStr(thisBoundLB) << "\t"
              << "BestLB = " << UtilDblToStr(m_nodeStats.objBest.first)
              << "\n";
             );
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "updateObjBoundLB()", m_param.LogDebugLevel, 2);
   return false;
}

//===========================================================================//
void DecompAlgoC::recomposeSolution(const double* solution,
                                    double*        rsolution)
{
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "recomposeSolution()", m_param.LogDebugLevel, 2);
   DecompConstraintSet*           modelCore   = m_modelCore.getModel();
   memcpy(rsolution, solution, modelCore->getNumCols() * sizeof(double));
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "recomposeSolution()", m_param.LogDebugLevel, 2);
}

//===========================================================================//
void DecompAlgoC::phaseDone()
{
   //1 = every iter
   //2 = only at end of node
#if 0
   DecompConstraintSet*           modelCore   = m_modelCore.getModel();

   if (m_param.CutDC == 2) {
      DecompCutList newCuts;
      printf("\n\n==================================== IN DECOMP\n\n");
      DecompAlgoD D(m_app, m_utilParam,
                    m_xhat, modelCore->getNumCols());
      //also might want to use the columns you get here for something...
      //heur for ubs, etc..
      //either returns a set of cuts or a decomposition, have
      //that wrap solve()?
      D.solveD(&newCuts);
      //a hidden advantage of decomp in BC?
      DecompSolution* bestSol = NULL;
      vector<DecompSolution*>::const_iterator it;
      double thisBound;
      double bestBoundUB = m_nodeStats.objBest.second;
      const vector<DecompSolution*>& xhatIPFeasD = D.getXhatIPFeas();

      for (it  = xhatIPFeasD.begin();
            it != xhatIPFeasD.end(); it++) {
         thisBound = (*it)->getQuality();
         UTIL_DEBUG(m_param.LogDebugLevel, 3,
                    (*m_osLog) << "From DECOMP, IP Feasible with Quality =";
                    (*m_osLog) << thisBound << endl;
                   );

         if ((*it)->getQuality() <= bestBoundUB) {
            bestBoundUB = (*it)->getQuality();
            bestSol     = (*it);
         }
      }

      //need to make copy of solution, since D.m_xhatIpFeas goes out of scope
      if (bestSol) {
         DecompSolution* bestSolCp = new DecompSolution(*bestSol);
         m_xhatIPFeas.push_back(bestSolCp);
         setObjBoundUB(bestSolCp->getQuality());
         m_xhatIPBest = bestSolCp;
         //m_xhatIPBest->print();
      }

      //this could also very likely return a new gUB -
      //  for the case when it does find a decomposition
      //  and luckily it is feasible to original?
      //STOP -- 6/6/08
      //if decomp is found, then can't use currently - just looking for
      //farkas -- if decomp is found this means that z_LP = z_DW for that
      //relaxation??
      printf("D.m_stopCriteria = %s\n",
             DecompAlgoStopStr[D.getStopCriteria()].c_str());
      //who deletes this memory? better to pass in newCuts..
      printf("\n\n====================================OUT DECOMP\n\n");
      //exit(1);
   }

#endif
}

//===========================================================================//
void DecompAlgoC::phaseUpdate(DecompPhase&   phase,
                              DecompStatus& status)
{
   bool         isCutPossible, mustSwitch, considerSwitch;
   DecompPhase  nextPhase  = PHASE_UNKNOWN;
   DecompStatus nextStatus = status;
   pair<double, double>& objBest    = m_nodeStats.objBest;
   int cutCallsTotal                = m_nodeStats.cutCallsTotal;
   int cutCallsRound                = m_nodeStats.cutCallsRound;
   int cutsThisCall                 = m_nodeStats.cutsThisCall;
   int cutsThisRound                = m_nodeStats.cutsThisRound;
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "phaseUpdate()", m_param.LogDebugLevel, 2);
   m_phaseLast = phase;
   UTIL_MSG(m_app->m_param.LogDebugLevel, 3,
            (*m_osLog) << "cutsThisRound      : " << cutsThisRound   << "\n";
            (*m_osLog) << "cutsThisCall       : " << cutsThisCall    << "\n";
            (*m_osLog) << "cutCallsTotal      : " << cutCallsTotal   << "\n";
            (*m_osLog) << "cutCallsRound      : " << cutCallsRound   << "\n";
            (*m_osLog) << "LimitTotalCutIters : " << m_param.LimitTotalCutIters << "\n";
            (*m_osLog) << "LimitRoundCutIters : " << m_param.LimitRoundCutIters << "\n";
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
   //--- if the lower bound meets the global ub, we are done
   //---
   //if(objBest.first >= (objBest.second - DecompEpsilon)){
   // nextPhase = PHASE_DONE;
   // goto PHASE_UPDATE_FINISH;
   //}

   //TODO: check infeasible case

   //---
   //--- if no cuts, then jump to finish
   //---
   if ((m_param.LimitTotalCutIters == 0) ||
         (m_param.LimitRoundCutIters == 0)) {
      nextPhase = PHASE_DONE;
      UTIL_DEBUG(m_param.LogDebugLevel, 4,
                 (*m_osLog) << "Done - no cuts allowed." << endl;);
      goto PHASE_UPDATE_FINISH;
   }

   //---
   //--- we have exceeded the cut iter limit we are done
   //---
   if (cutCallsTotal >= m_param.LimitTotalCutIters) {
      nextPhase = PHASE_DONE;
      UTIL_DEBUG(m_param.LogDebugLevel, 4,
                 (*m_osLog) << "Done total cut calls exceeds limit." << endl;);
      goto PHASE_UPDATE_FINISH;
   }

   if (status == STAT_INFEASIBLE) {
      nextPhase = PHASE_DONE;
      UTIL_DEBUG(m_param.LogDebugLevel, 4,
                 (*m_osLog) << "Done status INFEASIBLE." << endl;);
      goto PHASE_UPDATE_FINISH;
   }

   isCutPossible = (cutCallsTotal < m_param.LimitTotalCutIters);

   switch (phase) {
      /*case PHASE_INIT:
      {
      if(isCutPossible)
      nextPhase = PHASE_CUT;
      else
      nextPhase = PHASE_DONE;
      }
      break;*/
   case PHASE_CUT: {
      mustSwitch     = false;
      considerSwitch = false;

      if ((cutCallsTotal > 0) &&
            (!isCutPossible || (cutsThisCall == 0) || (cutsThisRound == 0))) {
         mustSwitch = true;
      }

      if (cutCallsRound >= m_param.LimitRoundCutIters) {
         considerSwitch = true;
      }

      //printf("isCutPossible =%d\n", isCutPossible);
      //printf("mustSwitch    =%d\n", mustSwitch);
      //printf("considerSwitch=%d\n", considerSwitch);

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
      assert(0);
   }

PHASE_UPDATE_FINISH:
   UTIL_MSG(m_param.LogDebugLevel, 3,
            (*m_osLog) << "PhaseOut: "   << DecompPhaseStr[nextPhase];
            (*m_osLog) << " StatusOut: " << DecompStatusStr[nextStatus];
            (*m_osLog) << endl;
           );
   phase  = nextPhase;
   status = nextStatus;
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "phaseUpdate()", m_param.LogDebugLevel, 2);
}
