//===========================================================================//
// This file is part of the DIP Solver Framework.                            //
//                                                                           //
// DIP is distributed under the Eclipse Public License as part of the        //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Authors: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)      //
//          Ted Ralphs, Lehigh University (ted@lehigh.edu)                   //
//          Jiadong Wang, Lehigh University (jiw408@lehigh.edu)              //
//                                                                           //
// Copyright (C) 2002-2019, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//                                                                           //
// Interface to Gurobi is Copyright 2015 Jazz Aviation LP                    //
//===========================================================================//

//===========================================================================//
#include "DecompApp.h"
#include "DecompAlgoPC.h"
#include "DecompSolverResult.h"
#include "DecompConstraintSet.h"
//===========================================================================//
#include "CoinWarmStartBasis.hpp"

using namespace std;

//#define   DO_INTERIOR //also in DecompAlgo

//===========================================================================//
void DecompAlgoPC::phaseInit(DecompPhase& phase)
{
   //---
   //--- solve the LP of the compact formulation
   //---  and current branching decisions
   //---
   //--- the main goal here is determine if a node is LP infeasible
   //---  due to the branching decisions - if we can determine this
   //---  here, we can just fathom the node and skip the expensive
   //---  PhaseI pricing to prove infeasible
   //---
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "phaseInit()", m_param.LogDebugLevel, 2);

   //---
   //--- set column bounds
   //---
   //TODO: reuse this memory - col size does not change
   //TODO: not getting cuts since only populate A'' in beginning
   //  by adding cuts this adds some of the integrality property
   //  otherwise the only thing causing infeasibility is the branching
   //TODO: to get cuts would want to populate m_auxSI with
   //  the current core model - not too bad since only once
   //  a node? but even without cuts, see if it helps
   if (m_auxSI) {
      UTIL_DEBUG(m_param.LogDebugLevel, 4,
                 (*m_osLog) << "Solve the LP of compact formulation." << endl;
                );
      int      c;
      int      numCols = m_modelCore.getModel()->getNumCols();
      int*     index   = new int[numCols];
      double* bounds  = new double[2 * numCols];

      if (!(index && bounds)) {
         throw UtilExceptionMemory("phaseInit", m_classTag);
      }

      for (c = 0; c < numCols; c++) {
         index[c]      = c;
         bounds[2 * c]   = m_colLBNode[c];
         bounds[2 * c + 1] = m_colUBNode[c];
      }

      m_auxSI->setColSetBounds(index, index + numCols, bounds);
      UTIL_DELARR(index);
      UTIL_DELARR(bounds);
      //---
      //--- solve LP relaxation
      //---
      m_auxSI->resolve();
      //---
      //--- two possible results that might help here:
      //---   (1) the LP is found infeasible                         -> fathom
      //---   (2) the LP bound is already greater than the global UB -> fathom
      //---
      UTIL_DEBUG(m_param.LogDebugLevel, 4,
                 (*m_osLog)
                 << "Iteration Count               : "
                 << m_auxSI->getIterationCount() << "\n"
                 << "isAbandoned()                 : "
                 << m_auxSI->isAbandoned() << "\n"
                 << "isProvenOptimal()             : "
                 << m_auxSI->isProvenOptimal() << "\n"
                 << "isProvenPrimalInfeasible()    : "
                 << m_auxSI->isProvenPrimalInfeasible() << "\n"
                 << "isProvenDualInfeasible()      : "
                 << m_auxSI->isProvenDualInfeasible() << "\n"
                 << "isPrimalObjectiveLimitReached : "
                 << m_auxSI->isDualObjectiveLimitReached() << "\n"
                 << "isDualObjectiveLimitReached   : "
                 << m_auxSI->isDualObjectiveLimitReached() << "\n"
                 << "isIterationLimitReached       : "
                 << m_auxSI->isIterationLimitReached() << "\n";
                );

      if (m_auxSI->isProvenPrimalInfeasible()) {
         UTIL_MSG(m_param.LogLevel, 3,
                  (*m_osLog) << "LP of Compact found Infeasible." << endl;
                 );
         phase  = PHASE_DONE;
      }
   }

   if (phase != PHASE_DONE)
      if (getNodeIndex() == 0 && !m_isStrongBranch) {
         phase = PHASE_PRICE1;
      }

   UTIL_MSG(m_param.LogLevel, 2,
            (*m_osLog) << "phase = " << DecompPhaseStr[phase] << endl;
           );
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "phaseInit()", m_param.LogDebugLevel, 2);
}

//===========================================================================//
void DecompAlgoPC::adjustMasterDualSolution()
{
   if (!m_param.DualStab) {
      return;
   }

   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "adjustMasterDualSolution()", m_param.LogDebugLevel, 2);
   //---
   //--- resize dual vectors
   //---
   int nRows = static_cast<int>(m_masterSI->getNumRows());
   m_dual.resize(nRows);
   m_dualRM.resize(nRows);
   m_dualST.resize(nRows);
   //---
   //--- calculate smoothed dual
   //---    pi_ST = alpha * pi_Bar + (1-alpha) * pi_RM
   //--- this is dual feasible because it is taking
   //---    a convex combination of previously dual feasible vectors
   //--- need to be careful here, as the init dual is 0, which might not
   //---    be dual feasible, therefore, in the first iteration, we need
   //---    to skip the smoothing and enforce that the first dual be set
   //---    to dualRM
   //---
   int            r;
   const double* u      = &m_dualSolution[0];
   double         alpha  = m_param.DualStabAlpha;
   double         alpha1 = 1.0 - alpha;
   copy(u, u + nRows, m_dualRM.begin()); //copy for sake of debugging

   //---
   //--- for both the first PhaseI and first PhaseII calls,
   //---   be sure to set the dual vector to dualRM as dual=0
   //---   might not be feasible
   //---
   if (m_param.LogDebugLevel >= 3) {
      (*m_osLog) << "m_firstPhase2Call = " << m_firstPhase2Call << endl;
   }

   if (((m_nodeStats.cutCallsTotal +
         m_nodeStats.priceCallsTotal) == 0) || m_firstPhase2Call) {
      if (m_param.LogDebugLevel >= 2) {
         (*m_osLog) << "Init dual to dualRM" << endl;
      }

      copy(m_dualRM.begin(), m_dualRM.end(), m_dual.begin());
   }

   if (m_firstPhase2Call) {
      m_app->initDualVector(m_dual);
   }

   for (r = 0; r < nRows; r++) {
      m_dualST[r] = (alpha * m_dual[r]) + (alpha1 * m_dualRM[r]);
   }

   //---
   //--- log for debugging
   //---
   if (m_param.LogDebugLevel >= 3) {
      const vector<string>& rowNames = m_masterSI->getRowNames();

      for (r = 0; r < m_masterSI->getNumRows(); r++) {
         if (!(UtilIsZero(m_dual[r]) &&
               UtilIsZero(m_dualRM[r]) && UtilIsZero(m_dualST[r]))) {
            if (r < static_cast<int>(rowNames.size())) {
               (*m_osLog) << "MASTER "
                          << DecompRowTypeStr[m_masterRowType[r]]
                          << " DUAL[ " << r << "->" << rowNames[r]
                          << "] = " << m_dual[r] << " RM = "
                          << m_dualRM[r] << " ST = " << m_dualST[r]
                          << endl;
            } else
               (*m_osLog) << "MASTER "
                          << DecompRowTypeStr[m_masterRowType[r]]
                          << " DUAL[ " << r
                          << "] = " << m_dual[r] << " RM = "
                          << m_dualRM[r] << " ST = " << m_dualST[r]
                          << endl;
         }
      }
   }

   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "adjustMasterDualSolution()", m_param.LogDebugLevel, 2);
}


//===========================================================================//
int DecompAlgoPC::adjustColumnsEffCnt()
{
   int            status         = DecompStatOk;
   int            colMasterIndex = -1;
   const double* redCost        = m_masterSI->getReducedCost();
   double         redCostI       = 0.0;
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "adjustColumnsEffCnt()", m_param.LogDebugLevel, 2);
   DecompVarList::iterator li;

   for (li = m_vars.begin(); li != m_vars.end(); li++) {
      colMasterIndex = (*li)->getColMasterIndex();
      redCostI       = redCost[colMasterIndex];
      assert(isMasterColStructural(colMasterIndex));

      if (redCostI > DecompEpsilon) {
         (*li)->decreaseEffCnt();
      } else {
         (*li)->increaseEffCnt();
      }

      UTIL_DEBUG(m_param.LogLevel, 4,
                 (*m_osLog) << "ColIndex= " << setw(5) << colMasterIndex
                 << " RedCost= " << UtilDblToStr(redCostI)
                 << " EffCnt= " << setw(5) << (*li)->getEffectiveness()
                 << endl;
                );
   }

   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "adjustColumnsEffCnt()", m_param.LogDebugLevel, 2);
   return status;
}

//===========================================================================//
int DecompAlgoPC::compressColumns()
{
   //---
   //--- periodically, get rid of ineffective columns
   //--- periodic:
   //---    every K iterations OR
   //---    numCols has inceased by X since last compression
   //---
   int status = DecompStatOk;
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "compressColumns()", m_param.LogDebugLevel, 2);
   m_stats.timerOther1.reset();
   int nHistorySize
   = static_cast<int>(m_nodeStats.objHistoryBound.size());

   if (nHistorySize > 0) {
      DecompObjBound& objBound
      = m_nodeStats.objHistoryBound[nHistorySize - 1];
      double masterUB  = objBound.thisBoundUB;
      double masterLB  = m_nodeStats.objBest.first;
      double masterGap = m_infinity;

      if (masterUB > -m_infinity &&
            masterUB <  m_infinity) {
         if (masterUB != 0.0) {
            masterGap = fabs(masterUB - masterLB) / masterUB;
         } else {
            masterGap = fabs(masterUB - masterLB);
         }
      }

      if (masterGap > m_param.CompressColumnsMasterGapStart) {
         return status;
      }
   } else {
      return status;
   }

   const int      CompressColsIterFreq      = m_param.CompressColumnsIterFreq;

   const double   CompressColsSizeMultLimit = m_param.CompressColumnsSizeMultLimit;

   const int      nMasterCols               = m_masterSI->getNumCols();

   const int      nMasterRows               = m_masterSI->getNumRows();

   int nColsSinceLast
   = nMasterCols - m_compressColsLastNumCols;

   int nIterSinceLast
   = m_nodeStats.priceCallsTotal - m_compressColsLastPrice;

   int nColsSinceLastLimit
   = static_cast<int>(ceil(m_compressColsLastNumCols *
                           CompressColsSizeMultLimit));

   UTIL_MSG(m_param.LogLevel, 4,
            (*m_osLog) << "nMasterCols              = "
            << nMasterCols << endl;
            (*m_osLog) << "m_compressColsLastNumCols= "
            << m_compressColsLastNumCols << endl;
            (*m_osLog) << "nColsSinceLast           = "
            << nColsSinceLast << endl;
            (*m_osLog) << "priceCallsTotal          = "
            << m_nodeStats.priceCallsTotal << endl;
            (*m_osLog) << "m_compressColsLastPrice  = "
            << m_compressColsLastPrice << endl;
            (*m_osLog) << "nItersSinceLast          = "
            << nIterSinceLast << endl;
           );

   if (nColsSinceLast < nColsSinceLastLimit &&
         nIterSinceLast < CompressColsIterFreq) {
      return status;
   }

   //TODO: reuse memory
   //TODO: using getBasics instead of getBasis since seems cheaper
   int    c;
   int*   basics  = new int[nMasterRows];
   bool* isBasic = new bool[nMasterCols];
   assert(basics && isBasic);
   UtilFillN(isBasic, nMasterCols, false);
   //---
   //--- COIN BUG: to use OSI::getBasics with CLP, you need to
   //---  enableSimplexInterface() - which has issues, so to get around
   //---  this, we will use the warm start object to get basis status
   //---  of variables
   //---
   //    m_masterSI->getBasics(basics);
   //    for(r = 0; r < nMasterRows; r++){
   //      c = basics[r];
   //      if(c < nMasterCols)
   //	   isBasic[c] = true;
   //}
#ifndef DO_INTERIOR
   bool            mustDeleteWS = false;
   CoinWarmStartBasis* warmStart
   = dynamic_cast<CoinWarmStartBasis*>(m_masterSI->getPointerToWarmStart(
                                          mustDeleteWS));

   for (c = 0; c < nMasterCols; c++) {
      if (warmStart->getStructStatus(c) == CoinWarmStartBasis::basic) {
         isBasic[c] = true;
      }
   }

   if (mustDeleteWS) {
      UTIL_DELPTR(warmStart);
   }

#endif
   //---
   //--- sanity check
   //---    m_vars should contain just the structural columns
   //---
   int nMasterColsStruct = 0;
   assert(nMasterCols == static_cast<int>(m_masterColType.size()));

   for (c = 0; c < nMasterCols; c++) {
      if (isMasterColStructural(c)) {
         nMasterColsStruct++;
      }
   }

   //assert(nMasterColsStruct == static_cast<int>(m_vars.size()));
   //several strategies here - we can sort by effCnt and
   // purge those with the worse (only those negative)
   //or we can just purge anything negative
   //or we can purge anything less than some threshold - but not sure
   //  how to set that threshold
   //decide what to purge based on effCnt but also be careful
   // and DO Not purge anything that is currently in the basis!
   // since eff count is based on > eps dual, you might have a degenerate
   // point that has 0 rc but is in basis - so need to check that
   int         shift = 0;
   int         colMasterIndex;
   vector<int> lpColsToDelete;
   vector<int> indexShift;
   UtilFillN(indexShift, nMasterCols, 0);
   UTIL_DEBUG(m_param.LogLevel, 5,
              (*m_osLog) << "VARS before compress:" << endl;
              printVars(m_osLog););
   DecompVarList::iterator li = m_vars.begin();
   int nCols       = 0;
   int nColsNoDel  = 0;
   int nColsBasic  = 0;
   int nColsEffPos = 0;

   while (li != m_vars.end()) {
      colMasterIndex             = (*li)->getColMasterIndex();
      indexShift[colMasterIndex] = shift;
      assert(isMasterColStructural(colMasterIndex));
      nCols++;

      //---
      //--- do not delete any columns that were marked "NoDelete"
      //---   these were degenerate points and deleting them can
      //---   cause cycling
      //---
      if (m_masterColType[colMasterIndex] == DecompCol_Structural_NoDelete
          || m_masterColType[colMasterIndex] == DecompCol_MasterOnly) {
         li++;
         nColsNoDel++;
         continue;
      }

      //---
      //--- do not delete any columns that are basic
      //--- do not delete any columns with non-negative effectiveness
      //---
      if (isBasic[colMasterIndex]) {
         nColsBasic++;
      }

      if ((*li)->getEffectiveness() >= 0) {
         nColsEffPos++;
      }

      if (isBasic[colMasterIndex] || ((*li)->getEffectiveness() >= 0)) {
         li++;
         continue;
      }

      UTIL_DEBUG(m_param.LogLevel, 4,
                 const double* masterSolution            = getMasterPrimalSolution();
                 (*m_osLog) << "CompressCol"
                 << " lpIndex= " << setw(5) << colMasterIndex
                 << " effCnt= "  << setw(2)  << (*li)->getEffectiveness()
                 << " currSol= " << setw(10)
                 << UtilDblToStr(masterSolution[colMasterIndex], 3) << endl;
                );
      //add this var to pool - THINK what exactly does that entail
      (*li)->resetEffectiveness();
      //this deletes the var object (we won't do this
      // once we move to pool)
      delete *li;
      li = m_vars.erase(li); //removes link in list
      lpColsToDelete.push_back(colMasterIndex);
      m_masterColType[colMasterIndex] = DecompCol_ToBeDeleted;
      shift++;
   }

   if (lpColsToDelete.size() > 0) {
      /*for(c = 0; c < m_masterSI->getNumCols(); c++){
      const string colN = m_masterSI->getColName(c);
      printf("Before Col[%4d] Name: %30s Type: %20s\n",
      c,
      colN.c_str(),
      DecompColTypeStr[m_masterColType[c]].c_str());
                }*/
      m_masterSI->deleteCols(static_cast<int>(lpColsToDelete.size()),
                             &lpColsToDelete[0]);
      m_cutpool.setRowsAreValid(false);
      UTIL_MSG(m_param.LogLevel, 3,
               (*m_osLog) << "Num Columns Deleted = "
               << lpColsToDelete.size()
               << " Cols = " << nCols
               << " NoDel = " << nColsNoDel
               << " Basic = " << nColsBasic
               << " EffPos = " << nColsEffPos
               << endl;
              );

      //---
      //--- now, we must update the mapping between LP index and
      //---  the index in the var list objects - but we might have
      //---  artificial columns lurking in between the LP columns
      //---
      //--- Example:
      //---   a=artificial
      //---   s=structural (either from original row, branch row or cut row)
      //---
      //---   lpColsToDelete = {6,7,15}
      //---     000000000011111111
      //---     012345678901234567
      //---     aaaassssaasssaasss
      //---
      //---     000000000011111111
      //---     012345678901234567
      //---     aaaas..saasssaa.ss
      //---     aaaassaasssaass
      //---   shift
      //---     00000..22222222.33
      //---
      //---
      //---
      //--- reset the master index in m_vars
      //---
      for (li = m_vars.begin(); li != m_vars.end(); li++) {
         colMasterIndex = (*li)->getColMasterIndex();
         (*li)->setColMasterIndex(colMasterIndex - indexShift[colMasterIndex]);
      }

      //---
      //--- delete the entries in vector m_masterColType
      //---   NOTE: this would be much faster if used list instead of vector
      //---
      vector<DecompColType>::iterator vi = m_masterColType.begin();

      while (vi != m_masterColType.end()) {
         if (*vi == DecompCol_ToBeDeleted) {
            vi = m_masterColType.erase(vi);
         } else {
            vi++;
         }
      }

      //---
      //--- sanity check
      //---    m_vars should contain just the structural columns
      //---
      int nMasterColsNew = m_masterSI->getNumCols();
      nMasterColsStruct  = 0;
      assert(nMasterColsNew == static_cast<int>(m_masterColType.size()));

      for (li = m_vars.begin(); li != m_vars.end(); li++) {
         colMasterIndex = (*li)->getColMasterIndex();
         assert(isMasterColStructural(colMasterIndex));
      }

      for (c = 0; c < nMasterColsNew; c++) {
         if (isMasterColStructural(c)) {
            nMasterColsStruct++;
         }
      }

      //assert(nMasterColsStruct == static_cast<int>(m_vars.size()));
      UTIL_DEBUG(m_param.LogLevel, 5,
                 (*m_osLog) << "VARS after compress:" << endl;
                 printVars(m_osLog););
      /*for(c = 0; c < m_masterSI->getNumCols(); c++){
      const string colN = m_masterSI->getColName(c);
      printf("After Col[%4d] Name: %30s Type: %20s\n",
      c,
      colN.c_str(),
      DecompColTypeStr[m_masterColType[c]].c_str());
                }*/
      //---
      //--- we deleted something, so reset the counters
      m_compressColsLastPrice   = m_nodeStats.priceCallsTotal;
      m_compressColsLastNumCols = m_masterSI->getNumCols();
      //---
      //--- if any vars were deleted, do a solution update to refresh
      //---
      status = solutionUpdate(m_phase, 99999, 99999);
   }

   m_stats.thisCompressCols.push_back(m_stats.timerOther1.getRealTime());
   UTIL_DELARR(basics);
   UTIL_DELARR(isBasic);
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "compressColumns()", m_param.LogDebugLevel, 2);
   return status;
}

//===========================================================================//
void DecompAlgoPC::phaseDone()
{
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "phaseDone()", m_param.LogDebugLevel, 2);

   if (m_param.SolveMasterAsMip                                &&
         getNodeIndex() % m_param.SolveMasterAsMipFreqNode == 0  &&
         m_stopCriteria != DecompStopTime) {
      solveMasterAsMIP();
   }

   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "phaseDone()", m_param.LogDebugLevel, 2);
}

//===========================================================================//
void DecompAlgoPC::solveMasterAsMIP()
{
   //---
   //--- if node was already found infeasible, just return
   //---
   if (m_status == STAT_INFEASIBLE) {
      return;
   }

   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "solveMasterAsMIP()", m_param.LogDebugLevel, 2);
   //---
   //--- no point in doing this if only one block, we check each
   //---  new column to see if it is feasible to original already
   //---
   assert(m_numConvexCon > 1);
   int  nMasterCols = m_masterSI->getNumCols();//lambda
   DecompConstraintSet* modelCore = m_modelCore.getModel();
   //---
   //--- set the master (generated) columns (lambda) to integer
   //--- set the master-onlys (that are integral) to integer
   //---
   int         colIndex;
   const char* intMarkerCore  = modelCore->getIntegerMark();

   for (colIndex = 0; colIndex < nMasterCols; colIndex++) {
      if (isMasterColStructural(colIndex)){
         m_masterSI->setInteger(colIndex);
      }
   }
   for (int i = 0; i < m_masterOnlyCols.size(); i++){
      if (intMarkerCore[m_masterOnlyCols[i]] == 'I'){
	 m_masterSI->setInteger(m_masterOnlyColsMap[m_masterOnlyCols[i]]);
      } 
   }
   if (m_param.LogDumpModel >= 2)
      printCurrentProblem(m_masterSI,
                          "masterProbRootIP",
                          m_nodeStats.nodeIndex,
                          m_nodeStats.cutCallsTotal,
                          m_nodeStats.priceCallsTotal);

   DecompSolverResult    result(m_infinity);

   if (m_param.DecompIPSolver == "SYMPHONY"){
      solveMasterAsMIPSym(&result);
   }else if (m_param.DecompIPSolver == "Cbc"){
      solveMasterAsMIPCbc(&result);
   }else if (m_param.DecompIPSolver == "CPLEX"){
      solveMasterAsMIPCpx(&result);
   }else if (m_param.DecompIPSolver == "Gurobi" ||
             m_param.DecompIPSolver == "Xpress"){
      solveMasterAsMIPOsi(&result);
   }else{
      throw UtilException("Unknown solver selected",
			  "solveMasterAsMIP", "DecompAlgoPC");
   }

   if (result.m_nSolutions) {
      double* rsolution = new double[modelCore->getNumCols()];

      if (!rsolution) {
         throw UtilExceptionMemory("solveMasterAsMIP", "DecompAlgoPC");
      }

      UTIL_MSG(m_param.LogLevel, 3,
               (*m_osLog) << "Solve as IP found a solution." << endl;);
      recomposeSolution(result.getSolution(0), rsolution);

      if (!isIPFeasible(rsolution))
         throw UtilException("Recomposed solution is not feasible",
                             "solveMasterAsMIP", "DecompAlgoPC");

      if (m_app->APPisUserFeasible(rsolution,
                                   modelCore->getNumCols(),
                                   m_param.TolZero)) {
         UTIL_MSG(m_param.LogLevel, 3,
                  (*m_osLog) << "Solution is app-feasible, nSolutions="
                  << (int)m_xhatIPFeas.size() << endl;);
         //check for dup sol - TODO: make func
         bool isDup = m_xhatIPFeas.size() > 0 ? true : false;
         vector<DecompSolution*>::iterator vit;

         for (vit  = m_xhatIPFeas.begin();
               vit != m_xhatIPFeas.end(); vit++) {
            const DecompSolution* xhatIPFeas = *vit;
            const double*          values
            = xhatIPFeas->getValues();

            for (int c = 0; c < modelCore->getNumCols(); c++) {
               if (!UtilIsZero(values[c] - rsolution[c])) {
                  isDup = false;
                  break;
               }
            }
         }

         if (isDup) {
            UTIL_MSG(m_param.LogLevel, 3,
                     (*m_osLog) << "Solution is a duplicate, not pushing."
                     << endl;);
         } else {
            DecompSolution* decompSol
            = new DecompSolution(modelCore->getNumCols(),
                                 rsolution,
                                 getOrigObjective());
            m_xhatIPFeas.push_back(decompSol);
            vector<DecompSolution*>::iterator vi;
            DecompSolution* viBest = NULL;
            double bestBoundUB = m_nodeStats.objBest.second;

            for (vi = m_xhatIPFeas.begin(); vi != m_xhatIPFeas.end(); vi++) {
               const DecompSolution* xhatIPFeas = *vi;

               if (xhatIPFeas->getQuality() <= bestBoundUB) {
                  bestBoundUB = xhatIPFeas->getQuality();
                  viBest = *vi;
               }
            }

            if (viBest) {
               //save the best
               setObjBoundIP(bestBoundUB);
               m_xhatIPBest = viBest;
            }
         }
      }

      if (m_param.LogDebugLevel >= 3) {
         int j;
         const vector<string>& colNames = modelCore->getColNames();

         for (j = 0; j < modelCore->getNumCols(); j++) {
            if (fabs(rsolution[j]) > DecompEpsilon) {
               if (j < static_cast<int>(colNames.size()))
                  printf("MASTER PRIM[%6d->%20s] = %12.10f\n",
                         j, colNames[j].c_str(), rsolution[j]);
               else
                  printf("MASTER PRIM[%6d] = %12.10f\n",
                         j, rsolution[j]);
            }
         }
      }

      UTIL_DELARR(rsolution);
   }

   //---
   //--- set the master columns back to continuous
   //---
   for (colIndex = 0; colIndex < nMasterCols; colIndex++) {
      if (isMasterColStructural(colIndex) ||
            isMasterColMasterOnly(colIndex)) {
         m_masterSI->setContinuous(colIndex);
      }
   }

   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "solveMasterAsMIP()", m_param.LogDebugLevel, 2);
}

//===========================================================================//
void DecompAlgoPC::solveMasterAsMIPSym(DecompSolverResult* result)
{
#ifdef DIP_HAS_SYMPHONY
   int colIndex;
   int numCols = m_masterSI->getNumCols();
   int nMasterCols = m_masterSI->getNumCols();//lambda
   int logIpLevel  = m_param.LogIpLevel;
   const char* intMarkerCore  = m_modelCore.getModel()->getIntegerMark();

   OsiSymSolverInterface* osiSym = new OsiSymSolverInterface();

#if 0
   OsiSolverInterface * m_masterClone = m_masterSI->clone();
 
   const CoinPackedMatrix* matrix_sym = (m_masterClone->getMatrixByRow());
   const double* col_lb = (m_masterClone->getColLower());
   const double* col_up = (m_masterClone->getColUpper());
   const double* row_lb = (m_masterClone->getRowLower());
   const double* row_up = (m_masterClone->getRowUpper());
   const double* obj_coef = (m_masterClone->getObjCoefficients());
   osiSym->assignProblem(const_cast<CoinPackedMatrix*&>(matrix_sym),
                          const_cast<double*&>(col_lb),
                          const_cast<double*&>(col_up), const_cast<double*&>(obj_coef),
                          const_cast<double*&>(row_lb), const_cast<double*&>(row_up));
   for (colIndex = 0; colIndex < nMasterCols; colIndex++) {
      if (isMasterColStructural(colIndex)){
         osiSym->setInteger(colIndex);
      }
   }
   for (int i = 0; i < m_masterOnlyCols.size(); i++){
      if (intMarkerCore[m_masterOnlyCols[i]] == 'I'){
	 osiSym->setInteger(m_masterOnlyColsMap[m_masterOnlyCols[i]]);
      } 
   }
#else

   osiSym->loadProblem(const_cast<CoinPackedMatrix&>(*m_masterSI->getMatrixByRow()),
		       const_cast<double*>(m_masterSI->getColLower()),
		       const_cast<double*>(m_masterSI->getColUpper()), 
		       const_cast<double*>(m_masterSI->getObjCoefficients()),
		       const_cast<double*>(m_masterSI->getRowLower()), 
		       const_cast<double*>(m_masterSI->getRowUpper()));

   for (colIndex = 0; colIndex < nMasterCols; colIndex++) {
      if (isMasterColStructural(colIndex)){
         osiSym->setInteger(colIndex);
      }
   }
   for (int i = 0; i < m_masterOnlyCols.size(); i++){
      if (intMarkerCore[m_masterOnlyCols[i]] == 'I'){
         osiSym->setInteger(m_masterOnlyColsMap[m_masterOnlyCols[i]]);
      }
   }
#endif
   assert(osiSym);
   sym_environment* env = osiSym->getSymphonyEnvironment();
   if (logIpLevel == 0){
     sym_set_int_param(env, "verbosity", -10);
   }
   else{
     sym_set_int_param(env, "verbosity", logIpLevel);
   }

   assert(env);
   osiSym->branchAndBound();
   int status = sym_get_status(env);

   if ((status == PREP_OPTIMAL_SOLUTION_FOUND) ||
         (status == TM_OPTIMAL_SOLUTION_FOUND)
         || (status == TM_TARGET_GAP_ACHIEVED)) {
      result->m_isOptimal = true;
      double* solution = new double[numCols];
      assert(solution);
      status = sym_get_col_solution(env, solution);
      result->m_nSolutions = 1;
      vector<double> solVec(solution, solution + numCols);
      result->m_solution.push_back(solVec);
      UTIL_DELARR(solution);

      if (status == FUNCTION_TERMINATED_ABNORMALLY)
         throw UtilException("sym_get_col_solution failure",
                             "solveMasterAsMIP", "DecompSubModel");
   } else {
      if (sym_is_proven_primal_infeasible(env)) {
         result->m_nSolutions = 0;
         result->m_isOptimal = true;
         //	 result->m_isCutoff = doCutoff;
      } else {
         //	 result->m_isCutoff = doCutoff;
         result->m_isOptimal = false ;
      }
   }

   if (status == (TM_ERROR__USER || TM_ERROR__COMM_ERROR
                  || TM_ERROR__NUMERICAL_INSTABILITY
                  || TM_ERROR__ILLEGAL_RETURN_CODE
                  || TM_ERROR__NO_BRANCHING_CANDIDATE)) {
      std::cerr << "Error: SYPHONMY IP solver status =  "
                << status << std::endl;
   }

   UTIL_DELPTR(osiSym);
#else
   throw UtilException("SYMPHONY selected as solver, but it's not available",
		       "solveMasterAsMIPSym", "DecompAlgoPC");
#endif
}

//===========================================================================//
void DecompAlgoPC::solveMasterAsMIPCbc(DecompSolverResult* result)
{
#ifdef DIP_HAS_CBC
   int  nMasterCols = m_masterSI->getNumCols();//lambda
   int  logIpLevel  = m_param.LogIpLevel;
   //TODO: what exactly does this do? make copy of entire model!?
   CbcModel cbc(*m_masterSI);
   cbc.setLogLevel(logIpLevel);
   cbc.setDblParam(CbcModel::CbcAllowableFractionGap,
                   m_param.SolveMasterAsMipLimitGap);
   cbc.setDblParam(CbcModel::CbcMaximumSeconds, m_param.SolveMasterAsMipTimeLimit);
   cbc.setDblParam(CbcModel::CbcCurrentCutoff, m_globalUB);
#if 0
   cbc.branchAndBound();
#else
   CbcMain0(cbc);
   //---
   //--- build argument list
   //---
   //TODO: time limit, cutoff,gap
   const char* argv[20];
   int    argc         = 0;
   string cbcExe       = "cbc";
   string cbcSolve     = "-solve";
   string cbcQuit      = "-quit";
   string cbcLog       = "-log";
   string cbcLogSet    = UtilIntToStr(logIpLevel);
   string cbcGap       = "-ratio";
   string cbcGapSet    = UtilDblToStr(m_param.SolveMasterAsMipLimitGap);
   string cbcTime      = "-seconds";
   string cbcTimeSet   = UtilDblToStr(m_param.SolveMasterAsMipTimeLimit);
   string cbcCutoff    = "-cutoff";
   string cbcCutoffSet = UtilDblToStr(m_globalUB, -1, 1.0e100);
   argv[argc++] = cbcExe.c_str();
   argv[argc++] = cbcLog.c_str();
   argv[argc++] = cbcLogSet.c_str();
   argv[argc++] = cbcGap.c_str();
   argv[argc++] = cbcGapSet.c_str();
   argv[argc++] = cbcTime.c_str();
   argv[argc++] = cbcTimeSet.c_str();
   argv[argc++] = cbcCutoff.c_str();
   argv[argc++] = cbcCutoffSet.c_str();
   argv[argc++] = cbcSolve.c_str();
   argv[argc++] = cbcQuit.c_str();
   //---
   //--- solve IP using argument list
   //---
   CbcMain1(argc, argv, cbc);
#endif
   //---
   //--- get solver status
   //---   comments based on Cbc2.3
   //---
   /** Final status of problem.
    *   -1  before branchAndBound
    *    0  finished - check isProvenOptimal or isProvenInfeasible
    *         to see if solution found (or check value of best solution)
    *    1  stopped - on maxnodes, maxsols, maxtime
    *    2  difficulties so run was abandoned
    *   (5  event user programmed event occurred)
   */
   const int statusSet[2] = {0, 1};
   result->m_solStatus    = cbc.status();

   if (!UtilIsInSet(result->m_solStatus, statusSet, 2)) {
      cerr << "Error: CBC IP solver status = " << result->m_solStatus << endl;
      //This shouldn't really cause an exception
      //throw UtilException("CBC solver status",
      //                    "solveMasterAsMIP", "DecompSubModel");
   }

   /** Secondary status of problem
    *   -1 unset (status_ will also be -1)
    *    0 search completed with solution
    *    1 linear relaxation not feasible (or worse than cutoff)
    *    2 stopped on gap
    *    3 stopped on nodes
    *    4 stopped on time
    *    5 stopped on user event
    *    6 stopped on solutions
    *    7 linear relaxation unbounded
    */
   const int statusSet2[4] = {0, 1, 2, 4};
   result->m_solStatus2 = cbc.secondaryStatus();

   //---
   //--- In root the subproblem should not be infeasible
   //---   unless due to cutoff. But, after branching it
   //---   can be infeasible.
   //---
   if (!UtilIsInSet(result->m_solStatus2, statusSet2, 4)) {
      cerr << "Warning: CBC IP solver 2nd status = "
           << result->m_solStatus2 << endl;
      //This shouldn't really cause an exception
      //throw UtilException("CBC solver 2nd status",
      //                    "solveMasterAsMIP", "DecompAlgoPC");
   }

   //---
   //--- update results object
   //---
   result->m_nSolutions = 0;
   result->m_isOptimal  = false;
   //TODO: can get multple solutions!
   //   how to retrieve?
   //TODO: look into setHotstartSolution... done automatically
   //   look at one call to the next
   //TODO: look into setNumberThreads
   //TODO: redo cpx in this same way - it could be stopping on time, not gap
   int nSolutions = cbc.getSolutionCount();
   result->m_nSolutions = nSolutions ? 1 : 0;

   if (cbc.isProvenOptimal() ||
         cbc.isProvenInfeasible()) {
      result->m_isOptimal  = true;
   }

   //---
   //--- get copy of solution
   //---
   result->m_objLB = cbc.getBestPossibleObjValue();

   if (nSolutions >= 1) {
      result->m_objUB = cbc.getObjValue();
      const double* solDbl = cbc.getColSolution();
      vector<double> solVec(solDbl, solDbl + nMasterCols);
      result->m_solution.push_back(solVec);
      assert(result->m_nSolutions ==
             static_cast<int>(result->m_solution.size()));
      //memcpy(result->m_solution,
      //     cbc.getColSolution(),
      //     nMasterCols * sizeof(double));
   }

#else
   throw UtilException("Cbc selected as solver, but it's not available",
		       "solveMasterAsMIPCbc", "DecompAlgoPC");
#endif
}

//===========================================================================//
void DecompAlgoPC::solveMasterAsMIPCpx(DecompSolverResult* result)
{
#ifdef DIP_HAS_CPX
   //---
   //--- get OsiCpx object from Osi object
   //--- get CPEXENVptr for use with internal methods
   //--- get CPXLPptr   for use with internal methods
   //---
   int  nMasterCols = m_masterSI->getNumCols();//lambda
   int  logIpLevel  = m_param.LogIpLevel;
   OsiCpxSolverInterface* osiCpx
   = dynamic_cast<OsiCpxSolverInterface*>(m_masterSI);
   CPXENVptr cpxEnv = osiCpx->getEnvironmentPtr();
   CPXLPptr  cpxLp  = osiCpx->getLpPtr();
   assert(cpxEnv && cpxLp);
   //---
   //--- set parameters
   //---
   int status = 0;

   if (logIpLevel) {
      status = CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_ON);

      if (status)
         throw UtilException("CPXsetintparam failure",
                             "solveMasterAsMIP", "DecompSubModel");

      status = CPXsetintparam(cpxEnv, CPX_PARAM_SIMDISPLAY, logIpLevel);

      if (status)
         throw UtilException("CPXsetintparam failure",
                             "solveMasterAsMIP", "DecompSubModel");
   } else {
      status = CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);

      if (status)
         throw UtilException("CPXsetintparam failure",
                             "solveMasterAsMIP", "DecompSubModel");
   }

   if (m_firstPhase2Call) {
      //---
      //--- if calling with first Phase2 call, it is meant to
      //---   "recombine" partial columns - i.e., if the user
      //---   produced a fully feasible solution that was then
      //---   separated into blocks - we want to be sure it
      //---   at least recombines it
      //--- so, make the stop on gap very small
      //---
      //--- TODO: we should get this incumbent in the system without
      //--- forcing the call to IP solver just to recombine
      //---
      status = CPXsetdblparam(cpxEnv, CPX_PARAM_EPGAP, 0.005); //0.5%
   } else {
      status = CPXsetdblparam(cpxEnv, CPX_PARAM_EPGAP,
                              m_param.SolveMasterAsMipLimitGap);
   }

   if (status)
      throw UtilException("CPXsetdblparam failure",
                          "solveMasterAsMIP", "DecompAlgoPC");

   status = CPXsetdblparam(cpxEnv, CPX_PARAM_TILIM,
                           m_param.SolveMasterAsMipTimeLimit);

   if (status)
      throw UtilException("CPXsetdblparam failure",
                          "solveMasterAsMIP", "DecompAlgoPC");

   status = CPXsetdblparam(cpxEnv, CPX_PARAM_CUTUP, m_globalUB);

   if (status)
      throw UtilException("CPXsetdblparam failure",
                          "solveMasterAsMIP", "DecompAlgoPC");

#if CPX_VERSION >= 1100
   status = CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, 1);

   if (status)
      throw UtilException("CPXsetintparam failure",
                          "solveMasterAsMIP", "DecompAlgoPC");

#endif
   //---
   //--- solve the MILP
   //---
   osiCpx->branchAndBound();
   //---
   //--- get solver status
   //---
   result->m_solStatus  = CPXgetstat(cpxEnv, cpxLp);
   result->m_solStatus2 = 0;

   //cout << "CPX IP solver status = " << result->m_solStatus << endl;
   //TEMP FIX?
   //THINK: if CPXMIP_INForUNBD, change to CPXMIP_INFEASIBLE,
   // I don't think there is anyway the price+branch heur could
   // be unbounded. But, what if the original full problem is unbounded?
   if (result->m_solStatus == CPXMIP_INForUNBD ) {
      result->m_solStatus = CPXMIP_INFEASIBLE;
   }

   const int statusSet[5] = {CPXMIP_OPTIMAL,
                             CPXMIP_OPTIMAL_TOL,
                             CPXMIP_INFEASIBLE,
                             CPXMIP_TIME_LIM_FEAS,
                             CPXMIP_TIME_LIM_INFEAS
                            };

   if (!UtilIsInSet(result->m_solStatus, statusSet, 5)) {
      cerr << "Error: CPX IP solver status = " << result->m_solStatus << endl;
      throw UtilException("CPX solver status",
                          "solveMasterAsMIP", "DecompAlgoPC");
   }

   //---
   //--- update results object
   //---
   result->m_nSolutions = 0;
   result->m_isOptimal  = false;

   if (result->m_solStatus == CPXMIP_OPTIMAL ||
         result->m_solStatus == CPXMIP_OPTIMAL_TOL) {
      result->m_nSolutions = 1;
      result->m_isOptimal  = true;
   } else {
      if (result->m_solStatus == CPXMIP_INFEASIBLE ||
            result->m_solStatus == CPXMIP_TIME_LIM_INFEAS) {
         result->m_nSolutions = 0;
         result->m_isOptimal  = true;
      }
      //STOP - could have stopped on time... not just gap... do
      //something like did in CBC
      else {
         //---
         //--- else it must have stopped on gap
         //---
         result->m_nSolutions = 1;
         result->m_isOptimal  = false;
      }
   }

   //---
   //--- get copy of solution
   //---
   status = CPXgetbestobjval(cpxEnv, cpxLp, &result->m_objLB);

   if (status)
      throw UtilException("CPXgetbestobjval failure",
                          "solveMasterAsMIP", "DecompAlgoPC");

   if (result->m_nSolutions >= 1) {
      status = CPXgetmipobjval(cpxEnv, cpxLp, &result->m_objUB);

      if (status)
         throw UtilException("CPXgetmipobjval failure",
                             "solveMasterAsMIP", "DecompAlgoPC");

      const double* solDbl = osiCpx->getColSolution();
      vector<double> solVec(solDbl, solDbl + nMasterCols);
      result->m_solution.push_back(solVec);
      assert(result->m_nSolutions ==
             static_cast<int>(result->m_solution.size()));
      //memcpy(result->m_solution,
      //     cbc.getColSolution(),
      //     nMasterCols * sizeof(double));
   }
   //---
   //--- set time back
   //---
   status = CPXsetdblparam(cpxEnv, CPX_PARAM_TILIM, m_infinity);

   if (status)
      throw UtilException("CPXsetdblparam failure",
                          "solveMasterAsMIP", "DecompAlgoPC");
#else
   throw UtilException("CPLEX selected as solver, but it's not available",
		       "solveMasterAsMIPCpx", "DecompAlgoPC");
#endif
}

//===========================================================================//
void DecompAlgoPC::solveMasterAsMIPOsi(DecompSolverResult* result)
{
   const int numCols = m_masterSI->getNumCols();

   m_masterSI->branchAndBound();

   result->m_isUnbounded = false;
   result->m_isOptimal   = false;
   result->m_isCutoff    = false;
   result->m_nSolutions  = 0;
   if (m_masterSI->isProvenOptimal()){
      const double *solution = m_masterSI->getColSolution();
      vector<double> solVec(solution, solution + numCols);
      result->m_solution.push_back(solVec);
      result->m_nSolutions++;
      result->m_isOptimal   = true;
   }else{
      result->m_isOptimal = true;
   }
}

//===========================================================================//
//because rowReform, this is very specific to PC
void DecompAlgoPC::addCutsToPool(const double*    x,
                                 DecompCutList& newCuts,
                                 int&            n_newCuts)
{
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "addCutsToPool()", m_param.LogDebugLevel, 2);
   int  r;
   int  cutIndex = 0;
   bool isDupCore;//also check relax?
   bool isDupPool;
   bool isViolated; //TODO: do something similiar to check for pos-rc vars
   bool addCut;
   DecompConstraintSet*           modelCore   = m_modelCore.getModel();
   DecompCutPool::iterator ci;
   DecompCutList::iterator li = newCuts.begin();

   while (li != newCuts.end()) {
      CoinPackedVector* row       = new CoinPackedVector();
      //---
      //--- create a row (in terms of original formulation, x), from a cut
      //---
      (*li)->expandCutToRow(row);
      //---
      //--- set the hash string (for quick duplicate checks)
      //---
      (*li)->setStringHash(row, m_infinity);
      //bool isOptViolated = false;
      //for(i = 0; i < m_optPoint.size(); i++){
      //isOptViolated = (*li)->calcViolation(row, &m_optPoint[i][0]);
      //if(isOptViolated){
      //    (*m_osLog) << "\n\nCUT VIOLATES OPT POINT";
      //    (*li)->print();
      //	 }
      // assert(!isOptViolated);
      //}
      //---
      //--- check the the cut is already in the model core
      //---   NOTE: if so this is an error (always?)
      //---
      addCut    = true;
      isDupCore = false;

      for (r = 0; r < modelCore->getNumRows(); r++) {
         //TODO: allow user to set hash
         // example: GSEC can be represented compactly with just S
         // or directly them they override an isSame( )
         if (modelCore->rowHash[r] == (*li)->getStrHash()) {
            (*m_osLog) << "CUT IS DUPLICATE with Core\n";
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
            double actViol;
            double relViol;
            double cutLB    = (*li)->getLowerBound();
            double cutUB    = (*li)->getUpperBound();
            double ax       = row->dotProduct(x);
            actViol = std::max<double>(cutLB - ax, ax - cutUB);
            actViol = std::max<double>(actViol, 0.0);

            if (UtilIsZero(ax)) {
               relViol = actViol;
            } else {
               relViol = actViol / std::fabs(ax);
            }

            //TODO: need status return not just assert
            //---
            //--- since it is already in LP core, the violation
            //---  should be very small
            //---
            if (relViol > 0.005) { //0.5% violated
               (*m_osLog) << "CUT actViol= " << actViol
                          << " relViol= "    << relViol << "\n";
               (*li)->print(m_osLog);
               assert(0);//0.1% violated
            }

            isDupCore = true;
            break;
         }
      }

      if (isDupCore) {
         addCut = false;
      } else {
         //---
         //--- is this cut already in pool
         //---  NOTE: this is not neccessarily an error, since
         //---   there could be a cut from a previous iteration
         //---   in the cut pool that was not entered because of
         //---   the limit on the number of cuts entered per iteration
         //---
         int cutIndexPool = 0;
         isDupPool = false;

         for (ci = m_cutpool.begin(); ci != m_cutpool.end(); ci++) {
            if ((*li)->getStrHash() == (*ci).getCutPtr()->getStrHash()) {
               UTIL_MSG(m_param.LogLevel, 3,
                        (*m_osLog) << "CUT "              << cutIndex
                        << " is Duplicate with Pool Cut " << cutIndexPool
                        << endl;
                        (*m_osLog) << "CUT           Hash = "
                        << (*li)->getStrHash() << endl;
                        (*m_osLog) << "CUT (in Pool) Hash = "
                        << (*ci).getCutPtr()->getStrHash() << endl;
                        (*li)->print();
                       );
               isDupPool = true;
               break;
            }

            cutIndexPool++;
         }

         if (isDupPool) {
            addCut = false;
         } else {
            isViolated = (*li)->calcViolation(row, x);//also sets it

            if (!isViolated) {
               //---
               //--- we are trying to add a cut that is NOT violated
               //---  NOTE: this is probably an error in the cut gen
               //---  THINK: are there cases where we want to add cuts
               //---    to pool even though we know they are not violated
               //---    at the current point? i.e., might be violated later?
               //---
               addCut = false;
               (*m_osLog) << "CUT "                         << cutIndex
                          << " is not violated! Not adding to pool.\n";
               (*m_osLog) << "CUT           Hash = "
                          << (*li)->getStrHash() << "\n";
               (*li)->print();
               assert(0);
            }
         }
      }

      if (addCut) {
         //---
         //--- create a row (in terms of reformulation, lambda), from row
         //---
         CoinPackedVector* rowReform
         = m_cutpool.createRowReform(modelCore->getNumCols(),
                                     row,
                                     m_vars);
	 int tempIndex(0); 
	 map<int, int >:: iterator mit;
	 for (int i = 0; i < row->getNumElements(); i++){
	   tempIndex = row->getIndices()[i];
	   mit =  m_masterOnlyColsMap.find(tempIndex);
	   if (mit != m_masterOnlyColsMap.end()){
	     rowReform->insert(mit->second, row->getElements()[i]);
	   } 
	 }

	
         if (!rowReform) {
            //TODO: need status return code for failure in -O
            (*m_osLog) << "ERROR in createRowReform\n";
            assert(0);
         } else {
            DecompWaitingRow waitingRow(*li, row, rowReform);
            //do this in a separate function so addCutsTo is not dependent
            //on passing in osolution for DecompVar
            //waitingRow.setViolation(x);//always on original solution!
            m_cutpool.push_back(waitingRow);
         }

         li++;
      } else {
         //---
         //--- cut is not being added to pool, delete memory
         //---
         UTIL_DELPTR(row);
         UTIL_DELPTR(*li);       //need to do?
         li = newCuts.erase(li); //does this call cut destructor?
         n_newCuts--;
      }

      cutIndex++;
   }

   CoinAssertDebug(n_newCuts >= 0);
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "addCutsToPool()", m_param.LogDebugLevel, 2);
}

//===========================================================================//
int DecompAlgoPC::addCutsFromPool()
{
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "addCutsFromPool()", m_param.LogDebugLevel, 2);
   //TODO: make this a parameter
   const int maxcuts_toadd = 100;//m_app->m_param.cut_maxcuts_periter;
   DecompConstraintSet*           modelCore   = m_modelCore.getModel();
   int n_newrows = CoinMin(static_cast<int>(m_cutpool.size()), maxcuts_toadd);
   int index = 0;
   //---
   //--- sort the cuts by violation
   //---  TODO: partial sort (limit by n_newrows)
   //---
   sort(m_cutpool.begin(), m_cutpool.end(), is_greater_thanD());
   //---
   //--- after sorting by violation, find the index that starts
   //---  where there are no violations (this can happen if pool
   //---  has leftover cuts from previous iterations due to limitation
   //---  on number of cuts entered per pass)
   //---
   DecompCutPool::iterator li;

   for (li = m_cutpool.begin(); li != m_cutpool.end(); li++) {
      if (m_param.LogDebugLevel >= 3) {
         (*m_osLog) << "CUT VIOLATION = " << (*li).getViolation() << endl;
      }

      if ((*li).getViolation() < DecompEpsilon) { //PARM
         break;
      }

      index++;
   }

   n_newrows = std::min<int>(n_newrows, index);

   if (n_newrows > 0) {
      m_varpool.setColsAreValid(false);
   }

   //TODO: look into coin build...
   double* rlb = new double[n_newrows];
   double* rub = new double[n_newrows];
   const CoinPackedVectorBase** rowReformBlock =
      new const CoinPackedVectorBase*[n_newrows];
   const CoinPackedVectorBase** rowBlock   =
      new const CoinPackedVectorBase*[n_newrows];
   //better design to have a "add row name"
   vector<string>& coreRowNames   = modelCore->getRowNamesMutable();
   vector<string>   colNames, rowNames;
   string           colName,  rowName;
   int              rowIndex, rowIndex0;
   int              colIndex, colIndex0;
   char             sense;
   double           rhs, range;
   index     = 0;
   rowIndex0 = m_masterSI->getNumRows();
   colIndex0 = m_masterSI->getNumCols();

   for (li = m_cutpool.begin(); li != m_cutpool.end(); li++) {
      if (index >= n_newrows) {
         break;
      }

      CoinPackedVector* rowReform = (*li).getRowReformPtr();
      CoinPackedVector* row       = (*li).getRowPtr();
      DecompCut*         cut       = (*li).getCutPtr();
      rlb[index]            = (*li).getLowerBound();
      rub[index]            = (*li).getUpperBound();
      rowReformBlock[index] = rowReform;
      rowBlock[index]       = row;
      rowIndex              = rowIndex0 + index;
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

   //---
   //--- add the new (lambda) rows to master
   //--- add the new (x)      rows to core
   //--- add the master row types
   //--- add the row names to core
   //--- add the row lb,ub to core
   //--- add a new artificial column for this cut (fix to 0)
   //---
   m_masterSI->addRows(n_newrows, rowReformBlock, rlb, rub);
   modelCore->M->appendRows(n_newrows, rowBlock);

   for (index = 0; index < n_newrows; index++) {
      if (rowNames.size()) {
         coreRowNames.push_back(rowNames[index]);
      }

      m_masterRowType.push_back(DecompRow_Cut);
      //TODO: make this a function
      UtilBoundToSense(rlb[index], rub[index],
                       m_infinity, sense, rhs, range);
      modelCore->rowLB.push_back(rlb[index]);
      modelCore->rowUB.push_back(rub[index]);
      modelCore->rowSense.push_back(sense);
      modelCore->rowRhs.push_back(rhs);
      rowIndex = rowIndex0 + index;
      colIndex = colIndex0 + index;

      switch (sense) {
      case 'L': {
         CoinPackedVector artCol;
         artCol.insert(rowIndex, -1.0);
         m_masterSI->addCol(artCol, 0.0, 0.0, 0.0);
         m_masterColType.push_back(DecompCol_ArtForCutL);
         m_masterArtCols.push_back(colIndex);
         colName = "sCL(c_" + UtilIntToStr(colIndex)
                   + "_" + UtilIntToStr(rowIndex) + ")";
         colNames.push_back(colName);
         colIndex++;
      }
      break;
      case 'G': {
         CoinPackedVector artCol;
         artCol.insert(rowIndex, 1.0);
         m_masterSI->addCol(artCol, 0.0, 0.0, 0.0);
         m_masterColType.push_back(DecompCol_ArtForCutG);
         m_masterArtCols.push_back(colIndex);
         colName = "sCG(c_" + UtilIntToStr(colIndex)
                   + "_" + UtilIntToStr(rowIndex) + ")";
         colNames.push_back(colName);
         colIndex++;
      }
      break;
      case 'E': {
         CoinPackedVector artColL;
         CoinPackedVector artColG;
         artColL.insert(rowIndex, -1.0);
         m_masterSI->addCol(artColL, 0.0, 0.0, 0.0);
         m_masterColType.push_back(DecompCol_ArtForCutL);
         m_masterArtCols.push_back(colIndex);
         colName = "sCL(c_" + UtilIntToStr(colIndex)
                   + "_" + UtilIntToStr(rowIndex) + ")";
         colNames.push_back(colName);
         artColG.insert(rowIndex,  1.0);
         m_masterSI->addCol(artColG, 0.0, 0.0, 0.0);
         m_masterColType.push_back(DecompCol_ArtForCutG);
         m_masterArtCols.push_back(colIndex);
         colName = "sCG(c_" + UtilIntToStr(colIndex)
                   + "_" + UtilIntToStr(rowIndex) + ")";
         colNames.push_back(colName);
         colIndex += 2;
      }
      break;
      default:
         assert(0);
      }

      rowIndex++;
   }

   //---
   //--- add the row names to master
   //--- add the row names to master
   //---
   if (rowNames.size() > 0)
      m_masterSI->setRowNames(rowNames, 0,
                              static_cast<int>(rowNames.size()), rowIndex0);

   if (colNames.size() > 0)
      m_masterSI->setColNames(colNames, 0,
                              static_cast<int>(colNames.size()), colIndex0);

   //---
   //--- clean up
   //---
   index = 0;

   for (li = m_cutpool.begin(); li != m_cutpool.end(); li++) {
      if (index >= n_newrows) {
         break;
      }

      (*li).deleteRowReform();
      (*li).deleteRow();
      (*li).clearCut();//need to do this?
      index++;
   }

   m_cutpool.erase(m_cutpool.begin(), li);
   UTIL_DELARR(rowReformBlock);
   UTIL_DELARR(rowBlock);
   UTIL_DELARR(rlb);
   UTIL_DELARR(rub);
#if 1
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 4,
              (*m_osLog) << "\nCUT POOL AFTER:\n";
              m_cutpool.print(m_osLog);
              (*m_osLog) << "\nCUTS AFTER:\n";
              printCuts(m_osLog);
             );
#endif
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "addCutsFromPool()", m_param.LogDebugLevel, 2);
   return n_newrows;
}


