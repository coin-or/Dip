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
#include "DecompAlgoD.h"
#include "DecompPortable.h"

#include "DecompCutOsi.h"

// ------------------------------------------------------------------------- //
const char* DecompAlgoD::classTag         = "\nD-ALGOD       : ";

//generateInitVars should be based on cost = -xhat? TODO

// ------------------------------------------------------------------------- //
void DecompAlgoD::createMasterProblem(DecompVarList& initVars)
{
   //---
   //--- Initialize the solver interface for the master problem.
   //--- PC: min 0(s lam)
   //---        s   lam    = xhat,
   //---     sum{s} lam_s  = 1   ,
   //---            lam_s >= 0   , s in F'[0]
   //---
   //CONSISTENCY???? look at when do TSP
   //NOTE: modelCore is NOT used in DECOMP!! weird....
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << classTag << " <---- createMasterProblem() ---- ";
             );
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,

   for (int i = 0; i < m_numOrigCols; i++) {
   if (!UtilIsZero(m_xhatD[i])) {
         (*m_osLog) << "\nxhatD[" << i << "] : " << m_xhatD[i];
      }
   }
             );
   //THINK? should initVars be in pool and do addVarsFromPool here?
   CoinPackedMatrix* M = new CoinPackedMatrix(true, 0, 0);
   M->setDimensions(m_numOrigCols + 1, 0);
   //m_nBaseCoreRows = m_nBaseCoreRows;
   //THINK? cutpool contains cut and expanded row - is that correct?
   //there should be a function that "addsVars" but from an argument sent
   //in... this time directly, later from the "pool"??? THINK??
   const int n_cols = static_cast<int>(initVars.size());
   double* colLB = new double[n_cols];
   double* colUB = new double[n_cols];
   double* obj   = new double[n_cols];
   double* denseCol = new double[m_numOrigCols + 1];
   int col_index     = 0;
   DecompVarList::iterator li;

   for (li = initVars.begin(); li != initVars.end(); li++) {
      UTIL_DEBUG(m_app->m_param.LogDebugLevel, 5,
                 (*li)->print(m_osLog);
                );
#if 0
      // ---
      // --- get dense column = A''s, append convexity constraint on end
      // ---
      m_modelCore.M->times((*li)->m_s, denseCol);
      denseCol[m_modelCore.getNumRows()] = 1.0;
      // ---
      // --- create a sparse column from the dense column
      // ---
      // THINK: do i need a DecompCol?
      // THINK: does this allocate memory for coinpackedvec twice?
      CoinPackedVector* sparseCol
      = UtilPackedVectorFromDense(m_modelCore.getNumRows() + 1,
                                  denseCol, m_app->m_param.TolZero);
      UTIL_DEBUG(m_app->m_param.LogDebugLevel, 5,
                 (*m_osLog) << "\nSparse Col: \n";
                 UtilPrintPackedVector(*sparseCol, m_osLog);
                );
      //
      // --- check for duplicates (against m_vars)
      // ?? or force initVars to be sent in with no dups?
      //
#endif
      CoinPackedVector* sparseCol = new CoinPackedVector((*li)->m_s);
      sparseCol->insert(m_numOrigCols, 1.0);//convexity
      //TODO: do in const blocks
      // ---
      // --- append the sparse column to the matrix
      // ---
      M->appendCol(*sparseCol);
      colLB[col_index] = 0.0; //THINK: (*li)->getLowerBound();
      colUB[col_index] = 1.0; //THINK: (*li)->getUpperBound();
      obj[col_index]   = 0.0; //(*li)->getOriginalCost(); //c.s
      col_index++;
      UTIL_DELPTR(sparseCol); //THINK
      //m_vars.push_back(*li);
   }

   //---
   //--- insert the initial set of variables into the master variable list
   //---
   //THINK: now doing in loop, so can check for dups
   m_vars.insert(m_vars.end(), initVars.begin(), initVars.end());
   //---
   //--- THINK: do we want to adjust m_modelCore directly here?
   //--- adjust row bounds for convexity constraint
   //---
   //WHY WOULD YOU DO THIS HERE? bad idea...
   //m_modelCore.rowLB.push_back(1.0);
   //m_modelCore.rowUB.push_back(1.0);
   vector<double> masterRowBound(m_xhatD, m_xhatD + m_numOrigCols);
   masterRowBound.push_back(1.0);
   printf("\nmasterRowBound.size() = %d", masterRowBound.size());
   printf("\nM->getNumRows()       = %d", M->getNumRows());
   //---
   //--- load the problem into master's solver interface
   //---
   assert(M->getNumRows() == static_cast<int>(masterRowBound.size()));
   assert(M->getNumRows() == static_cast<int>(masterRowBound.size()));
   m_masterSI->loadProblem(*M,  colLB, colUB, obj,
                           &masterRowBound[0],
                           &masterRowBound[0]);
   // ---
   // --- free local memory
   // ---
   //init_vars.clear(); //WHY?
   UTIL_DELPTR(M);
   UTIL_DELARR(denseCol);
   UTIL_DELARR(colLB);
   UTIL_DELARR(colUB);
   UTIL_DELARR(obj);
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << classTag << "  ---- createMasterProblem() ----> ";
             );
}


// ------------------------------------------------------------------------ //
decompPhase DecompAlgoD::phaseUpdate(const decompPhase   phase,
                                     const decompStat    stat,
                                     int&                n_newCuts,
                                     int&                n_newVars,
                                     int&                n_cutCalls,
                                     int&                n_priceCalls)
{
   decompPhase nextPhase = PHASE_UNKNOWN;
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << classTag << " <---- phaseUpdate() ---- ";
             );
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << "\nn_newCuts: "   << n_newCuts;
              (*m_osLog) << " n_newVars: "    << n_newVars;
              (*m_osLog) << " n_cutCalls: "   << n_cutCalls;
              (*m_osLog) << " n_priceCalls: " << n_priceCalls;
              (*m_osLog) << "\nPHASEIN : "    << decompPhaseStr[phase] << "\t";
              (*m_osLog) << "STAT IN: "       << decompStatStr[stat];
             );
   //stopping criteria - different override for each algo

   switch (phase) {
   case PHASE_INIT:

      if (stat == STAT_FEASIBLE) {
         nextPhase = PHASE_DONE;
      } else {
         nextPhase = PHASE_PRICE;
      }

      break;
   case PHASE_PRICE:

      //---
      //--- we are pricing
      //---
      if (stat == STAT_FEASIBLE) {
         nextPhase = PHASE_DONE;
      } else {
         //TODO: price iter limit
         //---
         //--- master is infeasible
         //---   if we find any neg redCost columns -> price again
         //---   if we find no  neg redCost columns -> produce farkas cut
         //---
         if (n_newVars > 0) {
            nextPhase = PHASE_PRICE;
         } else {
            //produce a farkas cut
            double* u = getDualRays(1)[0];
            CoinPackedVector cut;
            double lhs = 0.0;

            for (int i = 0; i < m_numOrigCols; i++) {
               cut.insert(i, u[i]);
               lhs += u[i] * m_xhatD[i];
            }

            double beta = -u[m_numOrigCols];
            OsiRowCut rc;
            rc.setRow(cut);
            rc.setLb(beta);
            rc.setUb(DecompInf);
            printf("\nFARKAS CUT lhs: %g, rhs: %g",
                   lhs, beta);
            assert((lhs - beta) < -0.00001);
            DecompCutOsi* decompCut = new DecompCutOsi(rc);
            decompCut->setStringHash();//TEST
            UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
                       decompCut->print(m_osLog);
                      );
            (*m_newCuts).push_back(decompCut);
            nextPhase = PHASE_DONE;
            //experiment... try lifting - random order?
            //to do so, we must have a way to tell subproblem solver
            //that some columns are fixed...
            //#define SMALL_IP
#ifdef SMALL_IP
            //- 0.22 x[0] + 0.06 x[1] <= -0.84
            //fix x[1] = 0, solve relaxed with c = (-0.22, 0.06)
            m_subprobSI->setColBounds(1, 0.0, 0.0);
            //what happens if this causes INF??
            double c[2];
            //(-u)x <= -beta=alpha
            //max -u -> min u
            c[0] = u[0];
            c[1] = u[1];
            DecompVarList potentialVars;
            m_app->solveRelaxed(c,
                                m_modelCore.objCoeff,
                                u[2],//alpha
                                2,
                                false,
                                false,
                                m_subprobSI,
                                potentialVars);
            //gamma =
            DecompVarList::iterator it;
            double gamma;

            for (it = potentialVars.begin(); it != potentialVars.end(); it++) {
               gamma = (*it)->getReducedCost(); //c*u  - alpha
            }

            gamma += u[2];
            gamma = -gamma; //max -u
            double liftCoeff = (-beta) - gamma;
            printf("\nliftCoeff = %g", liftCoeff);
            //set it back
            m_subprobSI->setColBounds(1, 0.0, 1.0);
#endif
         }
      }

      break;
   case PHASE_DONE:
      break;
   case PHASE_UNKNOWN:
   default:
      assert(0);//TODO: real error codes?
   }

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << "\nPHASEOUT: "
              << decompPhaseStr[nextPhase] << "\t";
             );
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << classTag << "  ---- phaseUpdate() ----> ";
             );
   return nextPhase;
}


// ------------------------------------------------------------------------ //
decompPhase DecompAlgoD::phaseUpdate(const decompPhase   phase,
                                     const decompStat    stat)
{
   decompPhase nextPhase = PHASE_UNKNOWN;
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << classTag << " <---- phaseUpdate() ---- ";
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
   //stopping criteria - different override for each algo

   switch (phase) {
   case PHASE_INIT:

      if (stat == STAT_FEASIBLE) {
         nextPhase = PHASE_DONE;
      } else {
         nextPhase = PHASE_PRICE;
      }

      break;
   case PHASE_PRICE:

      //---
      //--- we are pricing
      //---
      if (stat == STAT_FEASIBLE) {
         nextPhase = PHASE_DONE;
      } else {
         //TODO: price iter limit
         //---
         //--- master is infeasible
         //---   if we find any neg redCost columns -> price again
         //---   if we find no  neg redCost columns -> produce farkas cut
         //---
         if (m_varsThisCall > 0) {
            nextPhase = PHASE_PRICE;
         } else {
            //produce a farkas cut
            double* u = getDualRays(1)[0];
            CoinPackedVector cut;
            double lhs = 0.0;

            for (int i = 0; i < m_numOrigCols; i++) {
               cut.insert(i, u[i]);
               lhs += u[i] * m_xhatD[i];
            }

            double beta = -u[m_numOrigCols];
            OsiRowCut rc;
            rc.setRow(cut);
            rc.setLb(beta);
            rc.setUb(DecompInf);
            printf("\nFARKAS CUT lhs: %g, rhs: %g",
                   lhs, beta);
            assert((lhs - beta) < -0.00001);
            DecompCutOsi* decompCut = new DecompCutOsi(rc);
            decompCut->setStringHash();//TEST
            UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
                       decompCut->print(m_osLog);
                      );
            (*m_newCuts).push_back(decompCut);
            nextPhase = PHASE_DONE;
            //experiment... try lifting - random order?
            //to do so, we must have a way to tell subproblem solver
            //that some columns are fixed...
            //#define SMALL_IP
#ifdef SMALL_IP
            //- 0.22 x[0] + 0.06 x[1] <= -0.84
            //fix x[1] = 0, solve relaxed with c = (-0.22, 0.06)
            m_subprobSI->setColBounds(1, 0.0, 0.0);
            //what happens if this causes INF??
            double c[2];
            //(-u)x <= -beta=alpha
            //max -u -> min u
            c[0] = u[0];
            c[1] = u[1];
            DecompVarList potentialVars;
            m_app->solveRelaxed(c,
                                m_modelCore.objCoeff,
                                u[2],//alpha
                                2,
                                false,
                                false,
                                m_subprobSI,
                                potentialVars);
            //gamma =
            DecompVarList::iterator it;
            double gamma;

            for (it = potentialVars.begin(); it != potentialVars.end(); it++) {
               gamma = (*it)->getReducedCost(); //c*u  - alpha
            }

            gamma += u[2];
            gamma = -gamma; //max -u
            double liftCoeff = (-beta) - gamma;
            printf("\nliftCoeff = %g", liftCoeff);
            //set it back
            m_subprobSI->setColBounds(1, 0.0, 1.0);
#endif
         }
      }

      break;
   case PHASE_DONE:
      break;
   case PHASE_UNKNOWN:
   default:
      assert(0);//TODO: real error codes?
   }

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << "\nPHASEOUT: "
              << decompPhaseStr[nextPhase] << "\t";
             );
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << classTag << "  ---- phaseUpdate() ----> ";
             );
   return nextPhase;
}

/*---------------------------------------------------------------------------*/
//PC only? think about name?
void DecompAlgoD::recomposeSolution(const double* solution,
                                    double*        rsolution)
{
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << classTag << " <---- recomposeSolution()   ---- ";
             );
   UtilFillN(rsolution, m_modelCore->getNumCols(), 0.0);
   //m_setD.clear();
   int col_index = 0;
   DecompVarList::const_iterator li;

   for (li = m_vars.begin(); li != m_vars.end(); li++) {
      if (solution[col_index] > m_app->m_param.TolZero) {
         CoinPackedVector& v = (*li)->m_s;
         const int*     inds  = v.getIndices();
         const double* els   = v.getElements();

         for (int i = 0; i < v.getNumElements(); i++) {
            rsolution[inds[i]] += els[i] * solution[col_index];
         }

         //m_setD.push_back(*li);
      }

      col_index++;
   }

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << classTag << "  ---- recomposeSolution()   ---->";
             );
}
