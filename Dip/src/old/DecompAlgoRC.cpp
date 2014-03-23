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
#include "DecompAlgoRC.h"
#include "DecompPortable.h"
// ------------------------------------------------------------------------- //
const char* DecompAlgoRC::m_classTag         = "\nD-ALGORC      : ";

// ------------------------------------------------------------------------- //
void DecompAlgoRC::createMasterProblem(DecompVarList& initVars)
{
   //---
   //--- there is no master LP in RC, just initialize the dual vector
   //---
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << " <---- createMasterProblem() ---- ";
             );
   //THINK: is this the right place for this
   //TODO: give user option to feed in a good starting dual vector?
   //you might want to do that if switch between DW and RC... etc... THINK
   fill_n(back_inserter(m_u), m_modelCore->getNumRows(), 0.0);
   m_rc   = new double[m_modelCore->getNumCols()]; //better in constructor?
   //m_shat = new double[m_modelCore->getNumCols()];
   //m_vars here will contain all the shat's used - we can only
   //use one of them in subgradient, so just save the last one
   DecompVarList::iterator it = initVars.begin();
   assert(*it);
   m_vars.push_back(*it);
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*it)->print(m_osLog);
             );
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << "  ---- createMasterProblem() ----> ";
             );
}

// ------------------------------------------------------------------------- //
bool DecompAlgoRC::isDone()
{
   //iter count is checked by phaseUpdate
   //need to check step limit in here
   //need to check if ub-lb gap is small? isn't that always checked?
   printf("\nm_UB: %12.10f, m_LB: %12.10f", m_UB, m_LB);

   if ((m_step < 1.0e-3)           //step length too small
         || m_zeroSub                 //0 subgradient
         || UtilIsZero(m_UB - m_LB, 1.0e-3)) { //gap is small
      return true;
   }

   return false;
}

//do we need a var pool at all here?

// ------------------------------------------------------------------------- //
int DecompAlgoRC::addCutsFromPool()
{
   int n_newrows = DecompAlgo::addCutsFromPool();

   //faster way in STL?
   for (int r = 0; r < n_newrows; r++) {
      m_u.push_back(0.0);
   }

   //is this the right place to do this?
   //we were pricing out, then cutting, but when we price out,
   //step=0, so we need to start over with step size
   //best place for this would be in phaseUpdatE?
   if (n_newrows > 0) {
      m_step = 2.0;
   }

   return n_newrows;
}

// ------------------------------------------------------------------------- //
int DecompAlgoRC::generateVars(const decompStat   stat,
                               DecompVarList&     newVars,
                               double&            mostNegReducedCost)
{
   //really only returning one var here...
   //for RC, doesn't have to be negative??
   mostNegReducedCost = DecompInf;//bad name here
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << " <---- generateVars() ---- ";
             );
   //TODO: whenever a cut is added, if doing RC, you need to add an
   //element to u.... do we need to override gen cuts just for that?
   //gen cuts (x) is wrong here anyway... need gen cuts (s)
   assert(static_cast<int>(m_u.size()) == m_modelCore->getNumRows());
   //reduced cost = c - uA
   m_modelCore->M->transposeTimes(&m_u[0], m_rc);

   for (int c = 0; c < m_modelCore->getNumCols(); c++) {
      printf("\nRC[%d] -> c: %g - uA: %g = m_rc: %g",
             c, m_app->m_model.objCoeff[c], m_rc[c],
             m_app->m_model.objCoeff[c] - m_rc[c]);
      m_rc[c] = m_app->m_model.objCoeff[c] - m_rc[c];
   }

   double alpha = 0.0;
   DecompVarList potentialVars;
   //TODO: stat return, restrict how many? pass that in to user?
   //only take those with negative reduced cost?
   //check for dups here
   solveRelaxed(m_whichModel,
                m_rc, m_app->m_model.objCoeff, alpha,
                m_modelCore->getNumCols(), false, true, m_subprobSI[m_whichModel],
                potentialVars);//NO CHECK RC??
   DecompVarList::iterator it;
   double varRedCost;

   for (it = potentialVars.begin(); it != potentialVars.end(); it++) {
      varRedCost = (*it)->getReducedCost();
      newVars.push_back(*it);

      if (varRedCost < mostNegReducedCost) {
         mostNegReducedCost = varRedCost;
         m_shatVar = *it; //THINK
         m_shatVar->fillDenseArr(m_modelCore->getNumCols(), m_xhat);
      }

      //if(varRedCost < 0){ //TODO: strict, -dualTOL?
      //newVars.push_back(*it);
      //if(varRedCost < mostNegReducedCost)
      //mostNegReducedCost = varRedCost;
      //}
      //else{
      //UTIL_DELPTR(*it);
      //THINK: have to delete this memory?
      //}
   }

   potentialVars.clear(); //THINK? what does clear do exactly ?

   for (it = newVars.begin(); it != newVars.end(); it++)
      UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
                 (*it)->print(m_osLog);
                );

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << "  ---- generateVars() ----> ";
             );
   return static_cast<int>(newVars.size());
}

// ------------------------------------------------------------------------- //
decompStat DecompAlgoRC::solutionUpdate(const decompPhase phase,
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
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << " <---- solutionUpdate() ----\n ";
              (*m_osLog) << "\nm_iter " << m_iter;
             );
   m_UB = 1.05 * m_bestUpperBound; //TODO heuristics 1.05??
   //TODO: tols
   //how to allow hooks to other stabilization methods, bundle, etc... user
   //can simply derive from RC or base and recode this method...
   (*m_osLog) << "\nVARS m_vars:\n";
   printVars(m_osLog);
   //DecompVarList::reverse_iterator it = m_vars.rbegin();
   //DecompVar * shat = *it;
   DecompVar* shatVar = m_shatVar;
   assert(m_shatVar);
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << "\nshat: ";
              shatVar->print(m_osLog);
             );
   //make this part of class, else realloc every iter...
   //use vector, let it grow?
   int      r;
   //char     sense;
   //double   range;
   int      n_coreRows = m_modelCore->getNumRows();
   const double* rhs        = getRightHandSide();
   const char*    sense      = getRowSense();
   double* violation  = new double[n_coreRows];
   double* activity   = new double[n_coreRows];
   m_modelCore->M->times(shatVar->m_s, activity); //As
   assert(static_cast<int>(m_u.size()) == n_coreRows);
   // =, b  - Ax or Ax - b
   //>=, b  - Ax > 0 is a violation
   //<=, b  - Ax < 0 is a violation
   m_zeroSub = true;

   for (r = 0; r < n_coreRows; r++) {
      //UtilBoundToSense(m_modelCore->rowLB[r],
      //                 m_modelCore->rowUB[r],
      //                 DecompInf, sense, rhs[r], range);
      violation[r] = rhs[r] - activity[r];
      UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
                 (*m_osLog) << setprecision(8);
                 (*m_osLog) << "\nr: " << r << " vio: " << violation[r]
                 << " rhs: " << rhs[r] << " act: " << activity[r]
                 << " u: " << m_u[r] << " sense: " << sense[r];
                );
#if 1

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
         //use tol
         if (violation[r] < 0.0 && m_u[r] >= -1.0e-4 && m_u[r] <= 1.0e-4) {
            violation[r] = 0.0;
         }

         break;
         //case 'G':
      case 'L':

         //violation[i] = rhs[i] - activity[i];
         if (violation[r] > 0.0 && m_u[r] >= -1.0e-4 && m_u[r] <= 1.0e-4) {
            violation[r] = 0.0;
         }

         break;
      }

#endif
      UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
                 (*m_osLog) << " -> vio: " << violation[r];
                );

      //?? shouldn't it be if all are feasible?
      if (fabs(violation[r]) > 0.0001) {
         m_zeroSub = false;
      }
   }

   //when to half the step size?
#if 1
   //same stuff as m_tlb?? setTrueLowerBound?
   double bound    = shatVar->getReducedCost();//c - uA (is this set?)
   //double constant = calcConstant(m_u.size(), &m_u[0]);
   double constant = 0.0;

   for (r = 0; r < n_coreRows; r++) {
      constant += m_u[r] * rhs[r];
   }

   //needs rhs from OSI - this is why it was better to fake it
   //and have the information in OSI... but then carrying around
   //m_modelCore and OSI for no good reason... but... this is messier
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

   if (m_cntSameLB >= 10) {
      m_step /= 2.0;
      cout << "LB has not changed in " << m_cntSameLB
           << " iterations - halve step: " << m_step << endl;
      m_cntSameLB = 0;
   }

   printf("\nm_UB: %12.10f, m_LB: %12.10f", m_UB, m_LB);
   assert((m_UB - m_LB) > -0.0001);
   //if(count_sameLB >= m_app->m_param.RC_sameLBLimit){
   //step /= 2.0;
   //cout << "LB has not changed in " << count_sameLB
   //	   << " iterations - halve step: " << step << endl;
   //count_sameLB = 0;
   //}
#endif
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
              (*m_osLog) << "\nm_UB: " << m_UB << " m_LB: " << m_LB
              << " denom: " << denom << " m_step: " << m_step
              << " theta: " << theta;
             );

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
         //m_u[r] = min(0.0, m_u[r] + (theta * violation[r]));
         m_u[r] = max(0.0, m_u[r] - (theta * violation[r]));
         //u[i] = max(0.0, u[i] - (theta * violation[i]));
         //u[i] = max(0.0, u[i] + (theta * violation[i]));
         break;
      default:
         assert(0);
      }

      UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
                 (*m_osLog) << "\nr: " << r << " m_u: " << m_u[r];
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
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << "  ---- solutionUpdate() ----> ";
             );
   return STAT_FEASIBLE;
}
