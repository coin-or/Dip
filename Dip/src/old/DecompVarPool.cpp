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


#include "DecompVarPool.h"
#include "DecompConstraintSet.h"

// --------------------------------------------------------------------- //
bool DecompWaitingCol::setReducedCost(const double*       u,
                                      const decompStat    stat)
{
   double redCost;

   if (stat == STAT_FEASIBLE) {
      // ---
      // --- RC[s] = c[s] - u (A''s) - alpha
      // ---
      redCost = m_var->getOriginalCost() - m_col->dotProduct(u);
      //printf("\nHERE DecompWaitingCol::setReducedCost redCost: %g",
      //       redCost);
      m_var->setReducedCost(redCost);
      return redCost <= -0.0000000001;//m_app->m_param.dualTol;
   } else {
      // ---
      // --- RC[s] = u (A''s) + alpha -> dual ray
      // ---
      redCost = m_col->dotProduct(u);
      //redCost = -m_col->dotProduct(u);
      return redCost <= -0.0000000001;//m_app->m_param.dualTol;
   }
}


// --------------------------------------------------------------------- //
bool DecompVarPool::isDuplicate(const DecompWaitingCol& wcol)
{
   vector<DecompWaitingCol>::const_iterator vi;

   for (vi = begin(); vi != end(); vi++) {
      //TODO: this is very expensive
      //TODO: override DecompWaitingCol operator==
      //printf("\nHERE isDup");
      if ((*vi).getColPtr()->isEquivalent(*wcol.getColPtr())) {
         return true;
      }
   }

   return false;
}

#if 0
// --------------------------------------------------------------------- //
bool DecompVarPool::isDuplicate(const DecompVarList&     vars,
                                const DecompWaitingCol& wcol)
{
   DecompVarList::const_iterator vi;

   for (vi = vars.begin(); vi != vars.end(); vi++) {
      //TODO: this is very expensive
      //TODO: override DecompWaitingCol operator==
      //printf("\nHERE isDup");
      if ((*vi)->isEquivalent(*wcol.getVarPtr())) { //checks if s is equivalent
         return true;
      }
   }

   return false;
}
#else
// --------------------------------------------------------------------- //
bool DecompVarPool::isDuplicate(const DecompVarList&     vars,
                                const DecompWaitingCol& wcol)
{
   DecompVarList::const_iterator vi;

   for (vi = vars.begin(); vi != vars.end(); vi++) {
      if ((*vi)->getStrHash() == wcol.getVarPtr()->getStrHash()) {
         return true;
      }
   }

   return false;
}
#endif

/*-------------------------------------------------------------------------*/
bool DecompVarPool::setReducedCosts(const double*             u,
                                    const decompStat          stat,
                                    DecompVarPool::iterator   first,
                                    DecompVarPool::iterator   last)
{
   //printf("\nHERE DecompVarPool::setReducedCosts");
   bool found_negrc_var = false;

   for (DecompVarPool::iterator vi = first; vi != last; vi++) {
      // ---
      // --- calculate and set the reduced costs for the variables
      // --- which are pointed to in this pool, if any have rc < 0,
      // --- return true
      // ---
      found_negrc_var = (*vi).setReducedCost(u, stat) ? true : found_negrc_var;
   }

   return found_negrc_var;
}

// --------------------------------------------------------------------- //
//THINK: this is specific to PC and DC??
void DecompVarPool::reExpand(const DecompConstraintSet& modelCore,
                             const double                tolZero)
{
   //THIS IS WRONG...
   //in masterSI, we have
   //A'', convexity, cuts
   //in modelCore.M we have A'', cuts
   //the sparseCol that you come out with the end here has things in the wrong
   //order: //A'', cuts, convexity
   double* denseCol = new double[modelCore.getNumRows() + 1];
   vector<DecompWaitingCol>::iterator vi;

   for (vi = begin(); vi != end(); vi++) {
      // ---
      // --- get dense column = A''s, append convexity constraint on end
      // ---
      modelCore.M->times((*vi).getVarPtr()->m_s, denseCol);
      denseCol[modelCore.getNumRows()] = 1.0;
      // ---
      // --- create a sparse column from the dense column
      // ---
      CoinPackedVector* sparseCol
      = UtilPackedVectorFromDense(modelCore.getNumRows() + 1,
                                  denseCol, tolZero);
      (*vi).deleteCol();
      (*vi).setCol(sparseCol);
   }

   setColsAreValid(true);
   UTIL_DELARR(denseCol);
}

// --------------------------------------------------------------------- //
void DecompVarPool::print(ostream* os) const
{
   vector<DecompWaitingCol>::const_iterator vi;

   for (vi = begin(); vi != end(); vi++) {
      (*vi).getVarPtr()->print(os);
   }
}
