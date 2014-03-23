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


#include "DecompVar.h"
#include "DecompCutPool.h"
#include "DecompConstraintSet.h"

#if 0
// --------------------------------------------------------------------- //
bool DecompWaitingRow::calcViolation(const double* x)
{
   //always calculated wrt to original row!
   const double activity = m_row->dotProduct(x);
   double violation      = std::max<double>(m_cut->getLowerBound() - activity,
                           activity - m_cut->getUpperBound());
   violation = std::max<double>(0.0, violation);
   m_cut->setViolation(violation);
   return violation > 0.0000001;//param?
};
#endif


#if 0
// --------------------------------------------------------------------- //
bool DecompCutPool::isDuplicate(const DecompWaitingRow& wcol)
{
   return false;
   //THINK - same idea as cols for now?
#if 0
   vector<DecompWaitingRow>::const_iterator vi;

   for (vi = begin(); vi != end(); vi++) {
      //TODO: this is very expensive
      //TODO: override DecompWaitingRow operator==
      printf("\nHERE isDup");

      if ((*vi).getColPtr()->isEquivalent(*wcol.getColPtr())) {
         return true;
      }
   }

   return false;
#endif
}
#endif

/*-------------------------------------------------------------------------*/
bool DecompCutPool::calcViolations(const double*             x,
                                   DecompCutPool::iterator   first,
                                   DecompCutPool::iterator   last)
{
   bool found_violated_cut = false;

   for (DecompCutPool::iterator vi = first; vi != last; vi++) {
      // ---
      // --- calculate and set the violations for the cuts
      // --- which are pointed to in this pool, if any have vio > 0,
      // --- return true
      // ---
      found_violated_cut
      = (*vi).getCutPtr()->calcViolation((*vi).getRowPtr(), x) ?
        true : found_violated_cut;
   }

   return found_violated_cut;
}

/*-------------------------------------------------------------------------*/
void DecompCutPool::reExpand(const DecompVarList& vars,
                             const int             n_corecols)
{
   //---
   //--- For each waiting row in the cut pool, we need to reset
   //--- the row in the current master LP (in terms of reformulation)
   //--- to take into account any new columns.
   //---
   DecompCutPool::iterator vi;

   for (vi = begin(); vi != end(); vi++) {
      //only need to do this reformulation in PC...
      //make this re-expansion a function? - also called in addCutsToPool
      CoinPackedVector* rowReform = createRowReform(n_corecols,
                                    (*vi).getRowPtr(),
                                    vars);

      if (!rowReform) {
         assert(0);
         vi = erase(vi);//THINK...
      } else {
         (*vi).deleteRowReform();
         (*vi).setRowReform(rowReform);
      }

      //THINK: once we reexpand, do we need to reset violation?
      //no, this is done next.... this section just "re-expands"
   }

   setRowsAreValid(true);
}

/*------------------------------------------------------------------------*/
CoinPackedVector *
DecompCutPool::createRowReform(const int                  n_corecols,
                               const CoinPackedVector* row,
                               const DecompVarList&       vars)
{
   double coeff;
   //---
   //--- Create a dense row from the original sparse row (in terms of x).
   //---
   double* rowDense  = row->denseVector(n_corecols);
   //---
   //--- In order to expand to the reformulated row (in terms of lambda),
   //--- we need to substitute x = sum{s in F'} s lambda[s]
   //---
   //--- Example - Given a cut:
   //---   a[1]x[1] + a[2]x[2] >= b
   //---   a[1]x[1]             = a[1] (s1[1] lam[1] + s2[1] lam[2])
   //---              a[2]x[2]  = a[2] (s1[2] lam[1] + s2[2] lam[2])
   //--- So,   lam[1]'s coeff   = a[1] s1[1] + a[2] s1[2]
   //---       lam[2]'s coeff   = a[1] s2[1] + a[2] s2[2]
   //---
   int col_index                  = 0;
   CoinPackedVector* rowReform = new CoinPackedVector();
   DecompVarList::const_iterator vli;

   for (vli = vars.begin(); vli != vars.end(); vli++) {
      coeff = (*vli)->m_s.dotProduct(rowDense);

      if (fabs(coeff) > 1.0e-12) { //m_app->m_param.tolerance)
         rowReform->insert(col_index, coeff);
      }

      col_index++;
   }

   if (rowReform->getNumElements() == 0) {
      cerr << "row size 0, don't add it" << endl;
      UTIL_DELPTR(rowReform);
      rowReform = 0;
   }

   UTIL_DELARR(rowDense);
   //---
   //--- delete the temporary memory
   //---
   UTIL_DELARR(rowDense);
   return rowReform;
}

// --------------------------------------------------------------------- //
void DecompCutPool::print(ostream* os) const
{
   vector<DecompWaitingRow>::const_iterator vi;

   for (vi = begin(); vi != end(); vi++) {
      (*vi).getCutPtr()->print(os);
   }
}
