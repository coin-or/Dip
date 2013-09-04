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


#include "DecompVarPool.h"
#include "DecompConstraintSet.h"

using namespace std;

// --------------------------------------------------------------------- //
bool DecompWaitingCol::setReducedCost(const double*       u,
                                      const DecompStatus    stat)
{
   double redCost;

   if (stat == STAT_FEASIBLE) {
      // ---
      // --- RC[s] = c[s] - u (A''s) - alpha
      // ---
      redCost = m_var->getOriginalCost() - m_col->dotProduct(u);
      m_var->setReducedCost(redCost);
      return redCost <= -0.0000000001;//m_app->m_param.dualTol;
   } else {
      // ---
      // --- RC[s] = u (A''s) + alpha -> dual ray
      // ---
      redCost = -m_col->dotProduct(u);
      return redCost <= -0.0000000001;//m_app->m_param.dualTol;
   }
}


// --------------------------------------------------------------------- //
//use hash!
/*bool DecompVarPool::isDuplicate(const DecompWaitingCol & wcol){
   vector<DecompWaitingCol>::const_iterator vi;
   for(vi = begin(); vi != end(); vi++){
      //TODO: this is very expensive
      //TODO: override DecompWaitingCol operator==
      if((*vi).getColPtr()->isEquivalent(*wcol.getColPtr()))
         return true;
   }
   return false;
   }*/


// --------------------------------------------------------------------- //
bool DecompVarPool::isParallel(const DecompVarList&     vars,
                               const DecompWaitingCol& wcol,
                               const double             maxCosine)
{
   DecompVarList::const_iterator vi;
   int            j1, j2, index1, index2;
   double         cosine;
   DecompVar*     var    = wcol.getVarPtr();
   const int      block1 = var->getBlockId();
   const int      len1   = var->m_s.getNumElements();
   const int*     ind1   = var->m_s.getIndices();
   const double* els1   = var->m_s.getElements();
   const double   norm1  = var->getNorm();
   bool           isPara = false;

   if (len1 == 0) {
      return false;
   }

   for (vi = vars.begin(); vi != vars.end(); vi++) {
      //---
      //--- if different blocks, it doesn't matter if rest of var
      //---   is close to parallel
      //---
      const int      len2 = (*vi)->m_s.getNumElements();

      if ((*vi)->getBlockId() != block1 ||
            len2                == 0) {
         continue;
      }

      const int*     ind2 = (*vi)->m_s.getIndices();
      const double* els2 = (*vi)->m_s.getElements();
      const double   norm2 = (*vi)->getNorm();
      index1              = 0;
      index2              = 0;
      cosine              = 0.0;

      //---
      //--- calculate var1*var2 (both sparse)
      //---   var indices are assumed to be sorted increasing
      //---
      while (1) {
         j1 = ind1[index1];
         j2 = ind2[index2];

         if (j1 == j2) {
            cosine += els1[index1] * els2[index2];
            index1++;
            index2++;

            if (index2 >= len2 || index1 >= len1) {
               break;
            }
         } else if (j1 > j2) {
            index2++;

            if (index2 >= len2) {
               break;
            }
         } else {
            index1++;

            if (index1 >= len1) {
               break;
            }
         }
      }

      cosine /= norm1;
      cosine /= norm2;
      cosine = fabs(cosine);

      if (cosine > maxCosine) {
         isPara = true;
         printf("parallel: cosine=%g\n", cosine);
         break;
      }

      //printf("not parallel: cosine=%g\n", cosine);
   }

   return isPara;
}

// --------------------------------------------------------------------- //
bool DecompVarPool::isDuplicate(const DecompVarList&     vars,
                                const DecompWaitingCol& wcol)
{
   DecompVarList::const_iterator vi;
   DecompVar* var = wcol.getVarPtr();

   for (vi = vars.begin(); vi != vars.end(); vi++) {
      if (((*vi)->getBlockId() == var->getBlockId()) &&
            ((*vi)->getStrHash() == var->getStrHash())) {
         return true;
      }
   }

   return false;
}

// --------------------------------------------------------------------- //
bool DecompVarPool::isDuplicate(const DecompWaitingCol& wcol)
{
   vector<DecompWaitingCol>::const_iterator vi;
   DecompVar* var1 = wcol.getVarPtr();

   for (vi = begin(); vi != end(); vi++) {
      DecompVar* var2 = (*vi).getVarPtr();

      if ((var1->getBlockId() == var2->getBlockId()) &&
            (var1->getStrHash() == var2->getStrHash())) {
         return true;
      }
   }

   return false;
}

/*-------------------------------------------------------------------------*/
bool DecompVarPool::setReducedCosts(const double*             u,
                                    const DecompStatus          stat,
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
