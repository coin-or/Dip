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


//for testing, simulate a few simple steps of branch and bound
//the real control for this will be in ALPs or BiCeps
#include "DecompAlgo.h"
#include "DecompApp.h"
// --------------------------------------------------------------------- //
int DecompAlgo::chooseBranchVar(int&     branchedOnIndex,
                                double& branchedOnValue)
{
   //choose variables farthest from integer - based on x formulation
   vector<int>::iterator intIt;
   int    j;
   double x, dist, maxDist;
   //weird - integerVars are defined in Relax, not Core - ugh
   //should be for overall problem
   double obj = 0.0;
   printf("\nm_masterSI OBJ: %g",
          m_masterSI->getObjValue());
   const double* objCoeff = m_app->m_model.objCoeff; //only true for CPM
   maxDist         = 1.0e-6;//TODO: integer tolerance
   branchedOnIndex = -1;
   branchedOnValue =  0;

   for (intIt =  m_modelRelax->integerVars.begin();
         intIt != m_modelRelax->integerVars.end(); intIt++) {
      j   = *intIt;
      x   = m_xhat[j];
      obj += m_xhat[j] * objCoeff[j];
      dist = fabs(x - floor(x + 0.5));

      if (dist > 1.0e-10) {
         printf("\nx[%d]: %g, dist: %g", j, x, dist);
      }

      if (dist > maxDist) {
         maxDist         = dist;
         branchedOnIndex = j;
         branchedOnValue = x;
      }
   }

   printf("\nbranchedOnIndex = %d", branchedOnIndex);
   printf("\nbranchedOnValue = %g\n", branchedOnValue);
   return 0;
}

#if 0
// --------------------------------------------------------------------- //
//just to learn what is neccessary...
int DecompAlgo::branchAndBoundSimulate()
{
   int    branchedOnIndex;
   double branchedOnValue;
   //from root node
   chooseBranchVar(branchedOnIndex,
                   branchedOnValue);

   if (branchedOnIndex == -1) {
      return -1;
   }

   return 0;
}
#endif

// --------------------------------------------------------------------- //
int DecompAlgo::branch(int    branchedOnIndex,
                       double branchedOnValue)
{
   printf("\ndefault does nothing");
   //CPM only for now
   //branch down and branch up
   //will now be important to recognize infeasibility quickly...
   //we are going to branch, which means we are going to change the
   //column bounds of some x variable... how does effect the
   //branching for CPM is straightforward
   //branching for DW subproblem must take into account
   //and lambda too in master -> how enforce?
   //for PC, how does the branching generalize?
   //for 0-1, it is simple..
   //HOW DOES THIS WORK IN GENERAL?
   //branch on x[i] = 0
   //set lambda[s] = 0, for all s such that s[i] = 1
   //branch on x[i] = 0
   //set lambda[s] = 1, for all s such that s[i] = 0
   //branch method will be different for each algo - for sure
   return 0;
}
