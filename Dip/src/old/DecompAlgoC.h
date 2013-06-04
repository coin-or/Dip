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


#ifndef DECOMP_ALGOC_INCLUDED
#define DECOMP_ALGOC_INCLUDED

#include "DecompAlgo.h"

class DecompApp;
// --------------------------------------------------------------------- //
class DecompAlgoC : public DecompAlgo {
private:
   DecompAlgoC(const DecompAlgoC&);
   DecompAlgoC& operator=(const DecompAlgoC&);

private:
   static const char* m_classTag;

public:
   //inherited (from pure virtual) methods
   void createMasterProblem(DecompVarList& initVars);
   void recomposeSolution(const double* solution,
                          double*        rsolution);
   int  generateInitVars(DecompVarList& initVars) {
      printf("\ncut generateInitVars do nothing");
      return 0;
   }

   int generateVars(const decompStat   stat,
                    DecompVarList&     newVars,
                    double&            mostNegReducedCost) {
      assert(0);
      return 0;
   }
   decompPhase phaseUpdate(const decompPhase phase,
                           const decompStat  stat);
#if 0
   decompStat solutionUpdate(const decompPhase phase,
                             const int         maxInnerIter,
                             const int         maxOuterIter);

   decompPhase phaseUpdate(const decompPhase   phase,
                           const decompStat    stat,
                           int&                n_newCuts,
                           int&                n_newVars,
                           int&                n_cutCalls,
                           int&                n_priceCalls);
#endif

   //from virtual
   int branch(int    branchedOnIndex,
              double branchedOnValue);

public:
   //THINK: really DecompAlgoC should override param PriceIter = 0
   //because if user accidently does not do this, will cause errors
   DecompAlgoC(DecompApp* app)
      : DecompAlgo(CUT, app) {};
   ~DecompAlgoC() {};
};

#endif
