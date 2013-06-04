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

#ifndef DECOMP_ALGOPC_INCLUDED
#define DECOMP_ALGOPC_INCLUDED

#include "DecompAlgo.h"

class DecompApp;

//?? this should probably be a base class for different stablization methods
// --------------------------------------------------------------------- //
class DecompAlgoPC : public DecompAlgo {
private:
   DecompAlgoPC(const DecompAlgoPC&);
   DecompAlgoPC& operator=(const DecompAlgoPC&);

private:
   static const char* classTag;  //THINK - bad idea same name as base?

public:
   virtual void setMasterBounds(const double* lbs,
                                const double* ubs);

   //inherited (from pure virtual) methods
   virtual void createMasterProblem(DecompVarList& initVars);
   void createMasterStabilization();
   void addCutsToPool(const double*    x,
                      DecompCutList& newCuts,
                      int&            n_newCuts);
   int addCutsFromPool();

#if 0
   decompStat solutionUpdate(const decompPhase phase,
                             const int         maxInnerIter,
                             const int         maxOuterIter);
#endif

#if 0
   decompPhase phaseUpdate(const decompPhase   phase,
                           const decompStat    stat,
                           int&                n_newCuts,
                           int&                n_newVars,
                           int&                n_cutCalls,
                           int&                n_priceCalls);
#endif

   //inherited from virtual methods
   void recomposeSolution(const double* solution,
                          double*        rsolution);


public:
   DecompAlgoPC(DecompApp* app)
      : DecompAlgo(PRICE_AND_CUT, app) {};
   ~DecompAlgoPC() {};
};

#endif
