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

#ifndef DECOMP_ALGOD_INCLUDED
#define DECOMP_ALGOD_INCLUDED

#include "DecompAlgoPC.h"

class DecompApp;

//?? this should probably be a base class for different stablization methods
// --------------------------------------------------------------------- //
//derive from DecompAlgo or DecompAlgoPC?? THINK
class DecompAlgoD : public DecompAlgo {
private:
   DecompAlgoD(const DecompAlgoD&);
   DecompAlgoD& operator=(const DecompAlgoD&);

private:
   static const char* classTag;
   double*             m_xhatD;
   DecompCutList*      m_newCuts;

public:
   void solveD(DecompCutList* newCuts) {
      m_newCuts = newCuts;
      //need to change parameters to price, no cut
      m_app->m_param.LimitTotalCutIters = 0;
      m_app->m_param.LimitRoundCutIters = 0;
      m_app->m_param.LimitTotalPriceIters = 1000;
      m_app->m_param.LimitRoundPriceIters = 1000;
      DecompAlgo::solve();
   }

   //inherited (from pure virtual) methods
   void createMasterProblem(DecompVarList& initVars);

   decompPhase phaseUpdate(const decompPhase   phase,
                           const decompStat    stat,
                           int&                n_newCuts,
                           int&                n_newVars,
                           int&                n_cutCalls,
                           int&                n_priceCalls);
   decompPhase phaseUpdate(const decompPhase   phase,
                           const decompStat    stat);

   //=PO
   void recomposeSolution(const double* solution,
                          double*        rsolution);


public:
   //can pass this all in... to save time, and since we might
   //do francois idea of tightening P' as we go...
   //              DecompConstraintSet * m_modelRelax)
   DecompAlgoD(DecompApp*            app,
               double*               xhat,
               int                   numOrigCols) :
      DecompAlgo(DECOMP, app),
      m_newCuts(0) {
      //that would change them all
      //m_app->m_decompParam.CutIters = 0;
      m_xhatD = xhat;
      m_numOrigCols = numOrigCols;
   };
   ~DecompAlgoD() {};
};

#endif
