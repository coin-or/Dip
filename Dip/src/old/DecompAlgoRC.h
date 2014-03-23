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

/*-----------------------------------------------------------------------*/
/* Author: Matthew Galati (magh@lehigh.edu)                              */
/*                                                                       */
/* (c) Copyright 2004 Lehigh University. All Rights Reserved.            */
/*                                                                       */
/* This software is licensed under the Common Public License. Please see */
/* accompanying file for terms.                                          */
/*-----------------------------------------------------------------------*/

#ifndef DECOMP_ALGORC_INCLUDED
#define DECOMP_ALGORC_INCLUDED

#include "DecompAlgo.h"

class DecompApp;
// --------------------------------------------------------------------- //
class DecompAlgoRC : public DecompAlgo {
private:
   DecompAlgoRC(const DecompAlgoRC&);
   DecompAlgoRC& operator=(const DecompAlgoRC&);

private:
   static const char* m_classTag;

private:
   vector<double>   m_u;   //dual vector
   double*          m_rc;  //reduced cost
   double           m_UB;  //current best upper bound
   double           m_LB;  //current best lower bound

   int              m_cntSameLB;
   int              m_iter;
   double           m_step;
   bool             m_zeroSub;

   DecompVar*       m_shatVar;
   //double         * m_shat;

public:
   //inherited (from pure virtual) methods
   void createMasterProblem(DecompVarList& initVars);
   decompStat solutionUpdate(const decompPhase phase,
                             const int         maxInnerIter,
                             const int         maxOuterIter);
   //void addCutsToPool(const double  *  x,
   //                   DecompCutList & newCuts,
   //                   int           & n_newCuts) {assert(0);};
   int addCutsFromPool();
   int generateVars(const decompStat   stat,
                    DecompVarList&     newVars,
                    double&            mostNegReducedCost);

   bool isDone();

   const double* getRowPrice() const {
      return &m_u[0];
   }

public:
   DecompAlgoRC(DecompApp* app)
      : DecompAlgo(RELAX_AND_CUT, app),
        m_u(),
        m_rc(0),
        m_UB(DecompInf),
        m_LB(-DecompInf),
        m_cntSameLB(0),
        m_iter(0),
        m_step(2.0), //(0, 2] param?
        m_zeroSub(false),
        m_shatVar(0)
        //m_shat(0)
   {};
   ~DecompAlgoRC() {
      UTIL_DELARR(m_rc);
   };
};

#endif
