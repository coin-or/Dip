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



// --------------------------------------------------------------------- //
const char* DecompApp::m_classTag         = "\nD-APP         : ";

// --------------------------------------------------------------------- //
void DecompApp::startupLog()
{
   if (m_param.LogLevel > 0) {
      (*m_osLog)
            << "\n===========================================================\n"
            <<   "===========================================================\n"
            <<   "Welcome to DECOMP v" << DecompVersion << ".\n"
            <<   "  Copyright (c) 2004-2007 Lehigh University.\n"\
            <<   "===========================================================\n"
            <<   "===========================================================\n";
      m_param.dumpSettings(m_osLog);
   }
}

// --------------------------------------------------------------------- //
int DecompApp::createModel()
{
   UTIL_DEBUG(m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << " <---- createModel() ---- ";
             );
   //---
   //--- APP: create the user model A = A' union A''
   //---
   APPcreateModel(m_model.objCoeff, m_modelCore, m_modelRelax);
   //---
   //--- TODO: sanity checks that the user gave you all the model
   //--- information that is required, error codes if not
   //---
   //---
   //--- For each model core:
   //---   1.) set row senses and/or bounds (for relaxed too)
   //---   2.) create row hash
   //---   3.) set nBaseRows
   //---   4.) flip to row ordered, if neccessary (for relaxed too)
   //--- TODO: put timer on this ??
   //---
   DecompConstraintSet* modelCore = 0;
   map<int, DecompConstraintSet*>::iterator mdc;

   for (mdc = m_modelCore.begin(); mdc != m_modelCore.end(); mdc++) {
      modelCore = mdc->second;

      if (!modelCore->M) {
         continue;
      }

      if (modelCore->M->isColOrdered()) {
         modelCore->M->reverseOrdering();
      }

      modelCore->checkSenseAndBound();
      modelCore->createRowHash();
      modelCore->nBaseRows = modelCore->getNumRows();
   }

   DecompConstraintSet* modelRelax = 0;

   for (mdc = m_modelRelax.begin(); mdc != m_modelRelax.end(); mdc++) {
      modelRelax = mdc->second;

      if (!modelRelax->M) {
         continue;
      }

      if (modelRelax->M->isColOrdered()) {
         modelRelax->M->reverseOrdering();
      }

      modelRelax->checkSenseAndBound();
   }

   UTIL_DEBUG(m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << "  ---- createModel() ----> ";
             );
   return 0;//TODO: do return codes
}

// --------------------------------------------------------------------- //
int DecompApp::generateInitVars(DecompVarList& initVars,
                                int whichModel)
{
   //
   // ---
   // --- this function does nothing by default
   // ---
   //
   UTIL_DEBUG(m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << " <---- generateInitVars() ---- ";
             );
   UTIL_DEBUG(m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << "  ---- generateInitVars() ----> ";
             );
   return 0;
}

// --------------------------------------------------------------------- //
int DecompApp::generateCuts(const double*               x,
                            const DecompConstraintSet& modelCore,
                            const DecompConstraintSet& modelRelax,
                            DecompCutList&              newCuts)
{
   //
   // ---
   // --- this function does nothing by default
   // ---
   //
   UTIL_DEBUG(m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << " <---- generateCuts()     ---- ";
             );
   UTIL_DEBUG(m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << "  ---- generateCuts()     ----> ";
             );
   return 0;
}

/*-------------------------------------------------------------------------*/
void DecompApp::printOriginalColumn(const int index,
                                    ostream* os) const
{
   (*os) << index << " ";
}

/*-------------------------------------------------------------------------*/
void DecompApp::printOriginalSolution(const int      n_cols,
                                      const double* solution,
                                      ostream*       os) const
{
   (*os) << "\n";

   for (int i = 0; i < n_cols; i++) {
      if (!UtilIsZero(solution[i])) {
         printOriginalColumn(i, os);
         (*os) << " -> " << solution[i] << endl;
      }
   }
}



#if 0
// --------------------------------------------------------------------- //
void DecompApp::setOptimalPoint(vector< vector<double> > & optPoint)
{
   //
   // ---
   // --- this function does nothing by default
   // ---
   //
   UTIL_DEBUG(m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << " <---- setOptimalPoint()  ---- ";
             );
   UTIL_DEBUG(m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << "  ---- setOptimalPoint()  ----> ";
             );
}
#endif
