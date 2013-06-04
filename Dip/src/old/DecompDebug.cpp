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
#include "DecompAlgo.h"

#if 1
// --------------------------------------------------------------------- //
bool DecompAlgo::isDualRayInfProof(const double*            dualRay,
                                   const CoinPackedMatrix* rowMatrix,
                                   const double*            colLB,
                                   const double*            colUB,
                                   const double*            rowRhs,
                                   ostream*                 os)
{
   //---
   //--- Does dualRay provide a proof according to Farkas Lemma?
   //---    yA >= 0, yb < 0, or
   //---    yA <= 0, yb > 0 ??
   //---
   int      i;
   double   yb;
   bool     isProof = true;
   bool     ybPos   = true;
   double* yA = 0;
   const int m = rowMatrix->getNumRows();
   const int n = rowMatrix->getNumCols();
   //y^T b
   yb = 0.0;

   for (i = 0; i < m; i++) {
      yb += dualRay[i] * rowRhs[i];

      if (os) {
         (*os) << "\ni : " << i << " dualRay = " << dualRay[i]
               << " rowRhs = " << rowRhs[i] << " yb = " << yb;
      }
   }

   //TODO: tol
   if (yb > 1.0e-10) {
      ybPos = true;
   } else if (yb < -1.0e-10) {
      ybPos = false;
   } else {
      return isProof;
   }

   yA = new double[n];
   rowMatrix->transposeTimes(dualRay, yA);     //y^T A

   for (i = 0; i < n; i++) {
      if (os) {
         (*os) << "\nyA[" << i << "]:\t" << yA[i];
      }

      if (ybPos  && (yA[i] >  1.0e-6) ||
            !ybPos && (yA[i] < -1.0e-6)) {
         if (os) {
            (*os) << " -->isProof (false)";
         }

         isProof = false;
      }
   }

   UTIL_DELARR(yA);
   return isProof;
}

#else

// --------------------------------------------------------------------- //
bool DecompAlgo::isDualRayInfProof(const double*            dualRay,
                                   const CoinPackedMatrix* rowMatrix,
                                   const double*            colLB,
                                   const double*            colUB,
                                   const double*            rowRhs,
                                   ostream*                 os)
{
   //---
   //--- Assume:
   //---   Ax     >= b
   //---   y^T Ax >= y^T b, y >= 0 (for >=)
   //---
   //--- Let z[j] = u[j], if y^T A[j] > 0
   //---          = l[j], if y^T A[j] < 0
   //---          = arbitrary, otherwise
   //---
   //--- Then, WHY?
   //---   y^T b - y^T A z > 0 ==> contradiction
   //---
   //---  proof_p = y^T b - y^T A z
   //---
   int       i, j;
   double   yb, yAz;
   double* yA = 0;
   double*   z = 0;
   const int m = rowMatrix->getNumRows();
   const int n = rowMatrix->getNumCols();
   //TODO: check for out-of-mem conditions?
   yA = new double[n];
   rowMatrix->transposeTimes(dualRay, yA);     //y^T A
   //the ray is probably based on presolved problem
   //colUB is not really infinity, its 1.0
   z  = new double[n];

   for (j = 0; j < n; j++) {
      //if(yA[j] >= 0) z[j] = CoinMin(1.0e20, colUB[j]);
      //else           z[j] = colLB[j];
      if (yA[j] >= 0) {
         z[j] = CoinMin(1.0, colUB[j]);
      } else {
         z[j] = colLB[j];
      }
   }

   //y^T b
   yb = 0.0;

   for (i = 0; i < m; i++) {
      yb += dualRay[i] * rowRhs[i];
      (*os) << "\ni : " << i << " dualRay = " << dualRay[i]
            << " rowRhs = " << rowRhs[i] << " yb = " << yb;
   }

   //y^T A z
   yAz = 0.0;

   for (j = 0; j < n; j++) {
      yAz += yA[j] * z[j];
      (*os) << "\nj : " << j << " yA = " << yA[j]
            << " z = " << z[j] << " yAz = " << yAz;
   }

   (*os) << "\nyb - yAz = " << yb - yAz << endl;
   UTIL_DELARR(yA);
   UTIL_DELARR(z);

   //TODO: tol
   if (yb - yAz > 1.0e-3) {
      return true;
   } else {
      return false;
   }
}
#endif

// --------------------------------------------------------------------- //
void DecompAlgo::printBasisInfo(OsiSolverInterface* si,
                                ostream*             os)
{
   int      b, r, c;
   int*     basics   = 0;
   int*     rstat    = 0;
   int*     cstat    = 0;
   double* bInvRow  = 0;
   double* bInvARow = 0;
   const int n  = si->getNumCols();
   const int m  = si->getNumRows();
   char type[4] = {'F', 'B', 'U', 'L'};
   basics   = new int[m];
   bInvRow  = new double[m];
   bInvARow = new double[n];
   rstat    = new int[m];
   cstat    = new int[n];
   si->enableSimplexInterface(false);
   si->getBasics(basics);
   (*os) << "\n\nBasics: ";

   for (b = 0; b < m; b++) {
      (*os) << basics[b] << " ";
   }

   si->getBasisStatus(cstat, rstat);
   (*os) << "\ncstat: ";

   for (c = 0; c < n; c++) {
      (*os) << type[cstat[c]];
   }

   (*os) << "\n";
   (*os) << "rstat: ";

   for (r = 0; r < m; r++) {
      (*os) << type[rstat[r]];
   }

   (*os) << "\n";
   (*os) << "\nB-1:";

   for (r = 0; r < m; r++) {
      si->getBInvRow(r, bInvRow);
      (*os) << "\nB-1Row r: " << r << ": ";

      for (b = 0; b < m; b++) {
         (*os) << bInvRow[b] << " ";
      }
   }

   (*os) << "\nB-1A:";

   for (r = 0; r < m; r++) {
      si->getBInvARow(r, bInvARow);
      (*os) << "\nB-1ARow r: " << r << ": ";

      for (c = 0; c < n; c++) {
         (*os) << bInvARow[c] << " ";
      }
   }

   UTIL_DELARR(basics);
   UTIL_DELARR(bInvRow);
   UTIL_DELARR(bInvARow);
   UTIL_DELARR(rstat);
   UTIL_DELARR(cstat);
   si->disableSimplexInterface();
   //if you do this and want dual ray back, you need to resolve
   si->setHintParam(OsiDoPresolveInResolve, false, OsiHintDo);
   si->resolve();
   si->setHintParam(OsiDoPresolveInResolve, true, OsiHintDo);
}

// --------------------------------------------------------------------- //
void DecompAlgo::printCurrentProblem(const OsiSolverInterface* si,
                                     const string               baseName,
                                     const int                  nodeIndex,
                                     const int                  cutPass,
                                     const int                  pricePass,
                                     const bool                 printMps,
                                     const bool                 printLp)
{
   string filename = baseName
                     + ".n" + UtilIntToStr(nodeIndex)
                     + ".c" + UtilIntToStr(cutPass)
                     + ".p" + UtilIntToStr(pricePass);

   if (printMps) {
      si->writeMps(filename.c_str());
   }

   if (printLp)
      si->writeLp(filename.c_str(),
                  "lp", 1.0e-6, 10, 2, 1.0, true);
}

// --------------------------------------------------------------------- //
void DecompAlgo::printCurrentProblem(const OsiSolverInterface* si,
                                     const string               fileName,
                                     const bool                 printMps,
                                     const bool                 printLp)
{
   string filename = fileName;

   if (printMps) {
      si->writeMps(filename.c_str());
   }

   if (printLp)
      si->writeLp(filename.c_str(),
                  "lp", 1.0e-6, 10, 2, 1.0, true);
}

// --------------------------------------------------------------------- //
void DecompAlgo::printVars(ostream* os)
{
   DecompVarList::iterator it;
   int var_index = 0;

   for (it = m_vars.begin(); it != m_vars.end(); it++) {
      (*os) << "VAR " << var_index++ << " : ";
      (*it)->print(os, m_app);
      (*os) << endl;
   }

   (*os) << endl;
}

// --------------------------------------------------------------------- //
void DecompAlgo::printCuts(ostream* os)
{
   DecompCutList::iterator it;
   int cut_index = 0;

   for (it = m_cuts.begin(); it != m_cuts.end(); it++) {
      (*os) << "CUT " << cut_index++ << " : ";
      (*it)->print(os);
      (*os) << endl;
   }

   (*os) << endl;
}

// --------------------------------------------------------------------- //
void DecompAlgo::createFullMps(const string fileName)
{
   CoinAssert(m_algo == CUT);
   int n_integerVars = static_cast<int>(m_modelRelax->integerVars.size());
   m_masterSI->setInteger(&m_modelRelax->integerVars[0], n_integerVars);
   m_masterSI->writeMps(fileName.c_str());
   m_masterSI->setContinuous(&m_modelRelax->integerVars[0], n_integerVars);
}

// --------------------------------------------------------------------- //
void DecompAlgo::solveBruteForce()
{
   //
   // ---
   // --- Solve the original IP with a generic IP solver.
   // ---
   // ---  A simple sanity check for the case where it is possible to
   // ---  represent [A,b] in polynomial size, for example, SmallIP, 3AP.
   // ---
   //
   // --- welcome to decomp
   //
   //startupLog();
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << " <---- solveBruteForce() ---- ";
             );
   //
   //---  APP: create the user model A = A' union A''
   //
   int i;
   OsiIpSolverInterface* si = new OsiIpSolverInterface();
#ifdef __DECOMP_LP_CLP__
#ifdef __DECOMP_IP_CBC__
   si->getModelPtr()->messageHandler()->
   setLogLevel(m_app->m_param.LogLpLevel);
   si->getRealSolverPtr()->messageHandler()->
   setLogLevel(m_app->m_param.LogLpLevel);
#endif
#endif
   printf("\nSET MSG LEVEL = %d", m_app->m_param.LogLpLevel);
   si->messageHandler()->setLogLevel(m_app->m_param.LogLpLevel);
   //
   // ---
   // --- Construct the full model by combining:
   // ---   [A' ,b' ] modelRelax
   // ---   [A'',b''] modelCore
   // ---
   // THINK: what if A'',b'' = A,b, as in MMKP?
   si->loadProblem(*m_modelCore->M,
                   &m_modelCore->colLB[0],
                   &m_modelCore->colUB[0],
                   m_app->m_model.objCoeff,
                   &m_modelCore->rowLB[0],
                   &m_modelCore->rowUB[0]);

   //or, could use bottomAppendPackedMatrix
   for (i = 0; i < m_modelRelax->getNumRows(); i++) {
      si->addRow(m_modelRelax->M->getVector(i),
                 m_modelRelax->rowLB[i],
                 m_modelRelax->rowUB[i]);
   }

   int n_integerVars = static_cast<int>(m_modelRelax->integerVars.size());

   for (i = 0; i < n_integerVars; i++) {
      si->setInteger(m_modelRelax->integerVars[i]);
   }

   //si->writeMps("smallip");
   si->branchAndBound();
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              m_app->printOriginalSolution(si->getNumCols(),
                                           si->getColSolution());
             );
   UTIL_DELPTR(si);
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << "  ---- solveBruteForce()     ---->";
             );
}
