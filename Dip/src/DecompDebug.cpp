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

#include "DecompApp.h"
#include "DecompAlgo.h"
#include "DecompAlgoC.h"

using namespace std;

//===========================================================================//
bool DecompAlgo::checkPointFeasible(const DecompConstraintSet* model,
                                    const double*               x)
{
   //---
   //--- sanity check
   //---   Does the recomposed solution (x*) satisfy the core
   //---   constraints. If not, but in master solver OR in the
   //---   process of recomposed (the map).
   //---
   const  CoinPackedMatrix* M        = model->getMatrix();

   if (!M) {
      return true;
   }

   int    i;
   double actViol;
   double relViol;
   int    precision                   = 7;
   bool   isFeas                      = true;
   bool   hasColNames                 = false;
   bool   hasRowNames                 = false;
   const  int                nCols    = model->getNumCols();
   const  int                nRows    = model->getNumRows();
   const  double*            colLB    = model->getColLB();
   const  double*            colUB    = model->getColUB();
   const  double*            rowLB    = model->getRowLB();
   const  double*            rowUB    = model->getRowUB();
   const  vector<string>&    colNames = model->getColNames();
   const  vector<string>&    rowNames = model->getRowNames();
   double*                   ax       = new double[nRows];
   assert(M);
   assert(ax);

   if (colNames.size()) {
      hasColNames = true;
   }

   if (rowNames.size()) {
      hasRowNames = true;
   }

   //---
   //--- check column bounds
   //---
   for (i = 0; i < nCols; i++) {
      actViol = std::max<double>(colLB[i] - x[i], x[i] - colUB[i]);
      actViol = std::max<double>(actViol, 0.0);

      if (UtilIsZero(x[i], 1.0e-3)            ||
            (x[i] < 0 && UtilIsZero(colLB[i])) ||
            (x[i] > 0 && UtilIsZero(colUB[i]))) {
         relViol = actViol;
      } else {
         relViol = actViol / std::fabs(x[i]);
      }

      if (relViol > 0.0001) { //0.01% violated
         (*m_osLog) << "Point violates column " << i;

         if (hasColNames) {
            (*m_osLog) << " -> " << colNames[i];
         }

         (*m_osLog) << " LB= " << UtilDblToStr(colLB[i], precision)
                    << " x= "  << UtilDblToStr(x[i], precision)
                    << " UB= " << UtilDblToStr(colUB[i], precision)
                    << " RelViol= " << UtilDblToStr(relViol, precision)
                    << endl;

         //>1% violation is probably a bug, but <1% could be just
         //  round off error??? not sure about that
         if (relViol > 0.01) {
            isFeas = false;
         }
      }
   }

   //---
   //--- M * x = ax
   //---
   M->times(x, ax);

   //---
   //--- check row bounds
   //---
   for (i = 0; i < nRows; i++) {
      actViol = std::max<double>(rowLB[i] - ax[i], ax[i] - rowUB[i]);
      //printf("ax=%12.10f, actViol=%12.10f\n", ax[i], actViol);
      actViol = std::max<double>(actViol, 0.0);

      //printf("            actViol=%12.10f\n", actViol);
      if (m_param.LogDebugLevel >= 4) {
         CoinShallowPackedVector row = M->getVector(i);
         (*m_osLog) << "Row i: " << i;

         if (hasRowNames) {
            (*m_osLog) << " -> " << rowNames[i];
         }

         (*m_osLog) << " LB= "   << UtilDblToStr(rowLB[i], precision)
                    << " ax= "   << UtilDblToStr(ax[i], precision)
                    << " UB= "   << UtilDblToStr(rowUB[i], precision) << endl;
         //UtilPrintPackedVector(row);
      }

      if (UtilIsZero(ax[i], 1.0e-3)            ||
            (ax[i] < 0 && UtilIsZero(rowLB[i])) ||
            (ax[i] > 0 && UtilIsZero(rowUB[i]))) {
         relViol = actViol;
      } else {
         relViol = actViol / std::fabs(ax[i]);
      }

      if (relViol > 0.005) { //0.5% violated
         (*m_osLog) << "Point violates row " << i;

         if (hasRowNames) {
            (*m_osLog) << " -> " << rowNames[i];
         }

         (*m_osLog) << " LB= " << UtilDblToStr(rowLB[i], precision)
                    << " ax= "  << UtilDblToStr(ax[i], precision)
                    << " UB= " << UtilDblToStr(rowUB[i], precision)
                    << " RelViol= " << UtilDblToStr(relViol, precision)
                    << endl;

         //>5% violation is probably a bug, but <5% could be just
         //  round off error??? not sure about that
         if (relViol > 0.05) {
            isFeas = false;

            //---
            //--- if special case of relViol=actViol,
            //---   then check to see if possible round off issues
            //--- e.g., harp2 a[j]=1.0e9, actViol=1.0e3 is OK
            //---
            if (UtilIsZero(ax[i], 1.0e-3)            ||
                  (ax[i] < 0 && UtilIsZero(rowLB[i])) ||
                  (ax[i] > 0 && UtilIsZero(rowUB[i]))) {
               int                       k;
               CoinShallowPackedVector   row   = M->getVector(i);
               const int                 numNZ = row.getNumElements();
               const double*             els   = row.getElements();

               for (k = 0; k < numNZ; k++) {
                  if (fabs(els[k]) > 1.0e7) {
                     (*m_osLog) << "  row has a big coefficient "
                                << els[k] << endl;
                     isFeas = true;
                     break;
                  }
               }
            }
         }
      }
   }

   UTIL_DELARR(ax);
   return isFeas;
}

//===========================================================================//
void DecompAlgo::checkMasterDualObj()
{
   int    r;
   const int      nRows     = m_masterSI->getNumRows();
   const double* rowRhs    = m_masterSI->getRightHandSide();
   const double* dual      = m_masterSI->getRowPrice();
   const double   primalObj = m_masterSI->getObjValue();
   double dualObj = 0.0;

   for (r = 0; r < nRows; r++) {
      dualObj += dual[r] * rowRhs[r];
   }

   UTIL_DEBUG(m_param.LogDebugLevel, 4,
              (*m_osLog)
              << "checkMasterDualObj"
              << setw(10) << "primalObj="
              << setw(10) << UtilDblToStr(primalObj, 3)
              << setw(10) << "dualObj="
              << setw(10) << UtilDblToStr(dualObj, 3) << endl;
             );
   double actViol = std::fabs(primalObj - dualObj);
   double relViol = actViol;

   if (!UtilIsZero(dualObj, 1.0e-3)) {
      relViol = actViol / std::fabs(dualObj);
   }

   if (relViol > 1.0e-4) {
      cerr << "checkMasterDualObj"
           << setw(10) << "primalObj="
           << setw(10) << UtilDblToStr(primalObj, 3)
           << setw(10) << "dualObj="
           << setw(10) << UtilDblToStr(dualObj, 3) << endl;
      throw UtilException("primal and dual obj do not match",
                          "checkMasterDualObj", "DecompAlgo");
   }
}

//===========================================================================//
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
         (*os) << "i : " << i << " dualRay = " << dualRay[i]
               << " rowRhs = " << rowRhs[i] << " yb = " << yb << endl;
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
         (*os) << "yA[" << i << "]:\t" << yA[i];
      }

      //TODO: tol 1.0e-6 is too tight?
      if ((ybPos  && (yA[i] >  1.0e-2)) ||
            (!ybPos && (yA[i] < -1.0e-2))) {
         if (os) {
            (*os) << " -->isProof (false)" << endl;
         }

         isProof = false;
      } else if (os) {
         (*os) << endl;
      }
   }

   UTIL_DELARR(yA);
#if 0

   //sanity check
   if (!isProof)
      isProof
         = isDualRayInfProofCpx(dualRay, rowMatrix, colLB, colUB, rowRhs, os);

#endif
   return isProof;
}

//===========================================================================//
bool DecompAlgo::isDualRayInfProofCpx(const double*            dualRay,
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
   //---  proof_p = y^T b - y^T A z > 0
   //---
   //---  So, we want to maximize y^T A x to break the proof.
   //---
   int       i, j;
   double   yb, yAz;
   double* yA = 0;
   double*   z = 0;
   const int m = rowMatrix->getNumRows();
   const int n = rowMatrix->getNumCols();
   //TODO: check for out-of-mem conditions?
   yA = new double[n];
   UtilFillN(yA, n, 0.0);
   double* yA2 = new double[n];
   rowMatrix->transposeTimes(dualRay, yA2);     //y^T A

   for (i = 0; i < m; i++) {
      double yA_i = 0;
      CoinShallowPackedVector pv = rowMatrix->getVector(i);
      const int*     indI = pv.getIndices();
      const double* elsI = pv.getElements();
      const int      lenI = pv.getNumElements();

      for (int j = 0; j < lenI; j++) {
         yA_i += dualRay[indI[j]] * elsI[j];
         printf("i: %d, j: %d, indIj: %d, elsIj: %g ray: %g yA_i: %g\n",
                i, j, indI[j], elsI[j], dualRay[indI[j]], yA_i);
      }

      yA[i] = yA_i;

      if (!UtilIsZero(yA[i] - yA2[i])) {
         printf(" ---> yA: %g, yA2: %g\n", yA[i], yA2[i]);
      }

      fflush(stdout);
      CoinAssert(UtilIsZero(yA[i] - yA2[i]));
   }

   z  = new double[n];

   for (j = 0; j < n; j++) {
      if (yA[j] >= 0) {
         z[j] = CoinMin(1.0e20, colUB[j]);
      } else {
         z[j] = colLB[j];
      }
   }

   //y^T b
   yb = 0.0;

   for (i = 0; i < m; i++) {
      yb += dualRay[i] * rowRhs[i];

      if (os)
         (*os) << "\ni : " << i << " dualRay = " << dualRay[i]
               << " rowRhs = " << rowRhs[i] << " yb = " << yb;
   }

   //y^T A z
   yAz = 0.0;

   for (j = 0; j < n; j++) {
      yAz += yA[j] * z[j];

      if (os)
         (*os) << "\nj : " << j << " yA = " << yA[j]
               << " z = " << z[j] << " yAz = " << yAz;
   }

   if (os) {
      (*os) << "\nyb - yAz = " << yb - yAz << endl;
   }

   UTIL_DELARR(yA);
   UTIL_DELARR(z);

   //TODO: tol
   if (yb - yAz > 1.0e-3) {
      return true;
   } else {
      return false;
   }
}

//===========================================================================//
void DecompAlgo::printBasisInfo(OsiSolverInterface* si,
                                ostream*              os)
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
   //TODO: have to check sense?
   const double* rowRhs = si->getRightHandSide();
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
   //yb, where y is a row of B-1
   double yb = 0.0;
   (*os) << "\nB-1:";

   for (r = 0; r < m; r++) {
      yb = 0.0;
      si->getBInvRow(r, bInvRow);
      (*os) << "\nB-1Row r: " << r << ": ";

      for (b = 0; b < m; b++) {
         (*os) << bInvRow[b] << " ";
         //rowRhs is just orig row rhs? or change based on who is basic?
         yb += bInvRow[b] * rowRhs[b];
      }

      (*os) << " ---> yb: " << yb;
   }

   //all pos case? if yb < 0
   //all neg case? if yb > 0
   //  what if yb=0?
   (*os) << "\nB-1A:";
   bool allpos = true;
   bool allneg = true;

   for (r = 0; r < m; r++) {
      si->getBInvARow(r, bInvARow);
      (*os) << "\nB-1ARow r: " << r << ": ";
      allpos = true;
      allneg = true;

      for (c = 0; c < n; c++) {
         (*os) << bInvARow[c] << " ";

         if (bInvARow[c] < 0) {
            allpos = false;
         }

         if (bInvARow[c] > 0) {
            allneg = false;
         }
      }

      if (allpos) {
         (*os) << " ---> allpos";
      }

      if (allneg) {
         (*os) << " ---> allneg";
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

//===========================================================================//
void DecompAlgo::printCurrentProblemDual(OsiSolverInterface* si,
      const string         baseName,
      const int            nodeIndex,
      const int            cutPass,
      const int            pricePass)
{
   if (!si) {
      return;
   }

   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "printCurrentProblemDual()", m_param.LogDebugLevel, 2);
#ifdef __DECOMP_LP_CPX__
   OsiCpxSolverInterface* siCpx
      = dynamic_cast<OsiCpxSolverInterface*>(si);
   CPXENVptr env = siCpx->getEnvironmentPtr();
   CPXLPptr  lp  = siCpx->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL);
   string filename = DecompAlgoStr[m_algo] + "_" + baseName
                     + ".n" + UtilIntToStr(nodeIndex)
                     + ".c" + UtilIntToStr(cutPass)
                     + ".p" + UtilIntToStr(pricePass)
                     + ".dual.mps";
   double objShift;
   int status = CPXdualwrite(env, lp, filename.c_str(), &objShift);

   if (status)
      throw UtilException("CPXdualwrite failure",
                          "printCurrentProblemDual", "DecompAlgo");

   printf("objShift in dual = %g\n", objShift);
   UTIL_DEBUG(m_param.LogDebugLevel, 3,
              (*m_osLog) << "calling CPXdualwrite filename = "
              << filename << endl;
             );
#endif
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "printCurrentProblemDual()", m_param.LogDebugLevel, 2);
}

//===========================================================================//
void DecompAlgo::printCurrentProblem(const OsiSolverInterface* si,
                                     const string               baseName,
                                     const int                  nodeIndex,
                                     const int                  cutPass,
                                     const int                  pricePass,
                                     const int                  blockId,
                                     const bool                 printMps,
                                     const bool                 printLp)
{
   if (!si) {
      return;
   }

   string fileName = DecompAlgoStr[m_algo] + "_" + baseName
                     + ".n" + UtilIntToStr(nodeIndex)
                     + ".c" + UtilIntToStr(cutPass)
                     + ".p" + UtilIntToStr(pricePass);

   if (blockId != -1) {
      fileName += ".b" + UtilIntToStr(blockId);
   }

   printCurrentProblem(si, fileName, printMps, printLp);
}

//===========================================================================//
void DecompAlgo::printCurrentProblem(const OsiSolverInterface* si,
                                     const string               fileName,
                                     const bool                 printMps,
                                     const bool                 printLp)
{
   if (!si) {
      return;
   }

   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "printCurrentProblem()", m_param.LogDebugLevel, 2);
#ifdef __DECOMP_IP_CPX__
   //---
   //--- There is no derived OsiCpx::writeLp and the base writeLp does not
   //---   use names - for some reason (even though they are in Osi memory)
   //---
   //--- The characters [] are often used in names but not allowed by
   //---   CoinLp writer - so replace them here with ().
   //---
   int     i            = 0;
   int     nCols        = si->getNumCols();
   int     nRows        = si->getNumRows();
   char** colNamesChar = new char*[nCols];
   char** rowNamesChar = new char*[nRows + 1];

   for (i = 0; i < nCols; i++) {
      string colName  = si->getColName(i);
      replace(colName.begin(), colName.end(), '[', '(');
      replace(colName.begin(), colName.end(), ']', ')');
      colNamesChar[i] = new char[colName.size() + 1];
      copy(colName.begin(), colName.end(), colNamesChar[i]);
      colNamesChar[i][colName.size()] = '\0';
   }

   for (i = 0; i < nRows; i++) {
      string rowName  = si->getRowName(i);
      replace(rowName.begin(), rowName.end(), '[', '(');
      replace(rowName.begin(), rowName.end(), ']', ')');
      rowNamesChar[i] = new char[rowName.size() + 1];
      copy(rowName.begin(), rowName.end(), rowNamesChar[i]);
      rowNamesChar[i][rowName.size()] = '\0';
   }

   string objName = si->getObjName();
   //printf("objname=%s\n", objName.c_str());
   replace(objName.begin(), objName.end(), '[', '(');
   replace(objName.begin(), objName.end(), ']', ')');
   rowNamesChar[nRows] = new char[objName.size() + 1];
   copy(objName.begin(), objName.end(), rowNamesChar[nRows]);
   rowNamesChar[nRows][objName.size()] = '\0';
   //printf("nRows=%d objname=%s\n", nRows, rowNamesChar[nRows]);
#endif
   UTIL_DEBUG(m_param.LogDebugLevel, 3,

              if (printMps)
              (*m_osLog) << "calling writeMps fileName = "
              << fileName << endl;
              if (printLp)
                 (*m_osLog) << "calling writeLp  fileName = "
                 << fileName << endl;
                );

   if (printMps) {
#ifdef __DECOMP_IP_CPX__
      string fileNameMps = fileName + ".mps";
      si->writeMpsNative(fileNameMps.c_str(),
                         const_cast<const char**>(rowNamesChar),
                         const_cast<const char**>(colNamesChar), 1);
#else
      si->writeMps(fileName.c_str());
#endif
   }

   if (printLp) {
      double epsilon      = 1e-30;
      int    numberAcross = 5;
      int    decimals     = 10;
      string fileNameLp   = fileName + ".lp";
#ifdef __DECOMP_IP_CPX__
      si->writeLpNative(fileNameLp.c_str(),
                        rowNamesChar, colNamesChar,
                        epsilon, numberAcross, decimals);
#else
      //This works because the Osi object in this case is OsiClp
      // and Clp takes care of transferring the names.
      si->writeLp(fileName.c_str(), "lp",
                  epsilon, numberAcross, decimals);
#endif
   }

#ifdef __DECOMP_IP_CPX__

   for (i = 0; i < nCols; i++) {
      UTIL_DELARR(colNamesChar[i]);
   }

   for (i = 0; i < (nRows + 1); i++) {
      UTIL_DELARR(rowNamesChar[i]);
   }

   UTIL_DELARR(colNamesChar);
   UTIL_DELARR(rowNamesChar);
#endif
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "printCurrentProblem()", m_param.LogDebugLevel, 2);
}

/*
//===========================================================================//
void DecompAlgo::printCurrentProblem(const OsiSolverInterface * si,
                                     const string               fileName,
                                     const bool                 printMps,
                                     const bool                 printLp){
   string filename = fileName;
   if(printMps){
      //const char ** rowNames = NULL;
      //const char ** colNames = NULL;
      //TODO: col/row names have to explicitly pass in
#ifdef __DECOMP_IP_CPX__
      si->writeMpsNative(filename.c_str(), NULL, NULL, 1);
#else
      si->writeMps(filename.c_str());
#endif
   }
   if(printLp)
      si->writeLp(filename.c_str(), "lp", 1e-30, 5, 10);
}
*/

//===========================================================================//
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

//===========================================================================//
void DecompAlgo::createFullMps(const string fileName)
{
   CoinAssert(m_algo == CUT);
   DecompConstraintSet*           modelCore   = m_modelCore.getModel();
   int n_integerVars = static_cast<int>(modelCore->integerVars.size());
   m_masterSI->setInteger(&modelCore->integerVars[0], n_integerVars);
   m_masterSI->writeMps(fileName.c_str());
   m_masterSI->setContinuous(&modelCore->integerVars[0], n_integerVars);
}

//===========================================================================//
void DecompAlgo::printCuts(ostream* os)
{
   DecompCutList::iterator it;
   int cut_index = 0;

   for (it = m_cuts.begin(); it != m_cuts.end(); it++) {
      (*os) << "CUT " << cut_index++ << " : ";
      (*it)->print(os);
   }

   (*os) << endl;
}

//===========================================================================//
DecompSolverResult* DecompAlgoC::solveDirect(const DecompSolution* startSol)
{
   //---
   //--- Solve the original IP with a generic IP solver.
   //---
   //---  A simple sanity check for the case where it is possible to
   //---  represent [A,b] in polynomial size, for example, SmallIP, 3AP.
   //---
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "solveDirect()", m_param.LogDebugLevel, 2);
   DecompVarList dummy;
   int           i, nNodes;
   double        objLB      = -DecompInf;
   double        objUB      =  DecompInf;
   int           logIpLevel = m_param.LogLpLevel;
   DecompConstraintSet* modelCore = m_modelCore.getModel();
   int                   numInts   = modelCore->getNumInts();
   int                   numCols   = m_masterSI->getNumCols();
   double                timeLimit = m_param.LimitTime;
   //---
   //--- start timer
   //---
   UtilTimer timer;
   timer.start();
   //---
   //--- create a results object
   //---
   DecompSolverResult* result = new DecompSolverResult();
   //---
   //--- create the master problem
   //---
   createMasterProblem(dummy);
   //---
   //--- adjust ip solver log levels
   //---
   m_masterSI->messageHandler()->setLogLevel(logIpLevel);

   //---
   //--- set integer vars
   //---
   for (i = 0; i < numInts; i++) {
      m_masterSI->setInteger(modelCore->integerVars[i]);
   }

   //#define PERMUTE_STUFF
   /*#ifdef  PERMUTE_STUFF

   //---
   //--- randomly permute rows and cols for MIPLIB2010
   //---    delete random rows, append to end
   //---    delete random cols, append to end
   //---
   {
   int k, r, c, tmp;
   int nCols         = m_masterSI->getNumCols();
   int nRows         = m_masterSI->getNumRows();
   int rowsToPermute = static_cast<int>(nRows / 7);
   int colsToPermute = static_cast<int>(nCols / 7);
   vector<int>  newRowInd; //old row to new row index
   int rowToDelete[1];
   int colToDelete[1];

   srand(1);
   for(i = 0; i < nRows; i++){
   newRowInd.push_back(i);
   }
   for(i = 0; i < rowsToPermute; i++){
   r = UtilURand(0, nRows-1);
   //---
   //--- Example:
   //---     0,1,2,3,4,5,6 (r=2)
   //---  -> 0,1,3,4,5,6,2
   //---
   tmp = newRowInd[r];
   for(k = r; k < (nRows-1); k++){
   newRowInd[k] = newRowInd[k+1];
   }
   newRowInd[nRows-1] = tmp;

   const CoinPackedMatrix * M     = m_masterSI->getMatrixByRow();
   const double           * rowLB = m_masterSI->getRowLower();
   const double           * rowUB = m_masterSI->getRowUpper();
   const double             rLB   = rowLB[r];
   const double             rUB   = rowUB[r];
   CoinShallowPackedVector vecS = M->getVector(r);
   CoinPackedVector        vec(vecS);
   //printf("delete and move to end row r=%d %s\n",
   //     r, m_masterSI->getRowName(r).c_str());
   //vec.print();

   rowToDelete[0] = r;
   m_masterSI->deleteRows(1, rowToDelete);
   m_masterSI->addRow(vec, rLB, rUB);
   }

   for(i = 0; i < colsToPermute; i++){
   c = UtilURand(0, nCols-1);
   //---
   //--- Example:
   //---     0,1,2,3,4,5,6 (r=2)
   //---  -> 0,1,3,4,5,6,2
   //---
   const CoinPackedMatrix * M     = m_masterSI->getMatrixByCol();
   const double           * colLB = m_masterSI->getColLower();
   const double           * colUB = m_masterSI->getColUpper();
   const double           * objC  = m_masterSI->getObjCoefficients();
   const double             cLB   = colLB[c];
   const double             cUB   = colUB[c];
   const double             obj   = objC[c];
   const CoinShallowPackedVector vecS = M->getVector(c);
   CoinPackedVector              vec(vecS);

   /////////// THIS IS WRONG ///////////
   //TODO: copy integer info!

   colToDelete[0] = c;
   m_masterSI->deleteCols(1, colToDelete);
   m_masterSI->addCol(vec, cLB, cUB, obj);
   }

   printf("\n\nNew Row Map\n");
   for(i = 0; i < nRows; i++){
   printf("%10d%10d\n", i, newRowInd[i]);
   }
   }
   #endif   */

   //---
   //--- dump full milp
   //---
   if (m_param.LogDumpModel >= 2) {
      string fileName = "directMILP";
      printCurrentProblem(m_masterSI, fileName);
   }

#ifdef __DECOMP_IP_CBC__
   CbcModel cbc(*m_masterSI);
   CbcMain0(cbc);
   //---
   //--- build argument list
   //---
   const char* argv[20];
   int    argc         = 0;
   string cbcExe       = "cbc";
   string cbcSolve     = "-solve";
   string cbcQuit      = "-quit";
   string cbcLog       = "-log";
   string cbcLogSet    = UtilIntToStr(logIpLevel);
   string cbcTime      = "-seconds";
   string cbcTimeSet   = UtilIntToStr(timeLimit);
   argv[argc++] = cbcExe.c_str();
   argv[argc++] = cbcLog.c_str();
   argv[argc++] = cbcLogSet.c_str();
   argv[argc++] = cbcTime.c_str();
   argv[argc++] = cbcTimeSet.c_str();
   argv[argc++] = cbcSolve.c_str();
   argv[argc++] = cbcQuit.c_str();
   //---
   //--- solve IP using argument list
   //---
   CbcMain1(argc, argv, cbc);
   //---
   //--- get solver status
   //---   comments based on Cbc2.3
   //---
   const int statusSet[2] = {0, 1};
   int       solStatus    = cbc.status();
   int       solStatus2   = cbc.secondaryStatus();

   if (!UtilIsInSet(solStatus, statusSet, 2)) {
      cerr << "Error: CBC IP solver status = "
           << solStatus << endl;
      throw UtilException("CBC solver status", "solveDirect", m_classTag);
   }

   //---
   //--- get number of nodes
   //---
   nNodes = cbc.getNodeCount();
   //---
   //--- get objective and solution
   //---
   objLB = cbc.getBestPossibleObjValue();

   if (cbc.isProvenOptimal() || cbc.isSecondsLimitReached()) {
      objUB = cbc.getObjValue();

      if (result && cbc.getSolutionCount()) {
         const double* solDbl = cbc.getColSolution();
         vector<double> solVec(solDbl, solDbl + numCols);
         result->m_solution.push_back(solVec);
         result->m_nSolutions++;
         assert(result->m_nSolutions ==
                static_cast<int>(result->m_solution.size()));
         //copy(solution, solution+numCols, result->m_solution);
      }
   }

   //---
   //--- copy sol status into result
   //---
   if (result) {
      result->m_solStatus  = solStatus;
      result->m_solStatus2 = solStatus2;
   }

#endif
#ifdef __DECOMP_IP_CPX__
   OsiIpSolverInterface* masterSICpx =
      dynamic_cast<OsiIpSolverInterface*>(m_masterSI);
   CPXLPptr  cpxLp  = masterSICpx->getLpPtr();
   CPXENVptr cpxEnv = masterSICpx->getEnvironmentPtr();
   int       status = 0;
   masterSICpx->switchToMIP();//need?

   if (startSol) {
      int            nCols    = masterSICpx->getNumCols();
      int            beg[1]   = {0};
      int*           varInd   = new int[nCols];
      const double* solution = startSol->getValues();
      assert(nCols == startSol->getSize());
      UtilIotaN(varInd, nCols, 0);
      status = CPXaddmipstarts(cpxEnv, cpxLp,
                               1, nCols, beg, varInd, solution, NULL, NULL);

      if (status)
         throw UtilException("CPXaddmipstarts failure",
                             "solveDirect", "DecompAlgoC");

      UTIL_DELARR(varInd);
   }

   //---
   //--- set the time limit
   //---
   status = CPXsetdblparam(cpxEnv, CPX_PARAM_TILIM, timeLimit);

   if (status)
      throw UtilException("CPXsetdblparam failure",
                          "solveDirect", "DecompAlgoC");

   //---
   //--- solve the MILP
   //---
   m_masterSI->branchAndBound();
   //---
   //--- get solver status
   //---
   //---
   int solStatus = CPXgetstat(cpxEnv, cpxLp);

   if (result) {
      result->m_solStatus  = solStatus;
      result->m_solStatus2 = 0;
   }

   //---
   //--- get number of nodes
   //---
   nNodes  = CPXgetnodecnt(cpxEnv, cpxLp);
   //---
   //--- get objective and solution
   //---
   status = CPXgetbestobjval(cpxEnv, cpxLp, &objLB);

   if (status)
      throw UtilException("CPXgetbestobjval failure",
                          "solveDirect", "DecompAlgoC");

   //---
   //--- get objective and solution
   //---
   if (solStatus == CPXMIP_OPTIMAL     ||
         solStatus == CPXMIP_OPTIMAL_TOL ||
         solStatus == CPXMIP_TIME_LIM_FEAS) {
      status = CPXgetmipobjval(cpxEnv, cpxLp, &objUB);

      if (status)
         throw UtilException("CPXgetmipobjval failure",
                             "solveDirect", "DecompAlgoC");

      if (result) {
         const double* solDbl = m_masterSI->getColSolution();
         vector<double> solVec(solDbl, solDbl + numCols);
         result->m_solution.push_back(solVec);
         result->m_nSolutions++;
         assert(result->m_nSolutions ==
                static_cast<int>(result->m_solution.size()));
         //copy(solution, solution+numCols, result->m_solution);
      }
   }

   //---
   //--- copy sol status into result
   //---
   if (result) {
      result->m_solStatus  = solStatus;
      result->m_solStatus2 = 0;
   }

#endif

   //---
   //--- copy bounds into result
   //---
   if (result) {
      result->m_objUB = objUB;
      result->m_objLB = objLB;
   }

   //---
   //--- stop the timer, dump time to solve
   //---
   timer.stop();
   (*m_osLog) << "DIRECT SOLVE"
              << " Real=" << setw(10) << UtilDblToStr(timer.getRealTime(), 5)
              << " Cpu= " << setw(10) << UtilDblToStr(timer.getCpuTime() , 5)
              << " Nodes= " << setw(8) << nNodes
              << " objLB= " << setw(10) << UtilDblToStr(objLB, 3)
              << " objUB= " << setw(10) << UtilDblToStr(objUB, 3)
              << endl;
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "solveDirect()", m_param.LogDebugLevel, 2);
   return result;
}

