//===========================================================================//
// This file is part of the DIP Solver Framework.                            //
//                                                                           //
// DIP is distributed under the Eclipse Public License as part of the        //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Authors: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)      //
//          Ted Ralphs, Lehigh University (ted@lehigh.edu)                   //
//          Jiadong Wang, Lehigh University (jiw408@lehigh.edu)              //
//                                                                           //
// Copyright (C) 2002-2019, Lehigh University, Matthew Galati, Ted Ralphs    //
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
   //--- Need to deal with masterOnly variable

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
   const int      nRows     = m_masterSI->getNumRows();
   const double* rowRhs    = m_masterSI->getRightHandSide();
   const double* dual      = m_masterSI->getRowPrice();
   const double   primalObj = m_masterSI->getObjValue();
   double dualObj = 0.0;
   const int nCols = m_masterSI->getNumCols();
   const double* rc = m_masterSI->getReducedCost();
   const double* colLower = m_masterSI->getColLower();
   const double* colUpper = m_masterSI->getColUpper();
   //rStat might not be needed now, but will be needed
   // when we support ranged rows.
   int* rStat = new int[nRows];
   int* cStat = new int[nCols];
   m_masterSI->getBasisStatus(cStat, rStat);

   for (int c = 0; c < nCols; c++) {
      if (cStat[c] == 3) {
         dualObj += rc[c] * colLower[c];
      } else if (cStat[c] == 2 ) {
         dualObj += rc[c] * colUpper[c];
      }
   }

   for (int r = 0; r < nRows; r++) {
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

   UTIL_DELARR(rStat);
   UTIL_DELARR(cStat);
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
   if (m_param.DecompIPSolver == "CPLEX"){
#ifdef DIP_HAS_CPX
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
   }
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
   UTIL_DEBUG(m_param.LogDebugLevel, 3,
	      
	      if (printMps)
		 (*m_osLog) << "calling writeMps fileName = "
			    << fileName << endl;
	      if (printLp)
		 (*m_osLog) << "calling writeLp  fileName = "
			    << fileName << endl;
	      );
   

   if (m_param.DecompIPSolver == "CPLEX"){
#ifdef DIP_HAS_CPX
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
      if (printMps) {
	 string fileNameMps = fileName + ".mps";
	 si->writeMpsNative(fileNameMps.c_str(),
			    const_cast<const char**>(rowNamesChar),
			    const_cast<const char**>(colNamesChar), 1);
      }
      if (printLp) {
	 double epsilon      = 1e-30;
	 int    numberAcross = 5;
	 int    decimals     = 10;
	 string fileNameLp   = fileName + ".lp";
	 si->writeLpNative(fileNameLp.c_str(),
			   rowNamesChar, colNamesChar,
			   epsilon, numberAcross, decimals);
      }
      for (int i = 0; i < nCols; i++) {
	 UTIL_DELARR(colNamesChar[i]);
      }
      
      for (i = 0; i < (nRows + 1); i++) {
	 UTIL_DELARR(rowNamesChar[i]);
      }
      
      UTIL_DELARR(colNamesChar);
      UTIL_DELARR(rowNamesChar);
#endif
   }else{
      if (printMps) {
	 si->writeMps(fileName.c_str());
      }

      if (printLp) {
	 double epsilon      = 1e-30;
	 int    numberAcross = 5;
	 int    decimals     = 10;
	 string fileNameLp   = fileName + ".lp";
	 //This works because the Osi object in this case is OsiClp
	 // and Clp takes care of transferring the names.
	 si->writeLp(fileName.c_str(), "lp",
		     epsilon, numberAcross, decimals);
      }
   }

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
   if (m_param.DecompIPSolver == "CPLEX"){
      si->writeMpsNative(filename.c_str(), NULL, NULL, 1);
   }else{
      si->writeMps(filename.c_str());
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
      (*it)->print(m_infinity, os, m_app);
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
void DecompAlgo::checkDuals()
{
   //---
   //--- sanity check on duals returned
   //---   complementary slackness (c-uA)x = 0
   //---   also check that the given reduced cost matches the
   //---   hand calculation
   //---
   const double*            x     = m_masterSI->getColSolution();
   const double*            pi    = m_masterSI->getRowPrice();
   const int                nCols = m_masterSI->getNumCols();
   const CoinPackedMatrix* M     = m_masterSI->getMatrixByRow();
   double*                  uA    = new double[nCols];
   const double* objC = m_masterSI->getObjCoefficients();
   const double* rcLP = m_masterSI->getReducedCost();
   M->transposeTimes(pi, uA);
#ifndef DO_INTERIOR

   for (int i = 0; i < nCols; i++) {
      if (!UtilIsZero( x[i], 1.0e-5 ) &&
	  !UtilIsZero( (objC[i] - uA[i]) * x[i], 1.0e-4 ) ) {
	 printf("ERR in COMPL-SLACK i:%d objC:%15.10f uA:%15.10f x:%15.10f\n",
		i, objC[i], uA[i], x[i]);
	 fflush(stdout);
	 assert(0);
      }

      if (!UtilIsZero( (objC[i] - uA[i]) - rcLP[i], 1.0e-4 ) ) {
	 printf("ERR in RC i:%d objC:%15.10f uA:%15.10f RCLP:%15.10f\n",
		i, objC[i], uA[i], rcLP[i]);
	 fflush(stdout);
	 assert(0);
      }
   }

#endif
   UTIL_DELARR(uA);

}

//===========================================================================//
void DecompAlgo::checkReducedCost(const double *u, const double *u_adjusted)
{
   //---
   //--- sanity check - none of the columns currently in master
   //--- should have negative reduced cost
   //---   m_vars contains the variables (in x-space) that have
   //---   been pushed into the master LP (assumes no compression)
   //---

   DecompVarList::iterator it;
   int b, var_index = 0;
   double*       redCostX      = NULL;
   const double* objC = m_masterSI->getObjCoefficients();
   const double* rcLP = m_masterSI->getReducedCost();
   DecompConstraintSet* modelCore   = m_modelCore.getModel();
   const int     nCoreCols     = modelCore->getNumCols();
   double        alpha         = 0.0;
   int           nBaseCoreRows = modelCore->nBaseRows;
   const double* origObjective = getOrigObjective();

   for (it = m_vars.begin(); it != m_vars.end(); it++) {
      double redCost = 0.0;
      //m_s      is a sparse vector in x-space (the column)
      //redCostX is a dense  vector in x-space (the cost in subproblem)
      b       = (*it)->getBlockId();
      redCost = (*it)->m_s.dotProduct(redCostX);//??

      if ( (*it)->getVarType() == DecompVar_Point) {
	 alpha   = u[nBaseCoreRows + b];
      } else if ((*it)->getVarType() == DecompVar_Ray) {
	 alpha = 0;
      }

      assert(m_masterRowType[nBaseCoreRows + b] == DecompRow_Convex);
      assert(isMasterColStructural((*it)->getColMasterIndex()));
      UTIL_DEBUG(m_app->m_param.LogDebugLevel, 5,
		 (*m_osLog)
		 << "MasterColIndex = "
		 << setw(6) << (*it)->getColMasterIndex()
		 << "Block = "
		 << setw(3) << b
		 << "LPRedCost = " << setw(10)
		 << UtilDblToStr(rcLP[(*it)->getColMasterIndex()], 5)
		 << "CalcRedCost = " << setw(10)
		 << UtilDblToStr(redCost - alpha, 5)
		 << "ObjCost = " << setw(10)
		 << UtilDblToStr(objC[(*it)->getColMasterIndex()], 5)
		 << "Alpha = " << setw(10)
		 << UtilDblToStr(alpha, 5)
		 << endl; );
      //---
      //--- sanity check - none of the columns currently in master
      //---    should have negative reduced cost
      //--- unless they have been fixed to 0 by branching
      //---
      //const double * colLB = m_masterSI->getColLower();
      const double* colUB = m_masterSI->getColUpper();
      int            index = (*it)->getColMasterIndex();
      double         rcLPi = rcLP[index];

      if (rcLPi        < - m_param.RedCostEpsilon &&
	  colUB[index] > DecompEpsilon) {
	 (*m_osLog) << "VAR v-index:" << var_index++
		    << " m-index: " << (*it)->getColMasterIndex()
		    << " b-index: " << b
		    << " rcLP: "   << rcLPi
		    << endl;
	 (*it)->print(m_infinity, m_osLog, modelCore->colNames,
		      const_cast<double*>(redCostX));
	 (*m_osLog) << "******** ERROR ********" << endl;
	 assert(0);
      }

      //---
      //--- check that objective in LP and calculated objective match
      //---
      if (m_phase == PHASE_PRICE2) {
	 double objCalc = (*it)->m_s.dotProduct(origObjective);

	 if (!UtilIsZero(objCalc - objC[(*it)->getColMasterIndex()],
			 1.0e-3)) {
	    (*m_osLog) << "VAR v-index:" << var_index++
		       << " m-index: " << (*it)->getColMasterIndex()
		       << " b-index: " << b
		       << " objLP: "   << objC[(*it)->getColMasterIndex()]
		       << " objCalc: " << objCalc
		       << endl;
	    (*it)->print(m_infinity, m_osLog, modelCore->colNames,
			 const_cast<double*>(origObjective));
	    (*m_osLog) << "******** ERROR ********" << endl;
	    assert(0);
	 }
      }

      //---
      //--- check that LP reduced cost and calculated reduced cost
      //---  match up
      //--- in the case of using dual smoothing, we cannot do this check
      //---  since the calculated reduced cost is based on the smoothed
      //---  duals
      //---
      if (!m_param.DualStab &&
	  !UtilIsZero(rcLP[(*it)->getColMasterIndex()]
		      - (redCost - alpha), 1.0e-3)) {
	 //---
	 //--- this whole next section is an expansion of log
	 //---   when there is an issue found that the solver
	 //---   returns RC that doesn't match the one calculated
	 //---   based on core matrix and duals
	 //---
	 (*m_osLog) << "VAR v-index:" << var_index++
		    << " m-index: " << (*it)->getColMasterIndex()
		    << " b-index: " << b
		    << " rc: "      << redCost
		    << " alpha: "   << alpha
		    << " rc-a: "    << redCost - alpha
		    << " RCLP: "    << rcLP[(*it)->getColMasterIndex()]
		    << endl;
	 //---
	 //--- this, plus alpha shows the calculation of the red-cost in
	 //---   x-space, next, look at the same calculation in lambda-space
	 //---
	 (*it)->print(m_infinity, m_osLog, modelCore->colNames,
		      const_cast<double*>(redCostX));
	 (*m_osLog) << "******** ERROR ********" << endl;
	 //---
	 //--- the rows in lambda that show up should be ones where
	 //---   these components show up in original A''
	 //---
	 double* uA2    = new double[nCoreCols];
	 modelCore->M->transposeTimes(u_adjusted, uA2);
	 (*it)->print(m_infinity, m_osLog, modelCore->colNames,
		      const_cast<double*>(uA2));
	 UTIL_DELARR(uA2);
	 //---
	 //--- RC of a col of lambda-space is u.A[i]
	 //---
	 (*m_osLog) << " objLP: "
		    << UtilDblToStr(objC[(*it)->getColMasterIndex()], 4)
		    << endl;

	 //---
	 //--- recalc column of master A.s for this s
	 //---
	 if (m_algo != DECOMP) {
	    double* denseS = new double[modelCore->getNumCols()];
	    (*it)->fillDenseArr(modelCore->getNumCols(), denseS);
	    int r;
	    const CoinPackedMatrix* Mr = modelCore->getMatrix();

	    for (r = 0; r < modelCore->getNumRows(); r++) {
	       printf("\nROW %d\n", r);
	       CoinShallowPackedVector vec = Mr->getVector(r);
	       UtilPrintPackedVector(vec, m_osLog,
				     modelCore->getColNames(),
				     denseS);
	    }

	    UTIL_DELARR(denseS);
	 }

	 const CoinPackedMatrix* Mc = m_masterSI->getMatrixByCol();

            CoinShallowPackedVector vec
	       = Mc->getVector((*it)->getColMasterIndex());

            UtilPrintPackedVector(vec, m_osLog, m_masterSI->getColNames(), u);

            double uA = vec.dotProduct(u);

            (*m_osLog) << " objLP: "
                       << UtilDblToStr(objC[(*it)->getColMasterIndex()], 4)
                       << endl;

            (*m_osLog) << " uA   : "   << UtilDblToStr(uA, 4) << endl;

            (*m_osLog) << " RC   : "
                       << UtilDblToStr(objC[(*it)->getColMasterIndex()] - uA,
                                       4) << endl;

            (*m_osLog) << " RCLP : "
                       << UtilDblToStr(rcLP[(*it)->getColMasterIndex()], 4)
                       << endl;

            assert(0);

            (*m_osLog) << endl;
      } //END: if(!UtilIsZero(rcLP[(*it)->getColMasterIndex()] ...
   } //END: for(it = m_vars.begin(); it != m_vars.end(); it++)
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
   double        objLB      = -m_infinity;
   double        objUB      =  m_infinity;
   int           logIpLevel = m_param.LogIpLevel;
   DecompConstraintSet* modelCore = m_modelCore.getModel();
   int                   numInts   = modelCore->getNumInts();
   int                   numCols   = m_masterSI->getNumCols();
   double                timeLimit = m_param.TimeLimit;
   //---
   //--- start timer
   //---
   UtilTimer timer;
   timer.start();
   //---
   //--- create a results object
   //---
   DecompSolverResult* result = new DecompSolverResult(m_infinity);
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

   if (m_param.DecompIPSolver == "Cbc"){
#ifdef DIP_HAS_CBC
      CbcModel cbc(*m_masterSI);
      cbc.setLogLevel(logIpLevel);
      cbc.setDblParam(CbcModel::CbcMaximumSeconds, timeLimit);
      cbc.branchAndBound();
      const int statusSet[2] = {0, 1};
      int       solStatus    = cbc.status();
      int       solStatus2   = cbc.secondaryStatus();
      
      if (!UtilIsInSet(solStatus, statusSet, 2)) {
	 cerr << "Error: CBC IP solver status = "
	      << solStatus << endl;
	 throw UtilException("CBC solver status", "solveDirect", "solveDirect");
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
#else
      throw UtilException("Cbc selected as solver, but it's not available",
			  "solveDirect", "DecompDebug");
#endif
   }else if (m_param.DecompIPSolver == "CPLEX"){
#ifdef DIP_HAS_CPX
      OsiCpxSolverInterface* masterSICpx =
	 dynamic_cast<OsiCpxSolverInterface*>(m_masterSI);
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
      //---
      //--- set the thread limit, otherwise CPLEX will use all the resources
      //---
      status = CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, m_param.NumThreadsIPSolver);
      
      if (status)
	 throw UtilException("CPXsetdblparam failure",
			     "solveDirect", "DecompAlgoC");
      
      //---
      //--- solve the MILP
      //---
      UtilTimer timer1;
      timer1.start();
      m_masterSI->branchAndBound();
      timer1.stop();
      cout << "just after solving" << endl;
      cout << " Real=" << setw(10) << UtilDblToStr(timer1.getRealTime(), 5)
	   << " Cpu= " << setw(10) << UtilDblToStr(timer1.getCpuTime() , 5);
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
#else
      throw UtilException("CPLEX selected as solver, but it's not available",
			  "solveDirect", "DecompDebug");
#endif
   }else{
      throw UtilException("solveDirect not implemented for selected solver",
			  "solveDirect", "DecompDebug");
   }

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

