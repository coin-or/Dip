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
// Copyright (C) 2002-2009, Lehigh University, Matthew Galati, Ted Ralphs    //
// All Rights Reserved.                                                      //
//===========================================================================//

#include "DecompApp.h"
#include "DecompAlgo.h"
#include "DecompAlgoC.h"

//===========================================================================//
bool DecompAlgo::checkPointFeasible(const DecompConstraintSet * model,
                                    const double              * x){
   //---
   //--- sanity check
   //---   Does the recomposed solution (x*) satisfy the core 
   //---   constraints. If not, but in master solver OR in the
   //---   process of recomposed (the map).
   //---
   const  CoinPackedMatrix * M        = model->getMatrix();
   if(!M)
      return true;

   int    i;
   double actViol;
   double relViol;
   int    precision                   = 7;
   bool   isFeas                      = true;
   bool   hasColNames                 = false;
   bool   hasRowNames                 = false;
   const  int                nCols    = model->getNumCols();
   const  int                nRows    = model->getNumRows();
   const  double           * colLB    = model->getColLB();
   const  double           * colUB    = model->getColUB();
   const  double           * rowLB    = model->getRowLB();
   const  double           * rowUB    = model->getRowUB();
   const  vector<string>   & colNames = model->getColNames();
   const  vector<string>   & rowNames = model->getRowNames();
   double                  * ax       = new double[nRows];
   assert(M);
   assert(ax);

   if(colNames.size())
      hasColNames = true;
   if(rowNames.size())
      hasRowNames = true;

   //---
   //--- check column bounds
   //---
   for(i = 0; i < nCols; i++){
      actViol = std::max<double>(colLB[i] - x[i], x[i] - colUB[i]);
      actViol = std::max<double>(actViol, 0.0);
      if(UtilIsZero(x[i],1.0e-3)            ||
	 (x[i] < 0 && UtilIsZero(colLB[i])) ||
	 (x[i] > 0 && UtilIsZero(colUB[i])))	 
	 relViol = actViol;
      else
	 relViol = actViol / std::fabs(x[i]);      
      if(relViol > 0.0001){//0.01% violated
	 (*m_osLog) << "YYPoint violates column " << i;
	 if(hasColNames)
	    (*m_osLog) << " -> " << colNames[i];
	 (*m_osLog) << " LB= " << UtilDblToStr(colLB[i],precision)
		    << " x= "  << UtilDblToStr(x[i],precision)
		    << " UB= " << UtilDblToStr(colUB[i],precision) 
		    << " RelViol= " << UtilDblToStr(relViol,precision)
		    << endl;
	 //>1% violation is probably a bug, but <1% could be just
	 //  round off error??? not sure about that
	 if(relViol > 0.01)
	    isFeas = false;
      }
   }
   
   //---
   //--- M * x = ax
   //---
   M->times(x,ax);

   //---
   //--- check row bounds
   //---
   for(i = 0; i < nRows; i++){
      actViol = std::max<double>(rowLB[i] - ax[i], ax[i] - rowUB[i]);
      //printf("ax=%12.10f, actViol=%12.10f\n", ax[i], actViol);
      actViol = std::max<double>(actViol, 0.0);
      //printf("            actViol=%12.10f\n", actViol);
      if(m_param.LogDebugLevel >= 4){
	 CoinShallowPackedVector row = M->getVector(i);
         (*m_osLog) << "Row i: " << i;
         if(hasRowNames)
            (*m_osLog) << " -> " << rowNames[i];      
         (*m_osLog) << " LB= "   << UtilDblToStr(rowLB[i],precision)
                    << " ax= "   << UtilDblToStr(ax[i],precision)
                    << " UB= "   << UtilDblToStr(rowUB[i],precision) << endl;
	 //UtilPrintPackedVector(row);	 
      }
      if(UtilIsZero(ax[i],1.0e-3)            ||
	 (ax[i] < 0 && UtilIsZero(rowLB[i])) ||
	 (ax[i] > 0 && UtilIsZero(rowUB[i])))
	 relViol = actViol;
      else
	 relViol = actViol / std::fabs(ax[i]);      
      if(relViol > 0.005){//0.5% violated      
	 (*m_osLog) << "Point violates row " << i;
	 if(hasRowNames)
	    (*m_osLog) << " -> " << rowNames[i];
	 (*m_osLog) << " LB= " << UtilDblToStr(rowLB[i],precision)
		    << " ax= "  << UtilDblToStr(ax[i],precision)
		    << " UB= " << UtilDblToStr(rowUB[i],precision) 
		    << " RelViol= " << UtilDblToStr(relViol,precision)
		    << endl;
	 //>5% violation is probably a bug, but <5% could be just
	 //  round off error??? not sure about that
	 if(relViol > 0.05){
	    isFeas = false;
	 }
      }
   }
   UTIL_DELARR(ax);
   return isFeas;
}

//===========================================================================//
void DecompAlgo::checkMasterDualObj(){
				     
   int    r;
   const int      nRows     = m_masterSI->getNumRows();
   const double * rowRhs    = m_masterSI->getRightHandSide();
   const double * dual      = m_masterSI->getRowPrice();
   const double   primalObj = m_masterSI->getObjValue();
   double dualObj = 0.0;
   for(r = 0; r < nRows; r++){
      dualObj += dual[r] * rowRhs[r];
   }
   UTIL_DEBUG(m_param.LogDebugLevel, 4,
	      (*m_osLog)
	      << "checkMasterDualObj" 
	      << setw(10) << "primalObj="
	      << setw(10) << UtilDblToStr(primalObj,3)
	      << setw(10) << "dualObj="
	      << setw(10) << UtilDblToStr(dualObj, 3) << endl;
	      );
   double actViol = std::fabs(primalObj - dualObj);
   double relViol = actViol;
   if(!UtilIsZero(dualObj,1.0e-3)){
      relViol = actViol / std::fabs(dualObj);
   }
   if(relViol > 1.0e-4){
      cerr << "checkMasterDualObj" 
           << setw(10) << "primalObj="
           << setw(10) << UtilDblToStr(primalObj,3)
           << setw(10) << "dualObj="
           << setw(10) << UtilDblToStr(dualObj, 3) << endl;      
      throw UtilException("primal and dual obj do not match", 
                          "checkMasterDualObj", "DecompAlgo");
   }
}

//===========================================================================//
bool DecompAlgo::isDualRayInfProof(const double           * dualRay,
				    const CoinPackedMatrix * rowMatrix,
				    const double           * colLB,
				    const double           * colUB,
				    const double           * rowRhs,
				    ostream                * os){
   //---
   //--- Does dualRay provide a proof according to Farkas Lemma?
   //---    yA >= 0, yb < 0, or
   //---    yA <= 0, yb > 0 ??
   //---
   
   int      i;
   double   yb;
   bool     isProof = true;
   bool     ybPos   = true;

   double * yA = 0;  
   const int m = rowMatrix->getNumRows();
   const int n = rowMatrix->getNumCols();
      
   //y^T b
   yb = 0.0;
   for(i = 0; i < m; i++){
      yb += dualRay[i] * rowRhs[i];
      if(os){
         (*os) << "i : " << i << " dualRay = " << dualRay[i]
               << " rowRhs = " << rowRhs[i] << " yb = " << yb << endl;
      }
   }

   //TODO: tol
   if(yb > 1.0e-10)
      ybPos = true;
   else if(yb < -1.0e-10)
      ybPos = false;
   else
      return isProof;
   
   yA = new double[n];
   rowMatrix->transposeTimes(dualRay, yA);     //y^T A
   
   for(i = 0; i < n; i++){
      if(os){
         (*os) << "yA[" << i << "]:\t" << yA[i];
      }
      //TODO: tol 1.0e-6 is too tight?
      if((ybPos  && (yA[i] >  1.0e-2)) ||
	 (!ybPos && (yA[i] < -1.0e-2))){
	 if(os)
	    (*os) << " -->isProof (false)" << endl;
         isProof = false;
      }
      else
         if(os)
            (*os) << endl;
   }
   
   UTIL_DELARR(yA);


#if 0
   //sanity check
   if(!isProof)
      isProof 
	 = isDualRayInfProofCpx(dualRay, rowMatrix, colLB, colUB, rowRhs, os);
#endif
   return isProof;
}

//===========================================================================//
bool DecompAlgo::isDualRayInfProofCpx(const double           * dualRay,
				       const CoinPackedMatrix * rowMatrix,
				       const double           * colLB,
				       const double           * colUB,
				       const double           * rowRhs,
				       ostream                * os){
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

   double * yA = 0;
   double *  z = 0;
  
   const int m = rowMatrix->getNumRows();
   const int n = rowMatrix->getNumCols();

   //TODO: check for out-of-mem conditions?
   yA = new double[n];
   UtilFillN(yA, n, 0.0);
   
   double * yA2 = new double[n];
   rowMatrix->transposeTimes(dualRay, yA2);     //y^T A
   for(i = 0; i < m; i++) {
      double yA_i = 0;
      CoinShallowPackedVector pv = rowMatrix->getVector(i);
      const int    * indI = pv.getIndices();
      const double * elsI = pv.getElements();
      const int      lenI = pv.getNumElements();     
      for(int j = 0; j < lenI; j++){	
	 yA_i += dualRay[indI[j]] * elsI[j];
	 printf("i: %d, j: %d, indIj: %d, elsIj: %g ray: %g yA_i: %g\n",
		i, j, indI[j], elsI[j], dualRay[indI[j]], yA_i);
      }
      yA[i] = yA_i;
      if(!UtilIsZero(yA[i] - yA2[i])){
	 printf(" ---> yA: %g, yA2: %g\n", yA[i], yA2[i]);
      }
      fflush(stdout);
      CoinAssert(UtilIsZero(yA[i] - yA2[i]));
   }
   
   z  = new double[n];
   for(j = 0; j < n; j++){
      if(yA[j] >= 0) z[j] = CoinMin(1.0e20, colUB[j]);
      else           z[j] = colLB[j];
   }

   //y^T b
   yb = 0.0;
   for(i = 0; i < m; i++){
      yb += dualRay[i] * rowRhs[i];
      if(os)
	 (*os) << "\ni : " << i << " dualRay = " << dualRay[i]
	       << " rowRhs = " << rowRhs[i] << " yb = " << yb;
   }

   //y^T A z
   yAz = 0.0;
   for(j = 0; j < n; j++){
      yAz += yA[j] * z[j];
      if(os)
	 (*os) << "\nj : " << j << " yA = " << yA[j]
	       << " z = " << z[j] << " yAz = " << yAz;
   }

   if(os)
      (*os) << "\nyb - yAz = " << yb - yAz << endl;

   UTIL_DELARR(yA);
   UTIL_DELARR(z);

   //TODO: tol
   if(yb - yAz > 1.0e-3)
      return true;
   else
      return false;
}

//===========================================================================//
void DecompAlgo::printBasisInfo(OsiSolverInterface * si,
                                ostream             * os){
   int      b, r, c;
   int    * basics   = 0;
   int    * rstat    = 0;
   int    * cstat    = 0;
   double * bInvRow  = 0;
   double * bInvARow = 0;
   const int n  = si->getNumCols();
   const int m  = si->getNumRows();
   char type[4] = {'F','B','U','L'};
   //TODO: have to check sense?
   const double * rowRhs = si->getRightHandSide();
   
   basics   = new int[m];
   bInvRow  = new double[m];
   bInvARow = new double[n];
   rstat    = new int[m];
   cstat    = new int[n];
   
   si->enableSimplexInterface(false);
   si->getBasics(basics);
   (*os) << "\n\nBasics: ";
   for(b = 0; b < m; b++)
      (*os) << basics[b] << " ";
   si->getBasisStatus(cstat,rstat);
   
   (*os) << "\ncstat: ";
   for(c = 0; c < n; c++){
      (*os) << type[cstat[c]];
   }
   (*os) << "\n";
   (*os) << "rstat: ";
   for(r = 0; r < m; r++){
      (*os) << type[rstat[r]];
   }
   (*os) << "\n";

   //yb, where y is a row of B-1
   double yb = 0.0;

   (*os) << "\nB-1:";
   for(r = 0; r < m; r++){
      yb = 0.0;
      si->getBInvRow(r, bInvRow);
      (*os) << "\nB-1Row r: " << r << ": ";
      for(b = 0; b < m; b++){
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
   for(r = 0; r < m; r++){
      si->getBInvARow(r, bInvARow);
      (*os) << "\nB-1ARow r: " << r << ": ";
      allpos = true;
      allneg = true;
      for(c = 0; c < n; c++){
	 (*os) << bInvARow[c] << " ";
	 if(bInvARow[c] < 0)
	    allpos = false;
	 if(bInvARow[c] > 0)
	    allneg = false;
      }
      if(allpos)
	 (*os) << " ---> allpos";
      if(allneg)
	 (*os) << " ---> allneg";
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
void DecompAlgo::printCurrentProblem(const OsiSolverInterface * si,
                                     const string               baseName,
                                     const int                  nodeIndex,
                                     const int                  cutPass,
                                     const int                  pricePass,
                                     const int                  blockId,
                                     const bool                 printMps,
                                     const bool                 printLp){
  if(!si)
     return;
   string filename = DecompAlgoStr[m_algo] + "_" + baseName 
      + ".n" + UtilIntToStr(nodeIndex)
      + ".c" + UtilIntToStr(cutPass) 
      + ".p" + UtilIntToStr(pricePass);
   if(blockId != -1)
      filename += ".b" + UtilIntToStr(blockId);
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "printCurrentProblem()", m_param.LogDebugLevel, 2);
   
   UTIL_DEBUG(m_param.LogDebugLevel, 3,
	      if(printMps)
	      (*m_osLog) << "calling writeMps filename = " << filename << endl;
	      if(printLp)
	      (*m_osLog) << "calling writeLp  filename = " << filename << endl;
	      );
   if(printMps){
      si->writeMpsNative(filename.c_str(), NULL, NULL, 2);
      //si->writeMps(filename.c_str());
   }
   if(printLp)      
      si->writeLp(filename.c_str(), "lp", 1e-30, 5, 10);
   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "printCurrentProblem()", m_param.LogDebugLevel, 2);
}

//===========================================================================//
void DecompAlgo::printCurrentProblem(const OsiSolverInterface * si,
                                     const string               fileName,
                                     const bool                 printMps,
                                     const bool                 printLp){
   string filename = fileName; 
   if(printMps)
      si->writeMps(filename.c_str());
   if(printLp)
      si->writeLp(filename.c_str(), "lp", 1e-30, 5, 10);
}


//===========================================================================//
void DecompAlgo::printVars(ostream * os){
   
   DecompVarList::iterator it;
   int var_index = 0;
   for(it = m_vars.begin(); it != m_vars.end(); it++){
      (*os) << "VAR " << var_index++ << " : ";
      (*it)->print(os, m_app);
      (*os) << endl;
   }
   (*os) << endl;
}

//===========================================================================//
void DecompAlgo::createFullMps(const string fileName){
   CoinAssert(m_algo == CUT);

   DecompConstraintSet          * modelCore   = m_modelCore.getModel();   
   int n_integerVars = static_cast<int>(modelCore->integerVars.size());
   m_masterSI->setInteger(&modelCore->integerVars[0], n_integerVars);
   m_masterSI->writeMps(fileName.c_str());
   m_masterSI->setContinuous(&modelCore->integerVars[0], n_integerVars);
}

//===========================================================================//
void DecompAlgo::printCuts(ostream * os){
   DecompCutList::iterator it;
   int cut_index = 0;
   for(it = m_cuts.begin(); it != m_cuts.end(); it++){
      (*os) << "CUT " << cut_index++ << " : ";
      (*it)->print(os);
   }
   (*os) << endl;
}

//===========================================================================//
void DecompAlgoC::solveDirect(int                  timeLimit,
			      DecompSolverResult * result){

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

   DecompConstraintSet * modelCore = m_modelCore.getModel();
   int                   numInts   = modelCore->getNumInts();
   int                   numCols   = m_masterSI->getNumCols();


   //---
   //--- start timer
   //---
   UtilTimer timer; timer.start();
   
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
   for(i = 0; i < numInts; i++)
      m_masterSI->setInteger(modelCore->integerVars[i]);
   
   //---
   //--- dump full milp
   //---
   if(m_param.LogDumpModel >= 2){
      string fileName = "directMILP";
      printCurrentProblem(m_masterSI, fileName);
   }


#ifdef __DECOMP_IP_CBC__
   CbcModel cbc(*m_masterSI);
   CbcMain0(cbc);

   //---
   //--- build argument list
   //---
   const char * argv[20];
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
   const int statusSet[2] = {0,1};
   int       solStatus    = cbc.status();
   int       solStatus2   = cbc.secondaryStatus();
   if(!UtilIsInSet(solStatus, statusSet, 2)){
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
   if(cbc.isProvenOptimal() || cbc.isSecondsLimitReached()){      
      objUB = cbc.getObjValue();
      if(result && cbc.getSolutionCount()){
         const double * solution = cbc.getColSolution();
         copy(solution, solution+numCols, result->m_solution);
      }
   }

   //---
   //--- copy sol status into result
   //---
   if(result){
      result->m_solStatus  = solStatus;
      result->m_solStatus2 = solStatus2;
   }
#endif


#ifdef __DECOMP_IP_CPX__
   OsiIpSolverInterface * masterSICpx = 
      dynamic_cast<OsiIpSolverInterface*>(m_masterSI);
   CPXLPptr  cpxLp  = masterSICpx->getLpPtr();
   CPXENVptr cpxEnv = masterSICpx->getEnvironmentPtr();

   //---
   //--- set the time limit
   //---
   status = CPXsetdblparam(cpxEnv, CPX_PARAM_TILIM, timeLimit);
   if(status)
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
   if(result){
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
   if(status)
      throw UtilException("CPXgetbestobjval failure", 
                          "solveDirect", "DecompAlgoC");      
   
   //---
   //--- get objective and solution
   //---
   if(solStatus == CPXMIP_OPTIMAL     ||
      solStatus == CPXMIP_OPTIMAL_TOL ||
      solStatus == CPXMIP_TIME_LIM_FEAS){

      status = CPXgetmipobjval(cpxEnv, cpxLp, &objUB);
      if(status)
         throw UtilException("CPXgetmipobjval failure", 
                             "solveDirect", "DecompAlgoC");

      if(result){	 
         const double * solution = m_masterSI->getColSolution();
	 copy(solution, solution+numCols, result->m_solution);
      }
   }

   //---
   //--- copy sol status into result
   //---
   if(result){
      result->m_solStatus  = solStatus;
      result->m_solStatus2 = 0;
   }

#endif

   //---
   //--- copy bounds into result
   //---
   if(result){
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
}

