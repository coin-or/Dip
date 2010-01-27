//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)       //
//                                                                           //
// Copyright (C) 2002-2009, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

//===========================================================================//
#include "DecompModel.h"
#include "DecompSolverResult.h"
//===========================================================================//

//===========================================================================//
bool DecompAlgoModel::isPointFeasible(const double * x,
                                      const int      logLevel,
                                      const double   feasVarTol,
                                      const double   feasConTol){
                                         
   DecompConstraintSet    * model    = getModel();
   const CoinPackedMatrix * M        = model->getMatrix();
   const  vector<string>  & colNames = model->getColNames();
   const  vector<string>  & rowNames = model->getRowNames();
   if(!M)
      return true;

   int    c, r;
   int    precision   = 7;
   bool   isFeas      = true;
   bool   hasColNames = false;
   bool   hasRowNames = false;
   double ax          = 0.0;
   double clb         = 0.0;
   double cub         = 0.0;  
   double rlb         = 0.0;
   double rub         = 0.0;
   double actViol     = 0.0;
   double relViol     = 0.0;
   if(colNames.size())
      hasColNames = true;
   if(rowNames.size())
      hasRowNames = true;
   double feasVarTol10 = 10 * feasVarTol;
   double feasConTol10 = 10 * feasConTol;

   //---
   //--- do we satisfy all (active) column bounds
   //---
   vector<int>::const_iterator it;
   const vector<int> & activeColumns  = model->getActiveColumns();
   for(it = activeColumns.begin(); it != activeColumns.end(); it++){
      c        = *it;
      clb      = model->colLB[c];
      cub      = model->colUB[c];      
      actViol = std::max<double>(clb - x[c], x[c] - cub);
      actViol = std::max<double>(actViol, 0.0);
      if(UtilIsZero(x[c], feasVarTol) ||
	 (x[c] < 0 && UtilIsZero(clb)) ||
         (x[c] > 0 && UtilIsZero(cub)))
	 relViol = actViol;
      else
	 relViol = actViol / std::fabs(x[c]);      
      if(relViol > feasVarTol){
	 //Notify, but don't mark in feasible unless 10x worse.
	 UTIL_DEBUG(logLevel, 4,
                    cout << "Point violates column " << c;
                    if(hasColNames)
                       cout << " -> " << colNames[c];
                    cout << " LB= " << UtilDblToStr(clb,precision)
		    << " x= "  << UtilDblToStr(x[c],precision)
		    << " UB= " << UtilDblToStr(cub,precision) 
		    << " RelViol= " << UtilDblToStr(relViol,precision)
		    << endl;
		    );
	 if(relViol > feasVarTol10){
	    isFeas = false;
	    goto FUNC_EXIT;
	 }
      }      
   }
   
   //---
   //--- do we satisfy all row bounds
   //---   TODO: for core model, this includes branching rows and cuts
   //---         we can actually get away with just checking the 
   //---         original base rows
   //---
   for(r = 0; r < model->getNumRows(); r++){
      ax  = model->M->getVector(r).dotProduct(x);
      rlb = model->rowLB[r];
      rub = model->rowUB[r];
      actViol = std::max<double>(rlb - ax, ax - rub);     
      actViol = std::max<double>(actViol, 0.0);
      //printf("CORE r:%d rlb:%g ax:%g rub:%g actViol:%g\n", 
      //     r, rlb, ax, rub, actViol);
      if(UtilIsZero(ax, feasConTol)  ||
	 (ax < 0 && UtilIsZero(rlb)) ||
	 (ax > 0 && UtilIsZero(rub)))
	 relViol = actViol;
      else
	 relViol = actViol / std::fabs(ax);      
      if(relViol > feasConTol){
	 //Notify, but don't mark in feasible unless 10x worse.
	 UTIL_DEBUG(logLevel, 4,
		    cout 
		    << "Point violates row ax[" << r << "]: " 
		    << UtilDblToStr(ax)
		    << " LB: " << UtilDblToStr(rlb) 
		    << " UB: " << UtilDblToStr(rub) 
		    << " RelViol: " << UtilDblToStr(relViol)
		    << endl;
		    );	
	 if(relViol > feasConTol10){
	    isFeas = false;
	    goto FUNC_EXIT;
	 }
      }
   }

 FUNC_EXIT:
   UTIL_DEBUG(logLevel, 4,
	      cout << "isPointFeasible = " << isFeas << endl;
	      );
   return isFeas;
}

//===========================================================================//
void DecompAlgoModel::solveOsiAsIp(DecompSolverResult * result,
                                   DecompParam        & param,
                                   bool                 doExact,
                                   bool                 doCutoff,
                                   bool                 isRoot,
                                   double               cutoff){
   
   const int numCols    = m_osi->getNumCols();
   const int logIpLevel = param.LogLpLevel;

#ifdef __DECOMP_IP_CBC__
   //TODO: what exactly does this do? make copy of entire model!?
   CbcModel cbc(*m_osi);
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
   string cbcGap       = "-allowableGap";
   string cbcGapSet    = "0";
   string cbcCutoff    = "-cutoff";   
   string cbcCutoffSet = UtilDblToStr(cutoff);
   string cbcSLog      = "-slog";
   string cbcSLogSet   = "2";
   if(doExact)
      cbcGapSet = UtilDblToStr(param.SubProbGapLimitExact);
   else
      cbcGapSet = UtilDblToStr(param.SubProbGapLimitInexact);
   argv[argc++] = cbcExe.c_str();
   argv[argc++] = cbcLog.c_str();
   argv[argc++] = cbcLogSet.c_str();      
   //argv[argc++] = cbcSLog.c_str();    //for extra debugging
   //argv[argc++] = cbcSLogSet.c_str(); //for extra debugging
   argv[argc++] = cbcGap.c_str();   
   argv[argc++] = cbcGapSet.c_str();
   if(doCutoff){
      argv[argc++] = cbcCutoff.c_str();
      argv[argc++] = cbcCutoffSet.c_str();
   }
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
   /** Final status of problem.
    *   -1  before branchAndBound
    *    0  finished - check isProvenOptimal or isProvenInfeasible 
    *         to see if solution found (or check value of best solution)
    *    1  stopped - on maxnodes, maxsols, maxtime
    *    2  difficulties so run was abandoned
    *   (5  event user programmed event occurred)
   */   
   const int statusSet[2] = {0,1};
   result->m_solStatus    = cbc.status();
   if(!UtilIsInSet(result->m_solStatus, statusSet, 2)){
      cerr << "Error: CBC IP solver status = " << result->m_solStatus << endl;
      throw UtilException("CBC solver status", 
                          "solveOsiAsIp", "DecompAlgoModel");
   }

   /** Secondary status of problem
    *   -1 unset (status_ will also be -1)
    *    0 search completed with solution
    *    1 linear relaxation not feasible (or worse than cutoff)
    *    2 stopped on gap
    *    3 stopped on nodes
    *    4 stopped on time
    *    5 stopped on user event
    *    6 stopped on solutions
    *    7 linear relaxation unbounded
   */
   int       nSeta = 0;
   int       nSetb = 0;
   const int statusSet2a[2] = {0,2};   nSeta=2;
   const int statusSet2b[3] = {0,1,2}; nSetb=3;
   result->m_solStatus2 = cbc.secondaryStatus();

   //---
   //--- In root the subproblem should not be infeasible
   //---   unless due to cutoff. But, after branching it
   //---   can be infeasible.
   //---
   if(!doCutoff && isRoot){
      if(!UtilIsInSet(result->m_solStatus2, statusSet2a, nSeta)){
         cerr << "Error: CBC IP solver 2nd status = " 
              << result->m_solStatus2 << endl;
         throw UtilException("CBC solver 2nd status", 
                             "solveOsiAsIp", "DecompAlgoModel");
      }
   }
   else{
      if(!UtilIsInSet(result->m_solStatus2, statusSet2b, nSetb)){
         cerr << "Error: CBC IP solver 2nd status = " 
              << result->m_solStatus2 << endl;
         throw UtilException("CBC solver 2nd status", 
                             "solveOsiAsIp", "DecompAlgoModel");
      }
   }
   
   //---
   //--- update results object 
   //---
   result->m_nSolutions = 0;
   result->m_isOptimal  = false;
   result->m_isCutoff   = false;
   //printf("cbc.isProvenOptimal() = %d\n", cbc.isProvenOptimal());
   if(cbc.isProvenOptimal()){
      result->m_nSolutions = 1;
      result->m_isOptimal  = true;      
   }
   else{
      if(cbc.isProvenInfeasible()){
         result->m_nSolutions = 0;
         result->m_isCutoff   = doCutoff;
         result->m_isOptimal  = true;         
      }
      else{
         //---
         //--- else it must have stopped on gap
         //---
         result->m_nSolutions = 1;
         result->m_isCutoff   = doCutoff;
         result->m_isOptimal  = false;      
      }
   }
   
   //---
   //--- get copy of solution
   //---
   result->m_objLB = cbc.getBestPossibleObjValue();
   if(result->m_nSolutions >= 1){
      result->m_objUB = cbc.getObjValue();
      memcpy(result->m_solution, 
             cbc.getColSolution(), numCols * sizeof(double));
   }
#endif





#ifdef __DECOMP_IP_CPX__
   //---
   //--- get OsiCpx object from Osi object
   //--- get CPEXENVptr for use with internal methods
   //--- get CPXLPptr   for use with internal methods
   //---
   OsiCpxSolverInterface * osiCpx 
      = dynamic_cast<OsiCpxSolverInterface*>(m_osi);
   CPXENVptr cpxEnv = osiCpx->getEnvironmentPtr();
   CPXLPptr  cpxLp  = osiCpx->getLpPtr();
   assert(cpxEnv && cpxLp);
  
   //---
   //--- set parameters
   //---
   int status = 0;
   if(logIpLevel){
      status = CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_ON);
      if(status)
	 throw UtilException("CPXsetintparam failure", 
			     "solveOsiAsIp", "DecompAlgoModel");
      status = CPXsetintparam(cpxEnv, CPX_PARAM_SIMDISPLAY, logIpLevel);
      if(status)
	 throw UtilException("CPXsetintparam failure", 
			     "solveOsiAsIp", "DecompAlgoModel");
   }
   else{
      status = CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);
      if(status)
	 throw UtilException("CPXsetintparam failure", 
			     "solveOsiAsIp", "DecompAlgoModel");
   }
   
   if(doExact)
      status = CPXsetdblparam(cpxEnv, CPX_PARAM_EPGAP, 
			      param.SubProbGapLimitExact);
   else
      status = CPXsetdblparam(cpxEnv, CPX_PARAM_EPGAP, 
			      param.SubProbGapLimitInexact);
   if(status)
      throw UtilException("CPXsetdblparam failure", 
			  "solveOsiAsIp", "DecompAlgoModel");
   
   if(doCutoff)
      status = CPXsetdblparam(cpxEnv, CPX_PARAM_CUTUP, cutoff);
   else
      status = CPXsetdblparam(cpxEnv, CPX_PARAM_CUTUP, 1.0e+75);
   if(status)
      throw UtilException("CPXsetdblparam failure", 
			  "solveOsiAsIp", "DecompAlgoModel");

   //---
   //--- starting with CPX12, parallel MIP is on by default
   //---   we do not want that (usually)
   //--- TODO: make this a user option
   //---
#if CPX_VERSION >= 1100
   status = CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, 1);
   if(status)
      throw UtilException("CPXsetdblparam failure", 
			  "solveOsiAsIp", "DecompAlgoModel");
#endif
   //---
   //--- solve the MILP
   //---
   osiCpx->branchAndBound();	 
     
   //---
   //--- get solver status
   //---
   result->m_solStatus  = CPXgetstat(cpxEnv, cpxLp);
   result->m_solStatus2 = 0;

   const int statusSet1[2] = {CPXMIP_OPTIMAL,
			      CPXMIP_OPTIMAL_TOL};
   const int statusSet2[3] = {CPXMIP_OPTIMAL,
			      CPXMIP_OPTIMAL_TOL,
			      CPXMIP_INFEASIBLE};
   if(!UtilIsInSet(result->m_solStatus, statusSet2, 3)){
      cerr << "Error: CPX IP solver status = " << result->m_solStatus << endl;
      throw UtilException("CPX solver status", 
                          "solveOsiAsIp", "DecompAlgoModel");
   }

   //---
   //--- In root the subproblem should not be infeasible
   //---   unless due to cutoff. But, after branching it
   //---   can be infeasible.
   //---
   if(!doCutoff && isRoot){
      if(!UtilIsInSet(result->m_solStatus, statusSet1, 2)){
         cerr << "Error: CPX IP solver 2nd status = " 
              << result->m_solStatus << endl;
         throw UtilException("CPX solver status", 
                             "solveOsiAsIp", "DecompAlgoModel");
      }
   }
   else{
      if(!UtilIsInSet(result->m_solStatus, statusSet2, 3)){
         cerr << "Error: CPX IP solver 2nd status = " 
              << result->m_solStatus << endl;
         throw UtilException("CPX solver status", 
                             "solveOsiAsIp", "DecompAlgoModel");
      }
   }
   
   //---
   //--- update results object 
   //---
   result->m_nSolutions = 0;
   result->m_isOptimal  = false;
   result->m_isCutoff   = false;
   if(result->m_solStatus == CPXMIP_OPTIMAL ||
      result->m_solStatus == CPXMIP_OPTIMAL_TOL){
      result->m_nSolutions = 1;
      result->m_isOptimal  = true;      
   }
   else{
      if(result->m_solStatus == CPXMIP_INFEASIBLE){
         result->m_nSolutions = 0;
         result->m_isCutoff   = doCutoff;
         result->m_isOptimal  = true;         
      }
      else{
         //---
         //--- else it must have stopped on gap
         //---
         result->m_nSolutions = 1;
         result->m_isCutoff   = doCutoff;
         result->m_isOptimal  = false;      
      }
   }
   
   //---
   //--- get copy of solution
   //---
   status = CPXgetbestobjval(cpxEnv, cpxLp, &result->m_objLB);
   if(status)
      throw UtilException("CPXgetbestobjval failure", 
			  "solveOsiAsIp", "DecompAlgoModel");   
   if(result->m_nSolutions >= 1){
      status = CPXgetmipobjval(cpxEnv, cpxLp, &result->m_objUB);
      if(status)
	 throw UtilException("CPXgetmipobjval failure", 
			     "solveOsiAsIp", "DecompAlgoModel");
      memcpy(result->m_solution, 
             osiCpx->getColSolution(), numCols * sizeof(double));
   }
#endif
}
