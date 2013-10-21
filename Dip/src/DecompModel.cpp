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
// Copyright (C) 2002-2013, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

//===========================================================================//
#include "DecompModel.h"
#include "DecompSolverResult.h"
//===========================================================================//

using namespace std;

//===========================================================================//
bool DecompAlgoModel::isPointFeasible(const double* x,
                                      const bool     isXSparse,
                                      const int      logLevel,
                                      const double   feasVarTol,
                                      const double   feasConTol)
{
   DecompConstraintSet*     model    = getModel();

   if (!model) {
      return true;
   }

   const CoinPackedMatrix* M        = model->getMatrix();

   if (!M) {
      return true;
   }

   const  vector<string>&   colNames = model->getColNames();
   const  vector<string>&   rowNames = model->getRowNames();
   int    c, r, i;
   bool   isFeas      = true;
   bool   hasColNames = false;
   bool   hasRowNames = false;
   double xj          = 0.0;
   double ax          = 0.0;
   double clb         = 0.0;
   double cub         = 0.0;
   double rlb         = 0.0;
   double rub         = 0.0;
   double actViol     = 0.0;
   double relViol     = 0.0;

   if (colNames.size()) {
      hasColNames = true;
   }

   if (rowNames.size()) {
      hasRowNames = true;
   }

   double feasVarTol10 = 10 * feasVarTol;
   double feasConTol10 = 10 * feasConTol;
   //---
   //--- do we satisfy all (active) column bounds
   //---
   vector<int> ::const_iterator it;
   map<int, int>::const_iterator mcit;
   const vector<int>&   activeColumns  = model->getActiveColumns();
   bool                 isSparse       = model->isSparse();
   const map<int, int>& origToSparse   = model->getMapOrigToSparse();
   const map<int, int>& sparseToOrig   = model->getMapSparseToOrig();

   for (it = activeColumns.begin(); it != activeColumns.end(); it++) {
      if (isSparse) {
         mcit = origToSparse.find(*it);
         c    = mcit->second;
         xj   = isXSparse ? x[c] : x[*it];
      } else {
         c  = *it;
         xj = x[c];
         assert(!isXSparse);
      }

      clb      = model->colLB[c];
      cub      = model->colUB[c];
      UTIL_DEBUG(logLevel, 5,
                 int    precision   = 7;

      if (!UtilIsZero(xj)) {
      cout << "Point " << c;

      if (hasColNames) {
            cout << " -> " << colNames[c];
         }

         cout << " LB= " << UtilDblToStr(clb, precision)
              << " x= "  << UtilDblToStr(xj, precision)
              << " UB= " << UtilDblToStr(cub, precision)
              << endl;
      }
                );
      actViol = std::max<double>(clb - xj, xj - cub);
      actViol = std::max<double>(actViol, 0.0);

      if (UtilIsZero(xj, feasVarTol) ||
            (xj < 0 && UtilIsZero(clb)) ||
            (xj > 0 && UtilIsZero(cub))) {
         relViol = actViol;
      } else {
         relViol = actViol / std::fabs(xj);
      }

      if (relViol > feasVarTol) {
         //Notify, but don't mark in feasible unless 10x worse.
         UTIL_DEBUG(logLevel, 4,
                    int    precision   = 7;
                    cout << "Point violates column " << c;

                    if (hasColNames)
                    cout << " -> " << colNames[c];
                    cout << " LB= " << UtilDblToStr(clb, precision)
                    << " x= "  << UtilDblToStr(xj, precision)
                    << " UB= " << UtilDblToStr(cub, precision)
                    << " RelViol= " << UtilDblToStr(relViol, precision)
                    << endl;
                   );

         if (relViol > feasVarTol10) {
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
   for (r = 0; r < model->getNumRows(); r++) {
      if (isSparse) {
         if (isXSparse) {
            ax  = model->M->getVector(r).dotProduct(x);
         } else {
            //TODO: "sparse" is the wrong word here - should be using
            // the word projected for param, etc...
            //---
            //--- x is with respect to the original set of columns,
            //---   but the matrix is the sparse version
            //---
            const CoinShallowPackedVector v = model->M->getVector(r);
            const int      len = v.getNumElements();
            const int*     ind = v.getIndices();
            const double* els = v.getElements();
            ax = 0.0;

            for (i = 0; i < len; i++) {
               mcit = sparseToOrig.find(ind[i]);
               c    = mcit->second;
               ax  += x[c] * els[i];
            }
         }
      } else {
         ax  = model->M->getVector(r).dotProduct(x);
         assert(!isXSparse);
      }

      rlb = model->rowLB[r];
      rub = model->rowUB[r];
      actViol = std::max<double>(rlb - ax, ax - rub);
      actViol = std::max<double>(actViol, 0.0);

      //printf("CORE r:%d rlb:%g ax:%g rub:%g actViol:%g\n",
      //     r, rlb, ax, rub, actViol);
      if (UtilIsZero(ax, feasConTol)  ||
            (ax < 0 && UtilIsZero(rlb)) ||
            (ax > 0 && UtilIsZero(rub))) {
         relViol = actViol;
      } else {
         relViol = actViol / std::fabs(ax);
      }

      if (relViol > feasConTol) {
         //Notify, but don't mark in feasible unless 10x worse.
         UTIL_DEBUG(logLevel, 4,
                    cout
                    << "Point violates row ";

                    if (hasRowNames)
                    cout << " -> " << rowNames[r];
                    cout << " ax[" << r << "]: "
                    << UtilDblToStr(ax)
                    << " LB: " << UtilDblToStr(rlb)
                    << " UB: " << UtilDblToStr(rub)
                    << " RelViol: " << UtilDblToStr(relViol)
                    << endl;
                   );

         if (relViol > feasConTol10) {
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
void DecompAlgoModel::solveOsiAsIp(DecompSolverResult* result,
                                   DecompParam&         param,
                                   bool                 doExact,
                                   bool                 doCutoff,
                                   bool                 isRoot,
                                   double               cutoff)
{
   const int numCols    = m_osi->getNumCols();
   const int logIpLevel = param.LogLpLevel;
   double* solution = new double[numCols];
   assert(solution);
   //---
   //--- clear out any old solutions
   //---
   result->m_solution.clear();
#ifdef __DECOMP_IP_CBC__
   //TODO: what exactly does this do? make copy of entire model!?
   CbcModel cbc(*m_osi);
#if 0
   CbcMain0(cbc);
   //int i;
   //const double * colUB = cbc.getColUpper();
   //for(i = 0; i < cbc.getNumCols(); i++){
   //  printf("col %d -> ub: %g\n",
   //         i, colUB[i]);
   //}
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
   string cbcGap       = "-ratio";
   string cbcGapSet    = "0";
   string cbcTime      = "-seconds";
   string cbcTimeSet   = "0";
   string cbcCutoff    = "-cutoff";
   string cbcCutoffSet = UtilDblToStr(cutoff, -1, COIN_DBL_MAX);
   string cbcSLog      = "-slog";
   string cbcSLogSet   = "2";

   if (doExact) {
      cbcTimeSet = UtilDblToStr(param.SubProbTimeLimitExact, -1, COIN_DBL_MAX);
      cbcGapSet  = UtilDblToStr(param.SubProbGapLimitExact, -1, COIN_DBL_MAX);
   } else {
      cbcTimeSet = UtilDblToStr(param.SubProbTimeLimitInexact, -1, COIN_DBL_MAX);
      cbcGapSet  = UtilDblToStr(param.SubProbGapLimitInexact, -1, COIN_DBL_MAX);
   }

   bool   doTime       = false;
   double cbcMaxSecUB  = 1e100;

   if (doExact) {
      if (param.SubProbTimeLimitExact < cbcMaxSecUB) {
         doTime = true;
      }
   } else {
      if (param.SubProbTimeLimitInexact < cbcMaxSecUB) {
         doTime = true;
      }
   }

   argv[argc++] = cbcExe.c_str();
   argv[argc++] = cbcLog.c_str();
   argv[argc++] = cbcLogSet.c_str();
   //argv[argc++] = cbcSLog.c_str();    //for extra debugging
   //argv[argc++] = cbcSLogSet.c_str(); //for extra debugging
   argv[argc++] = cbcGap.c_str();
   argv[argc++] = cbcGapSet.c_str();

   if (doTime) {
      argv[argc++] = cbcTime.c_str();
      argv[argc++] = cbcTimeSet.c_str();
   }

   if (doCutoff) {
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
   #endif
   cbc.setLogLevel(0);
   cbc.branchAndBound();
   const int statusSet[2] = {0, 1};
   result->m_solStatus    = cbc.status();

   if (!UtilIsInSet(result->m_solStatus, statusSet, 2)) {
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
   const int statusSet2a[4] = {0, 2, 3, 4};
   nSeta = 4;
   const int statusSet2b[5] = {0, 1, 2, 4, 5};
   nSetb = 5;
   result->m_solStatus2 = cbc.secondaryStatus();

   //---
   //--- In root the subproblem should not be infeasible
   //---   unless due to cutoff. But, after branching it
   //---   can be infeasible.
   //---
   if (!doCutoff && isRoot) {
      if (!UtilIsInSet(result->m_solStatus2, statusSet2a, nSeta)) {
         cerr << "Error: CBC IP solver 2nd status = "
              << result->m_solStatus2 << endl;
         throw UtilException("CBC solver 2nd status",
                             "solveOsiAsIp", "DecompAlgoModel");
      }
   } else {
      if (!UtilIsInSet(result->m_solStatus2, statusSet2b, nSetb)) {
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
   if (cbc.isProvenOptimal()) {
      result->m_nSolutions = 1;
      result->m_isOptimal  = true;
   } else {
      if (cbc.isProvenInfeasible()) {
         result->m_nSolutions = 0;
         result->m_isCutoff   = doCutoff;
         result->m_isOptimal  = true;
      } else {
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

   if (result->m_nSolutions >= 1) {
      result->m_objUB = cbc.getObjValue();
      const double* solDbl = cbc.getColSolution();
      vector<double> solVec(solDbl, solDbl + numCols);
      result->m_solution.push_back(solVec);
      //memcpy(result->m_solution,
      //  cbc.getColSolution(), numCols * sizeof(double));
      assert(result->m_nSolutions ==
             static_cast<int>(result->m_solution.size()));
   }

#endif
#ifdef __DECOMP_IP_CPX__
   //---
   //--- get OsiCpx object from Osi object
   //--- get CPEXENVptr for use with internal methods
   //--- get CPXLPptr   for use with internal methods
   //---
   OsiCpxSolverInterface* osiCpx
      = dynamic_cast<OsiCpxSolverInterface*>(m_osi);
   CPXENVptr cpxEnv = osiCpx->getEnvironmentPtr();
   CPXLPptr  cpxLp  = osiCpx->getLpPtr();
   assert(cpxEnv && cpxLp);
   //---
   //--- set parameters
   //---
   int status = 0;

   if (logIpLevel) {
      status = CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_ON);

      if (status)
         throw UtilException("CPXsetintparam failure",
                             "solveOsiAsIp", "DecompAlgoModel");

      status = CPXsetintparam(cpxEnv, CPX_PARAM_SIMDISPLAY, logIpLevel);

      if (status)
         throw UtilException("CPXsetintparam failure",
                             "solveOsiAsIp", "DecompAlgoModel");
   } else {
      status = CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);

      if (status)
         throw UtilException("CPXsetintparam failure",
                             "solveOsiAsIp", "DecompAlgoModel");
   }

   if (doExact)
      status = CPXsetdblparam(cpxEnv, CPX_PARAM_EPGAP,
                              param.SubProbGapLimitExact);
   else
      status = CPXsetdblparam(cpxEnv, CPX_PARAM_EPGAP,
                              param.SubProbGapLimitInexact);

   if (status)
      throw UtilException("CPXsetdblparam failure",
                          "solveOsiAsIp", "DecompAlgoModel");

   if (doExact) {
      if (param.SubProbTimeLimitExact < COIN_DBL_MAX) {
         status = CPXsetdblparam(cpxEnv, CPX_PARAM_TILIM,
                                 param.SubProbTimeLimitExact);
      }
   } else {
      if (param.SubProbTimeLimitInexact < COIN_DBL_MAX) {
         status = CPXsetdblparam(cpxEnv, CPX_PARAM_TILIM,
                                 param.SubProbTimeLimitInexact);
      }
   }

   if (status)
      throw UtilException("CPXsetdblparam failure",
                          "solveOsiAsIp", "DecompAlgoModel");

   if (doCutoff) {
      status = CPXsetdblparam(cpxEnv, CPX_PARAM_CUTUP, cutoff);
   } else {
      status = CPXsetdblparam(cpxEnv, CPX_PARAM_CUTUP, 1.0e+75);
   }

   if (status)
      throw UtilException("CPXsetdblparam failure",
                          "solveOsiAsIp", "DecompAlgoModel");

   //---
   //--- starting with CPX12, parallel MIP is on by default
   //---   we do not want that (usually)
   //--- TODO: make this a user option
   //---
#if CPX_VERSION >= 1100
   status = CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, param.SubProbNumThreads);

   if (status)
      throw UtilException("CPXsetdblparam failure",
                          "solveOsiAsIp", "DecompAlgoModel");

   int startAlgo = 0;

   switch (param.SubProbSolverStartAlgo) {
   case DecompDualSimplex:
      startAlgo = CPX_ALG_DUAL;
      break;

   case DecompPrimSimplex:
      startAlgo = CPX_ALG_PRIMAL;
      break;

   case DecompBarrier:
      startAlgo = CPX_ALG_BARRIER;
      break;
   }

   status = CPXsetintparam(cpxEnv, CPX_PARAM_STARTALG, startAlgo);

   if (status)
      throw UtilException("CPXsetdblparam failure",
                          "solveOsiAsIp", "DecompAlgoModel");

   //---
   //--- check the mip starts solution pool, never let it get too
   //---   big, and refresh it periodically - assuming that the last
   //---   ones in the list are the last ones used - which would have the
   //---   best potential to help warm start
   //--- never let it get bigger than 10 solutions,
   //---   when refresh - keep only last 2
   //---
   int nMipStarts = CPXgetnummipstarts(cpxEnv, cpxLp);

   if (nMipStarts > 10) {
      status = CPXdelmipstarts(cpxEnv, cpxLp, 0, nMipStarts - 3);

      if (status)
         throw UtilException("CPXdelmipstarts failure",
                             "solveOsiAsIp", "DecompAlgoModel");
   }

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
   //printf("cplex status  = %d\n", result->m_solStatus);
   //printf("cplex status2 = %d\n", result->m_solStatus2);
   const int statusSet1[3] = {CPXMIP_OPTIMAL,
                              CPXMIP_OPTIMAL_TOL, //for stopping on gap
                              CPXMIP_TIME_LIM_FEAS
                             };
   const int statusSet2[4] = {CPXMIP_OPTIMAL,
                              CPXMIP_OPTIMAL_TOL, //for stopping on gap
                              CPXMIP_TIME_LIM_FEAS,
                              CPXMIP_INFEASIBLE
                             };

   if (!UtilIsInSet(result->m_solStatus, statusSet2, 4)) {
      cerr << "Error: CPX IP solver status = " << result->m_solStatus << endl;
      throw UtilException("CPX solver status",
                          "solveOsiAsIp", "DecompAlgoModel");
   }

   //---
   //--- In root the subproblem should not be infeasible
   //---   unless due to cutoff. But, after branching it
   //---   can be infeasible.
   //---
   if (!doCutoff && isRoot) {
      if (!UtilIsInSet(result->m_solStatus, statusSet1, 3)) {
         cerr << "Error: CPX IP solver 2nd status = "
              << result->m_solStatus << endl;
         throw UtilException("CPX solver status",
                             "solveOsiAsIp", "DecompAlgoModel");
      }
   } else {
      if (!UtilIsInSet(result->m_solStatus, statusSet2, 4)) {
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

   if (UtilIsInSet(result->m_solStatus, statusSet1, 3)) {
      int    i;
      int    nSols = CPXgetsolnpoolnumsolns(cpxEnv, cpxLp);
      double objVal;
      //printf("Number of solutions in solution pool = %d\n",
      //nSols);
      //TODO: currently just take up to the limit,
      //  but, should sort by objective and take n least?
      nSols = std::min<int>(nSols, param.SubProbNumSolLimit);

      for (i = 0; i < nSols; i++) {
         status = CPXgetsolnpoolobjval(cpxEnv, cpxLp, i, &objVal);

         if (status)
            throw UtilException("CPXgetsolnpoolobjval",
                                "solveOsiAsIp", "DecompAlgoModel");

         //printf("Sol %4d: Obj: %10g\n", i, objVal);
         status = CPXgetsolnpoolx(cpxEnv, cpxLp, i,
                                  solution, 0, numCols - 1);
         vector<double> solVec(solution, solution + numCols);
         result->m_solution.push_back(solVec);
         result->m_nSolutions++;
         assert(result->m_nSolutions ==
                static_cast<int>(result->m_solution.size()));
         //memcpy(result->m_solution,
         //	osiCpx->getColSolution(), numCols * sizeof(double));
      }

      result->m_nSolutions = nSols;
   }

   //printf("solStatus = %d\n", result->m_solStatus);
   if (result->m_solStatus == CPXMIP_OPTIMAL ||
         result->m_solStatus == CPXMIP_OPTIMAL_TOL) {
      result->m_isOptimal  = true;
   } else {
      if (result->m_solStatus == CPXMIP_INFEASIBLE) {
         result->m_nSolutions = 0;
         result->m_isCutoff   = doCutoff;
         result->m_isOptimal  = true;
      } else {
         //---
         //--- else it must have stopped on gap or time
         //---
         result->m_isCutoff   = doCutoff;
         result->m_isOptimal  = false;
      }
   }

   //---
   //--- get copy of solution
   //---
   status = CPXgetbestobjval(cpxEnv, cpxLp, &result->m_objLB);

   if (status)
      throw UtilException("CPXgetbestobjval failure",
                          "solveOsiAsIp", "DecompAlgoModel");

   if (result->m_nSolutions >= 1) {
      status = CPXgetmipobjval(cpxEnv, cpxLp, &result->m_objUB);

      if (status)
         throw UtilException("CPXgetmipobjval failure",
                             "solveOsiAsIp", "DecompAlgoModel");
   }

#endif
   UTIL_DELARR(solution);
}
