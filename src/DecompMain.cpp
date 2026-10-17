//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Authors: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)      //
//          Ted Ralphs, Lehigh University (ted@lehigh.edu)                   //
//          Jiadong Wang, Lehigh University (jiw408@lehigh.edu)              //
//                                                                           //
// Copyright (C) 2002-2019, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

//===========================================================================//
#include "UtilParameters.h"
//===========================================================================//
#include "DecompApp.h"
//===========================================================================//
#include "AlpsDecompModel.h"
//===========================================================================//
#include "DecompAlgoC.h"
#include "DecompAlgoPC.h"
#include "DecompAlgoRC.h"
//===========================================================================//
#include "UtilTimer.h"
#include <algorithm>
#ifdef _OPENMP
#include "omp.h"
#endif
using namespace std;

void blockNumberFinder(DecompParam utilParam,
                       std::vector<int>& blockNums,
                       const CoinPackedMatrix* matrix);

void DecompAuto(DecompApp milp,
                UtilParameters& utilParam,
                UtilTimer& timer,
                DecompMainParam& decompMainParam);

DecompSolverResult* solveDirect(const DecompApp& decompApp);

//===========================================================================//

int main(int argc, char** argv)
{
   try {

      UtilTimer timer;
      std::vector<int> blockNumCandidates;
      DecompMainParam decompMainParam;
      decompMainParam.timeSetupReal = 0.0;
      decompMainParam.timeSetupCpu  = 0.0;
      decompMainParam.timeSolveReal = 0.0;
      decompMainParam.timeSolveCpu  = 0.0;

      //---
      //--- Parse cpmmand line and store parameters
      //---
      UtilParameters utilParam(argc, argv);

      //---
      //--- start overall timer
      //---
      timer.start();

      //---
      //--- construct the instance
      //---
      DecompApp milp(utilParam);

      // get the current working Directory.
      char the_path[256];
      milp.m_param.CurrentWorkingDir = std::string(getcwd(the_path, 255));
      if (milp.m_param.LogDebugLevel >= 1) {
         std::cout << milp.m_param.CurrentWorkingDir << std::endl;
      }

      //---
      //--- Analyze the matrix
      //---
      const CoinPackedMatrix* m_matrix = milp.getMatrix();
      if (milp.m_param.BlockNumInput > 0) {
         milp.NumBlocks = milp.m_param.BlockNumInput;
         milp.m_param.Concurrent = false ;
         milp.m_param.NumBlocksCand = 0;
      }
      if (milp.m_param.NumBlocksCand == 0) {
	 blockNumberFinder(milp.m_param, blockNumCandidates, m_matrix);
      }
      // obtain the number of CPU (core)s on machines with operating
      // system Linux, Solaris, & AIX and Mac OS X
      // (for all OS releases >= 10.4, i.e., Tiger onwards)
      // other systems has different syntax to obtain the core number
#ifdef _OPENMP
      int numCPU = omp_get_num_procs();
#else
      int numCPU = 1;
#endif

      if (milp.m_param.LogDebugLevel > 1) {
         std::cout << "The number of cores is "
                   << numCPU << std::endl;
      }

      // the actual thread number is the minimum of
      // number of cores, total block numbers and the thread number
      // used in concurrent computations
      int numThreads = min(min(numCPU,
                               static_cast<int>(blockNumCandidates.size())),
                           milp.m_param.ConcurrentThreadsNum);
      std::vector<DecompApp> milpArray(static_cast<int>(numThreads + 1), milp);
      std::vector<DecompMainParam> decompMainParamArray(static_cast<int>
            (numThreads + 1),
            decompMainParam);
      std::vector<UtilTimer> timerArray(static_cast<int>(numThreads + 1),
                                        timer);
      std::vector<UtilParameters> utilParamArray(static_cast<int>
            (numThreads + 1),
            utilParam);

      if (milp.m_param.Concurrent == true ) {
         printf("===== START Concurrent Computations Process. =====\n");
#ifdef _OPENMP
#pragma omp parallel for
#endif

         for (int i = 0 ; i < (numThreads + 1); i++) {
            if (i == 0) {
               decompMainParamArray[i].doCut = true;
               decompMainParamArray[i].doPriceCut = false;
               decompMainParamArray[i].doDirect = true;
            } else {
               decompMainParamArray[i].doCut = false;
               decompMainParamArray[i].doPriceCut = true;
               decompMainParamArray[i].doDirect = false;
               milpArray[i].NumBlocks = blockNumCandidates[i - 1];
            }

            milpArray[i].m_threadIndex = i;
            DecompAuto(milpArray[i], utilParamArray[i],
                       timerArray[i], decompMainParamArray[i]);
         }
      } else {
         decompMainParam.doCut        = utilParam.GetSetting("doCut",        false);
         decompMainParam.doPriceCut   = utilParam.GetSetting("doPriceCut",   true);
         decompMainParam.doDirect     = utilParam.GetSetting("doDirect",     false);
         DecompAuto(milp, utilParam, timer, decompMainParam);
      }

      if (milp.m_param.Concurrent == true) {
         printf("===== FINISH Concurrent Computations Process. =====\n");
         printf("======== SUMMARY OF CONCURRENT COMPUTATIONS =======\n");
         cout << "Method" << setw(20) << "BlockNumber" << setw(20)
              << "WallClockTime" << setw(20) << "CPUTime" << setw(20)
              << "BestLB" << setw(25) << "BestUB" << endl;

         for (int i = 0 ; i < (numThreads + 1); i++) {
            if (i == 0) {
               cout << "B&C ";
            } else {
               cout << "B&P";
            }

            cout << setw(15);

            if (i == 0) {
               cout << "NA";
            } else {
               cout << milpArray[i].NumBlocks;
            }

            cout << setw(25) << setprecision(7)
                 << decompMainParamArray[i].timeSetupReal +
                 decompMainParamArray[i].timeSolveReal
                 << setw(23) << setprecision(7)
                 << decompMainParamArray[i].timeSetupCpu +
                 decompMainParamArray[i].timeSolveCpu
                 << setw(23) << setprecision(7)
                 << decompMainParamArray[i].bestLB
                 << setw(25) << setprecision(7)
                 << decompMainParamArray[i].bestUB
                 << endl;
         }
      }
   } catch (CoinError& ex) {
      cerr << "COIN Exception [ " << ex.message() << " ]"
           << " at " << ex.fileName()  << ":L" << ex.lineNumber()
           << " in " << ex.className() << "::" << ex.methodName() << endl;
      return 1;
   }

   return 0;
}

/*
 *   Determining the candidate block numbers based on the instance frequency
 *     table
 */
void blockNumberFinder(DecompParam decompParam,
                       std::vector<int>& blockNums,
                       const CoinPackedMatrix* matrix)
{
   if (decompParam.NumBlocksCand == 0) {
      return ;
   }

   const int* lengthRows = matrix->getVectorLengths();

   int numRows = matrix->getNumRows();

   // The following code creates a histogram table to store the
   // nonzero counts and number of rows
   std::map<int, int> histogram;

   for (int i = 0 ; i < numRows; ++i) {
      if (histogram.count(lengthRows[i]) > 0) {
         histogram[lengthRows[i]] += 1;
      } else {
         histogram.insert(std::pair<int, int>(lengthRows[i], 1));
      }
   }

   std::map<int, int>::iterator histIter;

   if (decompParam.LogDebugLevel >= 1) {
      std::ofstream histogramTable;
      std::string path1 = decompParam.CurrentWorkingDir + UtilDirSlash() +
                          "histogramTable.dat";
      histogramTable.open(path1.c_str());

      for (histIter = histogram.begin(); histIter != histogram.end();
            histIter++) {
         histogramTable << histIter->first << " " << histIter->second << "\n";
      }

      histogramTable.close();
   }

   // Aggregation steps, aggregate the entries in the histogram
   //     Number of nonzeros      Number of rows
   //     4                          8
   //     9                          8
   // After aggregation:
   //     9             8
   // Then put the number of rows into the candidates queue
   std::map<int, int>::iterator histIter2;
   std::map<int, int> histogram2;
   std::set<int> blocksNumTemp;

   for (histIter = histogram.begin(); histIter != histogram.end();
         ++histIter) {
      int keyvalue_pre = histIter->second ;
      int max = 0;
      std::map<int, int>::iterator histIterTemp = histIter;
      ++histIterTemp;

      for (histIter2 = histIterTemp; histIter2 != histogram.end();
            ++histIter2) {
         int keyvalue_curr = histIter2->second;
         // std::cout << " The current value of key is " << keyvalue_curr
         // <<std::endl;
         int max_inter = 0;

         if (keyvalue_pre == keyvalue_curr) {
            max_inter = (histIter->first > histIter2->first)
                        ?   histIter->first : histIter2->first;
         }

         max = (max_inter > max) ? max_inter : max;

         if (max && max != 1) {
            blocksNumTemp.insert(max);
            blockNums.push_back(max);
            histogram2.insert(std::pair<int, int>(histIter->second, max));
         }
      }

      if (histogram2.find(histIter->second) == histogram2.end() ) {
         histogram2.insert(std::pair<int, int>(histIter->second,
                                               histIter->first));
      }
   }

   if (decompParam.LogDebugLevel >= 1) {
      std::ofstream histogramTable1;
      std::string path2 = decompParam.CurrentWorkingDir + UtilDirSlash() +
                          "histogramTable1.dat";
      histogramTable1.open(path2.c_str());
      std::map<int, int >::iterator histIter3;

      for (histIter3 = histogram2.begin(); histIter3 != histogram2.end();
            ++histIter3) {
         histogramTable1 << histIter3->second << " " << histIter3->first << "\n";
      }

      histogramTable1.close();
   }

   int blockCands = std::min(decompParam.NumBlocksCand - static_cast<int>
                             (blocksNumTemp.size()),
                             static_cast<int>(histogram2.size()));

   if (blockCands > 0) {
      std::map<int, int >::iterator histIterLower  = histogram2.end();

      while (blockCands ) {
         --histIterLower;
         blockCands--;

         if (blocksNumTemp.find(histIterLower->second) == blocksNumTemp.end()) {
            blockNums.push_back(histIterLower->second);
         }
      }
   } else {
      int counter = decompParam.NumBlocksCand;
      std::set<int>:: iterator setIter = blocksNumTemp.begin();

      while (counter) {
         blockNums.push_back(*setIter);
         setIter++ ;
         --counter;
      }
   }
}

void DecompAuto(DecompApp milp,
                UtilParameters& utilParam,
                UtilTimer& timer,
                DecompMainParam& decompMainParam)
{
   if (milp.NumBlocks == 0 && milp.m_param.Concurrent) {
      milp.NumBlocks = 3;
   }

   if (milp.m_threadIndex == 0 && milp.m_param.Concurrent) {
      decompMainParam.timeSetupCpu  = timer.getCpuTime();
      decompMainParam.timeSetupReal = timer.getRealTime();
      //---
      //--- solve
      //---
      timer.start();
      DecompSolverResult* result = solveDirect(milp);
      timer.stop();
      decompMainParam.bestLB = result->m_objLB;
      decompMainParam.bestUB = result->m_objUB;
      decompMainParam.timeSolveCpu  = timer.getCpuTime();
      decompMainParam.timeSolveReal = timer.getRealTime();
      UTIL_DELPTR(result);
      return ;
   }

   //---
   //--- Initialize
   //---
   milp.initializeApp();

   //---
   //--- create the algorithm (a DecompAlgo)
   //---

   if ((decompMainParam.doCut + decompMainParam.doPriceCut) != 1)
      throw UtilException("doCut or doPriceCut must be set",
                          "main", "main");

   //---
   //--- create the algorithm object
   //---
   DecompAlgo* algo = NULL;
   if (decompMainParam.doCut) {
      algo = new DecompAlgoC(&milp, utilParam);
   }else{
      algo = new DecompAlgoPC(&milp, utilParam);
   }

   if (decompMainParam.doCut && decompMainParam.doDirect) {
      timer.stop();
      decompMainParam.timeSetupCpu  = timer.getCpuTime();
      decompMainParam.timeSetupReal = timer.getRealTime();
      //---
      //--- solve
      //---
      timer.start();
      DecompSolverResult* result = algo->solveDirect();
      timer.stop();
      decompMainParam.bestLB = result->m_objLB;
      decompMainParam.bestUB = result->m_objUB;
      decompMainParam.timeSolveCpu  = timer.getCpuTime();
      decompMainParam.timeSolveReal = timer.getRealTime();
      UTIL_DELPTR(result);
   } else {
      //---
      //--- create the driver AlpsDecomp model
      //---
      AlpsDecompModel alpsModel(utilParam, algo);
      timer.stop();
      decompMainParam.timeSetupCpu  = timer.getCpuTime();
      decompMainParam.timeSetupReal = timer.getRealTime();
      //---
      //--- solve
      //---
      timer.start();
      alpsModel.solve();
      timer.stop();

      if (milp.m_param.Concurrent == 1) {
         std::cout << "====== The thread number is " << milp.m_threadIndex
                   << "====" << std::endl;
         std::cout << "====== The block number is  " << milp.NumBlocks
                   << "====" << std::endl;
         std::cout << "====== Branch-and-Cut       " << decompMainParam.doCut
                   << "====" << std::endl;
         std::cout << "====== Branch-and-Price     " << decompMainParam.doPriceCut
                   << "====" << std::endl;
         std::cout << "                                          " << std::endl;
         std::cout << "                                          " << std::endl;
         std::cout << "                                          " << std::endl;
      }

      decompMainParam.timeSolveCpu  = timer.getCpuTime();
      decompMainParam.timeSolveReal = timer.getRealTime();
      //---
      //--- sanity check
      //---
      cout << setiosflags(ios::fixed | ios::showpoint);
      int statusCheck = alpsModel.getSolStatus();
      cout << "                                             " << endl;
      cout << "============== DECOMP Solution Info [Begin]: ============== \n";
      cout << "Status        = ";

      if ( !statusCheck) {
         cout << "Optimal" << endl;
      } else if (statusCheck == 1) {
         cout << "TimeLimit" << endl;
      } else if (statusCheck == 2) {
         cout << "NodeLimit" << endl;
      } else if (statusCheck == 3) {
         cout << "SolLimit" << endl;
      } else if (statusCheck == 4) {
         cout << "Feasible" << endl;
      } else if (statusCheck == 5) {
         cout << "Infeasible" << endl;
      } else if (statusCheck == 6) {
         cout << "NoMemory" << endl;
      } else if (statusCheck == 7) {
         cout << "Failed" << endl;
      } else if (statusCheck == 8) {
         cout << "Unbounded" << endl;
      } else {
         cout << "Unknown" << endl;
      }

      decompMainParam.bestLB = alpsModel.getGlobalLB();
      decompMainParam.bestUB = alpsModel.getGlobalUB();
      cout << "BestLB        = " << setw(10)
           << UtilDblToStr(alpsModel.getGlobalLB(), 5) << endl
           << "BestUB        = " << setw(10)
           << UtilDblToStr(alpsModel.getGlobalUB(), 5) << endl
           << "OptiGap       = " << setw(10)
           << UtilDblToStr(UtilCalculateGap(alpsModel.getGlobalLB(),
                                            alpsModel.getGlobalUB(),
					    milp.getDecompAlgo()->getInfinity()), 5)
           << endl
           << "Nodes         = "
           << alpsModel.getNumNodesProcessed() << endl
           << "SetupCPU      = " << decompMainParam.timeSetupCpu << endl
           << "SolveCPU      = " << decompMainParam.timeSolveCpu << endl
           << "TotalCPU      = " << decompMainParam.timeSetupCpu +
           decompMainParam.timeSolveCpu << endl
           << "SetupWallclock= " << decompMainParam.timeSetupReal << endl
           << "SolveWallclock= " << decompMainParam.timeSolveReal << endl
           << "TotalWallclock= " << decompMainParam.timeSetupReal
           + decompMainParam.timeSolveReal << endl   ;
      cout << "============== DECOMP Solution Info [END  ]: ============== \n";
      /* TODO: Add a global parameter to control the subproblem
               parallelization
      cout << "The parallel efficiency is "
           << timeSolveCpu/(milp.m_param.NumThreads*timeSolveReal)
           << endl;
      */
      //---
      //--- sanity check
      //---   if user defines bestLB==bestUB (i.e., known optimal)
      //---   and solved claims we have optimal, check that they match
      //---
      double epsilon  = 0.01; //1%
      double userLB   = milp.getBestKnownLB();
      double userUB   = milp.getBestKnownUB();
      double userDiff = fabs(userUB - userLB);

      if (alpsModel.getSolStatus() == AlpsExitStatusOptimal &&
            userDiff                  < epsilon) {
         double diff   = fabs(alpsModel.getGlobalUB() - userUB);
         double diffPer = userUB == 0 ? diff : diff / userUB;

         if (diffPer > epsilon) {
            cerr << setiosflags(ios::fixed | ios::showpoint);
            cerr << "ERROR. BestKnownLB/UB= "
                 << UtilDblToStr(userUB, 5)
                 << " but DIP claims GlobalUB= "
                 << UtilDblToStr(alpsModel.getGlobalUB(), 5)
                 << endl;
            throw UtilException("Invalid claim of optimal.",
                                "main", "MILPBlock");
         }
      }

      //---
      //--- get optimal solution
      //---

      if (milp.m_param.SolutionOutputToFile
            && alpsModel.getGlobalUB() < 1.e100) {
         const DecompSolution* solution = alpsModel.getBestSolution();
         const vector<string>& colNames = alpsModel.getColNames();
         string solutionFile;

         if (milp.m_param.SolutionOutputFileName == "") {
            string::size_type idx = milp.getInstanceName().rfind('/');
            string intanceNameWoDir;

            if (idx != string::npos) {
               intanceNameWoDir = milp.getInstanceName().substr(idx + 1);
            } else {
               intanceNameWoDir = milp.getInstanceName();
            }

            solutionFile = milp.m_param.CurrentWorkingDir + UtilDirSlash()
                           + intanceNameWoDir + ".sol";
         } else {
            solutionFile = milp.m_param.SolutionOutputFileName;
         }

         ofstream osSolution(solutionFile.c_str());
         osSolution.precision(16);
         const double* sol = solution->getValues();
         osSolution << "=obj=" << setw(10);
         osSolution.precision(8);
         osSolution << " " << alpsModel.getGlobalUB()
                    << std::endl;

         for (int i = 0; i < solution->getSize(); i++) {
            if (!UtilIsZero(sol[i])) {
               osSolution << colNames[i] << setw(10);
               osSolution.precision(8);
               osSolution << " " << sol[i] << std::endl;
            } else {
               osSolution << colNames[i] << setw(10);
               osSolution.precision(8);
               osSolution << " " << 0.0000000 << std::endl;
            }
         }

         osSolution.close();

         if (alpsModel.getSolStatus() == AlpsExitStatusOptimal) {
            std::cout << "Optimal Solution is " << std::endl;
            solution->print(colNames, 8);
            cout << " Optimal Solution can be found in the file "
                 << solutionFile  << endl;
         }
      }

      //---
      //--- free local memory
      //---
      delete algo;
   }
}

DecompSolverResult* solveDirect(const DecompApp& decompApp)
{
   //---
   //--- Solve the original IP with a generic IP solver
   //--- without going through the decomposition phase
   //--- this function is created such that the DIP can serves
   //--- as an interface to call standalone branch-and-cut solver
   //---

   OsiSolverInterface *m_problemSI;

   if (decompApp.m_param.DecompIPSolver == "SYMPHONY"){
#ifdef DIP_HAS_SYMPHONY
      m_problemSI = new OsiSymSolverInterface();
#else
      throw UtilException("SYMPHONY selected as solver, but it's not available",
			  "getOsiIpSolverInterface", "DecompAlgo");
#endif
   }else if (decompApp.m_param.DecompIPSolver == "Cbc"){
#ifdef DIP_HAS_CBC
      m_problemSI = new OsiCbcSolverInterface();
#else
      throw UtilException("Cbc selected as solver, but it's not available",
			  "getOsiIpSolverInterface", "DecompAlgo");
#endif
   }else if (decompApp.m_param.DecompIPSolver == "CPLEX"){
#ifdef DIP_HAS_CPX
      m_problemSI = new OsiCpxSolverInterface();
#else
      throw UtilException("CPLEX selected as solver, but it's not available",
			  "getOsiIpSolverInterface", "DecompAlgo");
#endif
   }else if (decompApp.m_param.DecompIPSolver == "Gurobi"){
#ifdef DIP_HAS_GRB
      m_problemSI = new OsiGrbSolverInterface();
#else
      throw UtilException("Gurobi selected as solver, but it's not available",
			  "getOsiIpSolverInterface", "DecompAlgo");
#endif
   }else if (decompApp.m_param.DecompIPSolver == "Xpress"){
#ifdef DIP_HAS_XPR
      m_problemSI = new OsiXprSolverInterface();
#else
      throw UtilException("Xpress selected as solver, but it's not available",
			  "getOsiIpSolverInterface", "DecompAlgo");
#endif
   }

   string fileName;

   if (decompApp.m_param.DataDir != "") {
      fileName = decompApp.m_param.DataDir + UtilDirSlash() +
                 decompApp.m_param.Instance;
   } else {
      fileName = decompApp.m_param.Instance;
   }

   std::cout << "The file name is " << fileName << std::endl;

   if (decompApp.m_param.Instance.empty()) {
      cerr << "================================================" << std::endl
           << "Usage:"
           << "./dip  --MILP:BlockFileFormat List" << std::endl
           << "       --MILP:Instance /FilePath/ABC.mps" << std::endl
           << "       --MILP:BlockFile /FilePath/ABC.block" << std::endl
           << "================================================" << std::endl
           << std::endl;
      exit(0);
   }

   m_problemSI->readMps(fileName.c_str());
   int numCols    = decompApp.m_mpsIO.getNumCols();
   int nNodes;
   double objLB   = -m_problemSI->getInfinity();
   double objUB   = m_problemSI->getInfinity();
   double timeLimit = decompApp.m_param.TimeLimit;
   UtilTimer timer;
   timer.start();
   DecompSolverResult* result = new DecompSolverResult(m_problemSI->getInfinity());
   if (decompApp.m_param.DecompIPSolver == "Cbc"){
#ifdef DIP_HAS_CBC
      CbcModel cbc(*m_problemSI);
      int logIpLevel = decompApp.m_param.LogIpLevel;
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
			  "solveDirect", "DecompMain");
#endif
   }else if (decompApp.m_param.DecompIPSolver == "CPLEX"){
#ifdef DIP_HAS_CPX
      OsiCpxSolverInterface* masterSICpx
	 = dynamic_cast<OsiCpxSolverInterface*>(m_problemSI);
      CPXLPptr  cpxLp  = masterSICpx->getLpPtr();
      CPXENVptr cpxEnv = masterSICpx->getEnvironmentPtr();
      int       status = 0;
      masterSICpx->switchToMIP();//need?
      //---
      //--- set the time limit
      //---
      status = CPXsetdblparam(cpxEnv, CPX_PARAM_TILIM, timeLimit);
      //---
      //--- set the thread limit, otherwise CPLEX will use all the resources
      //---
      status = CPXsetintparam(cpxEnv, CPX_PARAM_THREADS,
			      decompApp.m_param.NumThreadsIPSolver);
      
      if (status)
	 throw UtilException("CPXsetdblparam failure",
			     "solveDirect", "DecompAlgoC");
      
      //---
      //--- solve the MILP
      //---
      UtilTimer timer1;
      timer1.start();
      masterSICpx->branchAndBound();
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
	    const double* solDbl = masterSICpx->getColSolution();
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
			  "solveDirect", "DecompMain");
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
   cout << "DIRECT SOLVE"
        << " Real=" << setw(10) << UtilDblToStr(timer.getRealTime(), 5)
        << " Cpu= " << setw(10) << UtilDblToStr(timer.getCpuTime() , 5)
        << " Nodes= " << setw(8) << nNodes
        << " objLB= " << setw(10) << UtilDblToStr(objLB, 3)
        << " objUB= " << setw(10) << UtilDblToStr(objUB, 3)
        << endl;
   UTIL_DELPTR( m_problemSI);
   return result;
}
