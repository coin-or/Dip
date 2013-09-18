//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, Lehigh University                                 //
//                                                                           //
// Copyright (C) 2002-2013, Lehigh University, Matthew Galati, and Ted Ralphs//
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

void DecompAuto(DecompApp* milp,
                UtilParameters& utilParam,
                UtilTimer& timer,
                DecompMainParam& decompMainParam);

//===========================================================================//
int main(int argc, char** argv)
{
   try {
      //---
      //--- create the utility class for parsing parameters
      //---
      UtilParameters utilParam(argc, argv);
      DecompMainParam decompMainParam;
      UtilTimer timer;
      decompMainParam.timeSetupReal = 0.0;
      decompMainParam.timeSetupCpu  = 0.0;
      decompMainParam.timeSolveReal = 0.0;
      decompMainParam.timeSolveCpu  = 0.0;
      //---
      //--- start overall timer
      //---
      timer.start();
      DecompApp milp;
      //---
      //--- put the one of the functions in the constructor into the main
      //---
      std::vector<int> blockNumCandidates;
      milp.startupLog();
      // get the current working Directory.
      char the_path[256];
      std::string path(getcwd(the_path, 255));
      milp.m_param.CurrentWorkingDir = path;

      if (milp.m_param.LogDebugLevel >= 1) {
         std::cout << path << std::endl;
         std::cout << milp.m_param.CurrentWorkingDir << std::endl;
      }

      const CoinPackedMatrix* m_matrix = milp.readProblem(utilParam);
      blockNumberFinder(milp.m_param, blockNumCandidates, m_matrix);
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
      //      std::vector<DecompApp> milpArray(static_cast<int>(numThreads + 1), milp);
      DecompApp** milpArray = new DecompApp*[numThreads + 1];

      for (int i = 0; i < (numThreads + 1); i++) {
         milpArray[i] = new DecompApp();
      }

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
#ifdef _openmp
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
               milpArray[i]->m_param.NumBlocks = blockNumCandidates[i - 1];
            }

            milpArray[i]->m_param.ThreadIndex = i;
            DecompAuto(milpArray[i], utilParamArray[i],
                       timerArray[i], decompMainParamArray[i]);
         }
      } else {
         decompMainParam.doCut        = utilParam.GetSetting("doCut",        false);
         decompMainParam.doPriceCut   = utilParam.GetSetting("doPriceCut",   true);
         decompMainParam.doDirect     = utilParam.GetSetting("doDirect",     false);
         DecompAuto(&milp, utilParam, timer, decompMainParam);
      }

      for (int i = 0; i < (numThreads + 1); i++) {
         delete milpArray[i];
      }

      delete [] milpArray;

      if (milp.m_param.Concurrent == true) {
         printf("===== FINISH Concurrent Computations Process. =====\n");
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
 *   Determining the candidate block numbers based on the instance fequency
 *     table
 */
void blockNumberFinder(DecompParam utilParam,
                       std::vector<int>& blockNums,
                       const CoinPackedMatrix* matrix)
{
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

   if (utilParam.LogDebugLevel >= 1) {
      std::ofstream histogramTable;
      std::string path1 = utilParam.CurrentWorkingDir + UtilDirSlash() + "histogramTable.dat";
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

   if (utilParam.LogDebugLevel >= 1) {
      std::ofstream histogramTable1;
      std::string path2 = utilParam.CurrentWorkingDir + UtilDirSlash() + "histogramTable1.dat";
      histogramTable1.open(path2.c_str());
      std::map<int, int >::iterator histIter3;

      for (histIter3 = histogram2.begin(); histIter3 != histogram2.end();
            ++histIter3) {
         histogramTable1 << histIter3->second << " " << histIter3->first << "\n";
      }

      histogramTable1.close();
   }

   int blockCands = std::min(utilParam.NumBlocksCand - static_cast<int>(blocksNumTemp.size()),
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
      int counter = utilParam.NumBlocksCand;
      std::set<int>:: iterator setIter = blocksNumTemp.begin();

      while (counter) {
         blockNums.push_back(*setIter);
         setIter++ ;
         --counter;
      }
   }
}

void DecompAuto(DecompApp* milp,
                UtilParameters& utilParam,
                UtilTimer& timer,
                DecompMainParam& decompMainParam)
{
   //---
   //--- put the one of the functions in the constructor into the main
   //---
   milp->initializeApp(utilParam);
   //      milp.startupLog();
   //---
   //--- create the algorithm (a DecompAlgo)
   //---
   DecompAlgo* algo = NULL;

   if ((decompMainParam.doCut + decompMainParam.doPriceCut) != 1)
      throw UtilException("doCut or doPriceCut must be set",
                          "main", "main");

   //assert(doCut + doPriceCut == 1);

   //---
   //--- create the CPM algorithm object
   //---
   if (decompMainParam.doCut) {
      algo = new DecompAlgoC(milp, &utilParam);
   }

   //---
   //--- create the PC algorithm object
   //---
   if (decompMainParam.doPriceCut) {
      algo = new DecompAlgoPC(milp, &utilParam);
   }

   if (decompMainParam.doCut && decompMainParam.doDirect) {
      timer.stop();
      decompMainParam.timeSetupCpu  = timer.getCpuTime();
      decompMainParam.timeSetupReal = timer.getRealTime();
      //---
      //--- solve
      //---
      timer.start();
      algo->solveDirect();
      timer.stop();
      decompMainParam.timeSolveCpu  = timer.getCpuTime();
      decompMainParam.timeSolveReal = timer.getRealTime();
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

      if (milp->m_param.Concurrent == 1) {
         std::cout << "====== The thread number is " << milp->m_param.ThreadIndex
                   << "====" << std::endl;
         std::cout << "====== The block number is  " << milp->m_param.NumBlocks
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
      cout << "\n ============== DECOMP Solution Info [Begin]: ============= \n";
      cout << " Status        = ";

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

      cout << " BestLB        = " << setw(10)
           << UtilDblToStr(alpsModel.getGlobalLB(), 5) << endl
           << " BestUB        = " << setw(10)
           << UtilDblToStr(alpsModel.getGlobalUB(), 5) << endl
           << " Nodes         = "
           << alpsModel.getNumNodesProcessed() << endl
           << " SetupCPU      = " << decompMainParam.timeSetupCpu << endl
           << " SolveCPU      = " << decompMainParam.timeSolveCpu << endl
           << " TotalCPU      = " << decompMainParam.timeSetupCpu + decompMainParam.timeSolveCpu << endl
           << " SetupWallclock= " << decompMainParam.timeSetupReal << endl
           << " SolveWallclock= " << decompMainParam.timeSolveReal << endl
           << " TotalWallclock= " << decompMainParam.timeSetupReal + decompMainParam.timeSolveReal    ;
      cout << "\n ============== DECOMP Solution Info [END  ]: ============= \n";
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
      double userLB   = milp->getBestKnownLB();
      double userUB   = milp->getBestKnownUB();
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
      if (alpsModel.getSolStatus() == AlpsExitStatusOptimal) {
         string::size_type idx = milp->getInstanceName().rfind('/');
         string intanceNameWoDir;

         if (idx != string::npos) {
            intanceNameWoDir = milp->getInstanceName().substr(idx + 1);
         } else {
            intanceNameWoDir = milp->getInstanceName();
         }

         string solutionFile = milp->m_param.CurrentWorkingDir + UtilDirSlash()
                               + intanceNameWoDir + ".sol";
         ofstream osSolution(solutionFile.c_str());
         const DecompSolution* solution = alpsModel.getBestSolution();
         const vector<string>& colNames = alpsModel.getColNames();
         std::cout << "Optimal Solution is " << std::endl;
         solution->print(colNames, 8);
         cout << " Optimal Solution can be found in the file "
              << solutionFile  << endl;
         osSolution.close();
      }

      //---
      //--- free local memory
      //---
      delete algo;
   }
}
