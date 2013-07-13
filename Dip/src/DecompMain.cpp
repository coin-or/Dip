//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, Lehigh University                                 //
//                                                                           //
// Copyright (C) 2002-2011, Lehigh University, Matthew Galati, and Ted Ralphs//
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
#include "omp.h"

using namespace std; 

void preprocessApp(UtilParameters& utilParam,  
		   std::vector<int> &blockNums,  
		   DecompApp& milp);  

void* DecompAuto(DecompApp & milp, 
		  UtilParameters & utilParam, 
		  UtilTimer & timer, 
		  DecompMainParam & decompMainParam 
		  ); 

//===========================================================================//
int main(int argc, char ** argv){
   try{
      
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
      
      std::vector<int> blockNumCandidates;  
      //---
      //--- put the one of the functions in the constructor into the main
      //---
      milp.initializeApp(utilParam); 

      std::vector<int> blockNumCandidates;  

      milp.startupLog();
      preprocessApp(utilParam, blockNumCandidates, milp); 

      if(milp.m_param.Concurrent == 1 ) 
	{ 
	  
	  // obtain the number of CPU (core)s on machines with operating 
	  // system Linux, Solaris, & AIX and Mac OS X  
	  // (for all OS releases >= 10.4, i.e., Tiger onwards)  
	  // other systems has different syntax to obtain the core number 
	  
	  int numCPU = sysconf( _SC_NPROCESSORS_ONLN ); 
	  
	  std::cout << "The number of cores in this node is "  
		    << numCPU << std::endl; 
	  
	  int numThreads = min(numCPU, static_cast<int>(blockNumCandidates.size()));  
	  	  	  
	  printf("===== START Concurrent Computations Process. =====\n"); 

#pragma omp parallel for
	  for(int i = 0 ; i < (numThreads + 1); i++){ 
	    
	    std::vector<DecompApp> milpArray(static_cast<int>(numThreads + 1),milp); 
	    	    
	    std::vector<DecompMainParam> decompMainParamArray(static_cast<int>(numThreads + 1), 
							      decompMainParam); 	    	     

	    std::vector<UtilTimer> timerArray(static_cast<int>(numThreads + 1), 
					      timer); 

	    std::vector<UtilParameters> utilParamArray(static_cast<int>(numThreads + 1), 
						       utilParam); 

	    if (i == 0){ 
	      decompMainParamArray[i].doCut = true;  
	      decompMainParamArray[i].doPriceCut = false;  
	      decompMainParamArray[i].doDirect = true;         
	    } 
	    else{                
	      decompMainParamArray[i].doCut = false;  
	      decompMainParamArray[i].doPriceCut = true;  
	      decompMainParamArray[i].doDirect = false;                
	      milpArray[i].m_param.NumBlocks = blockNumCandidates[i-1];  
	    } 	    

	    milpArray[i].m_param.ThreadIndex = i; 

	    DecompAuto(milpArray[i], utilParamArray[i], 
		       timerArray[i], decompMainParamArray[i]);	  	    	    	  	  	    	    
	  } 
	} 
      else{ 		  
	  decompMainParam.doCut        = utilParam.GetSetting("doCut",        false); 	
	  decompMainParam.doPriceCut   = utilParam.GetSetting("doPriceCut",   true); 	
	  decompMainParam.doDirect     = utilParam.GetSetting("doDirect",     false); 
	  DecompAuto(milp, utilParam, timer, decompMainParam); 
      }
      printf("===== FINISH Concurrent Computations Process. =====\n");
   }
   catch(CoinError & ex){
      cerr << "COIN Exception [ " << ex.message() << " ]"
           << " at " << ex.fileName()  << ":L" << ex.lineNumber()
           << " in " << ex.className() << "::" << ex.methodName() << endl;
      return 1;
   }
   return 0;
} 

void* DecompSolve(DecompApp & milp,
			  UtilParameters & utilParam,
			  UtilTimer & timer,
			  DecompMainParam & decompMainParam){
  

       //---
       //--- put the one of the functions in the constructor into the main
      //---

      milp.initializeApp(utilParam); 


      //      milp.startupLog(); 

      //---
      //--- create the algorithm (a DecompAlgo)
      //---
      DecompAlgo * algo = NULL;
      if((decompMainParam.doCut + decompMainParam.doPriceCut) != 1)
         throw UtilException("doCut or doPriceCut must be set", 
                             "main", "main");
      //assert(doCut + doPriceCut == 1);

      //---
      //--- create the CPM algorithm object
      //---      
      if(decompMainParam.doCut)	 
         algo = new DecompAlgoC(&milp, &utilParam);
      
      //---
      //--- create the PC algorithm object
      //---
      if(decompMainParam.doPriceCut)
         algo = new DecompAlgoPC(&milp, &utilParam);
      
      if(decompMainParam.doCut && decompMainParam.doDirect){
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
      }
      else{
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


	 std::cout << "==========================================" << std::endl;
	 std::cout << "======   The thread number is  ============= " << std::endl;
	 std::cout << "============" << milp.m_param.ThreadIndex <<"============= " << std::endl;
	 std::cout << "======   The block number is ============ " << std::endl;	
	 std::cout << "===========" << milp.m_param.NumBlocks <<"============" << std::endl;
	 std::cout << "=========== Branch-and-Cut  " << decompMainParam.doCut
		   <<"============" << std::endl;
	 std::cout << "=========== Branch-and-Price  " << decompMainParam.doPriceCut 
		   <<"============" << std::endl;
	 
	 std::cout << "                                          " << std::endl;
	 std::cout << "                                          " << std::endl;
	 std::cout << "                                          " << std::endl;
	 std::cout << "                                          " << std::endl;


	 decompMainParam.timeSolveCpu  = timer.getCpuTime();
	 decompMainParam.timeSolveReal = timer.getRealTime();

	 //---
	 //--- sanity check
	 //---
	 cout << setiosflags(ios::fixed|ios::showpoint);
	 int statusCheck = alpsModel.getSolStatus(); 
	 cout << "                                             "<< endl;
	 cout << "\n ============== DECOMP Solution Info [Begin]: ============= \n";
	 cout << " Status        = ";
	 if ( !statusCheck)
	   cout << "Optimal" << endl;
	 else if (statusCheck == 1)
	   cout << "TimeLimit" << endl;
	 else if (statusCheck == 2)
	   cout << "NodeLimit" << endl;
	 else if (statusCheck == 3)
	   cout << "SolLimit" << endl;
	 else if (statusCheck == 4)
	   cout << "Feasible" << endl;
	 else if (statusCheck == 5)
	   cout << "Infeasible" << endl;
	 else if (statusCheck == 6)
	   cout << "NoMemory" << endl;
	 else if (statusCheck == 7)
	   cout << "Failed" << endl;
	 else if (statusCheck == 8)
	   cout << "Unbounded" << endl;
	 else 
	   cout << "Unknown" << endl;
	 
	 cout << " BestLB        = " << setw(10) 
	      << UtilDblToStr(alpsModel.getGlobalLB(),5) << endl
	      << " BestUB        = " << setw(10)
	      << UtilDblToStr(alpsModel.getGlobalUB(),5) << endl      
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
         double userLB   = milp.getBestKnownLB();
         double userUB   = milp.getBestKnownUB();
         double userDiff = fabs(userUB - userLB);
         if(alpsModel.getSolStatus() == AlpsExitStatusOptimal &&
            userDiff                  < epsilon){
            double diff   = fabs(alpsModel.getGlobalUB() - userUB);
	    double diffPer= userUB == 0 ? diff : diff / userUB;
            if(diffPer > epsilon){
	       cerr << setiosflags(ios::fixed|ios::showpoint);
               cerr << "ERROR. BestKnownLB/UB= " 
                    << UtilDblToStr(userUB,5) 
		    << " but DIP claims GlobalUB= " 
                    << UtilDblToStr(alpsModel.getGlobalUB(),5) 
		    << endl;
               throw UtilException("Invalid claim of optimal.",
                                   "main", "MILPBlock");
            }
         }
	 
         //---
         //--- get optimal solution
         //---     
         if(alpsModel.getSolStatus() == AlpsExitStatusOptimal){
	    string   solutionFile = milp.getInstanceName() + ".sol";
	    ofstream osSolution(solutionFile.c_str());
            const DecompSolution * solution = alpsModel.getBestSolution();
	    const vector<string> & colNames = alpsModel.getColNames();
            cout << " Optimal Solution can be found in the file "
		 << solutionFile  << endl;
            solution->print(colNames, 8, osSolution);
	    osSolution.close();
         }

	 //---
	 //--- free local memory
	 //---
	 delete algo;      
      }

}

   /**
    *  preprocessApp method takes the instances and performs the following
    *  1. MILP preprocessing ( tighting the variable bounds, remove redundant
    *                          constraints etc, on the TODO list)
    *
    *  2. determining the candidate block numbers based on the instance fequency
    *     table
    *
    */

void preprocessApp(UtilParameters& utilParam, 
		   std::vector<int> &blockNums,
		   DecompApp& milp)		 
{
   //---
   //--- get application parameters
   //---
  
   milp.m_param.getSettings(utilParam);

   if (milp.m_param.LogLevel >= 1) {
      milp.m_param.dumpSettings();
   }

   //---
   //--- read MILP instance (mps format)
   //---
   string fileName = milp.m_param.DataDir
                     + UtilDirSlash() + milp.m_param.Instance;

   if (milp.m_param.Instance.empty()) {
      cerr << "===========================================================" << std::endl
           << "Users need to provide correct "
           << "instance path and name" << std::endl
           << "                                                     " << std::endl
           << "Example: ./dip  --MILP:BlockFileFormat List" << std::endl
           << "                --MILP:Instance /FilePath/ABC.mps" << std::endl
           << "                --MILP:BlockFile /FilePath/ABC.block" << std::endl
           << "===========================================================" << std::endl
           << std::endl;
      throw UtilException("I/O Error.", "initializeApp", "DecompApp");
   }

   

   int rstatus = 0;
   bool foundFormat = false;

   if (milp.m_param.InstanceFormat == "") {
      string::size_type idx = fileName.rfind('.');

      if (idx != string::npos) {
         string extension = fileName.substr(idx + 1);

         if (extension == "MPS" || extension == "mps") {
	   milp.m_param.InstanceFormat = "MPS";
         } else if (extension == "LP" || extension == "lp") {
	   milp.m_param.InstanceFormat = "LP";
         }
      } else {
         cerr << "File format not specified and no file extension" << endl;
         throw UtilException("I/O Error.", "initializeApp", "DecompApp");
      }
   }

   if(milp.m_param.InstanceFormat =="MPS"){
     milp.m_mpsIO.messageHandler()->setLogLevel(milp.m_param.LogLpLevel);
   } else if (milp.m_param.InstanceFormat == "LP"){
     milp.m_lpIO.messageHandler()->setLogLevel(milp.m_param.LogLpLevel); 
   }

   if (milp.m_param.InstanceFormat == "MPS") {
      rstatus = milp.m_mpsIO.readMps(fileName.c_str());
      foundFormat = true;
   } else if (milp.m_param.InstanceFormat == "LP") {
      milp.m_lpIO.readLp(fileName.c_str());
      foundFormat = true;
   }

   if (!foundFormat) {
      cerr << "Error: Format = " << milp.m_param.InstanceFormat << " unknown."
           << endl;
      throw UtilException("I/O Error.", "initalizeApp", "DecompApp");
   }

   if (rstatus < 0) {
      cerr << "Error: Filename = " << fileName << " failed to open." << endl;
      throw UtilException("I/O Error.", "initalizeApp", "DecompApp");
   }

   /*
   if (milp.m_param.LogLevel >= 2)
     if(milp.m_param.InstanceFormat == "MPS"){
       (*m_osLog) << "Objective Offset = "      
		  << UtilDblToStr(milp.m_mpsIO.objectiveOffset()) << endl;
     } else if (milp.m_param.InstanceFormat == "LP"){
       (*m_osLog) << "Objective Offset = "      
		  << UtilDblToStr(milp.m_lpIO.objectiveOffset()) << endl;       
     }
   */
   //---
   //--- set best known lb/ub
   //---
   double offset = 0;

   if (milp.m_param.InstanceFormat == "MPS") {
      offset = milp.m_mpsIO.objectiveOffset();
   } else if (milp.m_param.InstanceFormat == "LP") {
      offset = milp.m_lpIO.objectiveOffset();
   }

   milp.setBestKnownLB(milp.m_param.BestKnownLB + offset);
   milp.setBestKnownUB(milp.m_param.BestKnownUB + offset);
   milp.preprocess();

   if (milp.m_param.Concurrent == 1) {

     if (milp.m_param.InstanceFormat == "MPS"){
       milp.m_matrix = milp.m_mpsIO.getMatrixByRow();
     } else if (milp.m_param.InstanceFormat == "LP"){
       milp.m_matrix = milp.m_lpIO.getMatrixByRow();
     }
     
      const int* lengthRows = milp.m_matrix->getVectorLengths();

      int numRows = 0;
      if (milp.m_param.InstanceFormat == "MPS"){ 
	numRows = milp.m_mpsIO.getNumRows();
      } else if (milp.m_param.InstanceFormat == "LP"){
	numRows = milp.m_lpIO.getNumRows(); 
      }
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

      if (milp.m_param.LogDebugLevel >= 1) {
         std::ofstream histogramTable;
         histogramTable.open("/home/jiw508/Dip-branch/build_test_main/bin/histogramTable.dat");

         for (histIter = histogram.begin(); histIter != histogram.end();
               ++histIter) {
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
      std::map<int, int>::iterator histIter2;
      std::map<int, int> histogram2(histogram);

      for (histIter = histogram.begin(); histIter != histogram.end();
            ++histIter) {
         int keyvalue_pre = histIter->second ;

         for (histIter2 = histogram2.begin(); histIter2 != histogram2.end();
               ++histIter2) {
            int keyvalue_curr = histIter2->second;
            // std::cout << " The current value of key is " << keyvalue_curr
            // <<std::endl;

            if (keyvalue_pre == keyvalue_curr) {
               if (histIter->first > histIter2->first) {
                  histogram2.erase(histIter2);
               }
            }

            // remove the entry with map values equal 1

            if (keyvalue_curr == 1) {
               histogram2.erase(histIter2);
            }
         }
      }

      if (milp.m_param.LogDebugLevel >= 1) {
         std::ofstream histogramTable1;
         histogramTable1.open("/home/jiw508/Dip-branch/build_test_ray/bin/histogramTable1.dat");

         for (histIter = histogram2.begin(); histIter != histogram2.end();
               ++histIter) {
            histogramTable1 << histIter->first << " " << histIter->second << "\n";
         }

         histogramTable1.close();
      }

      int blockCands = std::min(milp.m_param.NumBlocksCand, static_cast<int>(histogram2.size()));
      std::map<int, int>::iterator histIterLower  = histogram2.end();
      std::advance(histIterLower, - blockCands);
      histIter = histogram2.end();

      for (--histIter; histIter != histIterLower;
            --histIter) {
         blockNums.push_back(histIter->second);
      }
   }

   /*
   if(histogram2.size() == 1){
     histIter = histogram2.begin();
     numBlocksCand = (++histIter)->second;

     std::cout << "The block number is "
        << numBlocksCand << std::endl;

   }

   else if(histogram2.size()<=3){

     for(histIter = histogram2.begin(); histIter != histogram2.end();
    ++histIter)
       {
    int nB = histIter->first;
    histIter2 = ++histIter;
    if (nB < histIter2->first)
      {numBlocksCand = (++histIter)->second; }
    else
      {numBlocksCand = nB;}
       }
     std::cout << "The block number is "
          << numBlocksCand << std::endl;

     }

   else

     {
       int gcd = (++histogram2.begin())->second;

       for(histIter = histogram2.begin(); histIter != histogram2.end();
      ++histIter)
    {
      if(gcd !=1)
        gcd = UtilGcd(histIter->second, gcd);
      else
        {
          histIter2 = histIter;
          gcd = UtilGcd(histIter->second, (++histIter2)->second);
        }

    }

       if(gcd != 1){

    int mblockNum = (histogram2.find((histogram2.rbegin())->first))->second;

    std::cout << "The single one is " << mblockNum
   	   <<std::endl;

    std::cout << "The block number is "
   	   << gcd <<" " << gcd*2 << " "
   	   << gcd*4 <<" " << gcd*8 <<std::endl;
       }



     }

   */
   //   UTIL_DELARR(lengthRows);
}
