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


using namespace std; 

//===========================================================================//
int main(int argc, char ** argv){
   try{
      
      //---
      //--- create the utility class for parsing parameters
      //---
      UtilParameters utilParam(argc, argv);  

      bool doCut        = utilParam.GetSetting("doCut",        false);
      bool doPriceCut   = utilParam.GetSetting("doPriceCut",   true);
      bool doDirect     = utilParam.GetSetting("doDirect",     false);
      
      UtilTimer timer;
      double    timeSetupReal = 0.0;
      double    timeSetupCpu  = 0.0;
      double    timeSolveReal = 0.0;
      double    timeSolveCpu  = 0.0;      

      //---
      //--- start overall timer
      //---
      timer.start();
      
      //---
      //--- create the user application (a DecompApp)
      //---
      DecompApp milp(utilParam); 

      //---
      //--- put the one of the functions in the constructor into the main
      //---
      milp.initializeApp(utilParam); 

      milp.startupLog(); 

      //---
      //--- create the algorithm (a DecompAlgo)
      //---
      DecompAlgo * algo = NULL;
      if((doCut + doPriceCut) != 1)
         throw UtilException("doCut or doPriceCut must be set", 
                             "main", "main");
      //assert(doCut + doPriceCut == 1);

      //---
      //--- create the CPM algorithm object
      //---      
      if(doCut)	 
         algo = new DecompAlgoC(&milp, &utilParam);
      
      //---
      //--- create the PC algorithm object
      //---
      if(doPriceCut)
         algo = new DecompAlgoPC(&milp, &utilParam);
      
      if(doCut && doDirect){
	 timer.stop();
	 timeSetupCpu  = timer.getCpuTime();
	 timeSetupReal = timer.getRealTime();
	 
	 //---
	 //--- solve
	 //---
	 timer.start();      
	 algo->solveDirect();
	 timer.stop();
	 timeSolveCpu  = timer.getCpuTime();
	 timeSolveReal = timer.getRealTime();
      }
      else{
	 //---
	 //--- create the driver AlpsDecomp model
	 //---
	 AlpsDecompModel alpsModel(utilParam, algo);
	 
	 timer.stop();
	 timeSetupCpu  = timer.getCpuTime();
	 timeSetupReal = timer.getRealTime();
	 
	 //---
	 //--- solve
	 //---
	 timer.start();      
	 alpsModel.solve();
	 timer.stop();
	 timeSolveCpu  = timer.getCpuTime();
	 timeSolveReal = timer.getRealTime();

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
	      << " SetupCPU      = " << timeSetupCpu << endl
	      << " SolveCPU      = " << timeSolveCpu << endl
	      << " TotalCPU      = " << timeSetupCpu + timeSolveCpu << endl
	      << " SetupWallclock= " << timeSetupReal << endl
	      << " SolveWallclock= " << timeSolveReal << endl
	      << " TotalWallclock= " << timeSetupReal + timeSolveReal    ; 
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
   catch(CoinError & ex){
      cerr << "COIN Exception [ " << ex.message() << " ]"
           << " at " << ex.fileName()  << ":L" << ex.lineNumber()
           << " in " << ex.className() << "::" << ex.methodName() << endl;
      return 1;
   }
   return 0;
} 
