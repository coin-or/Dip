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
#include "ATM_DecompApp.h"
//===========================================================================//
#include "AlpsDecompModel.h"
//===========================================================================//
#include "DecompAlgoC.h"
#include "DecompAlgoPC.h"
#include "DecompAlgoRC.h"
//===========================================================================//
#include "UtilTimer.h"

//===========================================================================//
int main(int argc, char ** argv){
   try{

      //---
      //--- create the utility class for parsing parameters
      //---
      UtilParameters utilParam(argc, argv);  
            
      bool doGenRandom    = utilParam.GetSetting("doGenRandom",    false);
      int  randSeed       = utilParam.GetSetting("randSeed",       1    );
      int  randNumAtms    = utilParam.GetSetting("randNumAtms",    5    );
      int  randNumDates   = utilParam.GetSetting("randNumDates",   10   );
      
      bool doCut          = utilParam.GetSetting("doCut",          true);
      bool doPriceCut     = utilParam.GetSetting("doPriceCut",     false);
      bool doDirect       = utilParam.GetSetting("doDirect",       false);
      
      UtilTimer timer;
      double    timeSetupReal = 0.0;
      double    timeSetupCpu  = 0.0;
      double    timeSolveReal = 0.0;
      double    timeSolveCpu  = 0.0;

      //---
      //--- start overall timer
      //---
      timer.start();
      if(doGenRandom){
	 //---
	 //--- generate a random instance
	 //---
	 ATM_Instance instance;
	 instance.generateRandom(randNumAtms, randNumDates, randSeed);	 
      }
      else{
         //---
         //--- create the user application (a DecompApp)
         //---      
         ATM_DecompApp atm(utilParam); 
                  
         //---
         //--- create the algorithm (a DecompAlgo)
         //---
         DecompAlgo * algo = NULL;
         assert(doCut + doPriceCut == 1);

         //---
         //--- create the CPM algorithm object
         //---      
         if(doCut)
            algo = new DecompAlgoC(&atm, &utilParam);
         
         //---
         //--- create the PC algorithm object
         //---
         if(doPriceCut)
            algo = new DecompAlgoPC(&atm, &utilParam);
         
         
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
	    int             status = 0;
            AlpsDecompModel alpsModel(utilParam, algo);
	    
	    timer.stop();
            timeSetupCpu  = timer.getCpuTime();
            timeSetupReal = timer.getRealTime();
	    
            //---
            //--- solve
            //---
            timer.start();      
            status = alpsModel.solve();
            timer.stop();
            timeSolveCpu  = timer.getCpuTime();
            timeSolveReal = timer.getRealTime();

            //TODO: move doDirect solve into alpsModel so access
            //  solution the same way?

            //---
            //--- sanity check
            //---
            cout << setiosflags(ios::fixed|ios::showpoint);
	    cout << "Status= " << status 
		 << " BestLB= " << setw(10) 
                 << UtilDblToStr(alpsModel.getGlobalLB(),2)
                 << " BestUB= " << setw(10)
                 << UtilDblToStr(alpsModel.getGlobalUB(),2)        
                 << " Nodes= " << setw(6) 
                 << alpsModel.getNumNodesProcessed()
                 << " SetupCPU= "  << timeSetupCpu
                 << " SolveCPU= "  << timeSolveCpu 
                 << " TotalCPU= "  << timeSetupCpu + timeSolveCpu
                 << " SetupReal= " << timeSetupReal
                 << " SolveReal= " << timeSolveReal
                 << " TotalReal= " << timeSetupReal + timeSolveReal
                 << endl;      

	    //---
	    //--- now, initialize direct solve with best
	    //---   solution to PC
	    //--- TODO: only useful if stop early on time or nodes   
	    //--- TODO: cbc currently doesn't use warm-start, only cpx
	    //---
	    //DecompAlgo * algoC = new DecompAlgoC(&atm, &utilParam);
	    //algoC->solveDirect(algo->getXhatIPBest());
	    //delete algoC;
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
   }
   return 0;
}
      
