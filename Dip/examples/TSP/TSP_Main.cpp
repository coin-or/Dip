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
#include "TSP_DecompApp.h"
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
      
      bool   doCut      = utilParam.GetSetting("doCut",        true);
      bool   doPriceCut = utilParam.GetSetting("doPriceCut",   false);      

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
      TSP_DecompApp tsp(utilParam); 

      //---
      //--- create the algorithm (a DecompAlgo)
      //---
      DecompAlgo * algo = NULL;
      assert(doCut + doPriceCut == 1);

      //---
      //--- create the CPM algorithm object
      //---      
      if(doCut)	 
         algo = new DecompAlgoC(&tsp, &utilParam);

      //---
      //--- create the PC algorithm object
      //---
      if(doPriceCut)
         algo = new DecompAlgoPC(&tsp, &utilParam);

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

      //---
      //--- sanity check
      //---
      cout << setiosflags(ios::fixed|ios::showpoint);
      cout << "Status= "   << status  
           << " BestLB=  " << setw(10) 
           << UtilDblToStr(alpsModel.getGlobalLB(),5)
           << " BestUB= " << setw(10)
           << UtilDblToStr(alpsModel.getGlobalUB(),5)        
           << " Nodes= " << setw(6) 
           << alpsModel.getNumNodesProcessed()
           << " SetupCPU= "  << timeSetupCpu
           << " SolveCPU= "  << timeSolveCpu 
           << " TotalCPU= "  << timeSetupCpu + timeSolveCpu
           << " SetupReal= " << timeSetupReal
           << " SolveReal= " << timeSolveReal
           << " TotalReal= " << timeSetupReal + timeSolveReal
           << endl;      
	 
      if(status == AlpsExitStatusOptimal && tsp.getBestKnownUB() < 1.0e50){
         //---
         //--- the assumption here is that the BestKnownLB/UB is optimal
         //---
         double diff 
            = fabs(tsp.getBestKnownUB() - alpsModel.getGlobalUB());
         if(diff > 1.0e-4){
            cerr << "ERROR. BestKnownUB= " << tsp.getBestKnownUB()
                 << " but DECOMP claims GlobalUB= " 
                 << alpsModel.getGlobalUB() << endl;
            throw UtilException("Invalid claim of optimal.", 
                                "main", "DECOMP");
         }
      }
	 
      //---
      //--- free local memory
      //---
      delete algo;

   }
   catch(CoinError & ex){
      cerr << "COIN Exception [ " << ex.message() << " ]"
           << " at " << ex.fileName()  << ":L" << ex.lineNumber()
           << " in " << ex.className() << "::" << ex.methodName() << endl;
      return 1;
   }
   return 0;
}
      
