//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, Lehigh University                                 //
//                                                                           //
// Copyright (C) 2002-2009, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

//===========================================================================//
#include "UtilTimer.h"
#include "UtilParameters.h"
//===========================================================================//
#include "VRP_DecompApp.h"
//===========================================================================//
#include "AlpsDecompModel.h"
//===========================================================================//
#include "DecompAlgoC.h"
#include "DecompAlgoPC.h"
#include "DecompAlgoRC.h"
//===========================================================================//

//===========================================================================//
int main(int argc, char ** argv){
   try{
      //---
      //--- create the utility class for parsing parameters
      //---
      UtilParameters utilParam(argc, argv);  

      bool   doCut      = utilParam.GetSetting("doCut",         true);
      bool   doPriceCut = utilParam.GetSetting("doPriceCut",    false);
      string Instance   = utilParam.GetSetting("Instance", ".", "VRP");
      
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
      VRP_DecompApp vrp(utilParam); 
      
      //---
      //--- create the algorithm (a DecompAlgo)
      //---
      DecompAlgo * algo = NULL;
      assert(doCut + doPriceCut == 1);
      
      //---
      //--- create the CPM algorithm object
      //---      
      if(doCut)	 
         algo = new DecompAlgoC(&vrp, &utilParam);
      
      //---
      //--- create the PC algorithm object
      //---
      if(doPriceCut)
         algo = new DecompAlgoPC(&vrp, &utilParam);

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
      cout << "Instance = "  << Instance;
      if(alpsModel.getBestObj() > 1.0e16)
         cout << " Solution = 1.0e20";
      else
         cout << " Solution = " << alpsModel.getBestObj();
      cout << " [ "          << vrp.getBestKnownLB()
	   << " , "          << vrp.getBestKnownUB() << " ]"
	   << " SetupCPU = " << timeSetupCpu
	   << " SolveCPU = " << timeSolveCpu << endl;      
      //double diff = alpsModel.getBestObj() - vrp.getBestKnownUB();
      //CoinAssert(UtilIsZero(diff));
      //CoinAssert(alpsModel.getBestObj() >= vrp.getBestKnownUB());

      //---
      //--- free local memory
      //---
      delete algo;
   }
   catch(CoinError & ex){
      cerr << "COIN Exception [ " << ex.message() << " ]"
           << " at " << ex.fileName()  << ":L" << ex.lineNumber()
           << " in " << ex.className() << "::" << ex.methodName() << endl;
   }
   return 0;
}
      
