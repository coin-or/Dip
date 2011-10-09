//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, Lehigh University                                 //
//                                                                           //
// Copyright (C) 2002-2007, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

//===========================================================================//
#include "UtilParameters.h"
#include "MAD_DecompApp.h"

#include "AlpsDecompModel.h"
#include "AlpsKnowledgeBroker.h"
#include "DecompAlgoC2.h"
#include "DecompAlgoPC2.h"
#include "DecompAlgoRC.h"

#include "CoinError.hpp"
#include "AlpsTime.h"

//===========================================================================//
//#define CREATE_FULL

//===========================================================================//
int main(int argc, char ** argv){
   try{
      
      //---
      //--- create the utility class for parsing parameters
      //---
      UtilParameters utilParam(argc, argv);  
      
      bool useAlps    = utilParam.GetSetting("useAlps",      true); 
      
      bool doCut      = utilParam.GetSetting("doCut",        true);
      bool doPrice    = utilParam.GetSetting("doPrice",      false);
      bool doPriceCut = utilParam.GetSetting("doPriceCut",   false);
      bool doRelaxCut = utilParam.GetSetting("doRelaxCut",   false);

      string Instance = utilParam.GetSetting("Instance", ".", "MAD");

      AlpsTimer timer;
      double    timeSetupReal = 0.0;
      double    timeSetupCpu  = 0.0;
      double    timeSolveReal = 0.0;
      double    timeSolveCpu  = 0.0;

      timer.start();

      //---
      //--- create the user application (a DecompApp)
      //---
      MAD_DecompApp mad(utilParam); 
      mad.createModel();

#ifdef CREATE_FULL
      {
         //create full mps file for debugging
         string DataSubDir = utilParam.GetSetting("DataSubDir", ".",  "MAD");
         string::size_type pos = Instance.find_first_of(".");
         
         DecompAlgoC2 cutMps(&mad, &utilParam);
         cutMps.initSetup(&utilParam, "MAD");
         
         string fileName = DataSubDir + "_" + Instance.substr(0,pos);
         cutMps.createFullMps(fileName);
         exit(1);
      }
#endif



      
      //---
      //--- create the algorithm(s) (a DecompAlgo)
      //---
      DecompAlgoC2  * cut      = NULL;
      DecompAlgoPC2 * price    = NULL;
      DecompAlgoPC2 * pc       = NULL;
      DecompAlgoRC  * rc       = NULL;
      
      if(doCut){
         cut = new DecompAlgoC2(&mad, &utilParam);  
      }
      if(doPrice){
         price = new DecompAlgoPC2(&mad, &utilParam, "PRICE");
      }      
      if(doPriceCut){
         pc = new DecompAlgoPC2(&mad, &utilParam);
      }
      if(doRelaxCut){
         rc = new DecompAlgoRC(&mad, &utilParam);
      }

      if(useAlps){
         //---
         //--- create the driver AlpsDecomp model
         //---
         AlpsDecompModel alpsModel(utilParam);
         if(cut) 
            alpsModel.addDecompAlgo(cut);
         if(price) 
            alpsModel.addDecompAlgo(price);
         if(pc) 
            alpsModel.addDecompAlgo(pc);
         if(rc) 
            alpsModel.addDecompAlgo(rc);
         timer.stop();
         timeSetupCpu  = timer.getCpuTime();
         timeSetupReal = timer.getWallClock();
      
         timer.start();
         alpsModel.solve();
         timer.stop();
         timeSolveCpu  = timer.getCpuTime();
         timeSolveReal = timer.getWallClock();
      
         //---
         //--- sanity check
         //---
         cout << "Instance = "  << Instance
              << " NRows = "    << mad.getNOrigRows()
              << " Border = "   << mad.getNOrigRows() + alpsModel.getBestObj()
              << " Solution = " << alpsModel.getBestObj()
              << " [ " << mad.getBestKnownLB()
              << " , "   << mad.getBestKnownUB() << " ]"
              << " SetupCPU = " << timeSetupCpu
              << " SolveCPU = " << timeSolveCpu << endl;
         
         //double diff = alpsModel.getBestObj() - mad.getKnownOptimalBound();
         //CoinAssert(UtilIsZero(diff));

      }else{
         //---
         //--- just solve the bounding problem (root node)
         //---
      }

      if(cut)   delete cut;
      if(price) delete price;
      if(pc)    delete pc;
      if(rc)    delete rc;
   }
   catch(CoinError & ex){
      cerr << "COIN Exception:" << ex.message() << endl 
           << " from method "   << ex.methodName() << endl
           << " from class "    << ex.className() << endl; 
   }
}  
