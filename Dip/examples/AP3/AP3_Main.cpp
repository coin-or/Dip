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
#include "AP3_DecompApp.h"

#include "AlpsDecompModel.h"
#include "AlpsKnowledgeBroker.h"
#include "DecompAlgoC2.h"
#include "DecompAlgoPC2.h"
#include "DecompAlgoRC.h"

#include "CoinError.hpp"
#include "AlpsTime.h"

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

    bool doModelI   = utilParam.GetSetting("doModelI",     true);
    bool doModelJ   = utilParam.GetSetting("doModelJ",     false);
    bool doModelK   = utilParam.GetSetting("doModelK",     false);

    AlpsTimer timer;
    double    timeSetupReal = 0.0;
    double    timeSetupCpu  = 0.0;
    double    timeSolveReal = 0.0;
    double    timeSolveCpu  = 0.0;

    timer.start();

    //---
    //--- create the user application (a DecompApp)
    //---
    AP3_DecompApp ap3(utilParam); 
    ap3.createModel();
      
    //---
    //--- create the algorithm(s) (a DecompAlgo)
    //---
    DecompAlgoC2  * cut      = NULL;
    DecompAlgoPC2 * priceI   = NULL;
    DecompAlgoPC2 * priceJ   = NULL;
    DecompAlgoPC2 * priceK   = NULL;
    DecompAlgoPC2 * pcI      = NULL;
    DecompAlgoPC2 * pcJ      = NULL;
    DecompAlgoPC2 * pcK      = NULL;
    DecompAlgoRC  * rcI      = NULL;
    DecompAlgoRC  * rcJ      = NULL;
    DecompAlgoRC  * rcK      = NULL;

    if(doCut){	 
      cut = new DecompAlgoC2(&ap3, &utilParam);  
    }
    if(doPrice){
      CoinAssertHint(doModelI || doModelJ || doModelK,
		     "Error: must pick some base model to price");
      if(doModelI){
	 priceI = new DecompAlgoPC2(&ap3, &utilParam, "PRICE", true,
				    AP3_DecompApp::MODEL_I);
      }
      if(doModelJ){
	 priceJ = new DecompAlgoPC2(&ap3, &utilParam, "PRICE", true,
				    AP3_DecompApp::MODEL_J);
      }
      if(doModelK){
	 priceK = new DecompAlgoPC2(&ap3, &utilParam, "PRICE", true,
				    AP3_DecompApp::MODEL_K);
      }
    }
      
    if(doPriceCut){
      CoinAssertHint(doModelI || doModelJ || doModelK,
		     "Error: must pick some base model to price");
      if(doModelI){
	pcI = new DecompAlgoPC2(&ap3, &utilParam,
				AP3_DecompApp::MODEL_I);
      }
      if(doModelJ){
	pcJ = new DecompAlgoPC2(&ap3, &utilParam,
				AP3_DecompApp::MODEL_J);
      }
      if(doModelK){
	pcK = new DecompAlgoPC2(&ap3, &utilParam,
				AP3_DecompApp::MODEL_K);
      }
    }
    if(doRelaxCut){
      CoinAssertHint(doModelI || doModelJ || doModelK,
		     "Error: must pick some base model to price");
      if(doModelI){
	rcI = new DecompAlgoRC(&ap3, &utilParam,
			       AP3_DecompApp::MODEL_I);
      }
      if(doModelJ){
	rcJ = new DecompAlgoRC(&ap3, &utilParam,
			       AP3_DecompApp::MODEL_J);
      }
      if(doModelK){
	rcK = new DecompAlgoRC(&ap3, &utilParam,
			       AP3_DecompApp::MODEL_K);
      }
    }

    if(useAlps){
      //---
      //--- create the driver AlpsDecomp model
      //---
      AlpsDecompModel alpsModel(utilParam);
      if(cut) 
	alpsModel.addDecompAlgo(cut);
      if(priceI) 
	alpsModel.addDecompAlgo(priceI);
      if(priceJ) 
	alpsModel.addDecompAlgo(priceJ);
      if(priceK) 
	alpsModel.addDecompAlgo(priceK);
      if(pcI) 
	alpsModel.addDecompAlgo(pcI);
      if(pcJ) 
	alpsModel.addDecompAlgo(pcJ);         
      if(pcK) 
	alpsModel.addDecompAlgo(pcK);
      if(rcI) 
	alpsModel.addDecompAlgo(rcI);
      if(rcJ) 
	alpsModel.addDecompAlgo(rcJ);         
      if(rcK) 
	alpsModel.addDecompAlgo(rcK);         
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
      cout << "Instance = " << ap3.getInstanceName()
	   << " Solution = " << alpsModel.getBestObj()
	   << " SetupCPU = " << timeSetupCpu
	   << " SolveCPU = " << timeSolveCpu << endl;

      double diff = alpsModel.getBestObj() - ap3.getKnownOptimalBound();
      CoinAssert(UtilIsZero(diff));


    }else{
      //---
      //--- just solve the bounding problem (root node)
      //---
    }

    if(cut)    delete cut;
    if(priceI) delete priceI;
    if(priceJ) delete priceJ;
    if(priceK) delete priceK;
    if(pcI)    delete pcI;
    if(pcJ)    delete pcJ;
    if(pcK)    delete pcK;
    if(rcI)    delete rcI;
    if(rcJ)    delete rcJ;
    if(rcK)    delete rcK;
  }
  catch(CoinError & ex){
    cerr << "COIN Exception:" << ex.message() << endl 
	 << " from method "   << ex.methodName() << endl
	 << " from class "    << ex.className() << endl; 
  }
} 
