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
#include "UtilParameters.h"
//===========================================================================//
#include "SmallIP_DecompApp.h"
//===========================================================================//
#include "AlpsDecompModel.h"
//===========================================================================//
#include "DecompAlgoC.h"
#include "DecompAlgoPC.h"
#include "DecompAlgoRC.h"

//===========================================================================//
int main(int argc, char ** argv){
   try{
      
      //---
      //--- create the utility class for parsing parameters
      //---
      UtilParameters utilParam(argc, argv);  

      bool doCut          = utilParam.GetSetting("doCut",          true);
      bool doPriceCut     = utilParam.GetSetting("doPriceCut",     false);

      //---
      //--- create the user application (a DecompApp)
      //---      
      SmallIP_DecompApp sip(utilParam); 
      
      //---
      //--- create the algorithm (a DecompAlgo)
      //---
      DecompAlgo * algo = NULL;
      assert(doCut + doPriceCut == 1);

      //---
      //--- create the CPM algorithm object
      //---      
      if(doCut)
	 algo = new DecompAlgoC(&sip, &utilParam);  

      //---
      //--- create the PC algorithm object
      //---      
      if(doPriceCut)
	 algo = new DecompAlgoPC(&sip, &utilParam);
      
      //---
      //--- create the driver AlpsDecomp model
      //---
      AlpsDecompModel alpsModel(utilParam, algo);
      
      //---
      //--- solve
      //---      
      alpsModel.solve();

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
