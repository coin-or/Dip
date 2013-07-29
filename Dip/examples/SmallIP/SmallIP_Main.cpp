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
#if defined(VERSION1)
#include "SmallIP_DecompApp.h"
#elif defined (VERSION2)
#include "SmallIP_DecompApp2.h"
#endif
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

      bool doCut          = utilParam.GetSetting("doCut",          false);
      bool doPriceCut     = utilParam.GetSetting("doPriceCut",     false);
      bool doRelaxCut     = utilParam.GetSetting("doRelaxCut",     false);

      //---
      //--- create the user application (a DecompApp)
      //---      
      SmallIP_DecompApp sip(utilParam); 
      
      //---
      //--- create the algorithm (a DecompAlgo)
      //---
      DecompAlgo * algo = NULL;

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
      //--- create the PC algorithm object
      //---      
      if(doRelaxCut)
	 algo = new DecompAlgoRC(&sip, &utilParam);
      
      //---
      //--- create the driver AlpsDecomp model
      //---
      AlpsDecompModel alpsModel(utilParam, algo);
      
      //---
      //--- solve
      //---          
      alpsModel.solve();

      //---
      //--- sanity check that optimal solution is 3.0
      //---
      double epsilon   = 1.0e-5;
      double optimalUB = 3.0;      
      double diffUB    = alpsModel.getGlobalUB() - optimalUB;
      if(alpsModel.getSolStatus() != AlpsExitStatusOptimal ||
         fabs(diffUB)              > epsilon){
         throw UtilException("SmallIP bad solution.", "main", "SmallIP");
      }

      //---
      //--- get optimal solution
      //---      
      const DecompSolution * solution = alpsModel.getBestSolution();
      cout << "Optimal Solution" << endl;
      solution->print();
      if(fabs(solution->getQuality() - alpsModel.getGlobalUB()) > epsilon){
         throw UtilException("Best bound and best solution not matching.",
                             "main", "SmallIP");
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
