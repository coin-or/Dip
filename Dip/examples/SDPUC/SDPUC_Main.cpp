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
#include "SDPUC_DecompApp.h"
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
   try {

      //---
      //--- create the utility class for parsing parameters
      //---
      UtilParameters utilParam(argc, argv);          
      bool doCut          = utilParam.GetSetting("doCut",          true);
      bool doPriceCut     = utilParam.GetSetting("doPriceCut",     false);
      bool doDirect       = utilParam.GetSetting("doDirect",       false);
      int  timeLimit      = utilParam.GetSetting("timeLimit",      60);
      
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
      SDPUC_DecompApp mmkp(utilParam);
      //---
      //--- create the algorithm (a DecompAlgo)
      //---
      DecompAlgo * algo = NULL;
      assert(doCut + doPriceCut == 1);



      //---
      //--- create the CPM algorithm object
      //---      
      if(doCut)	 
         algo = new DecompAlgoC(&mmkp, &utilParam);
      //---
      //--- create the PC algorithm object
      //---
      if(doPriceCut)
         algo = new DecompAlgoPC(&mmkp, &utilParam);

   
      if(doCut && doDirect){
	 timer.stop();
	 timeSetupCpu  = timer.getCpuTime();
	 timeSetupReal = timer.getRealTime();

	 //--- create initial solution
	 //DecompSolution sol = SILCEP_DecompApp::createInitialSolution();

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
	 
	  //---
	 //--- Display solution
	 //---
   
	 SDPUC_Instance & inst = mmkp.getInstance();
	  int i = 0;
	 
	  int         numTimeperiods = inst.m_numTimeperiods;
	  int         numArcs        = inst.m_numArcs;
	  int         numNodes       = inst.m_numNodes;
	   int   numCols        =		 numArcs						//y1-vars
								 + 3 * numTimeperiods * numArcs	//y2-, z-, and x-vars
								 + numTimeperiods * numNodes;		//theta-vars
   int   numRows		  =		 numArcs * numTimeperiods;
   int	 col_yStartIndex  =		 0;
   int	 col_zStartIndex  =		 numArcs*(1+numTimeperiods);
   int	 col_xStartIndex  =		 numArcs * (1 + 2*numTimeperiods);
   int	 col_thetaStartIndex  =	 numArcs*(1 + 3*numTimeperiods) ;
	  double total_op_cost = 0;
	  double total_fix1_cost = 0;
	  double total_fix2_cost = 0;
	  const double * values = (alpsModel.getBestSolution()->getValues());
	  //cout << "__BEGIN SOL__\n"; 
	  for(i = 0; i < numCols; i++){
	  //cout << values[i] << "  ";
	  
	  if(i < numArcs) {
		 total_fix1_cost = total_fix1_cost + values[i] * mmkp.getObjective()[i];
	  }
	  else if(i < col_zStartIndex ) {
		 total_fix2_cost = total_fix2_cost + values[i] * mmkp.getObjective()[i];
	  }
	  else { //if(i > col_xStartIndex && i < col_zStartIndex ) {
		 total_op_cost = total_op_cost + values[i] * mmkp.getObjective()[i];
	  }
	  }	  

   cout << " - - - - \n";
   cout << "Total fix1 cost = " << total_fix1_cost << endl;
   cout << "Total fix2 cost = " << total_fix2_cost << endl;
   cout << "Total operational cost = " << total_op_cost << endl;
   cout << "Number of columns in Core : " << alpsModel.getNumCoreCols() << endl;
   cout << " - - - - " << endl;
   //---

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

