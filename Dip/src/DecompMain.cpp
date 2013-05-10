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

#include "DecompMainOneThread.h"

using namespace std; 

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
      
      //---
      //--- create the user application (a DecompApp)
      //---
      DecompApp milp(utilParam); 

	  
      std::vector<int> blockNumCandidates; 



      milp.startupLog();

      milp.preprocessApp(utilParam, blockNumCandidates); 
      
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

	  pthread_t* threads  = new pthread_t[numThreads +1] ; 
	  CoinAssertHint(threads, "Error: Out of Memory"); 


	  printf("===== START Concurrent Computations Process. =====\n");

	  for(int i = 0 ; i < (numThreads + 1); i++){
	    
	    std::vector<DecompApp> milpArray(static_cast<int>(numThreads + 1),milp);


	    std::vector<DecompMainParam> decompMainParamArray(static_cast<int>(numThreads + 1),
							    decompMainParam);


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

	    /*
	    pthread_create(&threads[i],
			   NULL ,
	        	   (void* ) DecompMainOneThread,			   
			   (void* ) &decompMainParamArray[i]);
	    //(vDecompMainOneThread(milpArray[i], utilParam, timer, decompMainParamArray[i]),
	    */
	   
	    
	    
	    std::cout << "==========================================" << std::endl;
	    std::cout << "======   The thread try is  ============= " << std::endl;
	    std::cout << "============" << i <<      "============= " << std::endl;
	    std::cout << "======   The block number is ============ " << std::endl;	
	    std::cout << "===========" << milpArray[i].m_param.NumBlocks <<"============" << std::endl;
	    std::cout << "=========== Branch-and-Cut  " << decompMainParamArray[i].doCut
		      <<"============" << std::endl;
	    std::cout << "=========== Branch-and-Price  " << decompMainParamArray[i].doPriceCut 
		      <<"============" << std::endl;

	    std::cout << "                                          " << std::endl;
	    std::cout << "                                          " << std::endl;
	    std::cout << "                                          " << std::endl;
	    std::cout << "                                          " << std::endl;
	    	    	  	  
	  }
	}
      else{


	decompMainParam.doCut        = utilParam.GetSetting("doCut",        false);

	decompMainParam.doPriceCut   = utilParam.GetSetting("doPriceCut",   true);

	decompMainParam.doDirect     = utilParam.GetSetting("doDirect",     false);
      
	DecompMainOneThread(milp, utilParam, timer, decompMainParam); 

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
