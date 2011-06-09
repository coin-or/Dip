//===========================================================================//
// This file is part of the DIP Solver Framework.                            //
//                                                                           //
// DIP is distributed under the Eclipse Public License as part of the        //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)       //
//                                                                           //
// Conceptual Design: Matthew Galati, SAS Institute Inc.                     //
//                    Ted Ralphs, Lehigh University                          //
//                                                                           //
// Copyright (C) 2002-2011, Lehigh University, Matthew Galati, Ted Ralphs    //
// All Rights Reserved.                                                      //
//===========================================================================//

#include "DecompApp.h"
#include "DecompVar.h"
#include "DecompConfig.h"

// --------------------------------------------------------------------- //
void DecompApp::startupLog(){
   if(m_param.LogLevel >= 0){
      (*m_osLog)
	 << "\n========================================================"
	 << "\n========================================================"
	 <<   "\nWelcome to the DIP Decomposition Framework"
	 <<   "\nCopyright 2002-2011 Lehigh University and others"
	 <<   "\nAll Rights Reserved"
	 <<   "\nDistributed under the Eclipse Public License 1.0"
	 <<   "\nVersion: " << DIP_VERSION
	 <<   "\nBuild Date: " << __DATE__
#ifdef DIP_SVN_REV
	 <<   "\nRevision Number: " << DIP_SVN_REV
#endif
	 << "\n========================================================"
	 << "\n========================================================"
	 << "\n";
      
   }
   if(m_param.LogLevel > 1){
      //m_param.dumpSettings(m_osLog);        
   }
}

// --------------------------------------------------------------------- //
int DecompApp::generateInitVars(DecompVarList & initVars){ 
   // ---
   // --- this function does nothing by default
   // ---
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "generateInitVars()", m_param.LogDebugLevel, 2);
   UtilPrintFuncEnd(m_osLog, m_classTag,
		      "generateInitVars()", m_param.LogDebugLevel, 2);
   return 0;
}

// --------------------------------------------------------------------- //
int DecompApp::generateCuts(const double  * x, 
			    DecompCutList & newCuts){

   // ---
   // --- this function does nothing by default
   // ---
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "generateCuts()", m_param.LogDebugLevel, 2);
   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "generateCuts()", m_param.LogDebugLevel, 2);
   return 0;
}

/*-------------------------------------------------------------------------*/
void DecompApp::printOriginalColumn(const int index, 
                                    ostream * os) const {
   (*os) << index << " ";
}

/*-------------------------------------------------------------------------*/
void DecompApp::printOriginalSolution(const int              n_cols, 
				      const vector<string> & colNames,
                                      const double         * solution,
                                      ostream              * os) const{
   int i;
   bool hasNames = false;

   //---
   //--- do we have column names?
   //---
   if(colNames.size() > 0)
      hasNames = true;

   (*os) << setiosflags(ios::fixed|ios::showpoint);
   for(i = 0; i < n_cols; i++){
      if(!UtilIsZero(solution[i])){
         printOriginalColumn(i,os);
         if(hasNames)
            (*os) << "\t" << colNames[i]
                  << "\t" << solution[i] << endl;
         else
            (*os) << "\t" << solution[i] << endl;
      }
   }
   (*os) << resetiosflags(ios::fixed|ios::showpoint|ios::scientific); 
}



#if 0
// --------------------------------------------------------------------- //
void DecompApp::setOptimalPoint(vector< vector<double> > & optPoint){

   //
   // ---
   // --- this function does nothing by default
   // ---
   //
   UTIL_DEBUG(m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << " <---- setOptimalPoint()  ----\n";
              );
   UTIL_DEBUG(m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << "  ---- setOptimalPoint()  ---->\n";
              );
}
#endif
