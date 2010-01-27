//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)       //
//                                                                           //
// Copyright (C) 2002-2007, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

#include "DecompApp.h"
#include "DecompVar.h"



// --------------------------------------------------------------------- //
void DecompApp::startupLog(){
   if(m_param.LogLevel > 0){
      (*m_osLog)
	 << "\n========================================================"
	 << "\n========================================================"
	 <<   "\nWelcome to DECOMP v" << DecompVersion << "."
	 <<   "\n  Copyright (c) 2004-2010 Lehigh University."
	 << "\n========================================================"
	 << "\n========================================================"
	 << "\n";
      
   }
   if(m_param.LogLevel > 1){
      //m_param.dumpSettings(m_osLog);        
   }
}

// --------------------------------------------------------------------- //
/*int DecompApp::createModel(){

   UTIL_MSG(m_param.LogLevel, 2,
	    (*m_osLog) << "Initial Model Setup\n";	    
	    );
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createModel()", m_param.LogDebugLevel, 2);
	    
   //---
   //--- APP: create the user model A = A' union A''
   //---
   APPcreateModel(m_model.objCoeff, m_modelCore, m_modelRelax);

   //---
   //--- TODO: sanity checks that the user gave you all the model
   //--- information that is required, error codes if not
   //---

   //---
   //--- For each model core:
   //---   1.) set row senses and/or bounds (for relaxed too)
   //---   2.) create row hash
   //---   3.) set nBaseRows
   //---   4.) flip to row ordered, if neccessary (for relaxed too)
   //--- TODO: put timer on this ??
   //---
   DecompConstraintSet * modelCore = 0;
   map<int, DecompConstraintSet*>::iterator mdc;
   for(mdc = m_modelCore.begin(); mdc != m_modelCore.end(); mdc++){
      modelCore = mdc->second;
      if(!modelCore->M)
         continue;
      if(modelCore->M->isColOrdered()){      
         modelCore->M->reverseOrdering();
      }
      modelCore->checkSenseAndBound();    
      modelCore->createRowHash();
      modelCore->nBaseRows = modelCore->getNumRows();

      //TODO: make this an option
      //---
      //--- if row/col names are not given, make up default ones
      //---
      int i, j;
      vector<string> & rowNames = modelCore->rowNames;
      vector<string> & colNames = modelCore->colNames;
	  if(rowNames.size() == 0){
	      for(i = 0; i < modelCore->getNumRows(); i++)
		     rowNames.push_back("core(" + UtilIntToStr(i) + ")");
		}
	  if(colNames.size() == 0){
		for(j = 0; j < modelCore->getNumCols(); j++)
			 colNames.push_back("x(" + UtilIntToStr(j) + ")");
		}
   }

   DecompConstraintSet * modelRelax = 0;
   map<int, vector< DecompConstraintSet* > >::iterator mdr;
   vector<DecompConstraintSet*>::iterator vi;
   for(mdr = m_modelRelax.begin(); mdr != m_modelRelax.end(); mdr++){
      for(vi = mdr->second.begin(); vi != mdr->second.end(); vi++){
	 modelRelax = *vi;
	 if(!modelRelax->M)
	    continue;
	 if(modelRelax->M->isColOrdered()){      
	    modelRelax->M->reverseOrdering();
	 }
	 modelRelax->checkSenseAndBound();    

         int i, j;
         vector<string> & rowNames = modelRelax->rowNames;
         vector<string> & colNames = modelRelax->colNames;
		 if(rowNames.size() == 0){
	         for(i = 0; i < modelRelax->getNumRows(); i++)
		        rowNames.push_back("core(" + UtilIntToStr(i) + ")");
		 }
		 if(colNames.size() == 0){
	         for(j = 0; j < modelRelax->getNumCols(); j++)
		        colNames.push_back("x(" + UtilIntToStr(j) + ")");
		 }
      }
   }

   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "createModel()", m_param.LogDebugLevel, 2);

   return 0;//TODO: do return codes
   }*/

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
int DecompApp::generateCuts(const double                       * x, 
			    DecompCutList                      & newCuts){

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
