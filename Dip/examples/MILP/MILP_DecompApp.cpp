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
#include "MILP_DecompApp.h"

//===========================================================================//
void MILP_DecompApp::initializeApp(UtilParameters & utilParam)  {
   
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "initializeApp()", m_appParam.LogLevel, 2);
   
   //---
   //--- get application parameters
   //---   
   m_appParam.getSettings(utilParam);   
   if(m_appParam.LogLevel >= 1)
      m_appParam.dumpSettings();

   //---
   //--- read MILP instance (mps format)
   //---
   string fileName = m_appParam.DataDir 
      + UtilDirSlash() + m_appParam.Instance;   

   m_mpsIO.messageHandler()->setLogLevel(m_param.LogLpLevel);

   int rstatus = m_mpsIO.readMps(fileName.c_str());   
   if(rstatus < 0){
      cerr << "Error: Filename = " << fileName << " failed to open." << endl;
      throw UtilException("I/O Error.", "initalizeApp", "MILP_DecompApp");
   }
   if(m_appParam.LogLevel >= 2)
      (*m_osLog) << "Objective Offset = " 
                 << UtilDblToStr(m_mpsIO.objectiveOffset()) << endl;

   //---
   //--- set best known lb/ub
   //---
   double offset = m_mpsIO.objectiveOffset();
   setBestKnownLB(m_appParam.BestKnownLB + offset);
   setBestKnownUB(m_appParam.BestKnownUB + offset);

   //---
   //--- create models
   //---
   createModels();

   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "initializeApp()", m_appParam.LogLevel, 2);
}

//===========================================================================//
void MILP_DecompApp::createModels(){

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createModels()", m_appParam.LogLevel, 2);
   

   //---
   //--- seed random number generator
   //---
   srand(m_appParam.RandomSeed);
   
   //---
   //--- how many rows to put into relaxation
   //---
   int       i, c, nRowsRelax, nRowsCore;
   const int nRows = m_mpsIO.getNumRows();
   const int nCols = m_mpsIO.getNumCols();
   nRowsRelax = static_cast<int>(ceil(nRows * m_appParam.RelaxPercent));
   nRowsRelax = min(nRows-1, nRowsRelax);
   nRowsCore  = nRows - nRowsRelax; 

   UTIL_MSG(m_appParam.LogLevel, 2,
            (*m_osLog) << "Instance    = " << m_appParam.Instance << endl;
            (*m_osLog) << " nRows      = " << nRows     << endl;
            (*m_osLog) << " nCols      = " << nCols     << endl;
            (*m_osLog) << " nRowsCore  = " << nRowsCore << endl;
            (*m_osLog) << " nRowsRelax = " << nRowsRelax 
            << " [ " << 100*nRowsRelax/nRows << " % ]" << endl;
            );


   //---
   //--- pick nRowsRelax random rows
   //---
   set<int>    relaxRows;
   while(static_cast<int>(relaxRows.size()) < nRowsRelax)
      relaxRows.insert(UtilURand(0,nRows-1));
   assert(static_cast<int>(relaxRows.size()) == nRowsRelax);

   //---
   //--- setup markers for core and relax rows
   //---
   int * rowsMarker = new int[nRows];
   int * rowsCore   = new int[nRowsCore];
   int * rowsRelax  = new int[nRowsRelax];
   UtilFillN(rowsMarker, nRows, 0);
   
   int nRowsCoreTmp  = 0;
   int nRowsRelaxTmp = 0;
   set<int>::iterator it;
   for(it = relaxRows.begin(); it != relaxRows.end(); it++)
      rowsMarker[*it] = 1;
   for(i = 0; i < nRows; i++){
      if(rowsMarker[i])
         rowsRelax[nRowsRelaxTmp++] = i;
      else
         rowsCore[nRowsCoreTmp++]   = i;
   }
   assert((nRowsRelaxTmp + nRowsCoreTmp) == nRows);

   UTIL_MSG(m_appParam.LogLevel, 3,
            (*m_osLog) << "Core  Rows:";
            for(i = 0; i < nRowsCore; i++)
               (*m_osLog) << rowsCore[i] << " ";
            (*m_osLog) << "\nRelax Rows:";
            for(i = 0; i < nRowsRelax; i++)
               (*m_osLog) << rowsRelax[i] << " ";
            (*m_osLog) << "\n";
            );
      
   //---
   //--- Construct the objective function.
   //---
   m_objective = new double[nCols];
   if(!m_objective)
      throw UtilExceptionMemory("createModels", "MMKP_DecompApp");
   memcpy(m_objective, 
          m_mpsIO.getObjCoefficients(), nCols * sizeof(double));
   setModelObjective(m_objective);

   //---
   //--- Construct the core matrix.
   //---
   m_modelRandCore.M = new CoinPackedMatrix(false, 0.0, 0.0);
   if(!m_modelRandCore.M)
      throw UtilExceptionMemory("createModels", "MILP_DecompApp");
   m_modelRandCore.reserve(nRowsCore, nCols);
   m_modelRandCore.M->submatrixOf(*m_mpsIO.getMatrixByRow(), 
                                   nRowsCore, rowsCore);   

   //---
   //--- Construct the relaxation matrix.
   //---
   m_modelRandRelax.M = new CoinPackedMatrix(false, 0.0, 0.0);
   if(!m_modelRandRelax.M)
      throw UtilExceptionMemory("createModels", "MILP_DecompApp");
   m_modelRandRelax.reserve(nRowsRelax, nCols);
   m_modelRandRelax.M->submatrixOf(*m_mpsIO.getMatrixByRow(), 
                                   nRowsRelax, rowsRelax);
         
   //---
   //--- set the row upper and lower bounds
   //--- set the col upper and lower bounds
   //---
   const double * rowLB = m_mpsIO.getRowLower();
   const double * rowUB = m_mpsIO.getRowUpper();
   const double * colLB = m_mpsIO.getColLower();
   const double * colUB = m_mpsIO.getColUpper();
   for(i = 0; i < nRowsCore; i++){
      m_modelRandCore.rowLB.push_back(rowLB[rowsCore[i]]);
      m_modelRandCore.rowUB.push_back(rowUB[rowsCore[i]]);
   }
   for(i = 0; i < nRowsRelax; i++){
      m_modelRandRelax.rowLB.push_back(rowLB[rowsRelax[i]]);
      m_modelRandRelax.rowUB.push_back(rowUB[rowsRelax[i]]);
   }
   copy(colLB, colLB + nCols, back_inserter( m_modelRandCore.colLB)  );
   copy(colUB, colUB + nCols, back_inserter( m_modelRandCore.colUB)  );
   copy(colLB, colLB + nCols, back_inserter( m_modelRandRelax.colLB) );
   copy(colUB, colUB + nCols, back_inserter( m_modelRandRelax.colUB) );

   //---
   //--- big fat hack... we don't deal with dual rays yet,
   //---  so, we assume subproblems are bounded
   //---
   //--- NOTE: might also need to tighten LBs
   //---
   for(c = 0; c < nCols; c++){
      //printf("c: %5d lb: %8.5f ub:%8.5f\n", c, colLB[c], colUB[c]);
      if(colUB[c] > 1.0e15){
	 printf("colUB[%d]: %g\n", c, colUB[c]);
	 m_modelRandRelax.colUB[c] = 1000;
	 m_modelRandCore.colUB[c]  = 1000;
      }
   }
 
   //---
   //--- set the indices of the integer variables of modelRelax
   //---
   const char * integerVars = m_mpsIO.integerColumns();
   if(integerVars){
      for(c = 0; c < nCols; c++){
         if(integerVars[c]){
            m_modelRandCore.integerVars.push_back(c);
            m_modelRandRelax.integerVars.push_back(c);         
         }
      }
   }

   //---
   //--- set core and relax systems in framework
   //---
   setModelCore (&m_modelRandCore,  "core");
   setModelRelax(&m_modelRandRelax, "relax");


   //---
   //--- free up local memory
   //---
   UTIL_DELARR(rowsMarker);
   UTIL_DELARR(rowsCore);
   UTIL_DELARR(rowsRelax);
   
   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "createModels()", m_appParam.LogLevel, 2);   
}


