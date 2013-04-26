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
#include "DecompVar.h"
#include "SmallIP_DecompApp.h"
//===========================================================================//
//parameters
const int LogLevel   = 1;
const int whichRelax = 1; //1 matches book chapter and thesis

//===========================================================================//
void SmallIP_DecompApp::createModels(){

   //---
   //--- This function does the work to create the different models 
   //---  that will be used. This memory is owned by the user. It will
   //---  be passed to the application interface and used by the algorithms.
   //---
   UtilPrintFuncBegin(m_osLog, m_classTag, "createModels()", LogLevel, 2);
   
   //--- 
   //--- Small Integer Program from: 
   //---   "Integer Programming: Theory and Practice", 
   //---   "Decomposition in Integer Linear Programming", Chapter 4
   //---
   //---  min  1 x[0]
   //--- 
   //---  s.t. 7 x[0] -  1 x[1] >= 13      (4.5)
   //---                 1 x[1] >= 1       (4.6)
   //---      -1 x[0] +  1 x[1] >= -3      (4.7)
   //---      -4 x[0] -  1 x[1] >= -27     (4.8)
   //---              -  1 x[1] >= -5      (4.9)
   //---     0.2 x[0] -  1 x[1] >= -4      (4.10)
   //---
   //---      -1 x[0] -  1 x[1] >= -8      (4.11)
   //---    -0.4 x[0] +  1 x[1] >= 0.3     (4.12)
   //---       1 x[0] +  1 x[1] >= 4.5     (4.13)
   //---       3 x[0] +  1 x[1] >= 9.5     (4.14)
   //---    0.25 x[0] -  1 x[1] >= -3      (4.15)
   //---         x[0], x[1]     >=0, <= 6  (4.16)
   //---         x[0], x[1]     integer    (4.17)
   //---
   //--- Model 1
   //---  Q'  = { x in R^2 | x satisfies (4.5-10, 16, 17} modelRelax
   //---  Q'' = { x in R^2 | x satisfies (4.11-16       } modelCore
   //---
   //--- Model 2
   //---  Q'  = { x in R^2 | x satisfies (4.11-17)      } modelRelax
   //---  Q'' = { x in R^2 | x satisfies (4.5-10, 16)   } modelCore
   //---

   //---
   //--- Construct the objective function (the original problem is 
   //---  a maximization, so we flip the sign to make it minimization).
   //---
   const int numCols = 2;
   m_objective = new double[numCols];
   if(!m_objective)
      throw UtilExceptionMemory("createModels", "SmallIP_DecompApp");      
   m_objective[0] = 1.0;
   m_objective[1] = 0.0;
   setModelObjective(m_objective);

   //---
   //--- build matrix part 1 (4.5 -10,16)
   //--- build matrix part 2 (4.11-15,16)
   //---
   const int numNzs1  = 10;
   const int numNzs2  = 10;
   const int numRows1 =  6;
   const int numRows2 =  5;
   bool      isRowOrdered  = false;
   
   int       rowIndices1[numNzs1] = {0,0,1,2,2,3,3,4,5,5};
   int       colIndices1[numNzs1] = {0,1,1,0,1,0,1,1,0,1};
   double    elements1  [numNzs1]   = { 7.0, -1.0,  1.0, -1.0,  1.0, 
					-4.0, -1.0, -1.0, 0.2, -1.0};   

   int       rowIndices2[numNzs2] = {0,0,1,1,2,2,3,3,4,4};
   int       colIndices2[numNzs2] = {0,1,0,1,0,1,0,1,0,1};
   double    elements2  [numNzs2] = {-1.0, -1.0, -0.4, 1.0,   1.0, 
				     1.0 ,  3.0,  1.0, 0.25, -1.0}; 

   m_modelPart1.M = new CoinPackedMatrix(isRowOrdered,
					 rowIndices1, colIndices1,
					 elements1,   numNzs1);   
   m_modelPart2.M = new CoinPackedMatrix(isRowOrdered,
					 rowIndices2, colIndices2,
					 elements2,   numNzs2);
   
   //---
   //--- set the row upper and lower bounds of part 1
   //---
   double   rowLB1[numRows1] = {13.0, 1.0, -3.0, -27.0, -5.0, -4.0};
   std::fill_n(back_inserter(m_modelPart1.rowUB), numRows1, DecompInf);
   std::copy  (rowLB1, rowLB1 + numRows1, back_inserter(m_modelPart1.rowLB));

   //---
   //--- set the column upper and lower bounds of part 1
   //---
   std::fill_n(back_inserter(m_modelPart1.colLB), numCols, 0.0);
   std::fill_n(back_inserter(m_modelPart1.colUB), numCols, 6.0);

   //---
   //--- set the integer variables for part 1
   //---
   m_modelPart1.integerVars.push_back(0);
   m_modelPart1.integerVars.push_back(1);

   //---
   //--- set the row upper and lower bounds of part 2
   //---
   double   rowLB2[numRows2] = {-8.0, 0.3, 4.5, 9.5, -3.0};
   std::fill_n(back_inserter(m_modelPart2.rowUB), numRows2, DecompInf);
   std::copy  (rowLB2, rowLB2 + numRows2, back_inserter(m_modelPart2.rowLB));

   //---
   //--- set the column upper and lower bounds of part 2
   //---
   std::fill_n(back_inserter(m_modelPart2.colLB), numCols, 0.0);
   std::fill_n(back_inserter(m_modelPart2.colUB), numCols, 6.0);

   //---
   //--- set the integer variables for part 2
   //---
   m_modelPart2.integerVars.push_back(0);
   m_modelPart2.integerVars.push_back(1);
   
   switch(whichRelax){
   case 1:
      {
	 //---
	 //--- Set Model 1
	 //---  Q'  = { x in R^2 | x satisfies (4.5-10, 16, 17} modelRelax
	 //---  Q'' = { x in R^2 | x satisfies (4.11-16       } modelCore
	 //---
	 setModelRelax(&m_modelPart1, "RELAX1");
	 setModelCore (&m_modelPart2, "CORE2");
      }
      break;
   case 2:
      {	 
	 //---
	 //--- Set Model 2
	 //---  Q'  = { x in R^2 | x satisfies (4.11-17)      } modelRelax
	 //---  Q'' = { x in R^2 | x satisfies (4.5-10, 16)   } modelCore
	 //---
	 setModelRelax(&m_modelPart2, "RELAX2");
	 setModelCore (&m_modelPart1, "CORE1");
      }
      break;
   }
   UtilPrintFuncEnd(m_osLog, m_classTag, "createModels()", LogLevel, 2);
}

//--------------------------------------------------------------------- //
int SmallIP_DecompApp::generateInitVars(DecompVarList & initVars){

   //---
   //--- generateInitVars is a virtual method and can be overriden
   //---   if the user has some idea how to initialize the list of 
   //---   initial variables (columns in the DW master)
   //---

   UtilPrintFuncBegin(m_osLog, m_classTag, "generateInitVars()", LogLevel, 2);

   //---
   //--- To follow the example in the chapter: (4,1) and (5,5)
   //---
   int    ind [2] = {0,1};
   double els1[2] = {4.0,1.0};
   double els2[2] = {5.0,5.0};
   initVars.push_back(new DecompVar(2, ind, els1, 4.0));
   initVars.push_back(new DecompVar(2, ind, els2, 5.0));
   
   UtilPrintFuncEnd(m_osLog, m_classTag, "generateInitVars()", LogLevel, 2);
   return static_cast<int>(initVars.size());
}

/*
//--------------------------------------------------------------------- //
int SmallIP_DecompApp::generateCuts(const double              * x, 
                                    const DecompConstraintSet & modelCore,
                                    const DecompConstraintSet & modelRelax,
                                    DecompCutList             & newCuts){


   //
   //---
   //--- For the sake of illustration, store the facets of P,
   //---   and send back any violated cuts.
   //---
   //--- Integer Hull = {(3,2),(3,3),(4,2),(4,),(5,3)}
   //---    0 +  1 x[0] - 1 x[1] >= 0
   //---   -3 +  1 x[0] + 0 x[1] >= 0
   //---   -2 +  0 x[0] + 1 x[1] >= 0
   //---    2 + -1 x[0] + 1 x[1] >= 0
   //---    8 + -1 x[0] - 1 x[1] >= 0
   //---

   const static double cuts[5][2] = {{ 1, -1},
				     { 1,  0},
				     { 0,  1},
				     {-1,  1},
				     {-1, -1}};
   const static double rhs[5]     = {0,3,2,-2,-8};

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "generateCuts()", m_param.LogDebugLevel, 2);


   int r;
   for(r = 0; r < 5; r++){
      if( ((x[0] * cuts[r][0]) + (x[1] * cuts[r][1])) < (rhs[r] - 0.00001)){
	 CoinPackedVector cut;
	 cut.insert(0, cuts[r][0]);
	 cut.insert(1, cuts[r][1]);
    
	 OsiRowCut rc;
	 rc.setRow(cut);      
	 rc.setLb(rhs[r]);
	 rc.setUb(DecompInf);
      
	 DecompCutOsi * decompCut = new DecompCutOsi(rc);
	 //the user should not have to do this hash - decompalgo should be doing this
	 decompCut->setStringHash();//TEST

	 UTIL_DEBUG(m_param.LogDebugLevel, 3,
		    (*m_osLog) << "Found violated cut:";
		    decompCut->print(m_osLog);
		    );

	 newCuts.push_back(decompCut);
      }
   }

   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "generateCuts()", m_param.LogDebugLevel, 2);

  
   return static_cast<int>(newCuts.size());
}
*/
