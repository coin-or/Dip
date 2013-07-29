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
#include "SmallIP_DecompApp2.h"
//===========================================================================//
//parameters
const int LogLevel   = 1;

//===========================================================================//
//---
//--- Version 2: To illustrate use of solveRelaxed function for a user-defined
//---   solver for the relaxation. Rather than using the built-in MILP solver.
//---

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
   
   //---
   //--- set the model core   
   //---  Q'  = { x in R^2 | x satisfies (4.5-10, 16, 17} modelRelax
   //---  Q'' = { x in R^2 | x satisfies (4.11-16       } modelCore
   //---
   setModelCore (&m_modelPart2, "CORE2");

   //---
   //--- set the model relax to NULL, in this case
   //---   solveRelaxed must be defined (for pricing algos)
   //---
   setModelRelax(NULL);

   //---
   //--- load the OSI object to be used in solve relaxed
   //---
   m_osi.messageHandler()->setLogLevel(0);      
   //m_osi.setHintParam(OsiDoReducePrint, true, OsiHintDo);
#ifdef __DECOMP_IP_CBC__   
   m_osi.getModelPtr()->setLogLevel(0);
#endif
   m_osi.loadProblem(*m_modelPart1.getMatrix(),
		     m_modelPart1.getColLB(),
		     m_modelPart1.getColUB(),
		     NULL,
		     m_modelPart1.getRowLB(),
		     m_modelPart1.getRowUB());
   m_osi.setInteger(m_modelPart1.getIntegerVars(), 
		     m_modelPart1.getNumInts());

   UtilPrintFuncEnd(m_osLog, m_classTag, "createModels()", LogLevel, 2);
}

//===========================================================================//
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

//===========================================================================//
DecompSolverStatus 
SmallIP_DecompApp::solveRelaxed(const int          whichBlock,
				const double     * redCostX,
				const double       convexDual,
				DecompVarList    & varList){
   
   //---
   //--- solveRelaxed is a virtual method and can be overriden
   //---   if the user wants to solve the subproblem using their own
   //---   solver, rather than the built-in solver
   //---
   //--- if the user does not define and set the the relaxation model
   //---   then this method must be defined
   //---
   UtilPrintFuncBegin(m_osLog, m_classTag, "solveRelaxed()", LogLevel, 2);

   DecompSolverStatus status = DecompSolStatNoSolution;

   //---
   //--- set the objective function of the subproblem to the current
   //---   reduced cost vector
   //---
   m_osi.setObjective(redCostX);
   
#ifdef __DECOMP_IP_CBC__   
   //---
   //--- because OsiCbc does not keep original column bounds
   //---  we must reset each time
   //--- inside DIP, we avoid this by using Cbc directly, not OsiCbc
   //---
   m_osi.setColLower(0, 0.0);
   m_osi.setColLower(1, 0.0);
   m_osi.setColUpper(0, 6.0);
   m_osi.setColUpper(1, 6.0);   
   m_osi.getModelPtr()->resetModel();
#endif

   //const double * colLB = m_osi.getColLower();
   //const double * colUB = m_osi.getColUpper();
   //for(int i = 0; i < m_osi.getNumCols(); i++)
   //printf("B i:%d lb:%g ub:%g\n", i, colLB[i], colUB[i]);

   //---
   //--- solve with OSI milp solver
   //---
   m_osi.branchAndBound();

   //for(int i = 0; i < m_osi.getNumCols(); i++)
   //printf("A i:%d lb:%g ub:%g\n", i, colLB[i], colUB[i]);

   //---
   //--- check that found optimal
   //---
   assert(!m_osi.isProvenPrimalInfeasible());
   assert(m_osi.isProvenOptimal());
   if(!m_osi.isProvenOptimal())
      return DecompSolStatNoSolution;
   else
      status = DecompSolStatOptimal;
   
   //TODO:
   //this is way too confusing for user to remember they need -alpha!
   //  let framework do that - also setting the block id - framework!

   //---
   //--- create a variable object from the optimal solution
   //---
   int            i;
   int            nOrigCols   = m_osi.getNumCols();
   double         varRedCost  = m_osi.getObjValue() - convexDual;
   double         varOrigCost = 0.0;
   const double * colSolution = m_osi.getColSolution(); 
   for(i = 0; i < m_osi.getNumCols(); i++)
      varOrigCost += colSolution[i] * m_objective[i];
   DecompVar * var = new DecompVar(nOrigCols,
				   colSolution,
				   varRedCost,
				   varOrigCost);   
   varList.push_back(var);      

   UtilPrintFuncEnd(m_osLog, m_classTag, "solveRelaxed()", LogLevel, 2);
   return status;
}
