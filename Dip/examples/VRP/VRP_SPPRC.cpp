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
#include "UtilMacros.h"
//===========================================================================//
#include "VRP_DecompApp.h"
//===========================================================================//
#include "DecompAlgo.h"

//===========================================================================//
void VRP_DecompApp::createModelESPPCC(DecompConstraintSet * model){
   //---
   //--- We create the ESPPCC model as an IP just for the sake of
   //---  debugging. The model is actually over a directed graph, so
   //---  the variables are not actually in the space of the original
   //---  model. So, we will not feed this model to the framework - but
   //---  we will use it for debuggin in the solveRelaxed function.
   //---
   //--- OR, we can feed it to framework to be called like others
   //---  but also have a function to map from this space to original space
   //---  that the user must fill-in.
   //---
   
   //---
   //--- TODO: enforce branching decisions
   //---

   //---
   //--- Elementary Shortest Path with Capacity Constraints
   //---  depot  = 0 (start) and n+1 (end)
   //---  C      = 1...n (customers)
   //---  N      = C union {0, n+1}
   //---  
   //--- Model with complete directed graph (|A| = |N|^2 edges).
   //--- For simplicity of index scheme, just use full graph and
   //--   fix edges to 0 that cannot exist (like self loops, etc).
   //---
   //--- sum{i in C, j in N} d[i]         x[i,j]   <= q
   //--- sum{j in N}                      x[0,j]    = 1
   //--- sum{i in N} x[i,h] - sum{j in N} x[h,j]    = 0, for h in C
   //--- sum{i in N}                      x[i,n+1]  = 1
   //---   x[i,j] = 0, for all {(i,j) in A : i=j or i=n+1 or j=0}
   //---   x[i,j] in {0,1} in A
   //---
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createModelESPPCC()", m_appParam.LogLevel, 2);

   UtilGraphLib & graphLib     = m_vrp.m_graphLib;
   const double   capacity     = graphLib.capacity;
   const int      numCustomers = graphLib.n_vertices - 1;
   const int      numVertices  = numCustomers + 2;
   const int      numCols      = numVertices * numVertices;
   const int      numRows      = numVertices + 1;
   int            colIndex, i, j, h;
   
   model->M = new CoinPackedMatrix(false, 0.0, 0.0);
   if(!model->M)
      throw UtilExceptionMemory("createModelPartMDKP", "MMKP_DecompApp");
   model->M->setDimensions(0, numCols);   
   model->reserve(numRows, numCols);

   //---
   //--- sum{i in C, j in N} d[i] x[i,j]   <= q
   //---
   CoinPackedVector rowCap;
   for(i = 1; i <= numCustomers; i++){
      colIndex = i * numVertices;
      for(j = 0; j < numVertices; j++){
         rowCap.insert(colIndex++, 1.0);
      }
   }
   model->appendRow(rowCap, -DecompInf, capacity);

   //---
   //--- sum{j in N} x[0,j]    = 1   
   //---
   CoinPackedVector rowFlowDepot1;
   colIndex = 0;
   for(j = 0; j < numVertices; j++)
      rowFlowDepot1.insert(colIndex++, 1.0);
   model->appendRow(rowFlowDepot1, 1.0, 1.0);
   
   //---
   //--- sum{i in N} x[i,n+1]  = 1
   //---
   CoinPackedVector rowFlowDepot2;
   colIndex = (numVertices-1) * numVertices;
   for(j = 0; j < numVertices; j++)
      rowFlowDepot2.insert(colIndex++, 1.0);
   model->appendRow(rowFlowDepot2, 1.0, 1.0);

   //---
   //--- sum{i in N} x[i,h] - sum{j in N} x[h,j] = 0, for h in C
   //---
   for(h = 1; h <= numCustomers; h++){
      CoinPackedVector row;
      for(j = 0; j < numVertices; j++){
	if(h==j)
	  continue;
         row.insert(diGraphIndex(j,h,numVertices),  1.0);
         row.insert(diGraphIndex(h,j,numVertices), -1.0);
      }
      model->appendRow(row, 0.0, 0.0);
   }

   //---   
   //--- set the col upper and lower bounds
   //---
   UtilFillN(model->colLB, numCols,  0.0);
   UtilFillN(model->colUB, numCols,  1.0);

   //---
   //--- x[i,j] = 0, for all {(i,j) in A : i=j or i=n+1 or j=0}
   //---
   colIndex = 0;
   for(i = 0; i < numVertices; i++){      
      for(j = 0; j < numVertices; j++){
         if(i == j                ||
            i == (numCustomers+1) ||
            j == 0){
            model->colUB[colIndex] = 0.0;
         }
         colIndex++;
      }
   }
   
   //---
   //--- set the indices of the integer variables of model
   //---
   UtilIotaN(model->integerVars, numCols, 0);

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createModelESPPCC()", m_appParam.LogLevel, 2);
}
