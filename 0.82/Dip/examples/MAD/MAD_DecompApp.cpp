//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, Lehigh University                                 //
//                                                                           //
// Copyright (C) 2002-2007, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

#include "DecompCutOsi.h"
#include "MAD_DecompApp.h"

#include "CoinMpsIO.hpp"
#include "CoinModel.hpp"

#define USE_QUALEX

//TODO: should be easy to write heuristics to get ip feasible points
//  from subproblem points - since capacity is only extra constraint
//  this should be point out as another advantage of Decomp methods... 

// --------------------------------------------------------------------- //
int MAD_DecompApp::isOrtho(const int    * rowInd1,
                           const int    * rowInd2,
                           const int      rowLen1,
                           const int      rowLen2){

   int i1, i2;
   
   /*
     ---
     --- CAREFUL: this is only correct if the input is sorted
     ---
   */
  
   i1 = 0;
   i2 = 0;

   while((i1 < rowLen1) && (i2 < rowLen2)){
      /*
	---
	--- if there is a row in common, then the rows are not orthogonal
	---
      */
      if(rowInd1[i1] == rowInd2[i2])
	 return 0;
    
      /*
	---
	--- advance the row pointer for the lesser valued row index
	---
      */
      if(rowInd1[i1] < rowInd2[i2])
	 ++i1;
      else
	 ++i2;
   }
   return 1;  
}

// --------------------------------------------------------------------- //
void MAD_DecompApp::initializeApp(UtilParameters & utilParam) {

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "initializeApp()", m_param.LogDebugLevel, 2);

  
   //---
   //--- get application parameters
   //---
   m_appParam.getSettings(utilParam);
   m_appParam.dumpSettings(m_osLog); //use message handler

   //---
   //--- read instance from lp file (from MADLIB)
   //---   http://elib.zib.de/pub/mp-testdata/madlib/index.html
   //---
   string lpFile  = m_appParam.DataDir    + UtilDirSlash();
   lpFile        += m_appParam.Instance;
   if(m_appParam.DataSubDir == "miplib"){
      lpFile     += ".p.lp";
   } else if(m_appParam.DataSubDir == "netlib"){
      lpFile     += ".ob4";
   }   
   m_instance.readLp(lpFile.c_str());
   
   m_nOrigRows = m_instance.getNumRows();   
   m_beta      = m_appParam.NumBlocks;

   //---
   //--- read best known lb/ub 
   //---
   string bestKnownFile  = m_appParam.DataDir + UtilDirSlash();
   bestKnownFile        += "madlib." + m_appParam.DataSubDir + ".opt";
   {
      ifstream is;
      string   instanceName;
      double   bestLB, bestUB;
      UtilOpenFile(is, bestKnownFile);
      while(!is.eof()){
         //---
         //--- these are the number of rows in the border (less is better)
         //---  
         is >> instanceName >> bestLB >> bestUB;
         //printf("Instance = %15s bestLB = %6g bestUB = %6g\n",
         //     instanceName.c_str(), bestLB, bestUB);
         
         instanceName = UtilStrTrim(instanceName);
         if(instanceName == m_appParam.Instance){
            //---
            //--- the paper solves z = max sum x,             
            //---    where x is an assignment to a block
            //---    so, nBorder = nRows - z
            //---    or, z       = nRows - nBorder
            //---
            //--- but we only do min, so we solve -z = min sum (-x)
            //----   so, -z      = nBorder - nRows
            //---    
            m_bestKnownLB = bestLB - m_nOrigRows;
            m_bestKnownUB = bestUB - m_nOrigRows;
            break;
         }
      }
   }


   //---
   //--- set capacity based on MADLIB study (if it is not set):
   //---  http://elib.zib.de/pub/mp-testdata/madlib/index.en.html
   //---
   if(m_appParam.Capacity != -1){
      m_kappa = m_appParam.Capacity;
   }
   else{
      m_kappa 
	 = static_cast<int>(ceil(static_cast<double>(m_nOrigRows)/m_beta));
      if(m_appParam.DataSubDir == "netlib" ||
	 m_appParam.DataSubDir == "equipart"){
	 m_kappa 
	    = static_cast<int>(ceil(static_cast<double>(m_nOrigRows)/m_beta));
      } else if(m_appParam.DataSubDir == "miplib" ||
		m_appParam.DataSubDir == "miplibT"){
	 m_kappa = static_cast<int>(ceil( 1.05 * m_nOrigRows / m_beta) ); 
      } else if(m_appParam.DataSubDir == "steiner"){
	 if(m_appParam.Instance[0] == 'g'){
	    m_kappa = 30;
	 } else if (m_appParam.Instance[0] == 'd'){
	    m_kappa = 50;
	 }
      }
   }

   UTIL_DEBUG(m_param.LogDebugLevel, 1,
	      (*m_osLog) 
              << "Instance = " << m_appParam.Instance << endl
              << "  nRows  = " << m_nOrigRows         << endl
              << "  bestLB = " << m_bestKnownLB       << endl
              << "  bestUB = " << m_bestKnownUB       << endl
              << "  Beta   = " << m_beta              << endl
              << "  Kappa  = " << m_kappa             << endl;
	      );
   
   int n_cols = m_nOrigRows * m_beta;
   
#ifdef __MAD_USE_CLIQUER__
   m_cliquer = new MAD_Cliquer(n_cols);
#endif 
#ifdef __MAD_USE_QUALEX__
   m_qualex  = new MAD_Qualex(n_cols);
#endif

   //TODO: not the best name - maybe column intersection graph?
   m_conflictGraph = graph_new(m_nOrigRows);
   m_auxMemPool.allocateMemory(m_nOrigRows, n_cols, m_beta);

   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "initializeApp()", m_param.LogDebugLevel, 2);
}

#if 1
// --------------------------------------------------------------------- //
void
MAD_DecompApp::APPcreateModel(double                        *& objCoeff,
                              map<int, DecompConstraintSet*> & modelCore,
                              map<int, vector<DecompConstraintSet* > > & modelRelax) {
   
   //---
   //--- createModel is a pure virtual method of DecompApp and must 
   //--- be derived by the application class to define the partitioning
   //--- of constraints into [A,b] = [A',b'] union [A'', b'']
   //---
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "APPcreateModel()", m_param.LogDebugLevel, 2);

   //---
   //--- x[i,b] in {0,1}, 1 = row i is assigned to block b
   //---   i in M = {1, ..., m} 
   //---   b in B = {1, ..., beta}
   //---   beta <= m (possibly an input = number of processors available)
   //---
   //--- max sum{i in M, b in B}  x[i,b]  <==>
   //--- min sum{i in M, b in B} -x[i,b]
   //--- s.t. 
   //--- (1) sum{b in B} x[i,b] <= 1, for i in M
   //--- (2) sum{i in M} x[i,b] <= k, for b in B
   //--- (3) x[i,b] + x[j,b']   <= 1, for b,b' in B, b != b',
   //---                                  i,j  in M, i != j, such that
   //---                                  a[i,k] != 0 != a[j,k], for some k
   //--- (4) x[i,b] in {0,1}
   //---
   //--- (1),(3), and (4) forms a (big) clique problem [modelRelax]
   //--- (2)              forms the core               [modelCore ]
   //---
   int n_cols  = m_nOrigRows * m_beta;

   //---
   //--- open memory for the objective coefficients of modelCore
   //---  who is responsible to open this memory? MAD_DecompApp
   //---  who is responsible to free this memory? DecompAlgo
   //---
   objCoeff    = new double[n_cols];
   CoinFillN(objCoeff, n_cols, -1.0);

   //---
   //--- set the constraint matrix of modelCore
   //---  who is responsible to open this memory? MAD_DecompApp
   //---  who is responsible to free this memory? DecompAlgo 
   //---
   //--- (2) sum{i in M} x[i,b] <= k, for b in B
   //---
   DecompConstraintSet * modelCoreCl  = new DecompConstraintSet();
   modelCoreCl->M = new CoinPackedMatrix(false, 0.0, 0.0);
   CoinAssertHint(modelCoreCl->M, "Error: Out of Memory");
   modelCoreCl->M->setDimensions(0, n_cols);
   modelCoreCl->M->reserve(m_beta, m_beta * m_nOrigRows);

   //TODO: for speed - do in blocks
   int b1, b2, i, j, b;
   for(b = 0; b < m_beta; b++){      
      CoinPackedVector row;
      for(i = 0; i < m_nOrigRows; i++){
	 row.insert(xIndex(i,b), 1.0);
      }
      modelCoreCl->M->appendRow(row);  
   }
   
   //---
   //--- set the row upper and lower bounds
   //--- set the col upper and lower bounds
   //---
   UtilFillN(modelCoreCl->rowLB,    m_beta,  -DecompInf);
   UtilFillN(modelCoreCl->rowUB,    m_beta,  static_cast<double>(m_kappa));
   UtilFillN(modelCoreCl->colLB,    n_cols,  0.0);
   UtilFillN(modelCoreCl->colUB,    n_cols,  1.0);



   //---
   //--- set the constraint matrix of modelRelax
   //---  who is responsible to open this memory? MAD_DecompApp
   //---  who is responsible to free this memory? DecompAlgo 
   //---
   //--- (1) sum{b in B} x[i,b] <= 1, for i in M
   //--- (3) x[i,b] + x[j,b']   <= 1, for b,b' in B, b != b',
   //---                                  i,j  in M, i != j, such that
   //---                                  a[i,k] != 0 != a[j,k], for some k
   //--- (4) x[i,b] in {0,1}
   //---
   DecompConstraintSet * modelRelaxCl  = new DecompConstraintSet();
   modelRelaxCl->M = new CoinPackedMatrix(false, 0.0, 0.0);
   CoinAssertHint(modelRelaxCl->M, "Error: Out of Memory");
   modelRelaxCl->M->setDimensions(0, n_cols);

   //how to efficiently construct the complement of conflict graph???
   int n_vertices = m_nOrigRows * m_beta;;
#ifdef __MAD_USE_CLIQUER__
   //construct complete graph -- then delete -- UGH!! TEMP   
   for(i = 0; i < n_vertices ; i++){
      for(j = (i+1); j < n_vertices; j++){
	 m_cliquer->addEdge(m_cliquer->m_g, i, j);
      }
   }
#endif
#ifdef __MAD_USE_QUALEX__
   for(i = 0; i < n_vertices ; i++){
      for(j = (i+1); j < n_vertices; j++){
	 m_qualex->addEdge(i, j);
      }
   }
#endif
   UTIL_DEBUG(m_param.LogDebugLevel, 4,
	      m_cliquer->printGraph(m_cliquer->m_g);
	      );

   //---
   //--- for speed, use CoinModel
   //---


   //TODO: consider not adding (1) to subproblem,
   //      at least for solving exact version....
   
   //---
   //--- (1) sum{b in B} x[i,b] <= 1, for i in M
   //--- 
   CoinModel cModel;  
   for(i = 0; i < m_nOrigRows; i++){
      CoinPackedVector row;
      for(b = 0; b < m_beta; b++){
	 row.insert(xIndex(i, b), 1.0);
      }
      cModel.addRow(row.getNumElements(), 
		    row.getIndices(), row.getElements());
      
#ifdef __MAD_USE_CLIQUER__
      for(b1 = 0; b1 < m_beta; b1++){
	 for(b2 = 0; b2 < m_beta; b2++){
	    if(b1 == b2)
	       continue;
	    m_cliquer->delEdge(m_cliquer->m_g,
			       xIndex(i, b1), xIndex(i, b2));
	 }
      }
#endif      
#ifdef __MAD_USE_QUALEX__
      for(b1 = 0; b1 < m_beta; b1++){
	 for(b2 = 0; b2 < m_beta; b2++){
	    if(b1 == b2)
	       continue;
	    m_qualex->removeEdge(xIndex(i, b1), xIndex(i, b2));
	 }
      }
#endif
   }
   UTIL_DEBUG(m_param.LogDebugLevel, 4,
	      m_cliquer->printGraph(m_cliquer->m_g);
	      );

   //---
   //--- In order to do this part fast, construct a row-ordered 
   //--- matrix which is sorted on column indices.
   //--- (3) x[i,b] + x[j,b']   <= 1, for b,b' in B, b != b',
   //---                                  i,j  in M, i != j, such that
   //---                                  a[i,k] != 0 != a[j,k], for some k
   //---
   CoinPackedMatrix Morder(*m_instance.getMatrixByRow());  
   Morder.orderMatrix();

   for(i = 0; i < m_nOrigRows; i++){
      for(j = (i+1); j < m_nOrigRows; j++){
	 CoinShallowPackedVector rowI = Morder.getVector(i);
	 CoinShallowPackedVector rowJ = Morder.getVector(j);
	 if(!isOrtho(rowI.getIndices(),
		     rowJ.getIndices(),
		     rowI.getNumElements(),
		     rowJ.getNumElements())){
	    //---
	    //--- add a constraint to the model for every pair of blocks
	    //---
	    for(b1 = 0; b1 < m_beta; b1++){
	       for(b2 = 0; b2 < m_beta; b2++){
		  if(b1 == b2)
		     continue;      
                  
		  CoinPackedVector row;
		  row.insert(xIndex(i, b1), 1.0);
		  row.insert(xIndex(j, b2), 1.0);
                  
		  cModel.addRow(row.getNumElements(),
				row.getIndices(), row.getElements());

#ifdef __MAD_USE_CLIQUER__
		  m_cliquer->delEdge(m_cliquer->m_g,
				     xIndex(i, b1), xIndex(j, b2));
		  m_cliquer->delEdge(m_cliquer->m_g,
				     xIndex(j, b1), xIndex(i, b2));
#endif

#ifdef __MAD_USE_QUALEX__
                  m_qualex->removeEdge(xIndex(i, b1), xIndex(j, b2));
                  m_qualex->removeEdge(xIndex(j, b1), xIndex(i, b2));
#endif                  
	       }
	    }
	    //if doing this way, then need to add an interface
	    //class for this graph object - fudge - do we use
	    //cliquer just for the graph object!?
	    //or should we have a graph object independent of all of
	    //of qualex/cliquer so we don't force user - boost?
	    //but then user must have boost/graph
	    GRAPH_ADD_EDGE(m_conflictGraph, i, j);
	 }
      }
   }
   UTIL_DEBUG(m_param.LogDebugLevel, 4,
	      m_cliquer->printGraph(m_cliquer->m_g);
	      );
   cModel.createPackedMatrix(*modelRelaxCl->M, cModel.associatedArray());

   //---
   //--- set the row upper and lower bounds
   //--- set the col upper and lower bounds
   //---
   UtilFillN(modelRelaxCl->rowLB, modelRelaxCl->M->getNumRows(),  -DecompInf);
   UtilFillN(modelRelaxCl->rowUB, modelRelaxCl->M->getNumRows(),  1.0);
   UtilFillN(modelRelaxCl->colLB, n_cols,                         0.0);
   UtilFillN(modelRelaxCl->colUB, n_cols,                         1.0);

   //---
   //--- set the indices of the integer variables of modelRelax
   //---
   UtilIotaN(modelRelaxCl->integerVars, n_cols, 0);
   
   //---
   //--- push core and relaxed into application object
   //---
   vector< DecompConstraintSet* > modelRelaxV;   
   modelRelaxV.push_back(modelRelaxCl);
   modelCore.insert(make_pair(MODEL_CLIQUE, modelCoreCl));
   modelRelax.insert(make_pair(MODEL_CLIQUE, modelRelaxV));


   UTIL_DEBUG(m_param.LogDebugLevel, 3,
	      (*m_osLog) << "\nCONFLICT GRAPH:\n";
	      graph_print(m_conflictGraph);
	      );
   
   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "APPcreateModel()", m_param.LogDebugLevel, 2);

}
#else

// --------------------------------------------------------------------- //
void
MAD_DecompApp::APPcreateModel(double                        *& objCoeff,
                              map<int, DecompConstraintSet*> & modelCore,
                              map<int, DecompConstraintSet*> & modelRelax) {
   
   //---
   //--- createModel is a pure virtual method of DecompApp and must 
   //--- be derived by the application class to define the partitioning
   //--- of constraints into [A,b] = [A',b'] union [A'', b'']
   //---

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "APPcreateModel()", m_param.LogDebugLevel, 2);
   
   //--- x[i,b] in {0,1}, 1 = row i is assigned to block b
   //---   i in M = {1, ..., m} 
   //---   b in B = {1, ..., beta}
   //---   beta <= m (possibly an input = number of processors available)
   //---
   //--- max sum{i in M, b in B}  x[i,b]  <==>
   //--- min sum{i in M, b in B} -x[i,b]
   //--- s.t. 
   //--- (1) sum{b in B} x[i,b] <= 1, for i in M
   //--- (2) sum{i in M} x[i,b] <= k, for b in B
   //--- (3) x[i,b] + x[j,b']   <= 1, for b,b' in B, b != b',
   //---                                  i,j  in M, i != j, such that
   //---                                  a[i,k] != 0 != a[j,k], for some k
   //--- (4) x[i,b] in {0,1}
   //---
   //--- (3),(4) form the clique problem [modelRelax]
   //--- (1),(2) form the core           [modelCore ]
   //---
   int n_cols  = m_nOrigRows * m_beta;

   //---
   //--- open memory for the objective coefficients of modelCore
   //---  who is responsible to open this memory? MAD_DecompApp
   //---  who is responsible to free this memory? DecompAlgo
   //---
   objCoeff    = new double[n_cols];
   CoinFillN(objCoeff, n_cols, -1.0);

   //---
   //--- set the constraint matrix of modelCore
   //---  who is responsible to open this memory? MAD_DecompApp
   //---  who is responsible to free this memory? DecompAlgo 
   //---
   //--- (1) sum{b in B} x[i,b] <= 1, for i in M
   //--- (2) sum{i in M} x[i,b] <= k, for b in B
   //---
   DecompConstraintSet * modelCoreCl  = new DecompConstraintSet();
   modelCoreCl->M = new CoinPackedMatrix(false, 0.0, 0.0);
   CoinAssertHint(modelCoreCl->M, "Error: Out of Memory");
   modelCoreCl->M->setDimensions(0, n_cols);
   modelCoreCl->M->reserve(m_nOrigRows + m_beta,
			   2 * m_beta * m_nOrigRows);

   //TODO: for speed - do in blocks
   int b1, b2, i, j, b;

   //TODO: make this optional to add to subproblem or not
   //---
   //--- (1) sum{b in B} x[i,b] <= 1, for i in M
   //---
   for(i = 0; i < m_nOrigRows; i++){
      CoinPackedVector row;
      for(b = 0; b < m_beta; b++){
	 row.insert(xIndex(i, b), 1.0);
      }
      modelCoreCl->M->appendRow(row);
   }
   
   //---
   //--- (2) sum{i in M} x[i,b] <= k, for b in B
   //---
   for(b = 0; b < m_beta; b++){      
      CoinPackedVector row;
      for(i = 0; i < m_nOrigRows; i++){
	 row.insert(xIndex(i,b), 1.0);
      }
      modelCoreCl->M->appendRow(row);  
   }
   
   //---
   //--- set the row upper and lower bounds
   //--- set the col upper and lower bounds
   //---
   double kappa = static_cast<double>(m_kappa);
   UtilFillN(modelCoreCl->rowLB,    m_nOrigRows,  -DecompInf);
   UtilFillN(modelCoreCl->rowUB,    m_nOrigRows,  1.0);
   UtilFillN(modelCoreCl->rowLB,    m_beta,       -DecompInf);
   UtilFillN(modelCoreCl->rowUB,    m_beta,       kappa);
   UtilFillN(modelCoreCl->colLB,    n_cols,       0.0);
   UtilFillN(modelCoreCl->colUB,    n_cols,       1.0);

   //---
   //--- set the constraint matrix of modelRelax
   //---  who is responsible to open this memory? MAD_DecompApp
   //---  who is responsible to free this memory? DecompAlgo 
   //---
   //--- (3) x[i,b] + x[j,b']   <= 1, for b,b' in B, b != b',
   //---                                  i,j  in M, i != j, such that
   //---                                  a[i,k] != 0 != a[j,k], for some k
   //--- (4) x[i,b] in {0,1}
   //---
   DecompConstraintSet * modelRelaxCl  = new DecompConstraintSet();
   modelRelaxCl->M = new CoinPackedMatrix(false, 0.0, 0.0);
   CoinAssertHint(modelRelaxCl->M, "Error: Out of Memory");
   modelRelaxCl->M->setDimensions(0, n_cols);

   //how to efficiently construct the complement of conflict graph???
   int n_vertices = m_nOrigRows * m_beta;;
#ifdef __MAD_USE_CLIQUER__
   //construct complete graph -- then delete -- UGH!! TEMP   
   for(i = 0; i < n_vertices ; i++){
      for(j = (i+1); j < n_vertices; j++){
	 m_cliquer->addEdge(m_cliquer->m_g, i, j);
      }
   }
   UTIL_DEBUG(m_param.LogDebugLevel, 4,
              (*m_osLog) << "m_cliquer complete graph:\n";
	      m_cliquer->printGraph(m_cliquer->m_g);
	      );
#endif
#ifdef __MAD_USE_QUALEX__
   for(i = 0; i < n_vertices ; i++){
      for(j = (i+1); j < n_vertices; j++){
	 m_qualex->addEdge(i, j);
      }
   }
   UTIL_DEBUG(m_param.LogDebugLevel, 4,
              (*m_osLog) << "m_qualex complete graph:\n";
	      m_qualex->printGraph(m_qualex->m_graphOrig);
	      );   
#endif

   //---
   //--- for speed, use CoinModel
   //---
   CoinModel cModel;  

   //---
   //--- In order to do this part fast, construct a row-ordered 
   //--- matrix which is sorted on column indices.
   //--- (3) x[i,b] + x[j,b']   <= 1, for b,b' in B, b != b',
   //---                                  i,j  in M, i != j, such that
   //---                                  a[i,k] != 0 != a[j,k], for some k
   //---
   //--- Example:
   //---  i:0  j:12, B={0,1,2}
   //---  i:12 i:0 , B={0,1,2}
   //---     no need to duplicate, can do just i < j
   //---    x[12,0] + x[0,1] <= 1
   //---    x[12,0] + x[0,2] <= 1
   //---    x[12,1] + x[0,0] <= 1
   //---    x[12,1] + x[0,2] <= 1
   //---    x[12,2] + x[0,0] <= 1
   //---    x[12,2] + x[0,1] <= 1
   //---

   CoinPackedMatrix Morder(*m_instance.getMatrixByRow());  
   Morder.orderMatrix();

   for(i = 0; i < m_nOrigRows; i++){
      for(j = (i+1); j < m_nOrigRows; j++){
	 CoinShallowPackedVector rowI = Morder.getVector(i);
	 CoinShallowPackedVector rowJ = Morder.getVector(j);
	 if(!isOrtho(rowI.getIndices(),
		     rowJ.getIndices(),
		     rowI.getNumElements(),
		     rowJ.getNumElements())){
            
	    //---
	    //--- add a constraint to the model for every pair of blocks
	    //---
	    for(b1 = 0; b1 < m_beta; b1++){
	       for(b2 = 0; b2 < m_beta; b2++){
		  if(b1 == b2)
		     continue;      
                  
		  CoinPackedVector row;
		  row.insert(xIndex(i, b1), 1.0);
		  row.insert(xIndex(j, b2), 1.0);

                  UTIL_DEBUG(m_param.LogDebugLevel, 5,
                             (*m_osLog) 
                             << "Row " << cModel.numberRows()
                             << " (" << i << "," << b1 << ")->" 
                             << xIndex(i, b1)
                             << " (" << j << "," << b2 << ")->" 
                             << xIndex(j, b2) << endl;
                             );
		  
                  cModel.addRow(row.getNumElements(),
				row.getIndices(), row.getElements());
                  
#ifdef __MAD_USE_CLIQUER__                  
		  m_cliquer->delEdge(m_cliquer->m_g,
				     xIndex(i, b1), xIndex(j, b2));
                  m_cliquer->delEdge(m_cliquer->m_g,
                                     xIndex(j, b1), xIndex(i, b2));
#endif
                  
#ifdef __MAD_USE_QUALEX__
                  m_qualex->removeEdge(xIndex(i, b1), xIndex(j, b2));
                  m_qualex->removeEdge(xIndex(j, b1), xIndex(i, b2));
#endif
                  
                  
                  
		  //if doing this way, then need to add an interface
		  //class for this graph object - fudge - do we use
		  //cliquer just for the graph object!?
		  //or should we have a graph object independent of all of
		  //of qualex/cliquer so we don't force user - boost?
		  //but then user must have boost/graph
		  GRAPH_ADD_EDGE(m_conflictGraph, i, j);
                  
	       }
	    }
	 }
      }
   }
#ifdef __MAD_USE_CLIQUER__
   UTIL_DEBUG(m_param.LogDebugLevel, 4,
              (*m_osLog) << "m_cliquer graph:\n";
	      m_cliquer->printGraph(m_cliquer->m_g);
	      );
#endif
#ifdef __MAD_USE_QUALEX__
   UTIL_DEBUG(m_param.LogDebugLevel, 4,
              (*m_osLog) << "m_qualex complete graph:\n";
	      m_qualex->printGraph(m_qualex->m_graphOrig);
	      );   
#endif
   
   cModel.createPackedMatrix(*modelRelaxCl->M, cModel.associatedArray());

   //---
   //--- set the row upper and lower bounds
   //--- set the col upper and lower bounds
   //---
   UtilFillN(modelRelaxCl->rowLB, modelRelaxCl->M->getNumRows(),  -DecompInf);
   UtilFillN(modelRelaxCl->rowUB, modelRelaxCl->M->getNumRows(),  1.0);
   UtilFillN(modelRelaxCl->colLB, n_cols,                         0.0);
   UtilFillN(modelRelaxCl->colUB, n_cols,                         1.0);

   //---
   //--- set the indices of the integer variables of modelRelax
   //---
   UtilIotaN(modelRelaxCl->integerVars, n_cols, 0);
   
   //---
   //--- push core and relaxed into application object
   //---
   modelCore.insert(make_pair(MODEL_CLIQUE, modelCoreCl));
   modelRelax.insert(make_pair(MODEL_CLIQUE, modelRelaxCl));


   UTIL_DEBUG(m_param.LogDebugLevel, 4,
	      (*m_osLog) << "\nCONFLICT GRAPH:\n";
	      graph_print(m_conflictGraph);
	      );


   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "APPcreateModel()", m_param.LogDebugLevel, 2);
}
#endif

//--------------------------------------------------------------------- //
//too many args - fix this
DecompStatus MAD_DecompApp::APPsolveRelaxed(const int             whichModel,
					    const double        * redCostX,
					    const double        * origCost,
					    const double          alpha,
					    const int             n_origCols,
					    const bool            checkRC,
					    const bool            checkDup,
					    bool                & isExact,
					    OsiSolverInterface  * m_subprobSI,
					    list<DecompVar*>    & vars){
   
   
   //TODO: are we using (1) or just (3) in clique definition?
   //---
   //--- This can be considered an independent set problem on 
   //--- the conflict graph defined by these constraints.
   //---
   //--- (1) sum{b in B} x[i,b] <= 1, for i in M
   //--- (3) x[i,b] + x[j,b']   <= 1, for b,b' in B, b != b',
   //---                                  i,j  in M, i != j, such that
   //---                                  a[i,k] != 0 != a[j,k], for some k
   //--- (4) x[i,b] in {0,1}
   //---
   //--- Note: Finding maximal cliques in a graph is equivalent to 
   //--- finding maximal independent sets in the complement of that graph.
   //---
   //--- We are using cliquer, so :
   //----  (1) we need to look at complement of conflict graph
   //---   (2) cliquer assumes maximization, but we want "least" reduced
   //---       cost - so, we need to flip the reduced cost.
   //---   (3) cliquer expects integral weights, so we need to scale
   //---   (4) cliquer only accepts positive vertex weights
   //---
   
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "APPsolveRelaxed()", m_param.LogDebugLevel, 2);
   
   
   CoinTimer thisTimer;
   thisTimer.restart();

   DecompStatus status = STAT_FEASIBLE;   
   isExact             = false;

   
   //---
   //--- try to generate a column with small reduced cost
   //--- that is feasible to original problem using primal heuristic
   //---
   int i;
   int nGreedyPts     = 0;
   int nGreedyPtsKept = 0;
   vector<GreedyPoint> greedyPoints;
   nGreedyPts = heuristicGreedy(greedyPoints,
				INCREASING, redCostX, origCost, redCostX);   
   for(i = 0; i < nGreedyPts; i++){
      GreedyPoint & gp = greedyPoints[i];
      
      //---
      //--- create a DecompVar from the greedySol
      //---
      UTIL_DEBUG(m_param.LogDebugLevel, 3,
		 (*m_osLog) 
		 << "Greedy DecompVar origCost = " << gp.solValueOrigCost
		 << " redCost = " << gp.solValueRedCost + alpha << endl;
		 printOriginalSolution(n_origCols, gp.solution, &cout);
		 );
      //TODO: might be more efficient if return sparse... greedySol 

      //TODO: -alpha???? check sign see what changed for MMKP

      if(!checkRC || ((gp.solValueRedCost + alpha) < -DecompEpsilon)){
	 vars.push_back(new DecompVar(n_origCols,
				      gp.solution,
				      gp.solValueRedCost + alpha,
				      gp.solValueOrigCost));
	 nGreedyPtsKept++;
      }			   
   }
     
   return STAT_FEASIBLE;




















   //stupid object if pass in graph each time
   int  n_verts = m_cliquer->getNumVertices(m_cliquer->m_g);

   UTIL_DEBUG(m_param.LogDebugLevel, 4,
	      printf("\nGRAPH m_g:\n");
	      m_cliquer->printGraph(m_cliquer->m_g);
	      );



   //if use IP solver over IS and min uA>=0, will pick all 0s
   //there is no incentive

   int b;
   //---
   //--- the goal is to find a subgraph with negative reduced cost
   //--- therefore, it will never make sense to include a node with
   //--- non-negative cost
   //---
   //how do we compile without cliquer?
   graph_t * m_gStar = m_cliquer->graphNew(n_verts);
   m_cliquer->copyGraphNonPos(m_gStar,
			      m_cliquer->m_g, redCostX); 
   printf("timer after copy graph = %12.10f\n",
	  thisTimer.timeElapsed());

   UTIL_DEBUG(m_param.LogDebugLevel, 4,
	      printf("\nGRAPH m_gStar:\n");
	      m_cliquer->printGraph(m_gStar);

	      for(i = 0; i < m_nOrigRows; i++){
		 for(b = 0; b < m_beta; b++){
		    cout << "x[ " << i << "," << b << " -> " << xIndex(i,b)
			 << " ] : " << redCostX[xIndex(i,b)] << endl;
		 }
	      }
	      );

   if(graph_edge_count(m_gStar) <= 0){
      m_cliquer->graphFree(m_gStar);
      return STAT_FEASIBLE; //think
   }


   
   
#ifdef USE_QUALEX
   //negative weights are bad -- qualex takes sqrt of weights
   //negative weights ok? or do we need an offset?
   //print to dimacs format
   //list<int> apInd;
   //m_cliquer->printGraphDimacs(m_gStar);

   double * redCostXNegOffset = m_auxMemPool.dblArrNCoreCols;
   memcpy(redCostXNegOffset, redCostX, n_origCols * sizeof(double));
   
   //---
   //--- flip reduced costs (max c == min -c)
   //---
   UtilNegateArr(n_origCols, redCostXNegOffset);

   //---
   //--- add a constant so that all vertex weights are positive, inc alpha
   //---
   double offset = 0.0;
   double minrc  = *min_element(redCostXNegOffset,
				redCostXNegOffset + n_origCols);
   if(minrc <= 0){
      offset = -minrc + 1;
      UtilAddOffsetArr(n_origCols, offset, redCostXNegOffset);
   }

   //---
   //--- if for initial vars, perturb the costs slightly,
   //---   finding max clique with all weights equal is harder
   //---
   if(!checkRC){
      //make srand a Util func - in case change random function?       
      UtilPerturbCost(n_origCols, n_origCols, 0.0, 1.0, redCostXNegOffset);
   }
  
   
   //m_cliquer->printWeightDimacs(m_gStar->n, redCostXNegOffset);   
   //m_cliquer->printGraphDimacs(m_cliquer->m_g);   
   //m_cliquer->printWeightDimacs(m_cliquer->m_g->n, redCostX);   



   UTIL_DEBUG(m_param.LogDebugLevel, 4,
	      for(i = 0; i < m_nOrigRows; i++){
		 for(b = 0; b < m_beta; b++){
		    cout << "x[ " << i << "," << b << " -> " << xIndex(i,b)
			 << " ] : " << redCostXNegOffset[xIndex(i,b)] << endl;
		 }
	      }
	      );



   m_qualex->setUpForSolve(redCostXNegOffset);

   //---
   //--- we need to remove non-negative vertices from the graph
   //--- data structure used by qualex - we are trying to find cliques
   //--- so the easiest thing to do, is to remove all edges to each
   //--- non-negative vertex
   //---
   Graph * m_graph = m_qualex->m_graph;
   int j;
#if 0
   for(i = 0; i < m_graph->n; i++){
      printf("\n%d: [%g] ", i, m_graph->weights[i]);
      vector<bool_vector>::iterator bi = m_graph->mates.begin() + i;
      for(j = 0; j < m_graph->n; j++){
	 if(bi->at(j))
	    printf(" %d", j);               
      }
   }
#endif

   m_qualex->removeNonNegVertices(redCostX);

#if 0
   for(i = 0; i < m_graph->n; i++){
      printf("\n%d: [%g] ", i, m_graph->weights[i]);
      vector<bool_vector>::iterator bi = m_graph->mates.begin() + i;
      for(j = 0; j < m_graph->n; j++){
	 if(bi->at(j))
	    printf(" %d", j);               
      }
   }
#endif

   printf("timer before greedy call = %12.10f\n",
	  thisTimer.timeElapsed());
   m_qualex->findMaxIndSetGreedy(redCostXNegOffset);
   printf("timer after greedy call = %12.10f\n",
	  thisTimer.timeElapsed());
   

   double                 varRedCost;
   double                 varOrigCost;
   varRedCost  = 0.0;
   list<int>::iterator it;
   for(it  = m_qualex->m_clique.begin(); 
       it != m_qualex->m_clique.end(); it++){         
      varRedCost  += redCostX[*it];
      UTIL_DEBUG(m_param.LogDebugLevel, 4,
		 printf("\nGreedy ind: %d redCostX: %g varRedCost: %g",
			*it, redCostX[*it], varRedCost);
		 );      
   }
   UTIL_DEBUG(m_param.LogDebugLevel, 4,                
	      printf("\nvarRC = %g, alpha = %g", varRedCost, alpha);
	      );
   if(!checkRC || ((varRedCost + alpha) < -DecompEpsilon)){
      //then no need to get into qualex
      UTIL_DEBUG(m_param.LogDebugLevel, 4,
		 (*m_osLog) 
		 << "\nGreedy is good enough rc = " 
		 << varRedCost+alpha << endl;
		 );
   }
   else{
      //status = STAT_UNKNOWN;
      //m_qualex->findMaxIndSetQualexMS();   //keeps failing use cbc
   }


   //---
   //--- store the solution as a DecompVar and push into list
   //---

   //if(status != STAT_UNKNOWN){
   if(!checkRC || ((varRedCost + alpha) < -DecompEpsilon)){
      
      CoinAssert(m_qualex->m_clique.size() > 0);
      varRedCost  = 0.0;
      varOrigCost = 0.0;
      vector<int>    apInd(m_qualex->m_clique.begin(),
                           m_qualex->m_clique.end());
      vector<double> apEls(m_qualex->m_clique.size(), 1.0);
      
      for(it  = m_qualex->m_clique.begin(); 
          it != m_qualex->m_clique.end(); it++){         
         varRedCost  += redCostX[*it];
         varOrigCost += origCost[*it];
         
         UTIL_DEBUG(m_param.LogDebugLevel, 4,
                    printf("\nind: %d redCostX: %g varRedCost: %g",
                           *it, redCostX[*it], varRedCost);
		 );      
      }
      UTIL_DEBUG(m_param.LogDebugLevel, 4,                
                 printf("\nvarRC = %g, alpha = %g", varRedCost, alpha);
	      );
      if(!checkRC || ((varRedCost + alpha) < -DecompEpsilon)){
         UTIL_DEBUG(m_param.LogDebugLevel, 4,
                    printf("\nPUSH var with RC = %g", varRedCost + alpha);
                    );
	 DecompVar * var = new DecompVar(apInd, apEls, 
					 varRedCost + alpha, varOrigCost); 
         vars.push_back(var);

	 //double tmpSol[100000];
	 //var->fillDenseArr(n_origCols, tmpSol); 
	 //printOriginalSolution(n_origCols, tmpSol, &cout);
      }
   } 
   printf("timer after setupvars = %12.10f\n",
	  thisTimer.timeElapsed());


     
#else //this is for cliquer
   
   int * redCostXInt = m_auxMemPool.intArrNCoreCols;
   int   alphaInt    = 0;

   //---
   //--- scale reduced costs (include alpha) to integers
   //---
   int scaleFactor = UtilScaleDblToIntArr(n_origCols,
					  redCostX, redCostXInt,
					  alpha,    &alphaInt);

   UTIL_DEBUG(m_param.LogDebugLevel, 4,
	      for(i = 0; i < m_nOrigRows; i++){
		 for(b = 0; b < m_beta; b++){
		    cout << "x[ " << i << "," << b << " -> " << xIndex(i,b)
			 << " ] : " << redCostXInt[xIndex(i,b)] << endl;
		 }
	      }
	      );


   //---
   //--- flip reduced costs (max c == min -c)
   //---
   UtilNegateArr(n_origCols, redCostXInt);

   UTIL_DEBUG(m_param.LogDebugLevel, 4,
	      for(i = 0; i < m_nOrigRows; i++){
		 for(b = 0; b < m_beta; b++){
		    cout << "x[ " << i << "," << b << " -> " << xIndex(i,b)
			 << " ] : " << redCostXInt[xIndex(i,b)] << endl;
		 }
	      }
	      );
   
   //---
   //--- add a constant so that all vertex weights are positive, inc alpha
   //---
   int offset = 0;
   int minrc  = *min_element(redCostXInt, redCostXInt + n_origCols);
   if(alphaInt < minrc)
      minrc = alphaInt;
   if(minrc <= 0){
      offset = -minrc + 1;
      UtilAddOffsetArr(n_origCols, offset, redCostXInt);
   }


   //the issue is "maximal"... you are finding the maximum clique ==
   //the maximum IS (in complemented graph), to find the minimum reduced cost
   //
   
   //---
   //--- set the vertex weights in the cliquer object
   //---
   m_cliquer->setVertWeight(m_gStar, redCostXInt);
   
   UTIL_DEBUG(m_param.LogDebugLevel, 4,
	      for(i = 0; i < m_nOrigRows; i++){
		 for(b = 0; b < m_beta; b++){
		    cout << "x[ " << i << "," << b << " -> " << xIndex(i,b)
			 << " ] : " << redCostX[xIndex(i,b)] << "\t" 
			 << redCostXInt[xIndex(i,b)] << endl;
		 }
	      }
	      m_cliquer->printGraph(m_gStar);
	      );

   if(checkRC){
      
      //---
      //--- we will be looking for cliques that have negative reduced cost
      //---    sum{rc[i] s[i]}     + alpha    < 0,
      //---
      //--- since cliquer is solving for maximal weighed clique,
      //--- we flipped the vertex weights 
      //---    sum{-rc[i] s[i]}    - alpha    > 0,
      //---    sum{-rcInt[i] s[i]} - alphaInt > 0,
      //---
      //--- so, the min_weight we will accept is alphaInt + epsilon
      //---

      //this won't work - might have to switch to find one... 
      //m_cliquer->cliqueFindAll(alphaInt + 2*offset, 0, 1);//THINK
     
      //this should find the max
      //m_cliquer->cliqueFindOne(0, 0, 1);
      m_cliquer->cliqueFindOne(m_gStar, 0, 0, 0);
   }
   else{
      //in general this will cosntruct all the maximal cliques - which is
      //bad for any large problems - need to limit, or just generate one
      //or use heuristic to generate a few...
      //if !checkRC, we can assume is coming from genInit, just gen one
      //m_cliquer->cliqueUnweightedFindAll(2, m_kappa, 0);

      //need to figure out a trivial ub on clique size based on costs sent
      //in? we need "size" < kappa
      //think - prune cliquer for restricted length? should be easy?
      //prune by weight and length - but then, does that solve all of
      //MAD with cliquer?
      
      //if we could generate all maximal ISs of size < kappa,
      //we'd be done - right?
      m_cliquer->cliqueFindOne(m_gStar, 0, 0, 0);
      //m_cliquer->cliqueFindAll(0, 0, 1);//THINK
   }

   UTIL_DEBUG(m_param.LogDebugLevel, 4,
	      m_cliquer->cliquePrintAll(m_gStar);
	      );

   //---
   //--- store the solution as a DecompVar and push into list
   //---
   double                 varRedCost;
   double                 varOrigCost;
   vector<int>::iterator  it;
   for(i = 0; i < m_cliquer->m_clique_count; i++) {      
      vector<int>    apInd;
      m_cliquer->cliquePopulateVector(i, apInd);
      
      varRedCost  = 0.0;
      varOrigCost = 0.0;
      vector<double> apEls(apInd.size(), 1.0);
      for(it = apInd.begin(); it != apInd.end(); it++){         
	 varRedCost  += redCostX[*it];
	 varOrigCost += origCost[*it];

	 UTIL_DEBUG(m_param.LogDebugLevel, 4,
		    printf("\nind: %d redCostX: %g varRedCost: %g",
			   *it, redCostX[*it], varRedCost);
		    );
         
      }
      UTIL_DEBUG(m_param.LogDebugLevel, 4,                
		 printf("\nvarRC = %g, alpha = %g", varRedCost, alpha);
		 );
      if(!checkRC || ((varRedCost + alpha) < -DecompEpsilon)){
	 UTIL_DEBUG(m_param.LogDebugLevel, 4,
		    printf("\nPUSH var with RC = %g\n", varRedCost + alpha);
		    );
	 vars.push_back(new DecompVar(apInd, apEls, 
				      varRedCost + alpha, varOrigCost));
      }
   }

   //---
   //--- now, we have to clean out the clique_list for the next pass
   //---
   m_cliquer->cliqueFreeMemory();
#endif

   
   m_cliquer->graphFree(m_gStar);

   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "APPsolveRelaxed()", m_param.LogDebugLevel, 2);


   
   return status; //think
}

//--------------------------------------------------------------------- //
int MAD_DecompApp::generateInitVars(DecompVarList & initVars){

   //---
   //--- x[i,b] in {0,1}, 1 = row i is assigned to block b
   //---   i in M = {1, ..., m} 
   //---   b in B = {1, ..., beta}
   //---   beta <= m (possibly an input = number of processors available)
   //---
   //--- min sum{i in M, b in B} -x[i,b]
   //--- s.t. 
   //--- (1) sum{b in B} x[i,b] <= 1, for i in M
   //--- (2) sum{i in M} x[i,b] <= k, for b in B
   //--- (3) x[i,b] + x[j,b']   <= 1, for b,b' in B, b != b',
   //---                                  i,j  in M, i != j, such that
   //---                                  a[i,k] != 0 != a[j,k], for some k
   //--- (4) x[i,b] in {0,1}
   //---

   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "generateInitVars()", m_param.LogDebugLevel, 2);

   //TODO: we should do some perturbation?
   const double * origCost = m_model.objCoeff; //all -1.0
   
   //note - if we do random perturbation and APPsolveRelaxed
   //starts calling heuristics, then just let base class do this work?
   
   //since all costs are -1.0, the sort would be arbitrary




   //---
   //--- try to generate a column with small reduced cost
   //--- that is feasible to original problem using primal heuristic
   //---
   int i;
   int nGreedyPts = 0;
   int nOrigCols  = m_nOrigRows * m_beta;   
   vector<GreedyPoint> greedyPoints;   
   nGreedyPts = heuristicGreedy(greedyPoints, INCREASING, origCost, origCost);
   for(i = 0; i < nGreedyPts; i++){
      GreedyPoint & gp = greedyPoints[i];
      
      //---
      //--- create a DecompVar from the greedySol
      //---
      UTIL_DEBUG(m_param.LogDebugLevel, 4,
		 (*m_osLog) 
		 << "Greedy DecompVar origCost = " << gp.solValueOrigCost;
		 printOriginalSolution(nOrigCols, gp.solution, &cout);
		 );
      //TODO: might be more efficient if return sparse... greedySol 
      initVars.push_back(new DecompVar(nOrigCols,
				       gp.solution,
				       0,
				       gp.solValueOrigCost));
   }

#if 0

   //---
   //--- residual capacity for blocks
   //---
   double * blockRes = m_auxMemPool.dblArrNBlocks;
   UtilFillN(blockRes, m_beta, static_cast<double>(m_kappa));

   //---
   //--- for marking if a row has already been assigned
   //---
   int    * isRowAssigned = m_auxMemPool.intArrNOrigRows;
   UtilFillN(isRowAssigned, m_nOrigRows, 0);
    
   //---
   //--- greedily assign rows to blocks
   //---    checking conflicts and capacity
   //---
   //--- for xhat[i,b],
   //---   if  block b has residual capacity 
   //---   and row i has not already been assigned
   //---   and for all b' != b
   //---       for all j in b', !isEdge(i,j)
   //---   then put i in b
   //---
   vector<int> * blocks = new vector<int>[m_beta];
   CoinAssertHint(blocks, "Error: Out of Memory");

   int  i, b, k, bp;
   int  n_cols  = m_nOrigRows * m_beta;   
   bool assignOk;
   vector<int>::iterator vit;
   pair<int,int> ib;

   for(i = 0; i < m_nOrigRows; i++)
      for(b = 0; b < m_beta; b++){

	 //for(k = 0; k < n_cols; k++){
	 //ib    = xIndexInv(k);
	 //i     = ib.first;
	 //b     = ib.second;      
	 if(isRowAssigned[i] || (blockRes[b] < DecompEpsilon))
	    continue;
      
      assignOk = true;
      for(bp = 0; bp < m_beta; bp++){
	 if(bp == b)
	    continue;
	 for(vit = blocks[bp].begin(); vit != blocks[bp].end(); vit++){
	    if(GRAPH_IS_EDGE_FAST(m_conflictGraph, i, *vit)){
	       assignOk = false;
	       break;
	    }
	 }
	 if(!assignOk)
	    break;
      }
      if(!assignOk)
	 continue;
         
      blockRes[b]      -= 1.0;
      isRowAssigned[i]  = 1;
      blocks[b].push_back(i);      
   }
   
   //---
   //--- place to store the greedy solution
   //---
   int      greedySize = 0;
   double * greedySol  = m_auxMemPool.dblArrNCoreCols;
   UtilFillN(greedySol, n_cols, 0.0);
   vector<int>    ind;
   vector<double> els;
   for(b = 0; b < m_beta; b++){
	 for(vit = blocks[b].begin(); vit != blocks[b].end(); vit++){
	    ind.push_back(xIndex(*vit,b));
	    els.push_back(1.0);
	    greedySol[xIndex(*vit,b)] = 1.0;
	    greedySize++;
	 }
      }

   

    
   initVars.push_back(new DecompVar(greedySize, &ind[0], &els[0],
				    0.0, -1.0 * greedySize));
   UTIL_DELARR(blocks);
#endif

    
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "generateInitVars()", m_param.LogDebugLevel, 2);

    

   return static_cast<int>(initVars.size());
}


//--------------------------------------------------------------------- //
void MAD_DecompApp::printOriginalColumn(const int   index, 
                                        ostream   * os) const {
   pair<int,int> p = xIndexInv(index);
   (*os) << "x[ " << p.first << " , " << p.second << " ]";
}

/*-------------------------------------------------------------------------*/
void MAD_DecompApp::printOriginalSolution(const int      n_cols, 
                                          const double * solution, 
                                          ostream      * os) const{
   int          i, j, b;
   int          border_size;
   double       xj;
   bool         isIntegral = true;

   DecompApp::printOriginalSolution(n_cols, solution, os);

   isIntegral = UtilIsIntegral(solution, n_cols);
   if(isIntegral){
      (*os) << "\nBlock Decomposition:";
      vector<bool> border(m_nOrigRows, true);      
      for(b = 0; b < m_beta; b++){
	 (*os) << "\nBLOCK " << b << ":\t";
	 for(i = 0; i < m_nOrigRows; i++){        
	    xj   = solution[xIndex(i,b)];
	    CoinAssertDebug(UtilIsIntegral(xj));
	    CoinAssertDebug(xj <  (1.0 + DecompEpsilon));
	    CoinAssertDebug(xj >  (    - DecompEpsilon));
	    if(xj > 0.5){
	       (*os) << i << " ";
	       border[i] = false;
	    }
	 }
      }    	
      border_size = count(border.begin(), border.end(), true);
      (*os) << "\nBORDER :\t";
      for(i = 0; i < m_nOrigRows; i++){
	 if(!border[i])
	    continue;
	 (*os) << i << " ";
      }
      (*os) << "\nBORDER Size =  " << border_size << "\n";


      const CoinPackedMatrix * M = m_instance.getMatrixByRow();
      for(b = 0; b < m_beta; b++){
	 (*os) << "\nBLOCK " << b << "\n";
	 for(i = 0; i < m_nOrigRows; i++){        
	    xj   = solution[xIndex(i,b)];
            if(xj > 0.5){
               CoinShallowPackedVector row = M->getVector(i);
               const int   rowLen  = row.getNumElements();
               const int * rowInd  = row.getIndices();
               (*os) << "Row i: " << i << "\t";
               for(j = 0; j < rowLen; j++){
                  (*os) << rowInd[j] << " ";
               }
               (*os) << endl;
            }
         }
      }
      (*os) << "\nBORDER\n";
      for(i = 0; i < m_nOrigRows; i++){
	 if(!border[i])
	    continue;
	 CoinShallowPackedVector row = M->getVector(i);
	 const int   rowLen  = row.getNumElements();
	 const int * rowInd  = row.getIndices();
	 (*os) << "Row i: " << i << "\t";
	 for(j = 0; j < rowLen; j++){
	    (*os) << rowInd[j] << " ";
	 }
	 (*os) << endl;
      }

      for(b = 0; b < m_beta; b++){
	 (*os) << "\nBLOCK " << b << "\n";
	 for(i = 0; i < m_nOrigRows; i++){        
	    xj   = solution[xIndex(i,b)];
            if(xj > 0.5){
               CoinShallowPackedVector row = M->getVector(i);
               const int   rowLen  = row.getNumElements();
               const int * rowInd  = row.getIndices();
               printRowMarks(rowInd, rowLen);
            }
         }
      }
      (*os) << "\nBORDER\n";
      for(i = 0; i < m_nOrigRows; i++){
	 if(!border[i])
	    continue;
	 CoinShallowPackedVector row = M->getVector(i);
	 const int   rowLen  = row.getNumElements();
	 const int * rowInd  = row.getIndices();
	 printRowMarks(rowInd, rowLen);
      }

   }
}

//TODO: visualization tool
//TODO: sanity check that really is feasible for MAD
