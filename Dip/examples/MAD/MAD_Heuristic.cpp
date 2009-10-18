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

#include "MAD_DecompApp.h"
#include "MAD_DecompSolution.h"


// --------------------------------------------------------------------- //
//name?
int MAD_DecompApp::APPheuristics(const double            * xhat,
				 const double            * origCost,
                                 vector<DecompSolution*> & xhatIPFeas){
   
   int nVars = 0;
   nVars += heuristicGreedy(DECREASING, xhat, origCost, xhatIPFeas);
   return nVars;
}

// --------------------------------------------------------------------- //
//TODO: this revisits an old question about vector<T> vs vector<T*>??
int MAD_DecompApp::heuristicGreedy(vector<GreedyPoint>     & greedyPoints,
				   const MAD_HeurSortDir     sortDir,
				   const double            * sortValues,
				   const double            * origCost,
				   const double            * redCost){

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "heuristicGreedy()", m_param.LogDebugLevel, 2);
   
   //---
   //--- One example of sortDir = DECREASING
   //---   order non-increasing based on xhat (current LP point)
   //---   so, we take the values that are at 1.0 first, etc... 
   //---
   //--- One example of sortDir = INCREASING
   //---   order non-dercreasing based on costs
   //---   so, we take the cheapest first
   //---
   pair<int,int> ib;
   int       k, i, b, bp;
   int       nGreedyPts        = 0;
   const int n_origCols        = m_nOrigRows * m_beta;
   pair<int,double> * colSort  = m_auxMemPool.pIntDblArrNCoreCols;
   for(i = 0; i < n_origCols; i++){
      colSort[i].first  = i;
      colSort[i].second = sortValues[i];
   }
   if(sortDir == INCREASING)
      sort(colSort, colSort + n_origCols, UtilIsLessThan<int,double>());
   else
      sort(colSort, colSort + n_origCols, UtilIsGreaterThan<int,double>());
   
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
   //TODO: use mem pool?
   vector<int> * blocks = new vector<int>[m_beta];
   CoinAssertHint(blocks, "Error: Out of Memory");

   bool assignOk;
   vector<int>::iterator vit;
   for(k = 0; k < n_origCols; k++){
      ib    = xIndexInv(colSort[k].first);
      i     = ib.first;
      b     = ib.second;      
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
   int      xInd;
   double   solValueOrigCost = 0.0;
   double   solValueRedCost  = 0.0; 
   double * greedySol    = m_auxMemPool.dblArrNCoreCols;
   UtilFillN(greedySol, n_origCols, 0.0);
   for(b = 0; b < m_beta; b++){
      for(vit = blocks[b].begin(); vit != blocks[b].end(); vit++){
	 xInd              = xIndex(*vit,b);
         greedySol[xInd]   = 1.0;
	 solValueOrigCost += origCost[xInd];
      }
   }
   if(redCost){
      for(b = 0; b < m_beta; b++){
	 for(vit = blocks[b].begin(); vit != blocks[b].end(); vit++){
	    solValueRedCost += redCost[xIndex(*vit,b)];
	 }
      }
   }
   
   //---
   //--- create a GreedyPoint and set solution ptr
   //---
   GreedyPoint gp;
   gp.solValueOrigCost = solValueOrigCost;
   gp.solValueRedCost  = solValueRedCost;
   gp.solution         = greedySol;
   greedyPoints.push_back(gp);
   
   nGreedyPts++;
   
   UTIL_DELARR(blocks);   
   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "heuristicGreedy()", m_param.LogDebugLevel, 2);
   
   return nGreedyPts;
}

// --------------------------------------------------------------------- //
int MAD_DecompApp::heuristicGreedy(const MAD_HeurSortDir     sortDir,
				   const double            * sortValues,
				   const double            * origCost,
                                   vector<DecompSolution*> & solVec){
   
   
   int       i;
   int       nGreedyPts = 0;
   const int n_origCols = m_nOrigRows * m_beta;
   
   vector<GreedyPoint> greedyPoints;  
   nGreedyPts = heuristicGreedy(greedyPoints, sortDir, sortValues, origCost);
   
   for(i = 0; i < nGreedyPts; i++){
      GreedyPoint & gp = greedyPoints[i];
      
      //---
      //--- create a DecompSolution from the greedySol
      //---
      //why not also use these as columns??
      UTIL_DEBUG(m_param.LogDebugLevel, 4,
		 (*m_osLog) 
		 << "Greedy DecompSolution cost = " << gp.solValueOrigCost
		 << endl;
		 printOriginalSolution(n_origCols, gp.solution, &cout);
		 );
      //TODO: for use in app solved - might not want this...  
      solVec.push_back(new MAD_DecompSolution(this,
					      n_origCols,
					      gp.solution,
					      gp.solValueOrigCost)
		       );
   }     
   return nGreedyPts;
}
