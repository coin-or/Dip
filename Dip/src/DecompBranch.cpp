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
// Copyright (C) 2002-2013, Lehigh University, Matthew Galati, Ted Ralphs    //
// All Rights Reserved.                                                      //
//===========================================================================//


#include "DecompAlgo.h"
#include "DecompApp.h"

// --------------------------------------------------------------------- //
bool DecompAlgo::
chooseBranchSet(std::vector< std::pair<int, double> >& downBranchLB,
                std::vector< std::pair<int, double> >& downBranchUB,
                std::vector< std::pair<int, double> >& upBranchLB,
                std::vector< std::pair<int, double> >& upBranchUB)
{
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "chooseBranchSet()", m_param.LogDebugLevel, 1);
   //---
   //--- Default branching in DIP is the most simple approach possible.
   //---   Choose variables farthest from integer - based on x formulation.
   //---
   std::vector<int>::iterator intIt;
   int    branchedOnIndex, j;
   double branchedOnValue, x, dist, maxDist;
   double obj              = 0.0;
   const double* objCoeff = getOrigObjective();
   DecompConstraintSet* modelCore = m_modelCore.getModel();
   maxDist         = DecompEpsilon;//TODO: parameter
	
	
   branchedOnIndex = -1;
   branchedOnValue =  0;
   CoinAssert(modelCore->integerVars.size() > 0);
   // const std::vector<std::string> & colNames = modelCore->getColNames();
      
   if(m_param.BranchPriorityMasterOnly == true 
      && !modelCore->masterOnlyCols.empty()){
          j = branchOnMasterOnly();
          std::cout << "The branch index is " << j << std::endl;
	if (j != -1){
            x = m_xhat[j];	
            obj +=m_xhat[j]*objCoeff[j];
            branchedOnIndex = j; 
            branchedOnValue = x; 
	    m_branchingImplementation = DecompBranchInMaster;      
     }
        else
          {
	   m_param.BranchPriorityMasterOnly = false; 
	}
      }
 
   else {
   for (intIt =  modelCore->integerVars.begin();
         intIt != modelCore->integerVars.end(); intIt++) {
      j   = *intIt;
      x   = m_xhat[j];
      obj += m_xhat[j] * objCoeff[j];
      dist = fabs(x - floor(x + 0.5));

      if (dist > maxDist) {
         maxDist         = dist;
         branchedOnIndex = j;
         branchedOnValue = x;
      }
    }
     if (std::find(m_branchedMasterOnly.begin(), m_branchedMasterOnly.end(), branchedOnIndex)
         != m_branchedMasterOnly.end()){
	  branchedOnIndex = -1; 
	}
   }
   std::map<int, int >:: iterator mit;

   if (branchedOnIndex != -1) {
      //---
      //--- Example x[0]=2.5:
      //---    x[0] <= 2 (down)
      //---    x[0] >= 3 (up  )
      //---
      mit = m_masterOnlyColsMap.find(branchedOnIndex);

      if (mit != m_masterOnlyColsMap.end()) {
         // it indicates the branched variable is a master-only variable
         // we need to set the branching method to branch in the master
         m_branchingImplementation = DecompBranchInMaster;
         std::cout << "branch in master " << branchedOnIndex << std::endl;
      }
      
      //std::cout << "The branching variable is " << branchedOnIndex
      //          << "  " << colNames[branchedOnIndex] << std::endl;
      downBranchUB.push_back(std::pair<int, double>(branchedOnIndex,
                             floor(branchedOnValue)));
      upBranchLB.push_back(std::pair<int, double>(branchedOnIndex,
                           ceil(branchedOnValue)));
      UTIL_MSG(m_param.LogDebugLevel, 3,
               int nColNames = static_cast<int>(modelCore->colNames.size());
               (*m_osLog) << "branchOnInd = " << branchedOnIndex << " -> ";

               if ( branchedOnIndex  < nColNames &&
                    branchedOnIndex >= 0)
               (*m_osLog) << modelCore->colNames[branchedOnIndex];
      else {
         m_app->printOriginalColumn(branchedOnIndex, m_osLog);
         }
      (*m_osLog) << "\tbranchOnVal = " << branchedOnValue  << "\n";
              );
      return true;
   } else {
      return false;
   }

   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "chooseBranchSet()", m_param.LogDebugLevel, 1);
}

