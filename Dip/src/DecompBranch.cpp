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
// Copyright (C) 2002-2009, Lehigh University, Matthew Galati, Ted Ralphs    //
// All Rights Reserved.                                                      //
//===========================================================================//


#include "DecompAlgo.h"
#include "DecompApp.h"
// --------------------------------------------------------------------- //
int DecompAlgo::chooseBranchVar(int    & branchedOnIndex,
				 double & branchedOnValue){
   //choose variables farthest from integer - based on x formulation
   vector<int>::iterator intIt;
   int    j;
   double x, dist, maxDist;
   
   //weird - integerVars are defined in Relax, not Core - ugh
   //should be for overall problem
  
   double obj = 0.0;
   const double * objCoeff = getOrigObjective();
   DecompConstraintSet          * modelCore   = m_modelCore.getModel();
  
   maxDist         = DecompEpsilon;//TODO: parameter
   branchedOnIndex = -1;
   branchedOnValue =  0;

   CoinAssert(modelCore->integerVars.size() > 0);
   for(intIt =  modelCore->integerVars.begin();
       intIt != modelCore->integerVars.end(); intIt++){
      j   = *intIt;
      x   = m_xhat[j];
      obj += m_xhat[j] * objCoeff[j];
      dist = fabs(x - floor(x+0.5));
      //printf("\nx[%d]: %g, dist: %g", j, x, dist);
      //if(dist > 1.0e-10)
      // printf("x[%d]: %g, dist: %g\n", j, x, dist);
      if(dist > maxDist){
	 maxDist         = dist;
	 branchedOnIndex = j;
	 branchedOnValue = x;
      }
   }
   
   //CoinAssert(branchedOnIndex >= 0);
   UTIL_MSG(m_param.LogDebugLevel, 3,
	    int nColNames = static_cast<int>(modelCore->colNames.size());
	    (*m_osLog) << "branchOnInd = " << branchedOnIndex << " -> ";
	    if( branchedOnIndex  < nColNames && 
		branchedOnIndex >= 0)
		(*m_osLog) << modelCore->colNames[branchedOnIndex];
	    else
	       m_app->printOriginalColumn(branchedOnIndex, m_osLog);
	    (*m_osLog) << "\tbranchOnVal = " << branchedOnValue  << "\n";
	    );
   return 0;  
}

