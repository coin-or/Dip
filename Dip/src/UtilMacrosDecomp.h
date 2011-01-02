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

#ifndef UTIL_MACROS_DECOMP_INCLUDED
#define UTIL_MACROS_DECOMP_INCLUDED

// =========================================================================
#define UTIL_USE_TIMERS

// =========================================================================
#include "Decomp.h"
#include "CoinError.hpp"
#include "CoinPackedVector.hpp"
#include "CoinHelperFunctions.hpp"

// =========================================================================
class DecompApp;

// =========================================================================
#ifdef UTIL_USE_TIMERS
#include "UtilTimer.h"
//NOTE: this is not thread safe!
static UtilTimer              globalTimer;
static map<string, UtilTimer> globalTimerFuncMap;
#endif

// =========================================================================
#define UtilException(msg,methodN,classN) \
   CoinError(msg,methodN,classN,__FILE__,__LINE__)
#define UtilExceptionMemory(methodN,classN) \
   UtilException("Out of memory",methodN,classN)


// =========================================================================
// Debug Macros                                                             
// =========================================================================
#ifdef UTIL_USE_TIMERS
// ------------------------------------------------------------------------- //
inline void UtilPrintFuncBegin(ostream      * os,
			       const string & classTag,
			       const string & funcName,
			       const int      logLevel,
			       const int      logLimit){

   const size_t nDashes = 30;
   string      funcKey       = classTag + funcName;
   UtilTimer & thisFuncTimer = globalTimerFuncMap[funcKey];
   thisFuncTimer.reset();
   if(logLevel >= logLimit){
      size_t i;
      string funcBegin = "<--- " + funcName + " ";
      for(i = funcBegin.size(); i < nDashes; i++){
	 funcBegin += "-";
      }
      (*os) << left << setw(9) << classTag << ": "
	    << setprecision(3) << setw(8) << globalTimer.getRealTime()
	    << " [CPU: " << setprecision(3) << setw(8) 
            << globalTimer.getCpuTime() << "] " << funcBegin << "\n";      
   }
}

// ------------------------------------------------------------------------- //
inline void UtilPrintFuncEnd(ostream      * os,
			     const string & classTag,
			     const string & funcName,
			     const int      logLevel,
			     const int      logLimit){
   
   const size_t nDashes = 30;
   string      funcKey       = classTag + funcName;
   UtilTimer & thisFuncTimer = globalTimerFuncMap[funcKey];
   if(logLevel >= logLimit){
      size_t i;
      string funcEnd = " --- " + funcName + " ";
      for(i = funcEnd.size(); i < nDashes; i++){
	 funcEnd += "-";
      }
      funcEnd += ">";
      (*os) << left << setw(9) << classTag << ": " 
	    << setprecision(3) << setw(8) << globalTimer.getRealTime() 
	    << " [CPU: " << setprecision(4) << setw(8) 
            << globalTimer.getCpuTime() << "] " << funcEnd << " funcT = "
	    << setprecision(3) << setw(8) << thisFuncTimer.getCpuTime()
	    << "\n";      
   }
}

#else

// ------------------------------------------------------------------------- //
inline void UtilPrintFuncBegin(ostream      * os,
			       const string & classTag,
			       const string & funcName,
			       const int      logLevel,
			       const int      logLimit){

   const int nDashes = 30;
   if(logLevel >= logLimit){
      int i;
      string funcBegin = "<--- " + funcName + " ";
      for(i = funcBegin.size(); i < nDashes; i++){
	 funcBegin += "-";
      }
      (*os) << left << setw(9) << classTag << ": " << funcBegin << "\n";      
   }
}

// ------------------------------------------------------------------------- //
inline void UtilPrintFuncEnd(ostream      * os,
			     const string & classTag,
			     const string & funcName,
			     const int      logLevel,
			     const int      logLimit){
   
   const int nDashes = 30;
   if(logLevel >= logLimit){
      int i;
      string funcEnd = " --- " + funcName + " ";
      for(i = funcEnd.size(); i < nDashes; i++){
	 funcEnd += "-";
      }
      funcEnd += ">";
      (*os) << left << setw(9) << classTag << ": " << funcEnd << "\n";      
   }
}

#endif

// =========================================================================
// COIN Macros       
// TODO: anything that depends on COIN should probably not be in util lib
// =========================================================================

// ------------------------------------------------------------------------- //
/**
 * Calculate gap: |(ub-lb)|/|lb|
 */
inline double UtilCalculateGap(const double boundLB, 
			       const double boundUB) {
   double gap = DecompInf;
   if(boundLB > -DecompInf && boundUB < DecompInf){
      if(boundLB != 0.0)
	 gap = fabs(boundUB-boundLB)/fabs(boundLB);
      else
	 gap = fabs(boundUB);
   }
   return gap;
}


// ------------------------------------------------------------------------- //
CoinPackedVector * UtilPackedVectorFromDense(const int      len,
                                             const double * dense,
                                             const double   etol);
void UtilPackedVectorFromDense(const int          len,
                               const double     * dense,
                               const double       etol,
                               CoinPackedVector & v);

void UtilPrintPackedVector(const CoinPackedVector & v,
                           ostream                * os  = &cout,
                           DecompApp              * app = 0);
void UtilPrintPackedVector(const CoinPackedVector & v,
                           ostream                * os,
                           const vector<string>         & colNames,
			   const double           * value = NULL);



#endif
