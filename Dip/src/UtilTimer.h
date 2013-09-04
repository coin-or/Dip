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
// Copyright (C) 2002-2013, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

#ifndef UTIL_TIMER_INCLUDED
#define UTIL_TIMER_INCLUDED

//===========================================================================//
#include "CoinTime.hpp"

//===========================================================================//
/* A timer used to record cpu and wallclock time. */
class UtilTimer {
private:
   /** Start, end markers. */
   double startCpu_;
   double finishCpu_;
   double startReal_;
   double finishReal_;

   /** Cpu time. */
   double cpu_;

   /** Real clock time. */
   double real_;

public:
   UtilTimer() {
      reset();
   }
   ~UtilTimer()  {}

   /** Reset. */
   inline void reset() {
      start();
      finishCpu_  = 0.0;
      finishReal_ = 0.0;
      cpu_        = 0.0;
      real_       = 0.0;
   }

   /** Start to count times. */
   inline void start() {
      startCpu_  = CoinCpuTime();
      startReal_ = CoinGetTimeOfDay();
   }

   /** Stop timer and computing times. */
   inline void stop() {
      finishCpu_  = CoinCpuTime();
      finishReal_ = CoinGetTimeOfDay();
      cpu_        = finishCpu_ - startCpu_;
      real_       = finishReal_ - startReal_;
   }

   /** Get cpu timee. */
   inline double getCpuTime() {
      finishCpu_ = CoinCpuTime();
      cpu_       = finishCpu_ - startCpu_;
      return cpu_;
   }

   /** Get cpu timee. */
   inline double getRealTime() {
      finishReal_ = CoinGetTimeOfDay();
      real_       = finishReal_ - startReal_;
      return real_;
   }

   /** Return whether the given amount of real time has
       elapsed since the timer was started */
   inline bool isPast(double limit) {
      return (getRealTime() > limit);
   }

};

#endif
