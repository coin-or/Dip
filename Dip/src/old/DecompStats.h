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


#ifndef DECOMP_STATS_INCLUDED
#define DECOMP_STATS_INCLUDED

#include "CoinTime.hpp"
#include "DecompPortable.h"

class DecompStats {

public:
   CoinTimer timerOverall;
   CoinTimer timerDecomp;
   CoinTimer timerOther1;
   CoinTimer timerOther2;

public:
   double totalOverall;

   double totalDecomp;
   double totalSolveRelax;
   double totalSolveRelaxApp;
   double totalSolUpdate;
   double totalGenCuts;
   double totalGenVars;

   double maxDecomp;
   double maxSolveRelax;
   double maxSolveRelaxApp;
   double maxSolUpdate;
   double maxGenCuts;
   double maxGenVars;

public:
   vector<double> thisDecomp;
   vector<double> thisSolveRelax;
   vector<double> thisSolveRelaxApp;
   vector<double> thisSolUpdate;
   vector<double> thisGenCuts;
   vector<double> thisGenVars;

public:
   void calculateStats();
   void printOverallStats (ostream* os = &cout); //ostream?
   void printDetailedStats(ostream* os = &cout); //ostream?

public:
   DecompStats() :

      timerOverall      (0),
      timerDecomp       (0),
      timerOther1       (0),
      timerOther2       (0),

      totalOverall      (0.0),

      totalDecomp       (0.0),
      totalSolveRelax   (0.0),
      totalSolveRelaxApp(0.0),
      totalSolUpdate    (0.0),
      totalGenCuts      (0.0),
      totalGenVars      (0.0),

      maxDecomp         (0.0),
      maxSolveRelax     (0.0),
      maxSolveRelaxApp  (0.0),
      maxSolUpdate      (0.0),
      maxGenCuts        (0.0),
      maxGenVars        (0.0)

   {
   }

   ~DecompStats() {}

};


#endif
