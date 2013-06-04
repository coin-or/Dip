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

// --------------------------------------------------------------------- //
#include "DecompStats.h"

#include <iomanip>
using namespace std;
// --------------------------------------------------------------------- //

// --------------------------------------------------------------------- //
void DecompStats::calculateStats ()
{
   //---
   //--- calculate stats totals
   //---
   totalDecomp          = accumulate(thisDecomp.begin(),
                                     thisDecomp.end(), 0.0);
   totalSolveRelax      = accumulate(thisSolveRelax.begin(),
                                     thisSolveRelax.end(), 0.0);
   totalSolveRelaxApp   = accumulate(thisSolveRelaxApp.begin(),
                                     thisSolveRelaxApp.end(), 0.0);
   totalSolUpdate       = accumulate(thisSolUpdate.begin(),
                                     thisSolUpdate.end(), 0.0);
   totalGenCuts         = accumulate(thisGenCuts.begin(),
                                     thisGenCuts.end(), 0.0);
   totalGenVars         = accumulate(thisGenVars.begin(),
                                     thisGenVars.end(), 0.0);
   //---
   //--- calculate stats max
   //---
   vector<double>::const_iterator it;

   if (thisDecomp.size() > 0) {
      it = max_element(thisDecomp.begin(), thisDecomp.end());
      maxDecomp = *it;
   }

   if (thisSolveRelax.size() > 0) {
      it = max_element(thisSolveRelax.begin(), thisSolveRelax.end());
      maxSolveRelax = *it;
   }

   if (thisSolveRelaxApp.size() > 0) {
      it = max_element(thisSolveRelaxApp.begin(), thisSolveRelaxApp.end());
      maxSolveRelaxApp = *it;
   }

   if (thisSolUpdate.size() > 0) {
      it = max_element(thisSolUpdate.begin(), thisSolUpdate.end());
      maxSolUpdate = *it;
   }

   if (thisGenCuts.size() > 0) {
      it = max_element(thisGenCuts.begin(), thisGenCuts.end());
      maxGenCuts = *it;
   }

   if (thisGenVars.size() > 0) {
      it = max_element(thisGenVars.begin(), thisGenVars.end());
      maxGenVars = *it;
   }
}

// --------------------------------------------------------------------- //
void DecompStats::printOverallStats (ostream* os)
{
   calculateStats();
   (*os).precision(2);
   (*os) << "\n========== DECOMP Statistics [BEGIN]: ========= ";
   (*os) << setw(40) << "\nTotal Overall         = "
         << setw(10) << totalOverall
         ;
   (*os) << setw(40) << "\nTotal Decomp          = "
         << setw(10) << totalDecomp
         << setw(10) << 100.0 * totalDecomp / totalOverall
         << setw(6)  << thisDecomp.size()
         << setw(6)  << maxDecomp
         ;
   (*os) << setw(40) << "\nTotal Solve Relax     = "
         << setw(10) << totalSolveRelax
         << setw(10) << 100.0 * totalSolveRelax / totalOverall
         << setw(6)  << thisSolveRelax.size()
         << setw(6)  << maxSolveRelax
         ;
   (*os) << setw(40) << "\nTotal Solve Relax App = "
         << setw(10) << totalSolveRelaxApp
         << setw(10) << 100.0 * totalSolveRelaxApp / totalOverall
         << setw(6)  << thisSolveRelaxApp.size()
         << setw(6)  << maxSolveRelaxApp
         ;
   (*os) << setw(40) << "\nTotal Solution Update = "
         << setw(10) << totalSolUpdate
         << setw(10) << 100.0 * totalSolUpdate / totalOverall
         << setw(6)  << thisSolUpdate.size()
         << setw(6)  << maxSolUpdate
         ;
   (*os) << setw(40) << "\nTotal Generate Cuts   = "
         << setw(10) << totalGenCuts
         << setw(10) << 100.0 * totalGenCuts / totalOverall
         << setw(6)  << thisGenCuts.size()
         << setw(6)  << maxGenCuts
         ;
   (*os) << setw(40) << "\nTotal Generate Vars   = "
         << setw(10) << totalGenVars
         << setw(10) << 100.0 * totalGenVars / totalOverall
         << setw(6)  << thisGenVars.size()
         << setw(6)  << maxGenVars
         ;
   (*os) << "\n========== DECOMP Statistics [END  ]: ========= \n";
}


// --------------------------------------------------------------------- //
void printDetailedStats(ostream* os = &cout)
{
}
