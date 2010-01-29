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

// --------------------------------------------------------------------- //
#include "UtilMacros.h"
#include "DecompStats.h"


#include <iomanip>
using namespace std;
// --------------------------------------------------------------------- //

// --------------------------------------------------------------------- //
void DecompStats::calculateStats (){
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
   totalCompressCols    = accumulate(thisCompressCols.begin(), 
                                     thisCompressCols.end(), 0.0);

   //---
   //--- calculate stats max
   //---   
   vector<double>::const_iterator it;   
   if(thisDecomp.size() > 0){
      it = max_element(thisDecomp.begin(), thisDecomp.end());
      maxDecomp = *it;
   }
   if(thisSolveRelax.size() > 0){
      it = max_element(thisSolveRelax.begin(), thisSolveRelax.end());
      maxSolveRelax = *it;
   }
   if(thisSolveRelaxApp.size() > 0){
      it = max_element(thisSolveRelaxApp.begin(), thisSolveRelaxApp.end());
      maxSolveRelaxApp = *it;
   }
   if(thisSolUpdate.size() > 0){
      it = max_element(thisSolUpdate.begin(), thisSolUpdate.end());
      maxSolUpdate = *it;
   }
   if(thisGenCuts.size() > 0){
      it = max_element(thisGenCuts.begin(), thisGenCuts.end());
      maxGenCuts = *it;
   }
   if(thisGenVars.size() > 0){
      it = max_element(thisGenVars.begin(), thisGenVars.end());
      maxGenVars = *it;
   }
   if(thisCompressCols.size() > 0){
      it = max_element(thisCompressCols.begin(), thisCompressCols.end());
      maxCompressCols = *it;
   }

}


// --------------------------------------------------------------------- //
void DecompNodeStats::printObjHistory(ostream * os) const {

   //---
   //--- combine the bounds into one vector and sort
   //---
   vector<DecompObjBound> objHistory(objHistoryLB);
   objHistory.insert(objHistory.begin(), 
		     objHistoryUB.begin(), objHistoryUB.end());

   sort(objHistory.begin(), objHistory.end());

   (*os) << setiosflags(ios::fixed|ios::showpoint);   
   (*os).precision(2);    
   (*os) << "\n========== OBJ History Node " 
	 << nodeIndex << " [BEGIN]: ========= " << endl;

   vector< DecompObjBound >::iterator it;
   for(it = objHistory.begin(); it != objHistory.end(); it++){
      if((*it).lbOrUb)
	 (*os) << setw(4) << "ub";
      else
	 (*os) << setw(4) << "lb";
      (*os) << setw(6)  << (*it).cutPass
	    << setw(6)  << (*it).pricePass
	    << setw(10) << UtilDblToStr((*it).timeStamp,5);
      if((*it).lbOrUb){
	 (*os) << setw(10) << "."
	       << setw(10) << "."
	       << setw(10) << UtilDblToStr((*it).thisBound,2)
	       << setw(10) << UtilDblToStr((*it).bestBound,2);
      }
      else{
	 (*os) << setw(10) << UtilDblToStr((*it).thisBound,2)
	       << setw(10) << UtilDblToStr((*it).bestBound,2)	    
	       << setw(10) << "."
	       << setw(10) << ".";	    
      }
      (*os) << endl;
   }
   (*os) << "\n========== OBJ History Node "
	 << nodeIndex << " [END]:   ========= " << endl;
}

// --------------------------------------------------------------------- //
void DecompNodeStats::printObjHistoryLB(ostream * os) const {
   (*os) << setiosflags(ios::fixed|ios::showpoint);   
   (*os).precision(2);    
   (*os) << "\n========== OBJ LB History Node " 
	 << nodeIndex << " [BEGIN]: ========= " << endl;   
   vector< DecompObjBound >::const_iterator it;
   for(it = objHistoryLB.begin(); it != objHistoryLB.end(); it++){
      (*os) << setw(6)  << (*it).cutPass
	    << setw(6)  << (*it).pricePass
	    << setw(10) << UtilDblToStr((*it).timeStamp,5)
	    << setw(10) << UtilDblToStr((*it).thisBound,2)
	    << setw(10) << UtilDblToStr((*it).bestBound,2)
	    << endl;
   }
   (*os) << "\n========== OBJ LB History Node "
	 << nodeIndex << " [END]:   ========= " << endl;
}

// --------------------------------------------------------------------- //
void DecompNodeStats::printObjHistoryUB(ostream * os) const {
   (*os) << setiosflags(ios::fixed|ios::showpoint);   
   (*os).precision(2);    
   (*os) << "\n========== OBJ UB History Node " 
	 << nodeIndex << " [BEGIN]: ========= " << endl;   
   vector< DecompObjBound >::const_iterator it;
   for(it = objHistoryUB.begin(); it != objHistoryUB.end(); it++){
      (*os) << setw(6)  << (*it).cutPass
	    << setw(6)  << (*it).pricePass
	    << setw(10) << UtilDblToStr((*it).timeStamp,5)
	    << setw(10) << UtilDblToStr((*it).thisBound,2)
	    << setw(10) << UtilDblToStr((*it).bestBound,2)
	    << endl;
   }
   (*os) << "\n========== OBJ UB History Node "
	 << nodeIndex << " [END]:   ========= " << endl;
}

// --------------------------------------------------------------------- //
void DecompStats::printOverallStats (ostream * os){
   
   calculateStats();

   (*os) << setiosflags(ios::fixed|ios::showpoint);   
   (*os).precision(2);    
   (*os) << "\n========== DECOMP Statistics [BEGIN]: ========= ";
   totalOverall = totalDecomp;
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
   (*os) << setw(40) << "\nTotal Compress Cols   = " 
         << setw(10) << totalCompressCols
         << setw(10) << 100.0 * totalCompressCols / totalOverall
         << setw(6)  << thisCompressCols.size()
         << setw(6)  << maxCompressCols
      ;
   (*os) << "\n========== DECOMP Statistics [END  ]: ========= \n";
}


// --------------------------------------------------------------------- //
void printDetailedStats(ostream * os = &cout){
}
