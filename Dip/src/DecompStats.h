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
// Conceptual Design: Matthew Galati, SAS Institute Inc.                     //
//                    Ted Ralphs, Lehigh University                          //
//                                                                           //
// Copyright (C) 2002-2009, Lehigh University, Matthew Galati, Ted Ralphs    //
// All Rights Reserved.                                                      //
//===========================================================================//

//===========================================================================//
#ifndef DecompStats_h_
#define DecompStats_h_

//===========================================================================//
#include "Decomp.h"
#include "UtilTimer.h"
//===========================================================================//


//===========================================================================//
typedef struct DecompObjBound DecompObjBound;
struct DecompObjBound
{
   int    lbOrUb;  //0=LB, 1=UB
   int    cutPass;
   int    pricePass;
   double timeStamp;
   double thisBound;
   double bestBound;

   bool operator<(const DecompObjBound & objBound) const {
      if(timeStamp < objBound.timeStamp)
	 return true;
      else 
	 return false;
   }
};

//===========================================================================//
class DecompNodeStats{
public:

   vector< DecompObjBound > objHistoryLB;
   vector< DecompObjBound > objHistoryUB;

   //vector< pair<double, double> > objHistory;
   pair<double, double>           objBest;

   int    nodeIndex;
   int    cutsThisRound;
   int    varsThisRound;
   int    cutsThisCall;
   int    varsThisCall;
   int    cutCallsTotal;
   int    priceCallsTotal;
   int    cutCallsRound;
   int    priceCallsRound;
   
public:
   void init(){
      objHistoryLB.clear();
      objHistoryUB.clear();
      objBest.first   = -DecompInf;
      objBest.second  =  DecompInf;
      nodeIndex       =  0;
      cutsThisRound   =  0;
      varsThisRound   =  0;
      cutsThisCall    =  0;
      varsThisCall    =  0;
      cutCallsTotal   =  0;
      priceCallsTotal =  0;
      cutCallsRound   =  0;
      priceCallsRound =  0;
   }

public:
   void printObjHistory(ostream * os = &cout) const;
   void printObjHistoryLB(ostream * os = &cout) const;
   void printObjHistoryUB(ostream * os = &cout) const;
   inline void resetCutRound() {
      cutCallsRound = 0;
      cutsThisRound = 0;
   }
   inline void resetPriceRound() {
      priceCallsRound = 0;
      varsThisRound   = 0;
   }

public:
   DecompNodeStats() :
      objHistoryLB(),
      objHistoryUB(),
      objBest     ()
   {
      init();
   }
};


//===========================================================================//
class DecompStats{

 public:
   UtilTimer timerOverall;
   UtilTimer timerDecomp;
   UtilTimer timerOther1;
   UtilTimer timerOther2;

 public:
   double totalOverall;

   double totalDecomp;
   double totalSolveRelax;
   double totalSolveRelaxApp;
   double totalSolUpdate;
   double totalGenCuts;
   double totalGenVars;
   double totalCompressCols;
   
   double maxDecomp;
   double maxSolveRelax;
   double maxSolveRelaxApp;
   double maxSolUpdate;
   double maxGenCuts;
   double maxGenVars;
   double maxCompressCols;

 public:
   vector<double> thisDecomp;
   vector<double> thisSolveRelax;
   vector<double> thisSolveRelaxApp;
   vector<double> thisSolUpdate;
   vector<double> thisGenCuts;
   vector<double> thisGenVars;
   vector<double> thisCompressCols;

 public:
   void calculateStats();
   void printOverallStats (ostream * os = &cout);//ostream?
   void printDetailedStats(ostream * os = &cout);//ostream?
   
 public:
   DecompStats() :
      
      totalOverall      (0.0),
      
      totalDecomp       (0.0),
      totalSolveRelax   (0.0),
      totalSolveRelaxApp(0.0),
      totalSolUpdate    (0.0),
      totalGenCuts      (0.0),
      totalGenVars      (0.0),
      totalCompressCols (0.0),

      maxDecomp         (0.0),
      maxSolveRelax     (0.0),
      maxSolveRelaxApp  (0.0),
      maxSolUpdate      (0.0),
      maxGenCuts        (0.0),
      maxGenVars        (0.0),
      maxCompressCols   (0.0)
      
      {
      }
   
   ~DecompStats() {}

};
//===========================================================================//

#endif

