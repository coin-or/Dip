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

//===========================================================================//
#ifndef DecompStats_h_
#define DecompStats_h_

//===========================================================================//
#include "Decomp.h"
#include "UtilTimer.h"
//===========================================================================//

//===========================================================================//
class DecompObjBound {
public:
   /**
    * The phase when bound was recorded.
    */
   int phase;
   /**
    * The cut pass when bound was recorded.
    */
   int cutPass;
   /**
    * The price pass when bound was recorded.
    */
   int pricePass;
   /**
    * The time stamp (from start) when bound was recorded.
    */
   double timeStamp;
   /**
    * The recorded continuous lower bound.
    */
   double thisBound;
   /**
    * The recorded continuous upper bound.
    */
   double thisBoundUB;
   /**
    * The best recorded continuous lower bound.
    *   global LB = max{active node lower bounds}
    */
   double bestBound;
   /**
    * The recorded integer upper bound.
    */
   double thisBoundIP;
   /**
    * The best recorded integer upper bound.
    *   global UB = min{node integer upper bounds}
    */
   double bestBoundIP;

   /**
    * Comparison operator for sorting on time.
    */
   bool operator<(const DecompObjBound& objBound) const {
      if (timeStamp < objBound.timeStamp) {
         return true;
      } else {
         return false;
      }
   }

public:
   DecompObjBound() :
      phase      (0),
      cutPass    (0),
      pricePass  (0),
      timeStamp  (0.0),
      thisBound  (-DecompInf),
      thisBoundUB( DecompInf),
      bestBound  (-DecompInf),
      thisBoundIP( DecompInf),
      bestBoundIP( DecompInf) {
   }

};

//===========================================================================//
class DecompNodeStats {
public:

   //---
   //--- Storage for the bound history for a node.
   //---    NOTE: we always assume a minimization problem
   //---

   /**
    * Storage of the bounds.
    *
    * For the continuous part:
    *   CPM  : Bounds on the objective of optimal master linear
    *          relaxation. Typically, this is an LP solved to optimality,
    *          so, LB = zCP = UB.
    *   PC/RC: Given bounds on the objective of optimal restricted master
    *          linear relaxation zPC_LB <= zPC* <= zPC_UB and a lower bound
    *          on the most negative reduced cost (RC_LB) extreme point (ray)
    *          from the subproblem polytope (for the associated master duals).
    *               LB = zPC_LB + RC_LB <= zPC* <= zPC_UB = UB
    */
   std::vector< DecompObjBound > objHistoryBound;

   /**
    * The global lower (.first) and upper (.second) bound.
    */
   std::pair<double, double> objBest;

   /**
    * The node index (in the branch-and-bound tree).
    */
   int    nodeIndex;

   /**
    * Number of cuts generated in this round of cut calls.
    */
   int    cutsThisRound;

   /**
    * Number of vars generated in this round of pricing calls.
    */
   int    varsThisRound;

   /**
    * Number of cuts generated in this particular cut call.
    */
   int    cutsThisCall;

   /**
    * Number of vars generated in this particular price call.
    */
   int    varsThisCall;

   /**
    * Number of cut calls in this node in total.
    */
   int    cutCallsTotal;

   /**
    * Number of price calls in this node in total.
    */
   int    priceCallsTotal;

   /**
    * Number of cut calls in this round.
    */
   int    cutCallsRound;

   /**
    * Number of price calls in this round.
    */
   int    priceCallsRound;

public:
   void init() {
      objHistoryBound.clear();
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
   void printObjHistoryBound  (std::ostream* os = &std::cout) const;
   inline void resetCutRound() {
      cutCallsRound = 0;
      cutsThisRound = 0;
   }
   inline void resetPriceRound() {
      priceCallsRound = 0;
      varsThisRound   = 0;
   }
   inline void resetBestLB() {
      objBest.first = -DecompInf;
   }
   inline DecompObjBound* getLastBound() {
      int nHistorySize = static_cast<int>(objHistoryBound.size());

      if (nHistorySize > 0) {
         return &(objHistoryBound[nHistorySize - 1]);
      } else {
         return 0;
      }
   }
   inline double getLastBoundThis() {
      double           thisBound = -DecompInf;
      DecompObjBound* lastBound = getLastBound();

      if (lastBound) {
         thisBound = lastBound->thisBound;
      }

      return thisBound;
   }

public:
   DecompNodeStats() :
      objHistoryBound(),
      objBest        () {
      init();
   }
};


//===========================================================================//
class DecompStats {

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
   double totalGenCutsApp;
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
   std::vector<double> thisDecomp;
   std::vector<double> thisSolveRelax;
   std::vector<double> thisSolveRelaxApp;
   std::vector<double> thisSolUpdate;
   std::vector<double> thisGenCuts;
   std::vector<double> thisGenCutsApp;
   std::vector<double> thisGenVars;
   std::vector<double> thisCompressCols;

public:
   void calculateStats();
   void printOverallStats (std::ostream* os = &std::cout); //ostream?
   void printDetailedStats(std::ostream* os = &std::cout); //ostream?

public:
   DecompStats() :

      totalOverall      (0.0),

      totalDecomp       (0.0),
      totalSolveRelax   (0.0),
      totalSolveRelaxApp(0.0),
      totalSolUpdate    (0.0),
      totalGenCuts      (0.0),
      totalGenCutsApp   (0.0),
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

