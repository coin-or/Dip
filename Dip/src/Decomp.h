//===========================================================================//
// This file is part of the DIP Solver Framework.                            //
//                                                                           //
// DIP is distributed under the Eclipse Public License as part of the        //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Authors: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)      //
//          Ted Ralphs, Lehigh University (ted@lehigh.edu)                   //
//          Jiadong Wang, Lehigh University (jiw408@lehigh.edu)              //
//                                                                           //
// Copyright (C) 2002-2019, Lehigh University, Matthew Galati, Ted Ralphs    //
// All Rights Reserved.                                                      //
//                                                                           //
// Interface to Gurobi is Copyright 2015 Jazz Aviation LP                    //
//===========================================================================//

//===========================================================================//
#ifndef Decomp_h_
#define Decomp_h_

//===========================================================================//
// Standard Headers                                                          //
//===========================================================================//
//---
//--- include the necessary standard libs
//---
#include <cstdio>
#include <cassert>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <algorithm>
#include <functional>
#include <string>
#include <map>
#include <limits>
#include <cmath>

#include "DecompConfig.h"

//===========================================================================//
// COIN Headers                                                              //
//===========================================================================//
//---
//--- include some standard COIN headers (depending on LP solver)
//---   depending on LP solver
//---
#include "CoinError.hpp"
#include "CoinFinite.hpp"
#include "CoinPackedVector.hpp"
#include "CoinPackedMatrix.hpp"

#ifdef DIP_HAS_CLP
#include "OsiClpSolverInterface.hpp"
#endif

#ifdef DIP_HAS_CPX
#include "cplex.h"
#include "OsiCpxSolverInterface.hpp"
#endif

#ifdef DIP_HAS_GRB
extern "C" {
#include "gurobi_c.h"
}
#include "OsiGrbSolverInterface.hpp"
#endif

#ifdef DIP_HAS_CBC
#include "OsiCbcSolverInterface.hpp"
#endif

#ifdef DIP_HAS_SYMPHONY
#include "symphony.h"
#include "OsiSymSolverInterface.hpp"
#endif

//===========================================================================//
// DECOMP Enums, Constants and Typedefs                                      //
//===========================================================================//

//===========================================================================//
//---
//--- DECOMP typedefs
//---
class DecompVar;
class DecompCut;
typedef std::list<DecompVar*> DecompVarList;
typedef std::list<DecompCut*> DecompCutList;

//===========================================================================//
//---
//--- DECOMP constants
//---
const double DecompBigNum   = 1.0e21;
const double DecompEpsilon  = 1.0e-6;
const double DecompZero     = 1.0e-14;

//===========================================================================//
//---
//--- DECOMP enums (for algorithms)
//---


struct DecompMainParam {
   bool doCut;
   bool doPriceCut;
   bool doDirect;
   double timeSetupCpu ;
   double timeSetupReal;
   double timeSolveCpu ;
   double timeSolveReal ;
   double bestLB;
   double bestUB;
};



enum DecompAlgoType {
   CUT,
   PRICE_AND_CUT,
   RELAX_AND_CUT,
   VOL_AND_CUT,
   DECOMP
};
const std::string DecompAlgoStr[5] = {
   "CUT",
   "PRICE_AND_CUT",
   "RELAX_AND_CUT",
   "VOL_AND_CUT",
   "DECOMP"
};

//---
//--- node stopping criteria
//---
enum DecompAlgoStop {
   DecompStopNo,
   DecompStopGap,
   DecompStopTailOff,
   DecompStopInfeasible,
   DecompStopBound,
   DecompStopTime,
   DecompStopIterLimit
};
const std::string DecompAlgoStopStr[7] = {
   "DecompStopNo",
   "DecompStopGap",
   "DecompStopTailOff",
   "DecompStopInfeasible",
   "DecompStopBound",
   "DecompStopTime",
   "DecompStopIterLimit"
};


//===========================================================================//
//---
//--- DECOMP enums (for phases)
//---
enum DecompPhase {
   PHASE_PRICE1,
   PHASE_PRICE2,
   PHASE_CUT,
   PHASE_DONE,
   PHASE_UNKNOWN
};
const std::string DecompPhaseStr[6] = {
   "PHASE_PRICE1",
   "PHASE_PRICE2",
   "PHASE_CUT",
   "PHASE_DONE",
   "PHASE_UNKNOWN"
};

//===========================================================================//
//---
//--- DECOMP enums (for status)
//---
enum DecompStatus {
   STAT_FEASIBLE,
   STAT_IP_FEASIBLE,
   STAT_INFEASIBLE,
   STAT_UNKNOWN
};
const std::string DecompStatusStr[3] = {
   "STAT_FEASIBLE",
   "STAT_INFEASIBLE",
   "STAT_UNKNOWN"
};

enum DecompPriceCutStrategy {
   Default,
   FavorPrice,
   FavorCut
};
const std::string DecompPriceCutStrategyStr[3] = {
   "Default",
   "Favor Price",
   "Favor Cut"
};

//===========================================================================//
enum DecompSolverStatus {
   DecompSolStatError = -1,
   DecompSolStatOptimal = 0,
   DecompSolStatFeasible = 1,
   DecompSolStatInfeasible = 2,
   DecompSolStatNoSolution = 3
};

//===========================================================================//
enum DecompGenericStatus {
   DecompStatOk          = 0,
   DecompStatError       = 1,
   DecompStatOutOfMemory = 2
};

enum DecompSolverType {
   DecompDualSimplex = 0,
   DecompPrimSimplex = 1,
   DecompBarrier     = 2
};

enum DecompRoundRobin {
   RoundRobinRotate    = 0,
   RoundRobinMostNegRC = 1
};

//===========================================================================//
enum DecompFunction {
   DecompFuncGeneric          = 0,
   DecompFuncGenerateInitVars = 1
};

enum DecompSubProbParallelType {
   SubProbScheduleStatic,
   SubProbScheduleDynamic,
   SubProbScheduleGuided,
   SubProbScheduleRuntime
};



//===========================================================================//
enum DecompRowType {
   //original row
   DecompRow_Original,
   //branching row
   DecompRow_Branch,
   //convexity constraint
   DecompRow_Convex,
   //row which is a cut
   DecompRow_Cut
};
const std::string DecompRowTypeStr[4] = {
   "DecompRow_Original",
   "DecompRow_Branch",
   "DecompRow_Convex",
   "DecompRow_Cut"
};

//===========================================================================//
//Corresponding to the class DecompVar
enum DecompVarType {
   // points generated from bounded subproblem
   DecompVar_Point = 0,
   // rays generated from unbounded subproblem
   DecompVar_Ray = 1
};



//===========================================================================//
enum DecompColType {
   //structural column
   DecompCol_Structural,
   //structural column (which should never be deleted)
   DecompCol_Structural_NoDelete,
   //master-only column
   DecompCol_MasterOnly,
   //artifical column for original row (L for <=)
   DecompCol_ArtForRowL,
   //artifical column for original row (G for >=)
   DecompCol_ArtForRowG,
   //artifical column for branching row (L for <=)
   DecompCol_ArtForBranchL,
   //artifical column for branching row (G for >=)
   DecompCol_ArtForBranchG,
   //artifical column for convexity row (L for <=)
   DecompCol_ArtForConvexL,
   //artifical column for convexity row (G for >=)
   DecompCol_ArtForConvexG,
   //artifical column for cut (L for <=)
   DecompCol_ArtForCutL,
   //artifical column for cutG(L for >=)
   DecompCol_ArtForCutG,
   //marker used for deletion
   DecompCol_ToBeDeleted
};
const std::string DecompColTypeStr[12] = {
   "DecompCol_Structural",
   "DecompCol_Structural_NoDelete",
   "DecompCol_MasterOnly",
   "DecompCol_ArtForRowL",
   "DecompCol_ArtForRowG",
   "DecompCol_ArtForBranchL",
   "DecompCol_ArtForBranchG",
   "DecompCol_ArtForConvexL",
   "DecompCol_ArtForConvexG",
   "DecompCol_ArtForCutL",
   "DecompCol_ArtForCutG",
   "DecompCol_ToBeDeleted"
};

enum DecompBranchingImplementation {
   DecompBranchInSubproblem,
   DecompBranchInMaster
};
/*
enum DecompNumericErrorType {



};
*/

//---
//--- COIN vectors can do some extra checking if this is true,
//---   but, it is expensive, so turn off when in optimized mode
//---
#ifdef NDEBUG
#define DECOMP_TEST_DUPINDEX false
#else
#define DECOMP_TEST_DUPINDEX true
#endif

#endif
