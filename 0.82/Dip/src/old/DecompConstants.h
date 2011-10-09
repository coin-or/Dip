//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)       //
//                                                                           //
// Copyright (C) 2002-2007, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

#ifndef DECOMP_CONSTANTS_INCLUDED
#define DECOMP_CONSTANTS_INCLUDED

#include "DecompTypes.h"

const char DecompVersion[10] = "0.1";
const double DecompEpsilon   = 1.0e-6;

#ifdef __DECOMP_LP_CLP__
#include "OsiClpSolverInterface.hpp"
const double DecompInf = OsiClpInfinity;
#endif

#ifdef __DECOMP_LP_CPX__
#include "OsiCpxSolverInterface.hpp"
const double DecompInf = CPX_INFBOUND;
#endif

enum decompAlgoType{CUT,
                    PRICE_AND_CUT,
                    RELAX_AND_CUT,
                    VOL_AND_CUT,
                    DECOMP};
enum decompPhase{PHASE_INIT,
                 PHASE_PRICE, 
                 PHASE_CUT, 
                 PHASE_DONE,
                 PHASE_UNKNOWN};

const char decompPhaseStr[5][20] = {"PHASE_INIT", 
                                    "PHASE_PRICE", 
                                    "PHASE_CUT", 
                                    "PHASE_DONE",
                                    "PHASE_UNKNOWN"};

enum decompStat{STAT_FEASIBLE, 
                STAT_INFEASIBLE,
                STAT_UNKNOWN};

const char decompStatStr[3][20] = {"STAT_FEASIBLE", 
                                   "STAT_INFEASIBLE",
                                   "STAT_UNKNOWN"};

#endif
