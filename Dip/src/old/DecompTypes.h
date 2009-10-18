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

#ifndef DECOMP_TYPES_INCLUDED
#define DECOMP_TYPES_INCLUDED

class DecompVar;
class DecompCut;

#include "DecompPortable.h"

#ifdef __DECOMP_LP_CLP__
#include "OsiClpSolverInterface.hpp"
typedef OsiClpSolverInterface OsiLpSolverInterface;
#endif

#ifdef __DECOMP_LP_CPX__
#include "OsiCpxSolverInterface.hpp"
typedef OsiCpxSolverInterface OsiLpSolverInterface;
#endif

#ifdef __DECOMP_IP_CBC__
#include "OsiCbcSolverInterface.hpp"
typedef OsiCbcSolverInterface OsiIpSolverInterface;
#endif

#ifdef __DECOMP_LP_CPX__
#include "OsiCpxSolverInterface.hpp"
typedef OsiCpxSolverInterface OsiIpSolverInterface;
#endif

typedef std::list<DecompVar*> DecompVarList;
typedef std::list<DecompCut*> DecompCutList;

#endif
