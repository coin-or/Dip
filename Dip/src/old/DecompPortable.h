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

#ifndef DECOMP_PORTABLE_INCLUDED
#define DECOMP_PORTABLE_INCLUDED

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
using namespace std;

#include "CoinPackedVector.hpp"
#include "CoinPackedMatrix.hpp"

#include "UtilMacros.h"

#include "DecompConfig.h"
#include "DecompConstants.h"
#include "DecompTypes.h"

#ifdef MODE_DEBUG_FULL
#define DECOMP_TEST_DUPINDEX true
#else
#define DECOMP_TEST_DUPINDEX false
#endif

#endif
