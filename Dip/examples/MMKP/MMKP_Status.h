//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, Lehigh University                                 //
//                                                                           //
// Copyright (C) 2002-2009, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

#ifndef MMKP_STATUS_INCLUDED
#define MMKP_STATUS_INCLUDED

//===========================================================================//
#include <string>
//===========================================================================//
//---
//--- return codes
//---
//===========================================================================//
enum MMKPStatus {
   MMKPStatusOk = 0,
   MMKPStatusError,      //general error
   MMKPStatusFileIO,     //error in i/o
   MMKPStatusOutOfMemory //out of memory
};
const std::string MMKPStatusStr[4] = {
   "MMKPStatusOk",
   "MMKPStatusError",
   "MMKPStatusFileIO",
   "MMKPStatusOutOfMemory"
};
#endif
