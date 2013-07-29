//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, Lehigh University                                 //
//                                                                           //
// Copyright (C) 2002-2013, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

#ifndef GAP_STATUS_INCLUDED
#define GAP_STATUS_INCLUDED

//===========================================================================//
#include <string>
//===========================================================================//
//---
//--- return codes
//---
//===========================================================================//
enum GAPStatus {
   GAPStatusOk = 0,
   GAPStatusError,      //general error
   GAPStatusFileIO,     //error in i/o
   GAPStatusOutOfMemory //out of memory
};
const std::string GAPStatusStr[4] = {
   "GAPStatusOk",
   "GAPStatusError",
   "GAPStatusFileIO",
   "GAPStatusOutOfMemory"
};
#endif
