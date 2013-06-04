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

#ifndef UtilMacrosAlps_h_
#define UtilMacrosAlps_h_

//===========================================================================//
#include "Alps.h"

//TODO: this is in BlisHelp, should perhaps be in a AlpsHelp level?

//===========================================================================//
class AlpsEncoded;
class CoinWarmStartBasis;

//===========================================================================//
/** Pack coin warm start into an encoded object. */
int UtilAlpsEncodeWarmStart(AlpsEncoded*               encoded,
                            const CoinWarmStartBasis* ws);

/** Unpack coin warm start from an encoded object. */
CoinWarmStartBasis* UtilAlpsDecodeWarmStart(AlpsEncoded&       encoded,
      AlpsReturnStatus* rc);

#endif
