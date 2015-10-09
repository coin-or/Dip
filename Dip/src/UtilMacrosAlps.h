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
// Copyright (C) 2002-2015, Lehigh University, Matthew Galati, Ted Ralphs    //
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
