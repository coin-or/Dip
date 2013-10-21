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

//copyright

//========================================================================== //
#include "UtilMacrosAlps.h"

#include "AlpsEncoded.h"
#include "CoinWarmStartBasis.hpp"

//========================================================================== //

//---
//--- helper functions that depend only on COIN and Alps
//---

//========================================================================== //
int UtilAlpsEncodeWarmStart(AlpsEncoded*               encoded,
                            const CoinWarmStartBasis* ws)
{
   int status = 0;
   int numCols = ws->getNumStructural();
   int numRows = ws->getNumArtificial();
   encoded->writeRep(numCols);
   encoded->writeRep(numRows);
   // Pack structural.
   int nint = (ws->getNumStructural() + 15) >> 4;
   encoded->writeRep(ws->getStructuralStatus(), nint * 4);
   // Pack artificial.
   nint = (ws->getNumArtificial() + 15) >> 4;
   encoded->writeRep(ws->getArtificialStatus(), nint * 4);
   return status;
}

//===========================================================================//
CoinWarmStartBasis* UtilAlpsDecodeWarmStart(AlpsEncoded&       encoded,
      AlpsReturnStatus* rc)
{
   //rc not used? not checked?
   int numCols;
   int numRows;
   encoded.readRep(numCols);
   encoded.readRep(numRows);
   int tempInt;
   // Structural
   int nint = (numCols + 15) >> 4;
   char* structuralStatus = new char[4 * nint];
   encoded.readRep(structuralStatus, tempInt);
   assert(tempInt == nint * 4);
   // Artificial
   nint = (numRows + 15) >> 4;
   char* artificialStatus = new char[4 * nint];
   encoded.readRep(artificialStatus, tempInt);
   assert(tempInt == nint * 4);
   CoinWarmStartBasis* ws = new CoinWarmStartBasis();

   if (!ws) {
      throw CoinError("Out of memory", "UtilAlpsDecodeWarmStart", "HELP");
   }

   ws->assignBasisStatus(numCols, numRows,
                         structuralStatus, artificialStatus);
   return ws;
}

