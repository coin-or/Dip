//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, Lehigh University                                 //
//                                                                           //
// Copyright (C) 2002-2007, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

#ifndef MAD_DECOMP_PARAM_INCLUDED
#define MAD_DECOMP_PARAM_INCLUDED

// --------------------------------------------------------------------- //
#include "UtilParameters.h"

using namespace std;

// --------------------------------------------------------------------- //
/*!
 * \class MAD_DecompParam
 * Storage for parameters for the Matrix Decomposition Problem (MAD).
 * 
 * \todo think about this design, register parameters, isoptional
 *       combos that violate, throw exceptions
 *
 */

// --------------------------------------------------------------------- //
class MAD_DecompParam{
 public:
   string DataDir;
   string DataSubDir;
   string Instance;
   int    NumBlocks;
   int    Capacity;

 public:
   void getSettings(UtilParameters & utilParam){
      static const char * common = "MAD";
      DataDir     = utilParam.GetSetting("DataDir",     ".",  common);
      DataSubDir  = utilParam.GetSetting("DataSubDir",  ".",  common);
      Instance    = utilParam.GetSetting("Instance",    ".",  common);
      NumBlocks   = utilParam.GetSetting("NumBlocks",   2,    common);
      Capacity    = utilParam.GetSetting("Capacity",   -1,    common);
   }

   void dumpSettings(ostream * os = &cout){
      static const char * common = "MAD";
      (*os) << "\n=====================================================\n"
            << "MAD_DECOMP PARAMETER SETTINGS \n";
      (*os) << common << ": DataDir    : " << DataDir     << endl;
      (*os) << common << ": DataSubDir : " << DataSubDir  << endl;
      (*os) << common << ": Instance   : " << Instance    << endl;
      (*os) << common << ": NumBlocks  : " << NumBlocks   << endl;
      (*os) << common << ": Capacity   : " << Capacity    << endl;
      (*os) <<   "=====================================================\n";
   }

 public:
   MAD_DecompParam():    
      DataDir   ("."),
      DataSubDir("."),
      Instance  ("."),
      NumBlocks (2  ),
      Capacity  (-1 )
   {
   };
   ~MAD_DecompParam() {};
};

#endif
