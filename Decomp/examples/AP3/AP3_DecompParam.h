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

#ifndef AP3_DECOMP_PARAM_INCLUDED
#define AP3_DECOMP_PARAM_INCLUDED

// --------------------------------------------------------------------- //
#include "UtilParameters.h"

// --------------------------------------------------------------------- //
/*!
 * \class AP3_DecompParam
 * Storage for parameters for the 3-Indexed Assignment Problem (AP3).
 * 
 * \todo think about this design, register parameters, isoptional
 *       combos that violate, throw exceptions
 *
 */

// --------------------------------------------------------------------- //
class AP3_DecompParam{
 public:
   string DataDir;
   string Instance;

 public:
   void getSettings(UtilParameters & utilParam){
      static const char * common = "AP3";
      DataDir  = utilParam.GetSetting("DataDir",  "",  common);
      Instance = utilParam.GetSetting("Instance", "",  common);    
   }

   void dumpSettings(ostream * os = &cout){
      static const char * common = "AP3";
      (*os) << "\n=====================================================\n"
            << "AP3_DECOMP PARAMETER SETTINGS \n";
      (*os) << common << ": DataDir   : " << DataDir  << endl;
      (*os) << common << ": Instance  : " << Instance << endl;
      (*os) <<   "=====================================================\n";
   }

 public:
   AP3_DecompParam():    
      DataDir(""),
      Instance("") {
   };
   ~AP3_DecompParam() {};
};

#endif
