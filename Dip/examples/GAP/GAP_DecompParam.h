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

#ifndef GAP_DECOMP_PARAM_INCLUDED
#define GAP_DECOMP_PARAM_INCLUDED

//===========================================================================//
#include "UtilParameters.h"

//===========================================================================//
/*!
 * \class GAP_DecompParam
 * Storage for parameters for the 
 *   Generalized Assignment Problem (GAP).
 * 
 */

//===========================================================================//
class GAP_DecompParam {
 public:
   int    LogLevel;
   string DataDir;
   string Instance;
   bool   UsePisinger;

 public:
   void getSettings(UtilParameters & utilParam){
      static const char * common = "GAP";
      LogLevel    = utilParam.GetSetting("LogLevel",       0, common);
      DataDir     = utilParam.GetSetting("DataDir",       "",  common);
      Instance    = utilParam.GetSetting("Instance",      "",  common);    
      UsePisinger = utilParam.GetSetting("UsePisinger", true,  common);
   }
   
   void dumpSettings(ostream * os = &cout){
      static const char * common = "GAP";
      (*os) << "\n=====================================================\n"
            << "GAP_DECOMP PARAMETER SETTINGS \n";
      (*os) << common << ": LogLevel    : " << LogLevel    << endl;
      (*os) << common << ": DataDir     : " << DataDir     << endl;
      (*os) << common << ": Instance    : " << Instance    << endl;
      (*os) << common << ": UsePisinger : " << UsePisinger << endl;
      (*os) <<   "=====================================================\n";
   }
   
 public:
   GAP_DecompParam():    
      LogLevel   (0 ),
      DataDir    (""),
      Instance   (""),
      UsePisinger(true)
         {};
   ~GAP_DecompParam() {};
};

#endif
