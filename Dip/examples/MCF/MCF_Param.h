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

#ifndef MCF_PARAM_INCLUDED
#define MCF_PARAM_INCLUDED

//===========================================================================//
#include "UtilMacros.h"
#include "UtilMacrosDecomp.h"
#include "UtilParameters.h"
//===========================================================================//

//===========================================================================//
/*!
 * \class MCF_DecompParam
 * Storage for parameters for the 
 *   Multi-Dimensional Multi-Choice Knapsack Problem (MCF).
 */

//===========================================================================//
class MCF_Param {
public:
   int    LogLevel;           //application log level
   int    UseSparse;          //use sparse version of relaxations
   string DataDir;            //data directory
   string Instance;           //name of instance

public:
   void getSettings(UtilParameters & utilParam){
      static const char * common = "MCF";
      LogLevel       = utilParam.GetSetting("LogLevel",        0,      common);
      UseSparse      = utilParam.GetSetting("UseSparse",       0,      common);
      DataDir        = utilParam.GetSetting("DataDir",        "",      common);
      Instance       = utilParam.GetSetting("Instance",       "",      common);
      if(!checkOptions())
         throw UtilException("Bad Parameter", "getSettings", "MCF_Param");
   }

   bool checkOptions(){
      bool optionsOk = true;
      return optionsOk;
   }

   void dumpSettings(ostream * os = &cout){
      static const char * common = "MCF";
      (*os) << "\n=====================================================\n"
            << "MCF_DECOMP PARAMETER SETTINGS \n";
      (*os) << common << ": LogLevel          : " << LogLevel      << endl;
      (*os) << common << ": UseSparse         : " << UseSparse     << endl;
      (*os) << common << ": DataDir           : " << DataDir       << endl;
      (*os) << common << ": Instance          : " << Instance      << endl;
      (*os) <<   "=====================================================\n";
   }
   
public:
   MCF_Param() :
      LogLevel          (0),
      UseSparse         (0),
      DataDir           (""),
      Instance          ("")
   {}
   ~MCF_Param() {};
};

#endif
