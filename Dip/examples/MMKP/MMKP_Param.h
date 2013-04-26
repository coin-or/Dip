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

#ifndef MMKP_PARAM_INCLUDED
#define MMKP_PARAM_INCLUDED

//===========================================================================//
#include "UtilMacros.h"
#include "UtilMacrosDecomp.h"
#include "UtilParameters.h"
//===========================================================================//

using namespace std;

//===========================================================================//
/*!
 * \class MMKP_DecompParam
 * Storage for parameters for the 
 *   Multi-Dimensional Multi-Choice Knapsack Problem (MMKP).
 */

//===========================================================================//
class MMKP_Param {
public:
   int    LogLevel;           //application log level
   string DataDir;            //data directory
   string Instance;           //name of instance
   string DataFormat;         //format of data
   string ModelNameCore;      //name of model core
   string ModelNameRelax;     //name of model relax
   string ModelNameRelaxNest; //name of nested model relax
   bool   UsePisinger;

public:
   void getSettings(UtilParameters & utilParam){
      static const char * common = "MMKP";
      LogLevel       = utilParam.GetSetting("LogLevel",        0,      common);
      DataDir        = utilParam.GetSetting("DataDir",        "",      common);
      Instance       = utilParam.GetSetting("Instance",       "",      common);
      DataFormat     = utilParam.GetSetting("DataFormat",     "",      common);
      ModelNameCore  = utilParam.GetSetting("ModelNameCore",  "MDKP0", common);
      ModelNameRelax = utilParam.GetSetting("ModelNameRelax", "MCKP0", common);
      ModelNameRelaxNest 
         = utilParam.GetSetting("ModelNameRelaxNest", "", common);
      UsePisinger    = utilParam.GetSetting("UsePisinger", true,  common);
      if(!checkOptions())
         throw UtilException("Bad Parameter", "getSettings", "MMKP_Param");
   }

   bool checkOptions(){
      bool optionsOk = true;
      if(!(ModelNameCore == "MDKP0" || 
	   ModelNameCore == "MCP"   ||
	   ModelNameCore == "MMKPHalf")){
         cerr << "Error: Parameter ModelNameCore = " << ModelNameCore
              << " is not a defined model choice." << endl;
         optionsOk = false;
      }
      if(!(ModelNameRelax == "MCKP0" || 
	   ModelNameRelax == "MDKP"  ||
	   ModelNameRelax == "MDKPHalf")){
         cerr << "Error: Parameter ModelNameRelax = " << ModelNameRelax
              << " is not a defined model choice." << endl;
         optionsOk = false;
      }
      if(!(ModelNameRelaxNest == "MC2KP0" || 
	   ModelNameRelaxNest == "MDKP"   ||
	   ModelNameRelaxNest == "MMKP"   ||
	   ModelNameRelaxNest == "")){
         cerr << "Error: Parameter ModelNameRelaxNest = " 
              << ModelNameRelaxNest
              << " is not a defined model choice." << endl;
         optionsOk = false;
      }
      //TODO: check for bad combos too
      return optionsOk;
   }

   void dumpSettings(ostream * os = &cout){
      static const char * common = "MMKP";
      (*os) << "\n=====================================================\n"
            << "MMKP_DECOMP PARAMETER SETTINGS \n";
      (*os) << common << ": LogLevel          : " << LogLevel      << endl;
      (*os) << common << ": DataDir           : " << DataDir       << endl;
      (*os) << common << ": Instance          : " << Instance      << endl;
      (*os) << common << ": DataFormat        : " << DataFormat    << endl;
      (*os) << common << ": ModelNameCore     : " << ModelNameCore     << endl;
      (*os) << common << ": ModelNameRelax    : " << ModelNameRelax    << endl;
      (*os) << common << ": ModelNameRelaxNest: " << ModelNameRelaxNest<< endl;
      (*os) << common << ": UsePisinger       : " << UsePisinger       << endl;
      (*os) <<   "=====================================================\n";
   }
   
public:
   MMKP_Param() :
      LogLevel          (0),
      DataDir           (""),
      Instance          (""),                   
      DataFormat        ("hifi"),//hifi, khan, or gsimon
      ModelNameCore     ("MDKP0"),
      ModelNameRelax    ("MCKP0"),
      ModelNameRelaxNest(""     ),
      UsePisinger       (true   ){}
   ~MMKP_Param() {};
};

#endif
