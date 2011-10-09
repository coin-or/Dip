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

#ifndef VRP_PARAM_INCLUDED
#define VRP_PARAM_INCLUDED

// --------------------------------------------------------------------- //
#include "UtilParameters.h"

// --------------------------------------------------------------------- //
/*!
 * \class VRP_Param
 * Storage for parameters for the Vehicle Routing Problem (VRP).
 * 
 */

// --------------------------------------------------------------------- //
class VRP_Param{  
public:
   int    LogLevel;           //application log level
   string DataDir;            //data directory
   string Instance;           //name of instance
   int    NumRoutes;          //number of routes
   string ModelNameCore;      //name of model core
   string ModelNameRelax;     //name of model relax
   string ModelNameRelaxNest; //name of nested model relax
   
public:
   void getSettings(UtilParameters & utilParam){
      static const char * common = "VRP";
      LogLevel     = utilParam.GetSetting("LogLevel",      0,  common);
      DataDir      = utilParam.GetSetting("DataDir",      "",  common);
      Instance     = utilParam.GetSetting("Instance",     "",  common);
      NumRoutes    = utilParam.GetSetting("NumRoutes",     1,  common);    
      ModelNameCore  
         = utilParam.GetSetting("ModelNameCore",  "", common);
      ModelNameRelax 
         = utilParam.GetSetting("ModelNameRelax", "", common);
      ModelNameRelaxNest 
         = utilParam.GetSetting("ModelNameRelaxNest", "", common);
   }

   void dumpSettings(ostream * os = &cout){
      static const char * common = "VRP";
      (*os) << "\n=====================================================\n"
            << "VRP_DECOMP PARAMETER SETTINGS \n";
      (*os) << common << ": LogLevel          : " << LogLevel          << endl;
      (*os) << common << ": DataDir           : " << DataDir           << endl;
      (*os) << common << ": Instance          : " << Instance          << endl;
      (*os) << common << ": NumRoutes         : " << NumRoutes         << endl;
      (*os) << common << ": ModelNameCore     : " << ModelNameCore     << endl;
      (*os) << common << ": ModelNameRelax    : " << ModelNameRelax    << endl;
      (*os) << common << ": ModelNameRelaxNest: " << ModelNameRelaxNest<< endl;
      (*os) <<   "=====================================================\n";
   }

public:
   VRP_Param():
      LogLevel          (0 ),
      DataDir           (""),
      Instance          (""),
      NumRoutes         (0 ),
      ModelNameCore     (""),
      ModelNameRelax    (""),
      ModelNameRelaxNest("")      
   {}
   ~VRP_Param() {};
};

#endif
