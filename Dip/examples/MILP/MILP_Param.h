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

#ifndef MILP_PARAM_INCLUDED
#define MILP_PARAM_INCLUDED

//===========================================================================//
#include "UtilParameters.h"

using namespace std;

//===========================================================================//
class MILP_Param{
public:
   int    LogLevel;   
   int    RandomSeed;
   double RelaxPercent;
   double BestKnownLB;
   double BestKnownUB;
   string DataDir;
   string Instance;   

public:
   void getSettings(UtilParameters & utilParam){
      static const char * common = "MILP";
      LogLevel     = utilParam.GetSetting("LogLevel",     0,       common);
      RandomSeed   = utilParam.GetSetting("RandomSeed",   1,       common);
      RelaxPercent = utilParam.GetSetting("RelaxPercent", 0.333,   common);
      BestKnownLB  = utilParam.GetSetting("BestKnownLB",  -1.e100, common);
      BestKnownUB  = utilParam.GetSetting("BestKnownUB",   1.e100, common);
      DataDir      = utilParam.GetSetting("DataDir",      "",      common);
      Instance     = utilParam.GetSetting("Instance",     "",      common);    
   }

   void dumpSettings(ostream * os = &cout){
      static const char * common = "MILP";
      (*os) << "\n=====================================================\n"
            << "MILP_DECOMP PARAMETER SETTINGS \n";
      (*os) << common << ": LogLevel    : " << LogLevel     << endl;
      (*os) << common << ": RandomSeed  : " << RandomSeed   << endl;
      (*os) << common << ": RelaxPercent: " << RelaxPercent << endl;
      (*os) << common << ": BestKnownLB : " << BestKnownLB  << endl;
      (*os) << common << ": BestKnownUB : " << BestKnownUB  << endl;
      (*os) << common << ": DataDir     : " << DataDir      << endl;
      (*os) << common << ": Instance    : " << Instance     << endl;
      (*os) << "\n=====================================================\n";
   }
   
public:
   MILP_Param():    
      LogLevel    (0     ),
      RandomSeed  (1     ),
      RelaxPercent(0.333 ),      
      BestKnownLB (-1e100),
      BestKnownUB ( 1e100),
      DataDir     (""    ),
      Instance    (""    ) {};
   ~MILP_Param() {};
};

#endif
