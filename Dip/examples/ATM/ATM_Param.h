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

#ifndef ATM_PARAM_INCLUDED
#define ATM_PARAM_INCLUDED

//===========================================================================//
#include "UtilParameters.h"

using namespace std;

//===========================================================================//
/*!
 * \class ATM_Param
 * Storage for parameters for the 
 *     ATM Cash Management Problem (ATM).
 */

//===========================================================================//
class ATM_Param {
public:
   int    LogLevel;          //application log level
   string DataDir;           //data directory
   string DataAtm;           //data file (atms)
   string DataDate;          //data file (dates)
   string DataAtmDate;       //data file (atms x dates)
   int    NumSteps;          //number of steps for discretization of NLP
   bool   UseTightModel;     //use tighter formulation of z=xy
   string ModelNameCore;     //name of model core
   string ModelNameRelax;    //name of model relax
   string ModelNameRelaxNest;//name of nested model relax

public:
   void getSettings(UtilParameters & utilParam){
      static const char * common = "ATM";
      LogLevel       = utilParam.GetSetting("LogLevel",         0, common);
      DataDir        = utilParam.GetSetting("DataDir",         "", common);
      DataAtm        = utilParam.GetSetting("DataAtm",         "", common);    
      DataDate       = utilParam.GetSetting("DataDate",        "", common);    
      DataAtmDate    = utilParam.GetSetting("DataAtmDate",     "", common);    
      NumSteps       = utilParam.GetSetting("NumSteps",        10, common);
      UseTightModel  = utilParam.GetSetting("UseTightModel", true, common);
      ModelNameCore  = utilParam.GetSetting("ModelNameCore",   "", common);
      ModelNameRelax = utilParam.GetSetting("ModelNameRelax",  "", common);
      ModelNameRelaxNest 
	 = utilParam.GetSetting("ModelNameRelaxNest", "", common);
   }
   
   void dumpSettings(ostream * os = &cout){
      static const char * common = "ATM";
      (*os) << "\n=====================================================\n"
            << "ATM_DECOMP PARAMETER SETTINGS \n";
      (*os) << common << ": LogLevel          : " << LogLevel          << endl;
      (*os) << common << ": DataDir           : " << DataDir           << endl;
      (*os) << common << ": DataAtm           : " << DataAtm           << endl;
      (*os) << common << ": DataDate          : " << DataDate          << endl;
      (*os) << common << ": DataAtmDate       : " << DataAtmDate       << endl;
      (*os) << common << ": NumSteps          : " << NumSteps          << endl;
      (*os) << common << ": UseTightModel     : " << UseTightModel     << endl;
      (*os) << common << ": ModelNameCore     : " << ModelNameCore     << endl;
      (*os) << common << ": ModelNameRelax    : " << ModelNameRelax    << endl;
      (*os) << common << ": ModelNameRelaxNest: " << ModelNameRelaxNest<< endl;
      (*os) <<   "=====================================================\n";
   }
   
public:
   ATM_Param():    
      LogLevel          (0),
      DataDir           (""),
      DataAtm           (""),
      DataDate          (""),
      DataAtmDate       (""),
      NumSteps          (0),
      UseTightModel     (true),
      ModelNameCore     (""),
      ModelNameRelax    (""),
      ModelNameRelaxNest(""){};
   ~ATM_Param() {};
};

#endif
