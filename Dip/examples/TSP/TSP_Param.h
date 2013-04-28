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

#ifndef TSP_PARAM_INCLUDED
#define TSP_PARAM_INCLUDED

// --------------------------------------------------------------------- //
#include "UtilParameters.h"

// --------------------------------------------------------------------- //
/*!
 * \class TSP_Param
 * Storage for parameters for the Traveling Salesman Problem (TSP).
 * 
 * \todo think about this design, register parameters, isoptional
 *       combos that violate, throw exceptions
 *
 */

// --------------------------------------------------------------------- //
class TSP_Param{  
public:
   int    LogLevel;
   string DataDir;
   string Instance;
   int    CutSubtoursX;
   //int    CutBlossomsX;
   //int    CutCombsX;
   string ModelNameCore;      //name of model core
   string ModelNameRelax;     //name of model relax

   
public:
   void getSettings(UtilParameters & utilParam){
      static const char * common = "TSP";
      LogLevel       = utilParam.GetSetting("LogLevel",        0, common);
      DataDir        = utilParam.GetSetting("DataDir",        "", common);
      Instance       = utilParam.GetSetting("Instance",       "", common);    
      CutSubtoursX   = utilParam.GetSetting("CutSubtoursX",    1, common);    
      //CutBlossomsX   = utilParam.GetSetting("CutBlossomsX",    1, common); 
      //CutCombsX      = utilParam.GetSetting("CutCombsX",       1, common);
      ModelNameCore  = utilParam.GetSetting("ModelNameCore",  "2MATCH", 
					    common);
      ModelNameRelax = utilParam.GetSetting("ModelNameRelax", "SUBTOUR", 
					    common);
   }

   void dumpSettings(ostream * os = &cout){
      static const char * common = "TSP";
      (*os) << "\n=====================================================\n"
            << "TSP_DECOMP PARAMETER SETTINGS \n";
      (*os) << common << ": LogLevel      : " << LogLevel       << endl;
      (*os) << common << ": DataDir       : " << DataDir        << endl;
      (*os) << common << ": Instance      : " << Instance       << endl;
      (*os) << common << ": CutSubtoursX  : " << CutSubtoursX   << endl;
      //(*os) << common << ": CutBlossomsX  : " << CutBlossomsX   << endl;
      //(*os) << common << ": CutCombsX     : " << CutCombsX      << endl;
      (*os) << common << ": ModelNameCore : " << ModelNameCore  << endl;
      (*os) << common << ": ModelNameRelax: " << ModelNameRelax << endl;
      (*os) <<   "=====================================================\n";
   }

public:
   TSP_Param():
      LogLevel      (0  ),
      DataDir       ("."),
      Instance      ("" ),
      CutSubtoursX  (1  ),
      //CutBlossomsX  (0  ),
      //CutCombsX     (0  ),   
      ModelNameCore ("2MATCH"),
      ModelNameRelax("SUBTOUR")

   {}
   ~TSP_Param() {};
};

#endif
