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

#ifndef TSP_DECOMP_PARAM_INCLUDED
#define TSP_DECOMP_PARAM_INCLUDED

// --------------------------------------------------------------------- //
#include "UtilParameters.h"

// --------------------------------------------------------------------- //
/*!
 * \class TSP_DecompParam
 * Storage for parameters for the Traveling Salesman Problem (TSP).
 * 
 * \todo think about this design, register parameters, isoptional
 *       combos that violate, throw exceptions
 *
 */

// --------------------------------------------------------------------- //
class TSP_DecompParam{  
public:
   string DataDir;
   string Instance;
   int    CutSubtoursX;
   int    CutBlossomsX;
   int    CutCombsX;
   string ModelNameCore;      //name of model core
   string ModelNameRelax;     //name of model relax

   
public:
   void getSettings(UtilParameters & utilParam){
      static const char * common = "TSP";
      DataDir        = utilParam.GetSetting("DataDir",        "", common);
      Instance       = utilParam.GetSetting("Instance",       "", common);    
      ModelNameCore  = utilParam.GetSetting("ModelNameCore",  "", common);
      ModelNameRelax = utilParam.GetSetting("ModelNameRelax", "", common);
      CutSubtoursX   = utilParam.GetSetting("CutSubtoursX",    1, common);    
      CutBlossomsX   = utilParam.GetSetting("CutBlossomsX",    1, common);    
      CutCombsX      = utilParam.GetSetting("CutCombsX",       1, common);    
   }

   void dumpSettings(ostream * os = &cout){
      static const char * common = "TSP";
      (*os) << "\n=====================================================\n"
            << "TSP_DECOMP PARAMETER SETTINGS \n";
      (*os) << common << ": DataDir       : " << DataDir        << endl;
      (*os) << common << ": Instance      : " << Instance       << endl;
      (*os) << common << ": ModelNameCore : " << ModelNameCore  << endl;
      (*os) << common << ": ModelNameRelax: " << ModelNameRelax << endl;
      (*os) << common << ": CutSubtoursX  : " << CutSubtoursX   << endl;
      (*os) << common << ": CutBlossomsX  : " << CutBlossomsX   << endl;
      (*os) << common << ": CutCombsX     : " << CutCombsX      << endl;
      (*os) <<   "=====================================================\n";
   }

public:
   TSP_DecompParam():
      DataDir       ("."),
      Instance      (""),
      ModelNameCore ("MDKP0"),
      ModelNameRelax("MCKP0"),
      CutSubtoursX  (1),
      CutBlossomsX  (0),
      CutCombsX     (0)    
   {}
   ~TSP_DecompParam() {};
};

#endif
