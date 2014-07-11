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

#ifndef DECOMP_PARAM_INCLUDED
#define DECOMP_PARAM_INCLUDED

#include "UtilParameters.h"
#include "DecompConstants.h"

// --------------------------------------------------------------------- //
class DecompParam {
private:
   DecompParam(const DecompParam&);
   DecompParam& operator=(const DecompParam&);

public:

   int    LogLevel;
   int    LogAppLevel;
   int    LogDebugLevel;
   int    LogLpLevel;    //name? inner solver
   unsigned int    LimitInitVars; //?? specific to PC? make own section ??
   double TolZero;
   int    LimitTotalCutIters;
   int    LimitTotalPriceIters;
   int    LimitRoundCutIters;
   int    LimitRoundPriceIters;
   double LimitTime;
   int    PriceMultiPoly;
   int    CutDC;
   int    CutCGL;
   //subsection?
   int    CutCglKnapC;
   int    CutCglFlowC;
   int    CutCglMir;
   int    CutCglClique;

public:
   void getSettings(UtilParameters& utilParam) {
      static const char* common = "DECOMP";
      LogLevel      = utilParam.GetSetting("LogLevel",      0,     common);
      LogAppLevel   = utilParam.GetSetting("LogAppLevel",   0,     common);
      LogDebugLevel = utilParam.GetSetting("LogDebugLevel", 0,     common);
      LogLpLevel    = utilParam.GetSetting("LogLpLevel",    0,     common);
      LimitInitVars = utilParam.GetSetting("LimitInitVars", 1,     common);
      TolZero       = utilParam.GetSetting("TolZero",
                                           DecompEpsilon, common);
      LimitTotalCutIters  = utilParam.GetSetting("LimitTotalCutIters",
                            2000, common);
      LimitTotalPriceIters = utilParam.GetSetting("LimitTotalPriceIters",
                             2000, common);
      LimitRoundCutIters  = utilParam.GetSetting("LimitRoundCutIters",
                            2000, common);
      LimitRoundPriceIters = utilParam.GetSetting("LimitRoundPriceIters",
                             2000, common);
      LimitTime           = utilParam.GetSetting("LimitTime",
                            600, common);
      //TODO: what if we want multi-poly on just 1st and 3rd - TODO
      PriceMultiPoly       = utilParam.GetSetting("PriceMultiPoly",
                             0,     common);
      CutDC                = utilParam.GetSetting("CutDC",
                             0,     common);
      CutCGL                = utilParam.GetSetting("CutCGL",
                              0,     common);
      CutCglKnapC                = utilParam.GetSetting("CutCglKnapC",
                                   0,     common);
      CutCglFlowC                = utilParam.GetSetting("CutCglFlowC",
                                   0,     common);
      CutCglMir                = utilParam.GetSetting("CutCglMir",
                                 0,     common);
      CutCglClique                = utilParam.GetSetting("CutCglClique",
                                    0,     common);
   }

   //this should be a parameter method, should parameter be an object?
   //have user register parameters, so can set usage too
   void dumpSettings(ostream* os = &cout) {
      static const char* common = "DECOMP";
      (*os) << "\n========================================================\n"
            << "DECOMP PARAMETER SETTINGS \n";
      (*os) << common << ": LogLevel      = " << LogLevel      << endl;
      (*os) << common << ": LogAppLevel   = " << LogAppLevel   << endl;
      (*os) << common << ": LogDebugLevel = " << LogDebugLevel << endl;
      (*os) << common << ": LogLpLevel    = " << LogLpLevel    << endl;
      (*os) << common << ": LimitInitVars = " << LimitInitVars << endl;
      (*os) << common << ": TolZero       = " << TolZero       << endl;
      (*os) << common << ": LimitTotalCutIters   = "
            << LimitTotalCutIters      << endl;
      (*os) << common << ": LimitTotalPriceIters = "
            << LimitTotalPriceIters    << endl;
      (*os) << common << ": LimitRoundCutIters   = "
            << LimitRoundCutIters      << endl;
      (*os) << common << ": LimitRoundPriceIters = "
            << LimitRoundPriceIters    << endl;
      (*os) << common << ": PriceMultiPoly= " << PriceMultiPoly  << endl;
      (*os) << common << ": CutDC         = " << CutDC           << endl;
      (*os) << common << ": CutCGL        = " << CutCGL          << endl;
      (*os) << common << ": CutCglKnapC   = " << CutCglKnapC     << endl;
      (*os) << common << ": CutCglFlowC   = " << CutCglFlowC     << endl;
      (*os) << common << ": CutCglMir     = " << CutCglMir       << endl;
      (*os) << common << ": CutCglClique  = " << CutCglClique    << endl;
      (*os) << "========================================================\n";
   }

public:
   DecompParam():
      LogLevel(0),
      LogAppLevel(0),
      LogDebugLevel(0),
      LogLpLevel(0),
      LimitInitVars(1),
      TolZero(DecompEpsilon),
      LimitTotalCutIters(2000),
      LimitTotalPriceIters(2000),
      LimitRoundCutIters(2000),
      LimitRoundPriceIters(2000),
      LimitTime(60),
      PriceMultiPoly(0),
      CutDC(0),
      CutCGL(0),
      CutCglKnapC(0),
      CutCglFlowC(0),
      CutCglMir(0),
      CutCglClique(0)
   {}

   ~DecompParam() {};
};

#endif
