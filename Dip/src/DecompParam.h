//===========================================================================//
// This file is part of the DIP Solver Framework.                            //
//                                                                           //
// DIP is distributed under the Eclipse Public License as part of the        //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)       //
//                                                                           //
// Conceptual Design: Matthew Galati, SAS Institute Inc.                     //
//                    Ted Ralphs, Lehigh University                          //
//                                                                           //
// Copyright (C) 2002-2009, Lehigh University, Matthew Galati, Ted Ralphs    //
// All Rights Reserved.                                                      //
//===========================================================================//

#ifndef DECOMP_PARAM_INCLUDED
#define DECOMP_PARAM_INCLUDED

//===========================================================================//
#include "Decomp.h"
#include "UtilMacros.h"
#include "UtilParameters.h"

//===========================================================================//
#define PARAM_getSetting(xstr, x) x = param.GetSetting(xstr, x, sec)

//===========================================================================//
class DecompParam{


   //----------------------------------------------------------------------//
   /**
    * @name Data.
    * @{
    */
   //----------------------------------------------------------------------//
public:   
   int    LogLevel;
   int    LogDebugLevel;
   int    LogLpLevel;    
   //=0 never
   //=1 only on error
   //=2 dump every model
   int    LogDumpModel;
   int    LimitInitVars; 

   int    DebugLevel;//=0 (default), =1 (extra checks on duals, etc)

   double TolZero;
   int    LimitTotalCutIters;
   int    LimitTotalPriceIters;
   int    LimitRoundCutIters;
   int    LimitRoundPriceIters;
   double LimitTime;


   //---
   //--- tailing off when average bound over TailoffLength iterations 
   //--- has changed less than TailoffPercent
   //---
   int    TailoffLength;
   double TailoffPercent;
   double MasterGapLimit;

   int    CompressColumns;
   int    CompressColumnsIterFreq;     //num iters between compress
   double CompressColumnsSizeMultLimit;//don't compress unless number of cols increased by this mult
   int    CutDC;
   int    CutCGL;
   
   int    CutCglKnapC;
   int    CutCglFlowC;
   int    CutCglMir;
   int    CutCglClique;
   int    CutCglOddHole;
   int    CutCglGomory;

   int    SubProbUseCutoff;

   double SubProbGapLimitExact;
   double SubProbGapLimitInexact;
   double SubProbTimeLimitExact;
   double SubProbTimeLimitInexact;
   int    SubProbNumThreads;
   
   //This option only works with Cpx:
   // DecompDualSimplex = 0,
   // DecompPrimSimplex = 1,
   // DecompBarrier     = 2
   int    SubProbSolverStartAlgo;

   //n = 0: do all blocks each time
   //n > 0: do all blocks every n*numBlocks iterations
   int    RoundRobinInterval;

   //TODO: named values in parameters?
   //0:RoundRobinRotate:    rotate through blocks in order 0...numBlocks-1
   //1:RoundRobinMostNegRC: choose the block with most neg reduced cost (in last iter)
   int    RoundRobinStrategy;

   //solve master as IP at end of each node (this should only be done
   //  if there are more than one blocks)
   //TODO: how often? after every pass?

   int    SolveMasterAsIp;         //{0,1}
   int    SolveMasterAsIpFreqNode; //solve every n nodes
   int    SolveMasterAsIpFreqPass; //solve every n passes (within one node)
   double SolveMasterAsIpLimitTime;
   double SolveMasterAsIpLimitGap;
   
   //0 = If a user function is defined, it will use the user function.
   //    If the user returns an exact solution, it will not run the built-in
   //    IP solve (default).
   //    If a user function is not defined, it will use the built-in IP solve.
   //1 = Use the built-in IP solve, even if there is a user defines a function.
   //2 = Calls the user defined function (if exists) and then calls built-in 
   //    IP solver (use this for debugging).
   int    SolveRelaxAsIp;

   int    InitVarsWithCutDC;
   int    InitVarsWithIP;
   int    InitVarsWithIPLimitTime;
   
   //solve compact formulation first before starting PhaseI
   //  hopefully identify infeasibiity in tree quicker
   int    InitCompactSolve;

   int    DualStab;
   double DualStabAlpha;

   int    BreakOutPartial;
   
   //when solving using IP solver, algorithm for initial relaxation
   //when solving using IP solver, algorithm for subproblems
   //  options= dual, primal, barrier
   //string IpAlgoStart;
   //string IpAlgoSub;

   int    BranchEnforceInSubProb;
   int    BranchEnforceInMaster;
   int    MasterConvexityLessThan; //0='E', 1='L'

   /**
    * @}
    */

   
   //-----------------------------------------------------------------------//
   /**
    * @name Helper functions.
    * @{
    */
   //-----------------------------------------------------------------------//
public:
   void getSettingsImpl(UtilParameters & param,
			const char     * sec){
      /** \todo: think about putting these into sections of structs */
      PARAM_getSetting("LogLevel",             LogLevel);      
      PARAM_getSetting("LogDebugLevel",        LogDebugLevel);
      PARAM_getSetting("LogLpLevel",           LogLpLevel);
      PARAM_getSetting("LogDumpModel",         LogDumpModel);
      PARAM_getSetting("LimitInitVars",        LimitInitVars);
      PARAM_getSetting("DebugLevel",           DebugLevel);
      PARAM_getSetting("TolZero",              TolZero);
      PARAM_getSetting("LimitTotalCutIters",   LimitTotalCutIters);
      PARAM_getSetting("LimitTotalPriceIters", LimitTotalPriceIters);
      PARAM_getSetting("LimitRoundCutIters",   LimitRoundCutIters);
      PARAM_getSetting("LimitRoundPriceIters", LimitRoundPriceIters);
      PARAM_getSetting("LimitTime",            LimitTime);
      
      PARAM_getSetting("TailoffLength",        TailoffLength);
      PARAM_getSetting("TailoffPercent",       TailoffPercent);       
      PARAM_getSetting("MasterGapLimit",       MasterGapLimit);
      PARAM_getSetting("CompressColumns",      CompressColumns);       
      PARAM_getSetting("CompressColumnsIterFreq",      CompressColumnsIterFreq);       
      PARAM_getSetting("CompressColumnsSizeMultLimit", CompressColumnsSizeMultLimit);       
      PARAM_getSetting("CutDC",                CutDC);
      PARAM_getSetting("CutCGL",               CutCGL);
      PARAM_getSetting("CutCglKnapC",          CutCglKnapC);
      PARAM_getSetting("CutCglFlowC",          CutCglFlowC);
      PARAM_getSetting("CutCglMir",            CutCglMir);
      PARAM_getSetting("CutCglClique",         CutCglClique);
      PARAM_getSetting("CutCglOddHole",        CutCglOddHole);
      PARAM_getSetting("CutCglGomory",         CutCglGomory);
      PARAM_getSetting("SubProbUseCutoff",     SubProbUseCutoff);
      PARAM_getSetting("SubProbGapLimitExact", SubProbGapLimitExact);
      PARAM_getSetting("SubProbGapLimitInexact",SubProbGapLimitInexact);
      PARAM_getSetting("SubProbTimeLimitExact",  SubProbTimeLimitExact);
      PARAM_getSetting("SubProbTimeLimitInexact",SubProbTimeLimitInexact);
      PARAM_getSetting("SubProbNumThreads",      SubProbNumThreads);
      PARAM_getSetting("SubProbSolverStartAlgo",SubProbSolverStartAlgo);
      PARAM_getSetting("RoundRobinInterval",   RoundRobinInterval);
      PARAM_getSetting("RoundRobinStrategy",   RoundRobinStrategy);
      PARAM_getSetting("SolveMasterAsIp",      SolveMasterAsIp);
      PARAM_getSetting("SolveMasterAsIpFreqNode",SolveMasterAsIpFreqNode);
      PARAM_getSetting("SolveMasterAsIpFreqPass",SolveMasterAsIpFreqPass);
      PARAM_getSetting("SolveMasterAsIpLimitTime", SolveMasterAsIpLimitTime);
      PARAM_getSetting("SolveMasterAsIpLimitGap",  SolveMasterAsIpLimitGap);
      PARAM_getSetting("SolveRelaxAsIp",       SolveRelaxAsIp);
      PARAM_getSetting("InitVarsWithCutDC",    InitVarsWithCutDC);
      PARAM_getSetting("InitVarsWithIP",       InitVarsWithIP);
      PARAM_getSetting("InitVarsWithIPLimitTime", InitVarsWithIPLimitTime);
      PARAM_getSetting("InitCompactSolve",     InitCompactSolve);
      PARAM_getSetting("DualStab",             DualStab);
      PARAM_getSetting("DualStabAlpha",        DualStabAlpha);
      PARAM_getSetting("BreakOutPartial",      BreakOutPartial);
      PARAM_getSetting("BranchEnforceInSubProb",  BranchEnforceInSubProb);
      PARAM_getSetting("BranchEnforceInMaster",   BranchEnforceInMaster);
      PARAM_getSetting("MasterConvexityLessThan", MasterConvexityLessThan);
   }

   inline void getSettings(UtilParameters & param){
      const string sec = "DECOMP";
      getSettingsImpl(param, sec.c_str());
   }
   
   inline void getSettings(UtilParameters & param,
			   const string   & sec){
      //---
      //--- first get any settings that apply across any DECOMP algo
      //---    from the [DECOMP] section
      //---
      getSettingsImpl(param, "DECOMP");
      
      //---
      //--- then get any settings that apply to this algo
      //---    from the [<sec>] section
      //---
      getSettingsImpl(param, sec.c_str());
      
   }

   /** \todo this should be derived from utilparam? */
   void dumpSettings(const string & sec,
		     ostream      * os = &cout){
      (*os) << "\n========================================================";
      (*os) << "\nDECOMP PARAMETER SETTINGS\n";
      UtilPrintParameter(os, sec, "LogLevel",            LogLevel);
      UtilPrintParameter(os, sec, "LogDebugLevel",       LogDebugLevel);
      UtilPrintParameter(os, sec, "LogLpLevel",          LogLpLevel);
      UtilPrintParameter(os, sec, "LogDumpModel",        LogDumpModel);
      UtilPrintParameter(os, sec, "LimitInitVars",       LimitInitVars);
      UtilPrintParameter(os, sec, "DebugLevel",          DebugLevel);
      UtilPrintParameter(os, sec, "TolZero",             TolZero);
      UtilPrintParameter(os, sec, "LimitTotalCutIters",  LimitTotalCutIters);
      UtilPrintParameter(os, sec, "LimitTotalPriceIters",LimitTotalPriceIters);
      UtilPrintParameter(os, sec, "LimitRoundCutIters",  LimitRoundCutIters);
      UtilPrintParameter(os, sec, "LimitRoundPriceIters",LimitRoundPriceIters);
      UtilPrintParameter(os, sec, "LimitTime",           LimitTime);
      
      UtilPrintParameter(os, sec, "TailoffLength",       TailoffLength);
      UtilPrintParameter(os, sec, "TailoffPercent",      TailoffPercent);      
      UtilPrintParameter(os, sec, "MasterGapLimit",      MasterGapLimit);      
      UtilPrintParameter(os, sec, "CompressColumns",     CompressColumns);
      UtilPrintParameter(os, sec, "CompressColumnsIterFreq",      CompressColumnsIterFreq);
      UtilPrintParameter(os, sec, "CompressColumnsSizeMultLimit", CompressColumnsSizeMultLimit);
      UtilPrintParameter(os, sec, "CutDC",               CutDC);
      UtilPrintParameter(os, sec, "CutCGL",              CutCGL);
      UtilPrintParameter(os, sec, "CutCglKnapC",         CutCglKnapC);
      UtilPrintParameter(os, sec, "CutCglFlowC",         CutCglFlowC);
      UtilPrintParameter(os, sec, "CutCglMir",           CutCglMir);
      UtilPrintParameter(os, sec, "CutCglClique",        CutCglClique);
      UtilPrintParameter(os, sec, "CutCglOddHole",       CutCglOddHole);
      UtilPrintParameter(os, sec, "CutCglGomory",        CutCglGomory);
      UtilPrintParameter(os, sec, "SubProbUseCutoff",    SubProbUseCutoff);
      UtilPrintParameter(os, sec, "SubProbGapLimitExact",   
                         SubProbGapLimitExact);
      UtilPrintParameter(os, sec, "SubProbGapLimitInexact", 
                         SubProbGapLimitInexact);
      UtilPrintParameter(os, sec, "SubProbTimeLimitExact",   
                         SubProbTimeLimitExact);
      UtilPrintParameter(os, sec, "SubProbTimeLimitInexact", 
                         SubProbTimeLimitInexact);
      UtilPrintParameter(os, sec, "SubProbNumThreads", SubProbNumThreads);
      UtilPrintParameter(os, sec, "SubProbSolverStartAlgo",
                         SubProbSolverStartAlgo);
      UtilPrintParameter(os, sec, "RoundRobinInterval",  RoundRobinInterval);
      UtilPrintParameter(os, sec, "RoundRobinStrategy",  RoundRobinStrategy);  
      UtilPrintParameter(os, sec, "SolveMasterAsIp",     SolveMasterAsIp);
      UtilPrintParameter(os, sec, "SolveMasterAsIpFreqNode",  SolveMasterAsIpFreqNode);
      UtilPrintParameter(os, sec, "SolveMasterAsIpFreqPass",  SolveMasterAsIpFreqPass);
      UtilPrintParameter(os, sec, "SolveMasterAsIpLimitTime", SolveMasterAsIpLimitTime);
      UtilPrintParameter(os, sec, "SolveMasterAsIpLimitGap",  SolveMasterAsIpLimitGap);
      UtilPrintParameter(os, sec, "SolveRelaxAsIp",     SolveRelaxAsIp);
      UtilPrintParameter(os, sec, "InitVarsWithCutDC",   InitVarsWithCutDC);
      UtilPrintParameter(os, sec, "InitVarsWithIP",   InitVarsWithIP);
      UtilPrintParameter(os, sec, "InitVarsWithIPLimitTime",   InitVarsWithIPLimitTime);
      UtilPrintParameter(os, sec, "InitCompactSolve",  InitCompactSolve);
      UtilPrintParameter(os, sec, "DualStab",          DualStab);
      UtilPrintParameter(os, sec, "DualStabAlpha",     DualStabAlpha);
      UtilPrintParameter(os, sec, "BreakOutPartial",   BreakOutPartial);
      UtilPrintParameter(os, sec, "BranchEnforceInSubProb",   BranchEnforceInSubProb);
      UtilPrintParameter(os, sec, "BranchEnforceInMaster",    BranchEnforceInMaster);
      UtilPrintParameter(os, sec, "MasterConvexityLessThan", MasterConvexityLessThan);
      (*os) << "========================================================\n";
   }

   void setDefaults(){
      LogLevel             = 0;
      LogDebugLevel        = 0;
      LogLpLevel           = 0;
      LogDumpModel         = 0;
      LimitInitVars        = 5;
      DebugLevel           = 0;
      TolZero              = DecompEpsilon;
      LimitTotalCutIters   = COIN_INT_MAX;
      LimitTotalPriceIters = COIN_INT_MAX;
      LimitRoundCutIters   = COIN_INT_MAX;
      LimitRoundPriceIters = COIN_INT_MAX;
      LimitTime            = COIN_DBL_MAX;     
      TailoffLength        = 10;
      TailoffPercent       = 0.10;
      MasterGapLimit       = 0.01;
      CompressColumns      = 1;
      CompressColumnsIterFreq      = 2;
      CompressColumnsSizeMultLimit = 1.20;
      CutDC                = 0;
      CutCGL               = 1;
      CutCglKnapC          = 1;
      CutCglFlowC          = 1;
      CutCglMir            = 1;
      CutCglClique         = 1;
      CutCglOddHole        = 1;
      CutCglGomory         = 1;
      SubProbUseCutoff     = 0;
      SubProbGapLimitExact   = 0.0001; // 0.01% gap
      SubProbGapLimitInexact = 0.1;    //10.00% gap
      SubProbTimeLimitExact   = COIN_DBL_MAX;
      SubProbTimeLimitInexact = COIN_DBL_MAX;
      SubProbNumThreads       = 1;
      SubProbSolverStartAlgo = DecompDualSimplex;
      RoundRobinInterval   = 0;
      RoundRobinStrategy   = RoundRobinRotate;
      SolveMasterAsIp          = 1;//TODO: turn off if one block
      SolveMasterAsIpFreqNode  = 1;
      SolveMasterAsIpFreqPass  = 1000;
      SolveMasterAsIpLimitTime = 30;
      SolveMasterAsIpLimitGap  = 0.05; //5% gap
      SolveRelaxAsIp           = 0;
      InitVarsWithCutDC        = 0;
      InitVarsWithIP           = 0;
      InitVarsWithIPLimitTime  = 10;
      InitCompactSolve         = 0;
      DualStab                 = 0;
      DualStabAlpha            = 0.10;
      BreakOutPartial          = 0;
      BranchEnforceInSubProb   = 0;
      BranchEnforceInMaster    = 1;
      MasterConvexityLessThan  = 0;
   }
   
   void dumpSettings(ostream * os = &cout){
      const string sec = "DECOMP";
      dumpSettings(sec, os);
   }
   /**
    * @}
    */

      
   //-----------------------------------------------------------------------//
   /**
    * @name Constructors and destructor.
    * @{
    */
   //-----------------------------------------------------------------------//
public:
   /**
    * Default constructors.
    */   
   DecompParam() {
      setDefaults();
   }
   
   /**
    * Destructor
    */   
   ~DecompParam() {}   
   /**
    * @}
    */
};

#endif
