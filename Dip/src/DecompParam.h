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
// Copyright (C) 2002-2013, Lehigh University, Matthew Galati, Ted Ralphs    //
// All Rights Reserved.                                                      //
//===========================================================================//

#ifndef DECOMP_PARAM_INCLUDED
#define DECOMP_PARAM_INCLUDED

//===========================================================================//
#include "Decomp.h"
#include "UtilMacros.h"
#include "UtilParameters.h"
#include "string"
//===========================================================================//
#define PARAM_getSetting(xstr, x) x = param.GetSetting(xstr, x, sec)

//===========================================================================//
class DecompParam {


   //----------------------------------------------------------------------//
   /**
    * @name Data.
    * @{
    */
   //----------------------------------------------------------------------//

public:
   int    LogLevel;
   int    LogDebugLevel;
   int    LogLpLevel;    //TODO: LpIpLevel separate

   //=0 never
   //=1 only on error
   //=2 dump every model

   int    LogDumpModel;

   /**
    * 0: print nothing
    * 1: print the node objective history
    */

   int    LogObjHistory;


   int    LimitInitVars;

   int    DebugLevel;//=0 (default), =1 (extra checks on duals, etc)

   double TolZero;
   int    LimitTotalCutIters;
   int    LimitTotalPriceIters;
   int    LimitRoundCutIters;
   int    LimitRoundPriceIters;
   double LimitTime;

   /**
    * Max number of nodes (copied from Alps parameters)
    */

   int    LimitNodes;

   //---
   //--- tailing off when average bound over TailoffLength iterations
   //--- has changed less than TailoffPercent
   //---

   int    TailoffLength;
   double TailoffPercent;
   double MasterGapLimit;

   //---
   //--- Strategy for switching from cutting to pricing
   //--- 0 = Default
   //--- 1 = Favor column generation
   //--- 2 = Favor cut generation

   int PCStrategy;

   int    CompressColumns;
   //num iters between compress
   int    CompressColumnsIterFreq;
   //don't compress unless number of cols increased by this mult
   double CompressColumnsSizeMultLimit;
   //do not start compression until master gap is within this limit
   double CompressColumnsMasterGapStart;
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
   int    SubProbNumSolLimit;

   //This option only works with Cpx:
   // DecompDualSimplex = 0,
   // DecompPrimSimplex = 1,
   // DecompBarrier     = 2

   int    SubProbSolverStartAlgo;

   //n = 0: do all blocks each time
   //n > 0: do all blocks every n iterations

   int    RoundRobinInterval;

   //TODO: named values in parameters?
   //NOT working
   //0:RoundRobinRotate:    rotate through blocks in order 0...numBlocks-1
   //1:RoundRobinMostNegRC: choose the block with most neg reduced cost
   //(in last iter)

   int    RoundRobinStrategy;

   //solve master as IP at end of each node (this should only be done
   //  if there are more than one blocks)
   //TODO: how often? after every pass?

   int    SolveMasterAsIp;         //{0,1}
   int    SolveMasterAsIpFreqNode; //solve every n nodes
   int    SolveMasterAsIpFreqPass; //solve every n passes (within one node)
   double SolveMasterAsIpLimitTime;
   double SolveMasterAsIpLimitGap;

   // DecompDualSimplex = 0,
   // DecompPrimSimplex = 1,
   // DecompBarrier     = 2

   int    SolveMasterUpdateAlgo;


   //0 = If a user function is defined, it will use the user function.
   //    If the user returns an exact solution, it will not run the built-in
   //    IP solve (default).
   //    If a user function is not defined, it will use the built-in IP solve.
   //1 = Use the built-in IP solve, even if there is a user defines a function.
   //2 = Calls the user defined function (if exists) and then calls built-in
   //    IP solver (use this for debugging).

   bool    SolveRelaxAsIp;

   int    InitVarsWithCutDC;
   int    InitVarsWithIP;
   int    InitVarsWithIPLimitTime;

   //solve compact formulation first before starting PhaseI
   //  hopefully identify infeasibiity in tree quicker

   int    InitCompactSolve;

   bool    DualStab;
   double DualStabAlpha;
   double DualStabAlphaOrig;

   bool    BreakOutPartial; //DISABLED for now

   //when solving using IP solver, algorithm for initial relaxation
   //when solving using IP solver, algorithm for subproblems
   //  options= dual, primal, barrier
   //string IpAlgoStart;
   //string IpAlgoSub;

   bool    BranchEnforceInSubProb;
   bool    BranchEnforceInMaster;
   int    MasterConvexityLessThan; //0='E', 1='L'
   double ParallelColsLimit;       //cosine of angle >, then consider parallel

   /**
    * Number of iterations to process in estimating bounds
    *  during strong branching.
    * CPM: this is simplex iterations of master
    * PC : this is outer price and cut iterations
    *         sets LimitTotalCutIters=LimitTotalPriceIters=BranchStrongIter
    *  THINK: or CPM could be cut passes... and solve master fully?
    *          which is expensive and clearly not standard strong branching
    */

   int    BranchStrongIter;

   /**
    * Number of threads to use in DIP.
    *
    * Currently, only used for solving the pricing problem for block
    * angular models. The subproblems (each block) are independent and
    * can be solved in parallel.
    */

   int NumThreads;

   /*
    * Check user columns for overlap. The default is true, but this
    * can be shut off to speed up the setup.
    */

   int DebugCheckBlocksColumns;

   /*
    * The block number for automatic decomposition
    */

   int NumBlocks;

   /*
    * The following parameters are extended from MILPBlock
    * applications but changed to MILP domain
    *
    */

   std::string DataDir;
   std::string Instance;
   std::string InstanceFormat;

   /*
   * The file defining which rows are in which blocks.
   */
   std::string BlockFile;

   /**
    * The format of BlockFile.
    *
    * (1) "List" or "LIST"
    * The block file defines those rows in each block.
    *   <block id>  <num rows in block>
    *   <row ids...>
    *   <block id>  <num rows in block>
    *   <row ids...>
    *
    * (2) "ZIBList" or "ZIBLIST"
    * The block file defines those rows in each block.
    *   NBLOCKS
    *   <numBlocks>
    *   BLOCK <block id>
    *   <row names...>
    *   BLOCK  <block id>
    *   <row names...>
    *
    * (3) "Pair" or "PAIR"
    * Each line is a block id to row id pair.
    *   <block id> <row id>
    *
    * (4) "PairName" or "PAIRNAME"
    * Each line is a block id to row name (matching mps) pair.
    *   <block id> <row name>
    */
   std::string BlockFileFormat;

   std::string PermuteFile;

   std::string InitSolutionFile;

   int UseNames; // col/row names for debugging
   int UseSparse; // create all blocks sparsely
   int FullModel; // create full model for CPM or direct
   double BestKnownLB;
   double BestKnownUB;
   double ColumnUB; // hack since missing extreme rays
   double ColumnLB; //hack since missing extreme rays

   int ObjectiveSense; //1=min, -1=max
   // variable indicates whether to use
   // multiple cores to compute concurrently

   bool Concurrent;

   // number of block candidates
   int NumBlocksCand;

   // time of concurrent CutOffTime to finalize
   // the choice of MILP solution method

   double ConcurrentCutOffTime;

   int ThreadIndex;

   std::string CurrentWorkingDir;

   bool SubProbParallel;

   int SubProbParallelType;

   int SubProbParallelChunksize;

   int ConcurrentThreadsNum;


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

   void getSettingsImpl(UtilParameters& param,
                        const char*      sec) {
      /** \todo: think about putting these into sections of structs */
      PARAM_getSetting("LogLevel",             LogLevel);
      PARAM_getSetting("LogDebugLevel",        LogDebugLevel);
      PARAM_getSetting("LogLpLevel",           LogLpLevel);
      PARAM_getSetting("LogDumpModel",         LogDumpModel);
      PARAM_getSetting("LogObjHistory",        LogObjHistory);
      PARAM_getSetting("LimitInitVars",        LimitInitVars);
      PARAM_getSetting("DebugLevel",           DebugLevel);
      PARAM_getSetting("TolZero",              TolZero);
      PARAM_getSetting("LimitTotalCutIters",   LimitTotalCutIters);
      PARAM_getSetting("LimitTotalPriceIters", LimitTotalPriceIters);
      PARAM_getSetting("LimitRoundCutIters",   LimitRoundCutIters);
      PARAM_getSetting("LimitRoundPriceIters", LimitRoundPriceIters);
      PARAM_getSetting("LimitTime",            LimitTime);
      PARAM_getSetting("LimitNodes",           LimitNodes);
      PARAM_getSetting("TailoffLength",        TailoffLength);
      PARAM_getSetting("TailoffPercent",       TailoffPercent);
      PARAM_getSetting("MasterGapLimit",       MasterGapLimit);
      PARAM_getSetting("PCStrategy",           PCStrategy);
      PARAM_getSetting("CompressColumns",      CompressColumns);
      PARAM_getSetting("CompressColumnsIterFreq",       CompressColumnsIterFreq);
      PARAM_getSetting("CompressColumnsSizeMultLimit",  CompressColumnsSizeMultLimit);
      PARAM_getSetting("CompressColumnsMasterGapStart", CompressColumnsMasterGapStart);
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
      PARAM_getSetting("SubProbGapLimitInexact", SubProbGapLimitInexact);
      PARAM_getSetting("SubProbTimeLimitExact",  SubProbTimeLimitExact);
      PARAM_getSetting("SubProbTimeLimitInexact", SubProbTimeLimitInexact);
      PARAM_getSetting("SubProbNumThreads",      SubProbNumThreads);
      PARAM_getSetting("SubProbNumSolLimit",     SubProbNumSolLimit);
      PARAM_getSetting("SubProbSolverStartAlgo", SubProbSolverStartAlgo);
      PARAM_getSetting("RoundRobinInterval",   RoundRobinInterval);
      PARAM_getSetting("RoundRobinStrategy",   RoundRobinStrategy);
      PARAM_getSetting("SolveMasterAsIp",      SolveMasterAsIp);
      PARAM_getSetting("SolveMasterAsIpFreqNode", SolveMasterAsIpFreqNode);
      PARAM_getSetting("SolveMasterAsIpFreqPass", SolveMasterAsIpFreqPass);
      PARAM_getSetting("SolveMasterAsIpLimitTime", SolveMasterAsIpLimitTime);
      PARAM_getSetting("SolveMasterAsIpLimitGap",  SolveMasterAsIpLimitGap);
      PARAM_getSetting("SolveMasterUpdateAlgo",    SolveMasterUpdateAlgo);
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
      PARAM_getSetting("ParallelColsLimit",       ParallelColsLimit);
      PARAM_getSetting("BranchStrongIter",        BranchStrongIter);
      PARAM_getSetting("NumThreads",              NumThreads);
      PARAM_getSetting("DebugCheckBlocksColumns", DebugCheckBlocksColumns);
      PARAM_getSetting("NumBlocks", NumBlocks);
      DataDir       = param.GetSetting("DataDir",       "",    "MILP");
      Instance      = param.GetSetting("Instance",      "",    "MILP");
      InstanceFormat = param.GetSetting("InstanceFormat", "",    "MILP");
      BlockFile     = param.GetSetting("BlockFile",     "",    "MILP");
      PermuteFile   = param.GetSetting("PermuteFile",   "",    "MILP");
      BlockFileFormat
         = param.GetSetting("BlockFileFormat",    "",    "MILP");
      InitSolutionFile
         = param.GetSetting("InitSolutionFile",   "",    "MILP");
      PARAM_getSetting("LogLevel", LogLevel);
      // PARAM_getSetting("DataDir",DataDir);
      //PARAM_getSetting("Instance",Instance);
      //PARAM_getSetting("BlockFile",BlockFile);
      //PARAM_getSetting("PermuteFile",PermuteFile);
      //PARAM_getSetting("BlockFileFormat",BlockFileFormat);
      //PARAM_getSetting("InitSolutionFile",InitSolutionFile);
      PARAM_getSetting("UseNames", UseNames);
      PARAM_getSetting("UseSparse", UseSparse);
      PARAM_getSetting("FullModel", FullModel);
      PARAM_getSetting("BestKnownLB", BestKnownLB);
      PARAM_getSetting("BestKnownUB", BestKnownUB);
      PARAM_getSetting("ColumnUB", ColumnUB);
      PARAM_getSetting("ColumnLB", ColumnLB);
      PARAM_getSetting("ObjectiveSense", ObjectiveSense);
      PARAM_getSetting("Concurrent", Concurrent);
      PARAM_getSetting("NumBlocksCand", NumBlocksCand);
      PARAM_getSetting("CconcurrentCutOffTime", ConcurrentCutOffTime);
      PARAM_getSetting("CurrentWorkingDir", CurrentWorkingDir);
      PARAM_getSetting("SubProbParallel", SubProbParallel);
      PARAM_getSetting("SubProbParallelType", SubProbParallelType);
      PARAM_getSetting("SubProbParallelChunksize", SubProbParallelChunksize);
      PARAM_getSetting("ConcurrentThreadsNum", ConcurrentThreadsNum);
      //---
      //--- store the original setting for DualStabAlpha
      //---
      DualStabAlphaOrig = DualStabAlpha;
   }

   inline void getSettings(UtilParameters& param) {
      const std::string sec = "DECOMP";
      getSettingsImpl(param, sec.c_str());
   }

   inline void getSettings(UtilParameters& param,
                           const std::string& sec) {
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
   void dumpSettings(const std::string& sec,
                     std::ostream*       os = &std::cout) {
      (*os) << "\n========================================================";
      (*os) << "\nDECOMP PARAMETER SETTINGS\n";
      UtilPrintParameter(os, sec, "LogLevel",            LogLevel);
      UtilPrintParameter(os, sec, "LogDebugLevel",       LogDebugLevel);
      UtilPrintParameter(os, sec, "LogLpLevel",          LogLpLevel);
      UtilPrintParameter(os, sec, "LogDumpModel",        LogDumpModel);
      UtilPrintParameter(os, sec, "LogObjHistory",       LogObjHistory);
      UtilPrintParameter(os, sec, "LimitInitVars",       LimitInitVars);
      UtilPrintParameter(os, sec, "DebugLevel",          DebugLevel);
      UtilPrintParameter(os, sec, "TolZero",             TolZero);
      UtilPrintParameter(os, sec, "LimitTotalCutIters",  LimitTotalCutIters);
      UtilPrintParameter(os, sec, "LimitTotalPriceIters", LimitTotalPriceIters);
      UtilPrintParameter(os, sec, "LimitRoundCutIters",  LimitRoundCutIters);
      UtilPrintParameter(os, sec, "LimitRoundPriceIters", LimitRoundPriceIters);
      UtilPrintParameter(os, sec, "LimitTime",           LimitTime);
      UtilPrintParameter(os, sec, "LimitNodes",          LimitNodes);
      UtilPrintParameter(os, sec, "TailoffLength",       TailoffLength);
      UtilPrintParameter(os, sec, "TailoffPercent",      TailoffPercent);
      UtilPrintParameter(os, sec, "MasterGapLimit",      MasterGapLimit);
      UtilPrintParameter(os, sec, "PCStrategy",          PCStrategy);
      UtilPrintParameter(os, sec, "CompressColumns",     CompressColumns);
      UtilPrintParameter(os, sec, "CompressColumnsIterFreq",
                         CompressColumnsIterFreq);
      UtilPrintParameter(os, sec, "CompressColumnsSizeMultLimit",
                         CompressColumnsSizeMultLimit);
      UtilPrintParameter(os, sec, "CompressColumnsMasterGapStart",
                         CompressColumnsMasterGapStart);
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
      UtilPrintParameter(os, sec, "SubProbNumThreads",  SubProbNumThreads);
      UtilPrintParameter(os, sec, "SubProbNumSolLimit", SubProbNumSolLimit);
      UtilPrintParameter(os, sec, "SubProbSolverStartAlgo",
                         SubProbSolverStartAlgo);
      UtilPrintParameter(os, sec, "RoundRobinInterval",  RoundRobinInterval);
      UtilPrintParameter(os, sec, "RoundRobinStrategy",  RoundRobinStrategy);
      UtilPrintParameter(os, sec, "SolveMasterAsIp",     SolveMasterAsIp);
      UtilPrintParameter(os, sec, "SolveMasterAsIpFreqNode",
                         SolveMasterAsIpFreqNode);
      UtilPrintParameter(os, sec, "SolveMasterAsIpFreqPass",
                         SolveMasterAsIpFreqPass);
      UtilPrintParameter(os, sec, "SolveMasterAsIpLimitTime",
                         SolveMasterAsIpLimitTime);
      UtilPrintParameter(os, sec, "SolveMasterAsIpLimitGap",
                         SolveMasterAsIpLimitGap);
      UtilPrintParameter(os, sec, "SolveMasterUpdateAlgo",
                         SolveMasterUpdateAlgo);
      UtilPrintParameter(os, sec, "SolveRelaxAsIp",     SolveRelaxAsIp);
      UtilPrintParameter(os, sec, "InitVarsWithCutDC",   InitVarsWithCutDC);
      UtilPrintParameter(os, sec, "InitVarsWithIP",   InitVarsWithIP);
      UtilPrintParameter(os, sec, "InitVarsWithIPLimitTime",
                         InitVarsWithIPLimitTime);
      UtilPrintParameter(os, sec, "InitCompactSolve",  InitCompactSolve);
      UtilPrintParameter(os, sec, "DualStab",          DualStab);
      UtilPrintParameter(os, sec, "DualStabAlpha",     DualStabAlpha);
      UtilPrintParameter(os, sec, "BreakOutPartial",   BreakOutPartial);
      UtilPrintParameter(os, sec, "BranchEnforceInSubProb",
                         BranchEnforceInSubProb);
      UtilPrintParameter(os, sec, "BranchEnforceInMaster",
                         BranchEnforceInMaster);
      UtilPrintParameter(os, sec, "MasterConvexityLessThan",
                         MasterConvexityLessThan);
      UtilPrintParameter(os, sec, "ParallelColsLimit", ParallelColsLimit);
      UtilPrintParameter(os, sec, "BranchStrongIter",  BranchStrongIter);
      UtilPrintParameter(os, sec, "NumThreads",        NumThreads);
      UtilPrintParameter(os, sec,
                         "DebugCheckBlocksColumns", DebugCheckBlocksColumns);
      UtilPrintParameter(os, sec, "NumBlocks",  NumBlocks);
      UtilPrintParameter(os, sec, "LogLevel",  LogLevel);
      UtilPrintParameter(os, sec, "DataDir",  DataDir);
      UtilPrintParameter(os, sec, "Instance",  Instance);
      UtilPrintParameter(os, sec, "InstanceFormat",  InstanceFormat);
      UtilPrintParameter(os, sec, "BlockFile",  BlockFile);
      UtilPrintParameter(os, sec, "PermuteFile",  PermuteFile);
      UtilPrintParameter(os, sec, "BlockFileFormat",  BlockFileFormat);
      UtilPrintParameter(os, sec, "InitSolutionFile",  InitSolutionFile);
      UtilPrintParameter(os, sec, "UseNames",  UseNames);
      UtilPrintParameter(os, sec, "UseSparse",  UseSparse);
      UtilPrintParameter(os, sec, "FullModel",  FullModel);
      UtilPrintParameter(os, sec, "BestKnownLB",  BestKnownLB);
      UtilPrintParameter(os, sec, "BestKnownUB",  BestKnownUB);
      UtilPrintParameter(os, sec, "ColumnUB",  ColumnUB);
      UtilPrintParameter(os, sec, "ColumnLB",  ColumnLB);
      UtilPrintParameter(os, sec, "ObjectiveSense",  ObjectiveSense);
      UtilPrintParameter(os, sec, "Concurrent", Concurrent);
      UtilPrintParameter(os, sec, "NumBlocksCand", NumBlocksCand);
      UtilPrintParameter(os, sec, "ConcurrentCutOffTime", ConcurrentCutOffTime);
      UtilPrintParameter(os, sec,  "ThreadIndex", ThreadIndex );
      UtilPrintParameter(os, sec,  "CurrentWorkingDir", CurrentWorkingDir);
      UtilPrintParameter(os, sec, "SubProbParallel", SubProbParallel);
      UtilPrintParameter(os, sec, "SubProbParallelType", SubProbParallelType);
      UtilPrintParameter(os, sec, "SubProbParallelChunksize", SubProbParallelChunksize);
      UtilPrintParameter(os, sec, "ConcurrentThreadsNum", ConcurrentThreadsNum);
      (*os) << "========================================================\n";
   }

   void setDefaults() {
      LogLevel             = 0;
      LogDebugLevel        = 0;
      LogLpLevel           = 0;
      LogDumpModel         = 0;
      LogObjHistory        = 0;
      LimitInitVars        = 5;
      DebugLevel           = 0;
      TolZero              = DecompEpsilon;
      LimitTotalCutIters   = COIN_INT_MAX;
      LimitTotalPriceIters = COIN_INT_MAX;
      LimitRoundCutIters   = COIN_INT_MAX;
      LimitRoundPriceIters = COIN_INT_MAX;
      LimitTime            = DecompBigNum;
      LimitNodes           = COIN_INT_MAX;
      TailoffLength        = 10;
      TailoffPercent       = 0.10;
      MasterGapLimit       = 1.0e-6;
      PCStrategy           = 0;
      CompressColumns      = 1;
      CompressColumnsIterFreq       = 2;
      CompressColumnsSizeMultLimit  = 1.20;
      CompressColumnsMasterGapStart = 0.20;
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
      SubProbTimeLimitExact   = DecompBigNum;
      SubProbTimeLimitInexact = DecompBigNum;
      SubProbNumThreads       = 1;
      SubProbNumSolLimit      = 1;
      SubProbSolverStartAlgo = DecompDualSimplex;
      RoundRobinInterval   = 0;
      RoundRobinStrategy   = RoundRobinRotate;
      SolveMasterAsIp          = 1;//TODO: turn off if one block
      SolveMasterAsIpFreqNode  = 1;
      SolveMasterAsIpFreqPass  = 1000;
      SolveMasterAsIpLimitTime = 30;
      SolveMasterAsIpLimitGap  = 0.05; //5% gap
      SolveRelaxAsIp           = false;
      SolveMasterUpdateAlgo    = DecompDualSimplex;
      InitVarsWithCutDC        = 0;
      InitVarsWithIP           = 0;
      InitVarsWithIPLimitTime  = 10;
      InitCompactSolve         = 0;
      DualStab                 = 0;
      DualStabAlpha            = 0.10;
      BreakOutPartial          = 0;
      BranchEnforceInSubProb   = 1;//usually much better if can
      BranchEnforceInMaster    = 0;
      MasterConvexityLessThan  = 0;
      ParallelColsLimit        = 1.0;
      BranchStrongIter         = 0;
      NumThreads               = 2;
      DebugCheckBlocksColumns  = true;
      NumBlocks                = 3;
      /*
       * parameters from MILPBlock and to be MILP
       */
      LogLevel                 = 0;
      DataDir                  = "";
      Instance                 = "";
      BlockFile                = "";
      BlockFileFormat          = "";
      PermuteFile              = "";
      InitSolutionFile         = "";
      UseNames                 = 0 ;
      UseSparse                = 1 ;
      FullModel                = 0 ;
      BestKnownLB              = -1.e100;
      BestKnownUB              = 1.e100;
      ColumnUB                 = 1.e20;
      ColumnLB                 = -1.e20;
      ObjectiveSense           = 1;
      Concurrent               = false;
      NumBlocksCand            = 10;
      ConcurrentCutOffTime     = 100;
      ThreadIndex              = 0;
      CurrentWorkingDir        = "";
      SubProbParallel          = false;
      SubProbParallelType      = SubProbScheduleDynamic;
      SubProbParallelChunksize = 1;
      ConcurrentThreadsNum     = 4;
   }

   void dumpSettings(std::ostream* os = &std::cout) {
      const std::string sec = "DECOMP";
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
