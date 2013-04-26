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

#ifndef MILPBLOCK_PARAM_INCLUDED
#define MILPBLOCK_PARAM_INCLUDED

//===========================================================================//
#include "UtilParameters.h"

using namespace std;

//===========================================================================//
class MILPBlock_Param{
public:
   int    LogLevel;   
   string DataDir;
   string Instance;

   /*
    * The file defining which rows are in which blocks.
    */
   string BlockFile;
   
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
    * (2) "Pair" or "PAIR"
    * Each line is a block id to row id pair.
    *   <block id> <row id> 
    *
    * (3) "PairName" or "PAIRNAME"
    * Each line is a block id to row name (matching mps) pair.
    *   <block id> <row name>     
    */
   string BlockFileFormat;

   string PermuteFile;
   string InitSolutionFile;
   int    UseNames;  //col/row names for debugging
   int    UseSparse; //create all blocks sparsely
   int    FullModel; //create full model for CPM or direct
   double BestKnownLB;
   double BestKnownUB;
   double ColumnUB; //hack since missing extreme rays
   double ColumnLB; //hack since missing extreme rays

   //TOOD: better solution for this
   int ObjectiveSense;   //1=min, -1=max

public:
   void getSettings(UtilParameters & utilParam){
      static const char * common = "MILPBlock";
      LogLevel     = utilParam.GetSetting("LogLevel",     0,     common);
      DataDir      = utilParam.GetSetting("DataDir",      "",    common);
      Instance     = utilParam.GetSetting("Instance",     "",    common);    
      BlockFile    = utilParam.GetSetting("BlockFile",    "",    common);
      PermuteFile  = utilParam.GetSetting("PermuteFile",  "",    common);
      BlockFileFormat 
         = utilParam.GetSetting("BlockFileFormat",    "",    common);    
      InitSolutionFile
         = utilParam.GetSetting("InitSolutionFile",    "",    common);    
      UseNames     = utilParam.GetSetting("UseNames",       0, common);
      UseSparse    = utilParam.GetSetting("UseSparse",      1, common);
      FullModel    = utilParam.GetSetting("FullModel",      0, common);
      BestKnownLB  = utilParam.GetSetting("BestKnownLB",  -1.e100, common);
      BestKnownUB  = utilParam.GetSetting("BestKnownUB",   1.e100, common);
      ColumnUB     = utilParam.GetSetting("ColumnUB",      1.e20,  common);
      ColumnLB     = utilParam.GetSetting("ColumnLB",     -1.e20,  common);
      ObjectiveSense= utilParam.GetSetting("ObjectiveSense",   1,  common);
   }

   void dumpSettings(ostream * os = &cout){
      static const char * common = "MILPBlock";
      (*os) << "\n=====================================================\n"
            << "MILPBlock_DECOMP PARAMETER SETTINGS \n";
      (*os) << common << ": LogLevel          : " << LogLevel         << endl;
      (*os) << common << ": DataDir           : " << DataDir          << endl;
      (*os) << common << ": Instance          : " << Instance         << endl;
      (*os) << common << ": BlockFile         : " << BlockFile        << endl;
      (*os) << common << ": PermuteFile       : " << PermuteFile      << endl;
      (*os) << common << ": BlockFileFormat   : " << BlockFileFormat  << endl;
      (*os) << common << ": InitSolutionFile  : " << InitSolutionFile << endl;
      (*os) << common << ": UseNames          : " << UseNames         << endl;
      (*os) << common << ": UseSparse         : " << UseSparse        << endl;
      (*os) << common << ": FullModel         : " << FullModel        << endl;
      (*os) << common << ": BestKnownLB       : " << BestKnownLB      << endl;
      (*os) << common << ": BestKnownUB       : " << BestKnownUB      << endl;
      (*os) << common << ": ColumnUB          : " << ColumnUB         << endl;
      (*os) << common << ": ColumnLB          : " << ColumnLB         << endl;
      (*os) << common << ": ObjectiveSense    : " << ObjectiveSense   << endl;
      (*os) << "\n=====================================================\n";
   }
   
public:
   MILPBlock_Param():    
      LogLevel        (0 ),
      DataDir         (""),
      Instance        (""),
      BlockFile       (""), 
      BlockFileFormat (""),
      PermuteFile     (""), 
      InitSolutionFile(""),
      UseNames        (0),
      UseSparse       (1),
      FullModel       (0),
      BestKnownLB     (-1.e100),
      BestKnownUB     ( 1.e100),
      ColumnUB        ( 1.e20),
      ColumnLB        (-1.e20),
      ObjectiveSense  (1)
{};
   ~MILPBlock_Param() {};
};

#endif
