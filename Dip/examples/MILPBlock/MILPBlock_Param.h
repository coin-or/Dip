//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, Lehigh University                                 //
//                                                                           //
// Copyright (C) 2002-2009, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

#ifndef MILPBLOCK_PARAM_INCLUDED
#define MILPBLOCK_PARAM_INCLUDED

//===========================================================================//
#include "UtilParameters.h"

//===========================================================================//
class MILPBlock_Param{
public:
   int    LogLevel;   
   string DataDir;
   string Instance;
   string BlockFile;
   string BlockFileFormat;
   string InitSolutionFile;
   double BestKnownLB;
   double BestKnownUB;
   double ColumnUB; //hack since missing extreme rays
   double ColumnLB; //hack since missing extreme rays
   //=1 if all master-only vars in one block
   //=0 if one block per master-only var
   int    MasterOnlyOneBlock;

public:
   void getSettings(UtilParameters & utilParam){
      static const char * common = "MILPBlock";
      LogLevel     = utilParam.GetSetting("LogLevel",     0,     common);
      DataDir      = utilParam.GetSetting("DataDir",      "",    common);
      Instance     = utilParam.GetSetting("Instance",     "",    common);    
      BlockFile    = utilParam.GetSetting("BlockFile",    "",    common);
      BlockFileFormat 
         = utilParam.GetSetting("BlockFileFormat",    "",    common);    
      InitSolutionFile
         = utilParam.GetSetting("InitSolutionFile",    "",    common);    
      BestKnownLB  = utilParam.GetSetting("BestKnownLB",  -1.e100, common);
      BestKnownUB  = utilParam.GetSetting("BestKnownUB",   1.e100, common);
      ColumnUB     = utilParam.GetSetting("ColumnUB",      1.e20,  common);
      ColumnLB     = utilParam.GetSetting("ColumnLB",     -1.e20,  common);
      MasterOnlyOneBlock = 
         utilParam.GetSetting("MasterOnlyOneBlock",            1,  common);
   }

   void dumpSettings(ostream * os = &cout){
      static const char * common = "MILPBlock";
      (*os) << "\n=====================================================\n"
            << "MILPBlock_DECOMP PARAMETER SETTINGS \n";
      (*os) << common << ": LogLevel          : " << LogLevel         << endl;
      (*os) << common << ": DataDir           : " << DataDir          << endl;
      (*os) << common << ": Instance          : " << Instance         << endl;
      (*os) << common << ": BlockFile         : " << BlockFile        << endl;
      (*os) << common << ": BlockFileFormat   : " << BlockFileFormat  << endl;
      (*os) << common << ": InitSolutionFile  : " << InitSolutionFile << endl;
      (*os) << common << ": BestKnownLB       : " << BestKnownLB      << endl;
      (*os) << common << ": BestKnownUB       : " << BestKnownUB      << endl;
      (*os) << common << ": ColumnUB          : " << ColumnUB         << endl;
      (*os) << common << ": ColumnLB          : " << ColumnLB         << endl;
      (*os) << common << ": MasterOnlyOneBlock: " 
            << MasterOnlyOneBlock << endl;
      (*os) << "\n=====================================================\n";
   }
   
public:
   MILPBlock_Param():    
      LogLevel       (0 ),
      DataDir        (""),
      Instance       (""),
      BlockFile      (""), 
      BlockFileFormat(""),
      InitSolutionFile(""),
      BestKnownLB    (-1.e100),
      BestKnownUB    ( 1.e100),
      ColumnUB       ( 1.e20),
      ColumnLB       (-1.e20),
      MasterOnlyOneBlock(1)
{};
   ~MILPBlock_Param() {};
};

#endif
