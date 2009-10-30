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
   double ColumnUB; //hack since missing extreme rays
   double ColumnLB; //hack since missing extreme rays

public:
   void getSettings(UtilParameters & utilParam){
      static const char * common = "MILPBlock";
      LogLevel     = utilParam.GetSetting("LogLevel",     0,     common);
      DataDir      = utilParam.GetSetting("DataDir",      "",    common);
      Instance     = utilParam.GetSetting("Instance",     "",    common);    
      BlockFile    = utilParam.GetSetting("BlockFile",    "",    common);    
      BlockFileFormat 
	 = utilParam.GetSetting("BlockFileFormat",    "",    common);    
      ColumnUB     = utilParam.GetSetting("ColumnUB",  1.0e20, common);
      ColumnLB     = utilParam.GetSetting("ColumnLB", -1.0e20, common);
   }

   void dumpSettings(ostream * os = &cout){
      static const char * common = "MILPBlock";
      (*os) << "\n=====================================================\n"
            << "MILPBlock_DECOMP PARAMETER SETTINGS \n";
      (*os) << common << ": LogLevel        : " << LogLevel        << endl;
      (*os) << common << ": DataDir         : " << DataDir         << endl;
      (*os) << common << ": Instance        : " << Instance        << endl;
      (*os) << common << ": BlockFile       : " << BlockFile       << endl;
      (*os) << common << ": BlockFileFormat : " << BlockFileFormat << endl;
      (*os) << common << ": ColumnUB        : " << ColumnUB        << endl;
      (*os) << common << ": ColumnLB        : " << ColumnLB        << endl;
      (*os) << "\n=====================================================\n";
   }
   
public:
   MILPBlock_Param():    
      LogLevel       (0 ),
      DataDir        (""),
      Instance       (""),
      BlockFile      (""), 
      BlockFileFormat(""),
      ColumnUB       ( 1.e20),
      ColumnLB       (-1.e20)
{};
   ~MILPBlock_Param() {};
};

#endif
