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

//===========================================================================//
#ifndef AlpsDecompParam_h_
#define AlpsDecompParam_h_

//===========================================================================//
#include "UtilParameters.h"

//===========================================================================//
/**
 * \class AlpsDecompParam
 * \brief
 * Parameters passed through to Alps. 
*/
//===========================================================================//

//===========================================================================//
class AlpsDecompParam {

   //----------------------------------------------------------------------//
   /**      
    * @name Data.
    * @{      
    */
   //----------------------------------------------------------------------//
   
public:
   /** The level of log file.
    *  - 0: no print to screen (Default)
    *  - 1: summary
    *  - 2: moderate
    *  - 3: verbose
    */
   int logFileLevel;
   
   /**
    * Print solution to screen and log if have a solution
    * and msgLevel and logFileLevel permits. Default: false.
    */
   bool printSolution;
   
   /**
    * Check memory. Default: false
    */
   bool checkMemory;
   
   /**
    * The level of printing messages on screen. Used to control master
    * and general messages.
    *  - 0: no print to screen
    *  - 1: summary
    *  - 2: moderate (Default)
    *  - 3: verbose
    */
   int msgLevel;
   
   /**
    * The max number of nodes can be processed. Default: ALPS_INT_MAX
   */
   int nodeLimit;
   
   /**
    * Node log interval. Default: 100
    */
   int nodeLogInterval;

   /**
    * The time limit (in seconds) of search. Default: ALPS_DBL_MAX
    */
   double timeLimit;
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
   void getSettings(UtilParameters & param){
      static const char * sec = "ALPS";
      logFileLevel    = param.GetSetting("logFileLevel",    0,            sec);
      printSolution   = param.GetSetting("printSolution",   false,        sec);
      checkMemory     = param.GetSetting("checkMemory",     false,        sec);
      msgLevel        = param.GetSetting("msgLevel",        2,            sec);
      nodeLimit       = param.GetSetting("nodeLimit",       ALPS_INT_MAX, sec);
      nodeLogInterval = param.GetSetting("nodeLogInterval", 10,           sec);
      timeLimit       = param.GetSetting("timeLimit",       ALPS_DBL_MAX, sec);
   }
   
   void dumpSettings(ostream * os = &cout){
      static const char * sec = "ALPS";
      (*os) << "\n========================================================\n"
            << "ALPS PARAMETER SETTINGS \n";
      (*os) << sec << ": logFileLevel    = " << logFileLevel    << endl;
      (*os) << sec << ": printSolution   = " << printSolution   << endl;
      (*os) << sec << ": checkMemory     = " << checkMemory     << endl;
      (*os) << sec << ": msgLevel        = " << msgLevel        << endl;
      (*os) << sec << ": nodeLimit       = " << nodeLimit       << endl;
      (*os) << sec << ": nodeLogInterval = " << nodeLogInterval << endl;
      (*os) << sec << ": timeLimit       = " << timeLimit       << endl;
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
   /**
    * Disable copy constructors.
    */
private:
   AlpsDecompParam(const AlpsDecompParam &);
   AlpsDecompParam & operator=(const AlpsDecompParam &);

public:   
   /**
    * Default constructors.
    */
   AlpsDecompParam() {}
   
   AlpsDecompParam(UtilParameters & utilParam) {
      getSettings(utilParam);
   }

   /**
    * Destructor
    */   
   ~AlpsDecompParam() {}
   /**
    * @}
    */
};

#endif
