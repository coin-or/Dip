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


   //-----------------------------------------------------------------------//
   /**
    * @name Helper functions.
    * @{
    */
   //-----------------------------------------------------------------------//
public:
   void getSettings(UtilParameters& param) {
      static const char* sec = "ALPS";
      logFileLevel    = param.GetSetting("logFileLevel",    0,            sec);
      printSolution   = param.GetSetting("printSolution",   false,        sec);
      checkMemory     = param.GetSetting("checkMemory",     false,        sec);
      msgLevel        = param.GetSetting("msgLevel",        2,            sec);
      nodeLimit       = param.GetSetting("nodeLimit",       ALPS_INT_MAX, sec);
      nodeLogInterval = param.GetSetting("nodeLogInterval", 10,           sec);

      if (msgLevel > 2) {
         dumpSettings();
      }
   }

   void dumpSettings(std::ostream* os = &std::cout) {
      static const char* sec = "ALPS";
      (*os) << "\n========================================================\n"
            << "ALPS PARAMETER SETTINGS \n";
      (*os) << sec << ": logFileLevel    = " << logFileLevel    << std::endl;
      (*os) << sec << ": printSolution   = " << printSolution   << std::endl;
      (*os) << sec << ": checkMemory     = " << checkMemory     << std::endl;
      (*os) << sec << ": msgLevel        = " << msgLevel        << std::endl;
      (*os) << sec << ": nodeLimit       = " << nodeLimit       << std::endl;
      (*os) << sec << ": nodeLogInterval = " << nodeLogInterval << std::endl;
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
   AlpsDecompParam() {}

   AlpsDecompParam(UtilParameters& utilParam) {
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
