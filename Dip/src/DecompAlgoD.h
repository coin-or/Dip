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
#ifndef DecompAlgoD_h_
#define DecompAlgoD_h_

//===========================================================================//
/**
 * \class DecompAlgoD
 * \brief Class for DECOMP algorithm Decomp.
 *
 */
//===========================================================================//
//THINK: derive from DecompAlgo or DecompAlgoPC?? THINK
//THINK: how can we reuse this object since call many times?
//   if init phase is feasible, we are done....


//===========================================================================//
#include "DecompAlgoPC.h"

//===========================================================================//
class DecompAlgoD : public DecompAlgoPC {
private:

   //----------------------------------------------------------------------//
   /**
    * @name Data.
    * @{
    */
   //----------------------------------------------------------------------//
   /**
    * Store the name of the class (for logging/debugging) - "who am I?"
    */
   std::string m_classTag;

   //TODO
   double*             m_xhatD;
   //TODO
   DecompCutList*      m_newCuts;
   //TODO
   int                 m_numOrigCols;




   /**
    * @}
    */

   //-----------------------------------------------------------------------//
   /**
    * @name Derived from virtual functions of DecompAlgoPC
    * @{
    */
   //-----------------------------------------------------------------------//
   /**
    * Create the master problem (all algorithms must define this function).
    */
   virtual void createMasterProblem(DecompVarList& initVars);
   virtual void masterMatrixAddArtCols(CoinPackedMatrix* masterM,
                                       double*            colLB,
                                       double*            colUB,
                                       double*            objCoeff,
                                       std::vector<std::string>&    colNames,
                                       int                startRow,
                                       int                endRow,
                                       char               origOrBranch);
   virtual void phaseUpdate(DecompPhase&   phase,
                            DecompStatus& status);
   virtual void phaseDone();

   /**
    * Set the current integer bound and update best/history.
    */
   virtual inline void setObjBoundIP(const double thisBound) {
      UtilPrintFuncBegin(m_osLog, m_classTag,
                         "setObjBoundIP()", m_param.LogDebugLevel, 2);

      if (thisBound < m_nodeStats.objBest.second) {
         UTIL_MSG(m_app->m_param.LogDebugLevel, 3,
                  (*m_osLog) << "New Global UB = "
                  << UtilDblToStr(thisBound) << std::endl;);
         //For DECOMP, don't update this object's global UB
         //  otherwise, we might stop to early, since it will
         //  compare to this object's lower bound.
         //m_nodeStats.objBest.second = thisBound;
      }

      UtilPrintFuncEnd(m_osLog, m_classTag,
                       "setObjBoundIP()", m_param.LogDebugLevel, 2);
   }



public:
   void solveD(DecompCutList* newCuts) {
      m_newCuts = newCuts;
      //need to change parameters to price, no cut
      m_param.LimitTotalCutIters   = 0;
      m_param.LimitRoundCutIters   = 0;
      m_param.LimitTotalPriceIters = 1000;
      m_param.LimitRoundPriceIters = 1000;
      m_param.SolveMasterAsIp      = 0;
      processNode(NULL);
   }


public:

   //-----------------------------------------------------------------------//
   /**
    * @name Constructors and destructor.
    * @{
    */
   //-----------------------------------------------------------------------//
private:
   /**
    * Disable copy constructors.
    */
   DecompAlgoD(const DecompAlgoD&);
   DecompAlgoD& operator=(const DecompAlgoD&);

public:

   /**
    * Default constructors.
    */
   DecompAlgoD(DecompApp*             app,
               UtilParameters*        utilParam,
               double*                xhat,
               int                    numOrigCols) //need to pass this? :
      :
      DecompAlgoPC(app, utilParam,
                   const_cast<std::string&>(DecompAlgoStr[DECOMP]), false),
      m_classTag   ("D-ALGOD"),
      m_xhatD      (xhat),
      m_newCuts    (0),
      m_numOrigCols(numOrigCols) { //need?
      std::string paramSection = DecompAlgoStr[DECOMP];
      m_algo              = DECOMP;
      initSetup(utilParam, paramSection);
   }

   //need this?
   DecompAlgoD(DecompApp*       app,
               UtilParameters* utilParam,
               std::string&          paramSection,
               double*          xhat,
               int              numOrigCols):
      DecompAlgoPC(app, utilParam, paramSection, false),
      m_classTag   ("D-ALGOD"),
      m_xhatD      (xhat),
      m_newCuts    (0),
      m_numOrigCols(numOrigCols) { //need?
      m_algo = DECOMP;
      initSetup(utilParam, paramSection);
   }

   /**
    * Destructor.
    */
   ~DecompAlgoD() {}
   /**
    * @}
    */
};

#endif
