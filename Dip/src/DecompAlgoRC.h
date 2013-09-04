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
#ifndef DecompAlgoRC_h_
#define DecompAlgoRC_h_

/** \todo Next: DecompAlgoVC - use Vol? or write from scratch? */
//===========================================================================//
#include "DecompAlgo.h"

//===========================================================================//
class DecompAlgoRC : public DecompAlgo {

   //----------------------------------------------------------------------//
   /**
    * @name Data.
    * @{
    */
   //----------------------------------------------------------------------//
private:
   /**
    * Store the name of the class (for logging/debugging) - "who am I?"
    */
   const std::string m_classTag;

private:
   std::vector<double>   m_u;   //dual vector
   double*          m_rc;  //reduced cost

   double           m_UB;
   double           m_LB;

   int              m_cntSameLB;
   int              m_iter;
   double           m_step;
   bool             m_zeroSub;

   DecompVar        m_shatVar;
   /**
    * @}
    */

   //-----------------------------------------------------------------------//
   /**
    * @name Derived from pure virtual functions of DecompAlgo.
    * @{
    */
   //-----------------------------------------------------------------------//
   /**
    * Create the master problem (all algorithms must define this function).
    */
   void createMasterProblem(DecompVarList& initVars);

   /**
    * Compose solution in x-space from current space.
    *  - PC: this recomposes x from lambda
    *  - C : this just copies over LP solution
    */
   //not pure
   void recomposeSolution(const double* solution,
                          double*        rsolution);


   /**
    * @}
    */



   //-----------------------------------------------------------------------//
   /**
    * @name Derived from virtual functions of DecompAlgo
    * @{
    */
   //-----------------------------------------------------------------------//
   DecompStatus solutionUpdate(const DecompPhase phase,
                               const int         maxInnerIter,
                               const int         maxOuterIter);
   int addCutsFromPool();
   int generateVars(const DecompStatus   stat,
                    DecompVarList&     newVars,
                    double&            mostNegReducedCost);
   bool updateObjBound(const double mostNegRC = -DecompBigNum);

   /**
    * Run the initial phase for processing node.
    */
   DecompPhase phaseInit();
   /**
    * Run the done phase for processing node.
    */
   void phaseDone();


   /**
    * @}
    */

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
   DecompAlgoRC(const DecompAlgoRC&);
   DecompAlgoRC& operator=(const DecompAlgoRC&);

public:
   /**
    * Default constructors.
    */
   DecompAlgoRC(DecompApp*             app,
                UtilParameters*        utilParam):
      DecompAlgo(RELAX_AND_CUT, app, utilParam),
      m_classTag("D-ALGORC"),
      m_u        (),
      m_rc       (NULL),
      m_UB       (DecompInf),
      m_LB       (-DecompInf),
      m_cntSameLB(0),
      m_iter     (0),
      m_step     (2.0), //(0, 2] param?
      m_zeroSub  (false),
      m_shatVar  () {
      std::string paramSection = DecompAlgoStr[RELAX_AND_CUT];
      initSetup(utilParam, paramSection);
   }

   DecompAlgoRC(DecompApp*       app,
                UtilParameters* utilParam,
                std::string&          paramSection) :
      DecompAlgo(RELAX_AND_CUT, app, utilParam),
      m_classTag ("D-ALGORC"),
      m_u        (),
      m_rc       (NULL),
      m_UB       (DecompInf),
      m_LB       (-DecompInf),
      m_cntSameLB(0),
      m_iter     (0),
      m_step     (2.0), //(0, 2] param?
      m_zeroSub  (false),
      m_shatVar  () {
      initSetup(utilParam, paramSection);
   }

   /**
    * Destructor.
    */
   ~DecompAlgoRC() {
      UTIL_DELARR(m_rc);
   }
   /**
    * @}
    */
   /**
    * @}
    */


public:
   bool isDone();
   //name - change to getDual?
   const double* getRowPrice() {
      return &m_u[0];
   }
   //user needs to do?
   //STOP
   void setInitObjUB(const double objUB) {
      m_UB = objUB;
   }


};
#endif
