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
#ifndef DecompAlgoC_h_
#define DecompAlgoC_h_

//===========================================================================//
/**
 * \class DecompAlgoC
 * \brief Class for DECOMP algorithm Cutting Plane Method.
 *
 */
//===========================================================================//

//===========================================================================//
#include "Decomp.h"
#include "DecompAlgo.h"

//===========================================================================//
class DecompAlgoC : public DecompAlgo {
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

   /**
    * Calculate the current objective LB, update the best, and
    * store in history.
    */
   bool updateObjBound(const double mostNegRC = -DecompBigNum);

   void phaseInit(DecompPhase& phase) {
      if (getNodeIndex() == 0) {
         phase = PHASE_CUT;
      }
   }
   void phaseDone();//chance to run DC

   /**
    * Update of the phase for process loop.
    */
   void phaseUpdate(DecompPhase&   phase,
                    DecompStatus& status);
   /**
    * Generate initial variables for master problem (PC/DC/RC).
    *   - in CPM, this does nothing
    */
   int generateInitVars(DecompVarList& initVars) {
      return 0;
   }

   void setMasterBounds(const double* lbs,
                        const double* ubs);
   void setSubProbBounds(const double* lbs,
                         const double* ubs) {};

public:
   virtual DecompSolverResult*
   solveDirect(const DecompSolution* startSol  = NULL);


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
   DecompAlgoC(const DecompAlgoC&);
   DecompAlgoC& operator=(const DecompAlgoC&);

public:
   /**
    * Default constructors.
    */
   DecompAlgoC(DecompApp*             app,
               UtilParameters*        utilParam):
      DecompAlgo(CUT, app, utilParam),
      m_classTag("D-ALGOC") {
      std::string paramSection = DecompAlgoStr[CUT];
      initSetup(utilParam, paramSection);
   }

   DecompAlgoC(DecompApp*       app,
               UtilParameters* utilParam,
               std::string&          paramSection):
      DecompAlgo(CUT, app, utilParam),
      m_classTag("D-ALGOC") {
      initSetup(utilParam, paramSection);
   }


   /**
    * Destructor.
    */
   ~DecompAlgoC() {}
   /**
    * @}
    */
};

#endif
