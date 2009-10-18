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
#ifndef DecompAlgoPC_h_
#define DecompAlgoPC_h_

//===========================================================================//
/**
 * \class DecompAlgoPC
 * \brief Class for DECOMP algorithm Price and Cut.
 *
 */
//===========================================================================//

//===========================================================================//
#include "DecompAlgo.h"

//===========================================================================//
class DecompAlgoPC : public DecompAlgo {
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
   string m_classTag;

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
   virtual void createMasterProblem(DecompVarList & initVars){
      DecompAlgo::createMasterProblem(initVars);
   }
   virtual int generateVarsFea(DecompVarList    & newVars, 
			       double           & mostNegReducedCost){
      return DecompAlgo::generateVarsFea(newVars, mostNegReducedCost);
   }
   virtual void phaseInit(DecompPhase & phase);
   
   /**
    * @}
    */
   
   //-----------------------------------------------------------------------//
   /**
    * @name Derived from virtual functions of DecompAlgo
    * @{
    */
   //-----------------------------------------------------------------------//
   //TODO
   void addCutsToPool(const double  *  x,
                      DecompCutList & newCuts,
                      int           & n_newCuts);

   //TODO
   void phaseDone();
   int  addCutsFromPool();
   void solutionUpdateAsIP();
   int  adjustColumnsEffCnt();
   int  compressColumns    ();


   
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
   DecompAlgoPC(DecompApp      * app,
                UtilParameters * utilParam,
                bool             doSetup    = true) :
      DecompAlgo(PRICE_AND_CUT, app, utilParam),
      m_classTag("D-ALGOPC") {
      if(doSetup){
         string paramSection = DecompAlgoStr[PRICE_AND_CUT];
         initSetup(utilParam, paramSection);
      }
   }
      
   DecompAlgoPC(DecompApp      * app,
                UtilParameters * utilParam,
                string         & paramSection,
                bool             doSetup    = true) :
      //is utilParam used in base class?
      DecompAlgo(PRICE_AND_CUT, app, utilParam),
      m_classTag("D-ALGOPC") {
      if(doSetup)
         initSetup(utilParam, paramSection);
   }
      
      
   /**
    * Destructor.
    */
   ~DecompAlgoPC(){}
   /**
    * @}
    */
};

#endif
