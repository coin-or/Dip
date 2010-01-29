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
// Conceptual Design: Matthew Galati, SAS Institute Inc.                     //
//                    Ted Ralphs, Lehigh University                          //
//                                                                           //
// Copyright (C) 2002-2009, Lehigh University, Matthew Galati, Ted Ralphs    //
// All Rights Reserved.                                                      //
//===========================================================================//

//#define DUAL_SMOOTHING

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
   vector<double> m_dual;    //duals from stabilized (if bound improved)
   vector<double> m_dualRM;  //duals from restricted master
   vector<double> m_dualST;  //duals from stabilized method

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


   virtual const double * getRowPrice()  {
      //---
      //--- return the duals to be used in pricing step
      //---
#ifdef DUAL_SMOOTHING
      //---
      //--- resize dual vectors
      //---
      size_t nRows = m_masterSI->getNumRows();
      m_dual.resize(nRows);
      m_dualRM.resize(nRows);
      m_dualST.resize(nRows);

      //---
      //--- calculate smoothed dual
      //---    pi_ST = alpha * pi_Bar + (1-alpha) * pi_RM
      //---
      int            r;
      const double * u      = m_masterSI->getRowPrice();
      double         alpha  = m_param.DualStabAlpha;
      double         alpha1 = 1.0 - alpha; 
      copy(u, u + nRows, m_dualRM.begin()); //copy for sake of debugging
      for(r = 0; r < nRows; r++){
         m_dualST[r] = (alpha * m_dual[r]) + (alpha1 * m_dualRM[r]);
      }      

      const vector<string> & rowNames = m_masterSI->getRowNames();
      for(r = 0; r < m_masterSI->getNumRows(); r++){
	 if(!(UtilIsZero(m_dual[r]) && UtilIsZero(m_dualRM[r]) && UtilIsZero(m_dualST[r]))){
	    if(r < static_cast<int>(rowNames.size())){
	       printf("MASTER %25s DUAL[%6d->%20s] = %12.10f RM = %12.10f ST = %12.10f\n",
		      DecompRowTypeStr[m_masterRowType[r]].c_str(),
		      r, rowNames[r].c_str(), m_dual[r], m_dualRM[r], m_dualST[r]);
	    }
	    else
	       printf("MASTER %25s DUAL[%6d] = %12.10f RM = %12.10f ST = %12.10f\n",
		      DecompRowTypeStr[m_masterRowType[r]].c_str(),
		      r, m_dual[r], m_dualRM[r], m_dualST[r]);
	 }
      }
	 
      return &m_dualST[0];
#else      
      return m_masterSI->getRowPrice();
#endif
   }


   virtual void setObjBoundLB(const double thisBound){
      if(thisBound > m_nodeStats.objBest.first){
#ifdef DUAL_SMOOTHING
	 printf("Bound improved %g to %g, update duals\n",
		m_nodeStats.objBest.first, thisBound);
	 copy(m_dualST.begin(), m_dualST.end(), m_dual.begin());
#endif
      }
      DecompAlgo::setObjBoundLB(thisBound);
   }
   
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

      //---
      //--- do any parameter overrides of the defaults here
      //---    by default turn off gomory cuts for PC
      //---
      m_param.CutCglGomory = 0;

      //---
      //--- run init setup
      //---
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

      //---
      //--- do any parameter overrides of the defaults here
      //---    by default turn off gomory cuts for PC
      //---
      m_param.CutCglGomory = 0;

      //---
      //--- run init setup
      //---
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
