//===========================================================================//
// This file is part of the DIP Solver Framework.                            //
//                                                                           //
// DIP is distributed under the Eclipse Public License as part of the        //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Authors: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)      //
//          Ted Ralphs, Lehigh University (ted@lehigh.edu)                   //
//          Jiadong Wang, Lehigh University (jiw408@lehigh.edu)              //
//                                                                           //
// Copyright (C) 2002-2015, Lehigh University, Matthew Galati, Ted Ralphs    //
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
   std::string m_classTag;

   /**
    * Dual vector.
    *
    * The stabilized dual (copied from m_dualST) when the bound has improved.
    * This dual vector is the one used in reduced-cost calculations when
    * using a stabilized dual method (m_param.DualStab > 0).
    */
   std::vector<double> m_dual;

   /**
    * Dual vector from restricted master.
    *
    * A copy of the dual vector from the restricted master.
    */
   std::vector<double> m_dualRM;

   /**
    * Dual vector stabilized.
    *
    * The stabilized dual from dual stabilization method.
    */
   std::vector<double> m_dualST;

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
   virtual void createMasterProblem(DecompVarList& initVars) {
      DecompAlgo::createMasterProblem(initVars);
   }
   virtual int generateVars(DecompVarList&     newVars,
			    double&            mostNegReducedCost) {
      return DecompAlgo::generateVars(newVars, mostNegReducedCost);
   }
   virtual void phaseInit(DecompPhase& phase);

   /**
    * Get current dual solution for master problem.
    *
    * When using dual stabilization, this comes from the stabilized dual
    * vector (m_dualST). Otherwise, it comes from m_dualSolution (which
    * comes directly from the LP solver).
    *
    */
   virtual const double* getMasterDualSolution() const {
      //---
      //--- return the duals to be used in pricing step
      //---
      if (m_param.DualStab) {
         return &m_dualST[0];
      } else {
         return &m_dualSolution[0];
      }
   }

   /**
    * Adjust the current dual solution for master problem.
    *
    * When using dual stabilization, this adjusts based on Wengtes smoothing.
    */
   virtual void adjustMasterDualSolution();



   /**
    *
    */
   virtual void setObjBound(const double thisBound,
                            const double thisBoundUB) {
      UtilPrintFuncBegin(m_osLog, m_classTag,
                         "setObjBound()", m_param.LogDebugLevel, 2);

      if (m_param.DualStab) {
         if (thisBound > (m_nodeStats.objBest.first + DecompEpsilon)) {
            //(*m_osLog) << "Bound improved " << m_nodeStats.objBest.first
            //       << " to " << thisBound << " , update duals" << endl;
            copy(m_dualST.begin(), m_dualST.end(), m_dual.begin());
         }
      }

      DecompAlgo::setObjBound(thisBound, thisBoundUB);
      UtilPrintFuncEnd(m_osLog, m_classTag,
                       "setObjBound()", m_param.LogDebugLevel, 2);
   }

   /**
    * Set the current integer bound and update best/history.
    */
   virtual inline void setObjBoundIP(const double thisBound) {
      DecompAlgo::setObjBoundIP(thisBound);
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
   void addCutsToPool(const double*    x,
                      DecompCutList& newCuts,
                      int&            n_newCuts);

   //TODO
   void phaseDone();
   int  addCutsFromPool();
   void solveMasterAsMIP();
   void solveMasterAsMIPSym(DecompSolverResult* result);
   void solveMasterAsMIPCbc(DecompSolverResult* result);
   void solveMasterAsMIPCpx(DecompSolverResult* result);
   void solveMasterAsMIPGrb(DecompSolverResult* result);
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

   std::vector<double>& getDualBest() {
      return m_dual;
   }
   std::vector<double>& getDualRMP() {
      return m_dualRM;
   }


   /**
    * Default constructors.
    */
   DecompAlgoPC(DecompApp*       app,
                UtilParameters& utilParam,
                bool             doSetup    = true,
		const DecompAlgoType   algo = PRICE_AND_CUT) :
      DecompAlgo(algo, app, utilParam),
      m_classTag("D-ALGOPC") {
      //---
      //--- do any parameter overrides of the defaults here
      //---    by default turn off gomory cuts for PC
      //---
      m_param.CutCglGomory = 0;
   }

   /**
    * Destructor.
    */
   ~DecompAlgoPC() {}
   /**
    * @}
    */
};

#endif
