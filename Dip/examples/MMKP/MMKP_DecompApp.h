//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, Lehigh University                                 //
//                                                                           //
// Copyright (C) 2002-2013, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

#ifndef MMKP_DECOMPAPP_INCLUDED
#define MMKP_DECOMPAPP_INCLUDED

//===========================================================================//
#include "DecompApp.h"

//===========================================================================//
#include "MMKP_Instance.h"
#include "MMKP_MCKnap.h"
#include "MMKP_Param.h"
#include "MMKP_MemPool.h"
//===========================================================================//

//===========================================================================//
/*!
 * \class MMKP_DecompApp
 * A DecompApp for solving the 
 *     Multi-Dimensional Mulit-Choice Knapsack Problem (MMKP).
 * 
 * \see
 * DecompApp
 *
 */

//===========================================================================//
class MMKP_DecompApp : public DecompApp {
private:
   /** Class id tag (for log / debugging). */
   const string m_classTag;

   /** Application specific parameters. */
   MMKP_Param m_appParam;  

   /** MMKP problem instance data */
   MMKP_Instance m_instance;

   /** MMKP_MCKnap object for each knapsack row. */
   vector<MMKP_MCKnap*> m_mcknap;

   /** The model objective coefficients (original space). */
   double * m_objective;

   /** The various model constraint systems used for different algos. */
   vector<DecompConstraintSet*> m_models;
   
   /** Auxiliary memory storage. */
   MMKP_MemPool m_auxMemPool;
 
public:
   /* @name Inherited (from virtual) methods. */

   /** Solve the relaxed problem. */
   virtual DecompSolverStatus solveRelaxed(const int          whichBlock,
                                           const double     * redCostX,
                                           const double       convexDual,
                                           DecompVarList    & varList);

   /** Print an original column (format for this app). */
   virtual void printOriginalColumn(const int   index, 
                                    ostream   * os = &cout) const;

public:
   /** @name Helper functions (public). */   

   /** Initialize application. */
   void initializeApp(UtilParameters & utilParam);

   /* Create models. */
   void createModels();

   /* Create MCKP model. */
   void createModelPartMCP(DecompConstraintSet * model);
   void createModelPartMCKP(DecompConstraintSet * model,
                            int                   whichKnap = 0);
   void createModelPartMC2KP(DecompConstraintSet * model,
                             int                   whichKnap1,
                             int                   whichKnap2);
   void createModelPartMCKP(DecompConstraintSet * model,
                            vector<int>         & whichKnaps);
   void createModelPartMMKPHalf(DecompConstraintSet * model);
   void createModelPartMMKP(DecompConstraintSet * model);
   
   /* Create MDKP model. */
   void createModelPartMDKPCompl(DecompConstraintSet * model,
                                 int                   whichKnap = 0);
   void createModelPartMDKP(DecompConstraintSet * model);
   void createModelPartMDKPHalf(DecompConstraintSet * model);
   void createModelPartMDKP(DecompConstraintSet * model,
                            vector<int>         & whichKnaps);

public:
   /** @name Constructor and Destructor */
   
   /** Default constructor. Takes an instance of UtilParameters */
   MMKP_DecompApp(UtilParameters & utilParam) :
      DecompApp   (utilParam),
      m_classTag  ("MMKP-APP"),
      m_mcknap    (),
      m_objective (NULL),
      m_auxMemPool() 
   {
      initializeApp(utilParam);                    
   }
   
   virtual ~MMKP_DecompApp() {
      UtilDeleteVectorPtr(m_mcknap);
      UtilDeleteVectorPtr(m_models);
      UTIL_DELARR(m_objective);
   };
};

#endif
