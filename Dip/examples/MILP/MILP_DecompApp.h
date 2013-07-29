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

#ifndef MILP_DECOMPAPP_INCLUDED
#define MILP_DECOMPAPP_INCLUDED

//===========================================================================//
#include "DecompApp.h"
#include "MILP_Param.h"
//===========================================================================//
#include "CoinMpsIO.hpp"
//===========================================================================//

/*!
 * \class MILP_DecompApp
 * A DecompApp to illustrate a basic usage of Decomp.
 * 
 * \see
 * DecompApp
 */
//===========================================================================//
class MILP_DecompApp : public DecompApp{
private:
   /** Class id tag (for log / debugging). */
   const string m_classTag;
   
   /** MPS object for reading MILP instances. */
   CoinMpsIO m_mpsIO;
   
   /** Application specific parameters. */
   MILP_Param m_appParam;

   /** The model objective coefficients (original space). */
   double * m_objective;
   
   /** The model constraint systems used for different algos. */
   DecompConstraintSet m_modelRandCore;
   DecompConstraintSet m_modelRandRelax;
      
private:
   /** @name Helper functions (private). */

   /** Initialize application. */
   void initializeApp(UtilParameters & utilParam);

   /** Create model part. */
   void createModels();
   
public:
   /** @name Constructor and Destructor */
   MILP_DecompApp(UtilParameters & utilParam) : 
      DecompApp  (utilParam),
      m_classTag ("MILP-APP"),
      m_objective(NULL)
   {
      initializeApp(utilParam); //can there be a default?
   }
   
   virtual ~MILP_DecompApp() {
      UTIL_DELARR(m_objective);
   }
};

#endif
