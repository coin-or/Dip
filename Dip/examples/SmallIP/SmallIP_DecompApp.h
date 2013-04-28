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

#ifndef SMALLIP_DECOMPAPP_INCLUDED
#define SMALLIP_DECOMPAPP_INCLUDED

//===========================================================================//
#include "DecompApp.h"

using namespace std;

//===========================================================================//
/*!
 * \class SmallIP_DecompApp
 * A DecompApp to illustrate a basic usage of Decomp.
 * 
 * \see
 * DecompApp
 */

//===========================================================================//
class SmallIP_DecompApp : public DecompApp {
private:
   /** Class id tag (for log / debugging). */
   const string m_classTag;

   /** The model objective coefficients (original space). */
   double * m_objective;

   /** The various model constraint systems used for different algos. */
   DecompConstraintSet m_modelPart1;
   DecompConstraintSet m_modelPart2;

public:
   /** @name Helper functions (public). */   

   /* Create models. */
   void createModels();
   
public:
   /* @name Inherited (from virtual) methods. */
   virtual int generateInitVars(DecompVarList & initVars);  

public:
   SmallIP_DecompApp(UtilParameters & utilParam) : 
      DecompApp  (utilParam),
      m_classTag ("SMALL-APP"),
      m_objective(NULL)
   {
      createModels();
   }
  
   virtual ~SmallIP_DecompApp() {
      UTIL_DELARR(m_objective);
   };
};

#endif
