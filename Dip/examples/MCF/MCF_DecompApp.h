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

#ifndef MCF_DECOMPAPP_INCLUDED
#define MCF_DECOMPAPP_INCLUDED

//===========================================================================//
#include "DecompApp.h"

//===========================================================================//
#include "MCF_Instance.h"
#include "MCF_Param.h"
//===========================================================================//

//===========================================================================//
/*!
 * \class MCF_DecompApp
 * A DecompApp for solving the 
 *     (Integer) Multi-Commodity Flow Problem (MCF)
 * 
 * \see
 * DecompApp
 *
 */

//===========================================================================//
class MCF_DecompApp : public DecompApp {
private:
   /** Class id tag (for log / debugging). */
   const string m_classTag;

   /** Application specific parameters. */
   MCF_Param m_appParam;  

   /** MCF problem instance data */
   MCF_Instance m_instance;

   /** The model objective coefficients (original space). */
   double * m_objective;

   /** Model constraint systems. */
   vector<DecompConstraintSet*> m_models;

public:
   /** @name Helper functions (public). */   

   /** Initialize application. */
   void initializeApp(UtilParameters & utilParam);

   /* Create models. */
   void createModels();
   void createModelCore(DecompConstraintSet * model);
   void createModelRelax(DecompConstraintSet * model,
                         int                   commId);
   void createModelRelaxSparse(DecompConstraintSet * model,
                               int                   commId);

public:
   /** @name Constructor and Destructor */
   
   /** Default constructor. Takes an instance of UtilParameters */
   MCF_DecompApp(UtilParameters & utilParam) :
      DecompApp   (utilParam),
      m_classTag  ("MCF-APP"),
      m_objective (NULL)
   {
      initializeApp(utilParam);                    
   }
   
   virtual ~MCF_DecompApp() {
      UTIL_DELARR(m_objective);
      UtilDeleteVectorPtr(m_models);
   };
};

#endif
