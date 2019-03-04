//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Authors: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)      //
//          Ted Ralphs, Lehigh University (ted@lehigh.edu)                   //
//          Jiadong Wang, Lehigh University (jiw408@lehigh.edu)              //
//                                                                           //
// Copyright (C) 2002-2019, Lehigh University, Matthew Galati, and Ted Ralphs//
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
   double* objective;

   /** Model constraint systems. */
   vector<DecompConstraintSet*> m_models;

   DecompConstraintSet* modelRelax;
   DecompConstraintSet* modelCore;

public:
   /** @name Helper functions (public). */


   /** Initialize application. */
   void initializeApp();

   /* Create models. */
   void createModels();
   void createModelCore(DecompConstraintSet* model);
   void createModelRelax(DecompConstraintSet* model,
                         int                   commId);
   void createModelRelaxSparse(DecompConstraintSet* model,
                               int                   commId);

public:
   /** @name Constructor and Destructor */

   /** Default constructor. Takes an instance of UtilParameters */
   MCF_DecompApp(UtilParameters& utilParam) :
      DecompApp   (utilParam),
      m_classTag  ("MCF-APP"),
      objective   (   NULL  ),
      modelRelax  (   NULL  ),
      modelCore   (   NULL  )
   {
      //---
      //--- get application parameters
      //---
      m_appParam.getSettings(utilParam);
      
      if (m_appParam.LogLevel >= 1) {
	 m_appParam.dumpSettings(m_osLog);
      }
      
      initializeApp();
   }

   virtual ~MCF_DecompApp() {
      UTIL_DELARR(objective);
      UtilDeleteVectorPtr(m_models);
      UTIL_DELPTR(modelCore);
   };
};

#endif
