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

#ifndef MMKPMODEL_INCLUDED
#define MMKPMODEL_INCLUDED

//===========================================================================//
#include "UtilMacros.h"
//===========================================================================//
#include "DecompConstraintSet.h"
//===========================================================================//
#include "MMKP_Param.h"
#include "MMKP_Instance.h"
//===========================================================================//

//===========================================================================//
/*!
 * \class MMKPModel
 * A class for storage of the different models used in:
 *     Multi-Dimensional Mulit-Choice Knapsack Problem (MMKP).
 */

//===========================================================================//
class MMKP_Model {
private:
   /** Class id tag (for log / debugging). */
   const string m_classTag;

   
public:
   /* Read data instance. */
   int readInstance();


public:
   /** @name Access methods */
   inline const MMKP_Instance & getInstance() const {
      return m_instance;
   }
   inline const MMKP_Param & getParam() const {
      return m_appParam;
   }
   inline const double * getObjective() const {
      return m_objective;
   }
   inline DecompConstraintSet * getModelCore() const {
      return getModel(m_appParam.ModelNameCore);
   }
   inline DecompConstraintSet * getModelRelax() const {
      return getModel(m_appParam.ModelNameRelax);
   }

   inline DecompConstraintSet * getModel(string modelName) const {
      map<string, DecompConstraintSet*>::const_iterator it;
      it = m_models.find(modelName);
      if(it == m_models.end()){
         cout << "Error: model with name " << modelName << " not defined."
              << endl;
         assert(it != m_models.end());
         return NULL;
      }
      return it->second;
   }
   
public:
   /** @name Constructor and Destructor */

   /** Default constructor. Takes an instance of UtilParameters */
   MMKP_Model(UtilParameters & utilParam) : 
      m_classTag ("MMKP-MOD"),
      m_objective(NULL      ) {
      //---
      //--- get application parameters
      //---
      m_appParam.getSettings(utilParam);
      if(m_appParam.LogLevel >= 1)
         m_appParam.dumpSettings(); 
   }
   
   ~MMKP_Model() {
      UTIL_DELARR(m_objective);
      UtilDeleteMapPtr(m_models);
   };
};

#endif
