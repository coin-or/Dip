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

#ifndef GAP_DECOMPAPP3_INCLUDED
#define GAP_DECOMPAPP3_INCLUDED

//===========================================================================//
//---
//--- Version 3:
//---   for relaxation solver, use built-in MILP solver and dense format
//--- Version 4:
//---   for relaxation solver, use built-in MILP solver and sparse format
//--- 
//===========================================================================//


//===========================================================================//
#include "DecompApp.h"

//===========================================================================//
#include "GAP_Instance.h"
#include "GAP_DecompParam.h"

//===========================================================================//
/*!
 * \class GAP_DecompApp
 * A DecompApp for solving the 
 *     Generalized Assignment Problem (GAP).
 * 
 * \see
 * DecompApp
 *
 */

//===========================================================================//
class GAP_DecompApp : public DecompApp {
private:
   /** Class id tag (for log / debugging). */
   const string m_classTag;

   /** GAP problem instance data */
   GAP_Instance m_instance;

   /** Application specific parameters. */
   GAP_DecompParam m_appParam;  

   /** The model objective coefficients (original space). */
   double * m_objective;
   
   /** Store pointers to the various model constraint systems,
       so their memory can be deleted. */
   vector<DecompConstraintSet*> m_models;
   

public:
   /* @name Inherited (from virtual) methods. */
   /** Print an original column (format for this app). */
   void printOriginalColumn(const int   index, 
                            ostream   * os = &cout) const;

public:
   /** @name Helper functions (public). */   

   /** Guts of constructor. */
   void initializeApp(UtilParameters & utilParam);

   /** Helper methods for indexing. */
   inline const int getOffsetI(const int i) const {
      return i * m_instance.getNTasks();
   }
   inline const int getIndexIJ(const int i,
                               const int j) const {
      return (i * m_instance.getNTasks()) + j;
   }
   
   inline pair<int,int> getIndexInv(const int index) const {      
      return make_pair(index / m_instance.getNTasks(), 
                       index % m_instance.getNTasks());
   }

   /** Creation of the various model constraint systems. */ 
   int createModels();
   int createModelPartAP(DecompConstraintSet * model);
   int createModelPartKP(DecompConstraintSet * model);
   int createModelPartKP(DecompConstraintSet * model, 
                         int                   whichKnap);
   int createModelPartKP(DecompConstraintSet * model, 
                         vector<int>         & whichKnaps);

public:
   /** Some access methods to private data. */
   inline const GAP_Instance & getInstance() const {
      return m_instance;
   }
   inline const GAP_DecompParam & getParam() const {
      return m_appParam;
   }
   inline const double * getObjective() const {
      return m_objective;
   }
 
public:
   /** @name Constructor and Destructor */

   /** Default constructor. Takes an instance of UtilParameters */
   GAP_DecompApp(UtilParameters & utilParam) : 
      DecompApp   (utilParam),
      m_classTag  ("GAP-APP"),
      m_objective (NULL)
   {
      initializeApp(utilParam);
   }
   
   virtual ~GAP_DecompApp() {
      UTIL_DELARR(m_objective);
      UtilDeleteVectorPtr(m_models);
   };
};

#endif
