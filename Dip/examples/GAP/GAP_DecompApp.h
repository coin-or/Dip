//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, Lehigh University                                 //
//                                                                           //
// Copyright (C) 2002-2013 Lehigh University, Matthew Galati, and Ted Ralphs //
// All Rights Reserved.                                                      //
//===========================================================================//

#ifndef GAP_DECOMPAPP_INCLUDED
#define GAP_DECOMPAPP_INCLUDED

// --------------------------------------------------------------------- //
#include "DecompApp.h"

// --------------------------------------------------------------------- //
#include "GAP_Instance.h"
#include "GAP_KnapPisinger.h"
#include "GAP_DecompParam.h"

// --------------------------------------------------------------------- //
/*!
 * \class GAP_DecompApp
 * A DecompApp for solving the 
 *     Generalized Assignment Problem (GAP).
 * 
 * \see
 * DecompApp
 *
 */

// --------------------------------------------------------------------- //
class GAP_DecompApp : public DecompApp {
private:
   /** Class id tag (for log / debugging). */
   const string m_classTag;

   /** GAP problem instance data */
   GAP_Instance m_instance;

   /** Application specific parameters. */
   GAP_DecompParam m_appParam;  

   /** GAP_Knapsack object for each knapsack row. */
   vector<GAP_KnapPisinger*> m_knap;
   
   /** The model objective coefficients (original space). */
   double * m_objective;
   
   /** The various model constraint systems used for different 
       algorithms, keyed by a unique string (model name). */
   map<string, DecompConstraintSet*> m_models;
   

public:
   /* @name Inherited (from virtual) methods. */
   /** Solve the relaxed problem. */
   DecompSolverStatus solveRelaxed(const int             whichBlock,
				   const double        * redCostX,
				   const double          convexDual,
				   list<DecompVar*>    & vars);   

   /** Print an original column (format for this app). */
   void printOriginalColumn(const int   index, 
                            ostream   * os = &cout) const;

public:
   /** @name Helper functions (public). */   

   /** Guts of constructor. */
   void initializeApp(UtilParameters & utilParam);

   /** TODO comment */
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

   /** TODO comment */
   int createModels();
   int createModelPartAP(DecompConstraintSet * model);

public:
   inline const GAP_Instance & getInstance() const {
      return m_instance;
   }
   inline const GAP_DecompParam & getParam() const {
      return m_appParam;
   }
   inline const double * getObjective() const {
      return m_objective;
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
   GAP_DecompApp(UtilParameters & utilParam) : 
      DecompApp   (utilParam),
      m_classTag  ("GAP-APP"),
      m_objective (NULL)
   {
      initializeApp(utilParam);
   }
   
   virtual ~GAP_DecompApp() {
      UtilDeleteVectorPtr(m_knap);
      UTIL_DELARR(m_objective);
      UtilDeleteMapPtr(m_models);
   };
};

#endif
