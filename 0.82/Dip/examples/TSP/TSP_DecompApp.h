//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, Lehigh University                                 //
//                                                                           //
// Copyright (C) 2002-2007, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

#ifndef TSP_DECOMPAPP_INCLUDED
#define TSP_DECOMPAPP_INCLUDED

//===========================================================================//
#include "Decomp.h"
#include "DecompApp.h"
#include "TSP_Instance.h"
#include "TSP_Param.h"
//===========================================================================//

/*!
 * \class TSP_DecompApp
 * A DecompApp for solving the Traveling Salesman Problem.
 * 
 * \see
 * DecompApp
 *
 */

// --------------------------------------------------------------------- //
class TSP_DecompApp : public DecompApp {
private:
   /** Class id tag (for log / debugging). */
   const string m_classTag;

   /** Application specific parameters. */
   TSP_Param m_appParam;  

   /** Storage of TSP instance. */
   TSP_Instance m_tsp;

   /** The model objective coefficients (original space). */
   double * m_objective;

   /** The various model constraint systems used for different algos. */
   vector<DecompConstraintSet*> m_models;
   
public:
   /* @name Inherited (from virtual) methods. */

   /** Solve the relaxed problem. */
   virtual DecompSolverStatus solveRelaxed(const int          whichBlock,
					   const double     * redCostX,
					   const double       convexDual,
					   DecompVarList    & varList);
   
   virtual int generateCuts(const double  * x, 
			    DecompCutList & newCuts);
   
   virtual bool APPisUserFeasible(const double * x, 
				  const int      n_cols,
				  const double   tolZero);
   virtual void printOriginalColumn(const int   index, 
				    ostream   * os) const;
   
      
private: 
   /** @name Helper functions (private). */   

   /** Guts of constructor. */
   void initializeApp(UtilParameters & utilParam);

   /** Create models. */
   void createModels();
   
   /** Create the two-matching constraints. */
   void createModel2Match(DecompConstraintSet * modelCS);
   
   /** Create trivial subtour elimination constraints. */
   void createModelTrivialSEC(DecompConstraintSet * modelCS);

   int generateCutsSubtour(DecompCutList & newCuts);

   void solveOneTree(const double               * cost, 
                     const double                 alpha,
                     vector< pair<int,double> > & edge_cost,
                     DecompVarList              & vars,
                     Graph                      & g);
   
public:
   /** @name Constructor and Destructor */

   /** Default constructor. Takes an instance of UtilParameters */
   TSP_DecompApp(UtilParameters & utilParam) : 
      DecompApp(utilParam),
      m_classTag("TSP-APP"),
      m_objective(NULL)
   {
      initializeApp(utilParam);         
   }
   
   virtual ~TSP_DecompApp() {
      UtilDeleteVectorPtr(m_models);
      UTIL_DELARR(m_objective);
   };
};

#endif



