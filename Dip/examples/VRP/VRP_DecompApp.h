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

#ifndef VRP_DECOMPAPP_INCLUDED
#define VRP_DECOMPAPP_INCLUDED

// --------------------------------------------------------------------- //
#include "Decomp.h"
#include "DecompApp.h"
#include "VRP_Boost.h"
#include "VRP_CVRPsep.h"
#include "VRP_Instance.h"
#include "VRP_Param.h"
// --------------------------------------------------------------------- //
#define VRP_DECOMPAPP_USECONCORDE
#ifdef VRP_DECOMPAPP_USECONCORDE
#include "VRP_Concorde.h"
#endif
// --------------------------------------------------------------------- //

/*!
 * \class VRP_DecompApp
 * A DecompApp for solving the Traveling Salesman Problem.
 * 
 * \see
 * DecompApp
 *
 */

// --------------------------------------------------------------------- //
class VRP_DecompApp : public DecompApp {
private:
   /** Class id tag (for log / debugging). */
   const string m_classTag;

   /** Application specific parameters. */
   VRP_Param m_appParam;  

   /** Storage of TSP instance. */
   VRP_Instance    m_vrp;

   /** Interface class for CVRPSEP methods. */
   VRP_CVRPsep     m_cvrpSep;

   /** Interface class for Boost methods. */
   VRP_Boost       m_boost;

   /** Interface class for Concorde methods. */
#ifdef VRP_DECOMPAPP_USECONCORDE
   VRP_Concorde    m_concorde;
#endif

   /** The model objective coefficients (original space). */
   double * m_objective;

   /** The various model constraint systems used for different algos. */
   vector<DecompConstraintSet*>   m_models;
   DecompConstraintSet          * m_modelESPPRC;

public:
   /* @name Inherited (from virtual) methods. */

   /** Solve the relaxed problem. */
   virtual DecompSolverStatus solveRelaxed(const int          whichBlock,
                                           const double     * redCostX,
                                           const double       convexDual,
                                           DecompVarList    & varList);


   virtual int generateCuts(const double              * x, 
                            DecompCutList             & newCuts);

   //TODO: change this name... 
   virtual bool APPisUserFeasible(const double * x, 
                                  const int      n_cols,
                                  const double   tolZero);
   virtual void printOriginalColumn(const int   index, 
                                    ostream   * os) const;
   
      
   
public:
   /** @name Helper functions (public). */   

   void initializeApp(UtilParameters & utilParam);

   /* Create models. */
   void createModels();

   /* Create two-degree model. */
   void createModelTwoDegree(DecompConstraintSet * model);

   /* Create ESPPCC model. */
   void createModelESPPCC(DecompConstraintSet * model);

   const int diGraphIndex(int i, int j, int numVertices){
      return i*numVertices + j;
   }
   const pair<int,int> diGraphBothEnds(int index, int numVertices){
      int i = index / numVertices;
      return make_pair(i, index % numVertices);
   }

public:
   /** @name Constructor and Destructor */

   /** Default constructor. Takes an instance of UtilParameters */
   VRP_DecompApp(UtilParameters & utilParam) : 
      DecompApp(utilParam),
      m_classTag("VRP-APP"),
	m_objective(NULL),
	m_modelESPPRC(NULL)
   {
      initializeApp(utilParam);         
   }
   
   virtual ~VRP_DecompApp() {
      UtilDeleteVectorPtr(m_models);
      UTIL_DELARR(m_objective);
   };
};

#endif
