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

#ifndef MAD_DECOMPAPP_INCLUDED
#define MAD_DECOMPAPP_INCLUDED

#define __MAD_USE_QUALEX__
#define __MAD_USE_CLIQUER__

// --------------------------------------------------------------------- //
#include "DecompApp.h"
#include "MAD_MemPool.h"
#include "MAD_DecompParam.h"
#include "MAD_DecompDebug.h"
#include "CoinLpIO.hpp"
// --------------------------------------------------------------------- //

#ifdef __MAD_USE_CLIQUER__
#include "MAD_Cliquer.h"
#define cliquer_t MAD_Cliquer
#else
#define cliquer_t void
#endif

#ifdef __MAD_USE_QUALEX__
#include "MAD_Qualex.h"
#define MAD_Qualex MAD_Qualex
#else
#define MAX_Qualex void
#endif

// --------------------------------------------------------------------- //
typedef enum MAD_HeurSortDir{
   INCREASING,
   DECREASING
};

// --------------------------------------------------------------------- //
typedef struct GreedyPoint GreedyPoint;
struct GreedyPoint {
   double * solution;
   double   solValueOrigCost;    
   double   solValueRedCost; 
};


// --------------------------------------------------------------------- //

/*!
 * \class MAD_DecompApp
 * A DecompApp for solving the Matrix Decomposition Problem (MAD).
 * 
 * \see
 * DecompApp
 *
 */

// --------------------------------------------------------------------- //
class MAD_DecompApp : public DecompApp {
 private:
   /** Class id tag (for log / debugging). */
   const string m_classTag;

 private:
   /** MAD problem instance data */
   CoinLpIO  m_instance;
   int       m_nOrigRows;
   int       m_beta;
   int       m_kappa;

   /** Storage for cliquer graph */
   cliquer_t * m_cliquer;

   /** Storage for MAD_Qualex interface */
   MAD_Qualex * m_qualex;

   //THINK about data structures! if have cliquer around vs qualex
   //or just use generic boost? in the case that we have neither?
   graph_t   * m_conflictGraph;

   /** Auxiliary memory storage. */
   MAD_MemPool m_auxMemPool;

 
 protected:
   /** Application specific parameters. */
   MAD_DecompParam m_appParam;  
   
 public:
   /** Model types */
   enum ModelType {
      /** Relaxation is a Maximum Weighted Clique Problem */
      MODEL_CLIQUE
   };
   
 public:
   /* @name Inherited (from pure virtual) methods. */

   /** Create the application model(s). */
   //TODO: model object? vs objCoeff?
   void APPcreateModel(double                        *& objCoeff,
                       map<int, DecompConstraintSet*> & modelCore,
                       map<int, vector<DecompConstraintSet* > > & modelRelax);

 public:
   /* @name Inherited (from virtual) methods. */

   /** Solve the relaxed problem. */
   //TOOD: too messy?
   DecompStatus APPsolveRelaxed(const int             whichModel,
				const double        * redCostX,
				const double        * origCost,
				const double          alpha,
				const int             n_origCols,
				const bool            checkRC,
				const bool            checkDup,
				bool                & isExact,
				OsiSolverInterface  * m_subprobSI,
				list<DecompVar*>    & vars);
   
   /** TODO */
   int APPheuristics(const double            * xhat,
		     const double            * origCost,
		     vector<DecompSolution*> & xhatIPFeas);   

    /** TODO */
    int generateInitVars(DecompVarList & initVars);

   /** Print an original column (format for this app). */
   void printOriginalColumn(const int   index, 
                            ostream   * os = &cout) const;
   
   /** Print the full original column (format for this app). */
   void printOriginalSolution(const int      n_cols, 
                              const double * solution, 
                              ostream      * os) const;

 public:
   /** @name Helper functions (public). */   
   
   /** Global index for column x[i,b]. */
   inline const int xIndex(const int i, 
                           const int b) const{ 
      return i * m_beta + b;
   }

   /** Return the indices (i,b) given the index. */
   inline pair<int,int> xIndexInv(const int index) const{      
     return make_pair(index / m_beta, index % m_beta);
   }

   /** Are the two vectors orthogonal? */
   int isOrtho(const int    * rowInd1,
               const int    * rowInd2,
               const int      rowLen1,
               const int      rowLen2);

   /** Guts of constructor. */
   void initializeApp(UtilParameters & utilParam);

   /** TODO */
   int heuristicGreedy(vector<GreedyPoint>     & greedyPoints,
		       const MAD_HeurSortDir     sortDir,
		       const double            * sortValues,
		       const double            * origCost,
		       const double            * redCost = NULL);
   int heuristicGreedy(const MAD_HeurSortDir     sortDir,
		       const double            * sortValues,
		       const double            * origCost,
		       vector<DecompSolution*> & solVec);
      



   void printRowMarks(const int * rowInd,
                      const int   rowLen) const;


   /** Access method for member data. */
   inline const int getNOrigRows() const { return m_nOrigRows; }
   inline const int getBeta()      const { return m_beta;      }
   inline const int getKappa()     const { return m_kappa;     }


   
 private:
   /** @name Copy Constructors
    *
    * Disable the default copy constructors.
    *
    */
   MAD_DecompApp(const MAD_DecompApp &);
   MAD_DecompApp & operator=(const MAD_DecompApp &);

 public:
   /** @name Constructor and Destructor */

   /** Default constructor. Takes an instance of UtilParameters */
   MAD_DecompApp(UtilParameters & utilParam) : 
      DecompApp(utilParam), 
      m_classTag("MAD-APP"),
      m_beta(0),
      m_kappa(0),
      m_cliquer(0),
      m_conflictGraph(0),
      m_auxMemPool(),
//      m_bestKnownLB(-1.0e20),
//      m_bestKnownUB( 1.0e20),
      m_appParam()
      {
         initializeApp(utilParam);
      }
   
   virtual ~MAD_DecompApp() {  
      UTIL_DELPTR(m_cliquer);
      UTIL_DELPTR(m_qualex);
      graph_free(m_conflictGraph);
   };
};

#endif
