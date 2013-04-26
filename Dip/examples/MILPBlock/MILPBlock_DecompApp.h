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

#ifndef MILPBLOCK_DECOMPAPP_INCLUDED
#define MILPBLOCK_DECOMPAPP_INCLUDED

//===========================================================================//
#include "DecompApp.h"
#include "MILPBlock_Param.h"
//===========================================================================//
#include "CoinMpsIO.hpp"
//===========================================================================//

/*!
 * \class MILPBlock_DecompApp
 * A DecompApp to illustrate a basic usage of Decomp.
 * 
 * \see
 * DecompApp
 */
//===========================================================================//
class MILPBlock_DecompApp : public DecompApp{
private:
   /** Class id tag (for log / debugging). */
   const string m_classTag;
   
   /** MPS object for reading MILPBlock instances. */
   CoinMpsIO m_mpsIO;
   
   /** Application specific parameters. */
   MILPBlock_Param m_appParam;

   /** The model objective coefficients (original space). */
   double * m_objective;
   
   /** The model constraint systems used for different algos. */
   DecompConstraintSet *          m_modelC;
   map<int, DecompConstraintSet*> m_modelR;

   /** Definition of blocks (by rows). */
   map<int, vector<int> > m_blocks;

private:

   /* @name Inherited (from virtual) methods. */

   /** Generate init columns. */
   virtual int generateInitVars(DecompVarList & initVars);

   /** @name Helper functions (private). */

   /** Initialize application. */
   void initializeApp(UtilParameters & utilParam);

   /** Create model parts. */
   void                  createModels();
   DecompConstraintSet * createModelPart(const int   nRowsPart,
                                         const int * rowsPart);
   void createModelPart(DecompConstraintSet * model,
			const int             nRowsPart,
			const int           * rowsPart);
   void createModelPartSparse(DecompConstraintSet * model,
			      const int             nRowsPart,
			      const int           * rowsPart);   
   void                  createModelMasterOnlys(vector<int> & masterOnlyCols);
   void                  readInitSolutionFile(DecompVarList & initVars);


   /** Read block file. */
   void readBlockFile();

   /** Find the active columns for some block. */
   void findActiveColumns(const vector<int> & rowsPart,
                          set<int>          & activeColsSet);
   
public:
   /** User access methods. */
   
   /** Get instance name. */
   const string getInstanceName(){
      return m_appParam.Instance;
   }
  
public:
   /** @name Constructor and Destructor */
   MILPBlock_DecompApp(UtilParameters & utilParam) : 
      DecompApp  (utilParam),
      m_classTag ("MILPB-APP"),
      m_objective(NULL)
   {
      initializeApp(utilParam); //can there be a default?
   }
   
   virtual ~MILPBlock_DecompApp() {
      UTIL_DELARR(m_objective);
      UTIL_DELPTR(m_modelC);
      UtilDeleteMapPtr(m_modelR);
   }
};

#endif
