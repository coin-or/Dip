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

#ifndef ATM_DECOMPAPP_INCLUDED
#define ATM_DECOMPAPP_INCLUDED

//===========================================================================//
#include "DecompApp.h"
//===========================================================================//
#include "ATM_Param.h"
#include "ATM_Instance.h"


//===========================================================================//
/*!
 * \class ATM_DecompApp
 * A DecompApp for solving the 
 *     ATM Cash Management Problem (ATM).
 * 
 * \see
 * DecompApp
 *
 */

//===========================================================================//
class ATM_DecompApp : public DecompApp {
private:
   /** Class id tag (for log / debugging). */
   const string m_classTag;

   /** ATM problem instance data */
   ATM_Instance m_instance;
 
   /** Application specific parameters. */
   ATM_Param m_appParam;  

   /** The model objective coefficients (original space). */
   double * m_objective;

   /** The various model constraint systems used for different algos. */
   vector<DecompConstraintSet*> m_models;
   
public:
   /* @name Inherited (from virtual) methods. */
   virtual bool APPisUserFeasible(const double * x,
                                  const int      nCols,
                                  const double   tolZero);

public:
   /** @name Helper functions (public). */   

   /** Initialize application. */
   void initializeApp(UtilParameters & utilParam);

   /* Create models. */
   void createModels();


   //create model
   int createConZtoX(DecompConstraintSet * model,
		     const int             atmIndex);
   int createConPickOne(DecompConstraintSet * model,
			const int             atmIndex);
   int createConCount(DecompConstraintSet * model,
		      const int             atmIndex);

   DecompConstraintSet * createModelCore1(bool includeCount = true);
   DecompConstraintSet * createModelRelax1(const int a,
					   bool      includeCount = true);

   DecompConstraintSet * createModelCore2();
   DecompConstraintSet * createModelRelax2(const int d);

   DecompConstraintSet * createModelCoreCount();
   DecompConstraintSet * createModelRelaxCount();
   
   //counts
   inline const int numCoreCols() const {
      //---    for a in A, d in D:
      //---       f+[a,d], f-[a,d] >= 0
      //---                f-[a,d] <= w[a,d]
      //---       v[a,d] in {0,1}
      //---    for a in A:
      //---       x2[a] in [0,1], x3[a] >= 0
      //---    for a in A, t in T
      //---       x1[a,t] in {0,1}, z[a,t] >= 0
      int nAtms  = m_instance.getNAtms();
      int nPairs = m_instance.getNPairs();
      int nCols  = (2 * (getNAtmsSteps() + nPairs + nAtms)) + nPairs;
      return nCols;
   }
   
   inline const int getNAtmsSteps() const {
      if(m_appParam.UseTightModel)
	 return m_instance.getNAtms() * (m_appParam.NumSteps+1);
      else
	 return m_instance.getNAtms() * m_appParam.NumSteps;
   }
   
   //---
   //--- Columns
   //---   x1[a,t] (binary)
   //---   z [a,t]
   //---   f+[a,d]
   //---   f-[a,d]
   //---   x2[a]
   //---   x3[a]
   //---   v[a,d] (binary)
   //---
   void createModelColumns(DecompConstraintSet * model,
			   const int             atmIndex  = -1,
                           const int             dateIndex = -1);

   inline const int getColOffset_x1() const {
      return 0;
   }
   inline const int getColOffset_z() const {
      return getColOffset_x1() + getNAtmsSteps();
   }
   inline const int getColOffset_fp() const {
      return getColOffset_z() + getNAtmsSteps();
   }
   inline const int getColOffset_fm() const {
      return getColOffset_fp() + m_instance.getNPairs();
   }
   inline const int getColOffset_x2() const {
      return getColOffset_fm() + m_instance.getNPairs();
   }
   inline const int getColOffset_x3() const {
      return getColOffset_x2() + m_instance.getNAtms();
   }
   inline const int getColOffset_v() const {
      return getColOffset_x3() + m_instance.getNAtms();
   }

   inline const int colIndex_x1(const int a, const int t) const {
      int tLen = m_appParam.NumSteps;
      if(m_appParam.UseTightModel)
	 tLen++;
      return getColOffset_x1() + ((a * tLen) + t);
   }
   inline const int colIndex_z(const int a, const int t) const {
      int tLen = m_appParam.NumSteps;
      if(m_appParam.UseTightModel)
	 tLen++;
      return getColOffset_z() + ((a * tLen) + t);
   }
   inline const int colIndex_fp(const int pairIndex) const {
      return getColOffset_fp() + pairIndex;         
   }
   inline const int colIndex_fm(const int pairIndex) const {
      return getColOffset_fm() + pairIndex;
   }
   inline const int colIndex_x2(const int a) const {
      return getColOffset_x2() + a;
   }
   inline const int colIndex_x3(const int a) const {
      return getColOffset_x3() + a;
   }
   inline const int colIndex_v(const int pairIndex) const {
      return getColOffset_v() + pairIndex;         
   }

   //column names
   void addColumnNamesA(DecompConstraintSet * model,
			const string          prefix,
			const int             offset);
   void addColumnNamesAT(DecompConstraintSet * model,
			 const string          prefix,
			 const int             offset);
   void addColumnNamesAD(DecompConstraintSet * model,
			 const string          prefix,
			 const int             offset);
   
 
public:
   /** @name Constructor and Destructor */

   /** Default constructor. Takes an instance of UtilParameters */
   ATM_DecompApp(UtilParameters & utilParam) : 
      DecompApp   (utilParam),
      m_classTag  ("ATM-APP"),
      m_objective (NULL)
   {
      initializeApp(utilParam);
   };
   
   virtual ~ATM_DecompApp() {
      UTIL_DELARR(m_objective);
      UtilDeleteVectorPtr(m_models);
   };
};

#endif
