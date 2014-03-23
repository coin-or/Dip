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

//this is a lambda_s variable
//assumption throughout is lambda is {0,1}

#ifndef DECOMP_VAR_INCLUDED
#define DECOMP_VAR_INCLUDED

#include "DecompPortable.h"
#include "UtilHash.h"

#include <iostream>
using namespace std;

class DecompApp;

// --------------------------------------------------------------------- //
class DecompVar {
private:
   DecompVar(const DecompVar&);
   DecompVar& operator=(const DecompVar&);

public:
   //THINK: or user overriden way to store it (like a tree)
   //and a function which says expandVarToCol - just as in cuts
   CoinPackedVector m_s;//this is the var in terms of x-space

private:
   //TODO: lb, ub, "type"?
   double           m_origCost;
   double           m_redCost; //(c - uA'')s - alpha
   int              m_effCnt;  //effectiveness counter
   string           m_strHash;

public:
   inline double getOriginalCost()  const {
      return m_origCost;
   }
   inline double getReducedCost()   const {
      return m_redCost;
   }
   inline double getEffectiveness() const {
      return m_effCnt;
   }
   inline double getLowerBound()    const {
      return 0.0;         //TODO
   }
   inline double getUpperBound()    const {
      return DecompInf;   //TODO
   }
   inline string getStrHash()       const {
      return m_strHash;
   }

   inline void   setReducedCost(const double redCost) {
      m_redCost = redCost;
   }
   bool   isEquivalent(const DecompVar& dvar) {
      return m_s.isEquivalent(dvar.m_s);
   }

   void fillDenseArr(int      len,
                     double* arr);

public:
   virtual void  print(ostream*    os  = &cout,
                       DecompApp* app = 0) const;

public:
   DecompVar(const vector<int>    & ind,
             const vector<double> & els,
             const double           redCost,
             const double           origCost) :
      m_s       (),
      m_origCost(origCost),
      m_redCost (redCost),
      m_effCnt  (0),
      m_strHash () {
      if (ind.size() > 0) {
         m_s.setVector(static_cast<int>(ind.size()),
                       &ind[0], &els[0], DECOMP_TEST_DUPINDEX);
         m_strHash = UtilCreateStringHash(static_cast<int>(ind.size()),
                                          &ind[0], &els[0]);
      }
   }

   DecompVar(const int              len,
             const int*             ind,
             const double*          els,
             const double           redCost,
             const double           origCost) :
      m_s       (),
      m_origCost(origCost),
      m_redCost (redCost),
      m_effCnt  (0),
      m_strHash () {
      if (len > 0) {
         m_s.setVector(len, ind, els, DECOMP_TEST_DUPINDEX);
         m_strHash = UtilCreateStringHash(len, ind, els);
      }
   }

   DecompVar(const int              denseLen,
             const double*          denseArray,
             const double           redCost,
             const double           origCost) :
      m_s       (DECOMP_TEST_DUPINDEX),
      m_origCost(origCost),
      m_redCost (redCost),
      m_effCnt  (0),
      m_strHash () {
      UtilPackedVectorFromDense(denseLen, denseArray, DecompEpsilon, m_s);

      if (m_s.getNumElements() > 0) {
         m_strHash = UtilCreateStringHash(denseLen, denseArray);
      }
   }

   virtual ~DecompVar() {};
};

#endif
