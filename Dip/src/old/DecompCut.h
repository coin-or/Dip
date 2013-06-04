//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)       //
//                                                                           //
// Copyright (C) 2002-2007, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//


#ifndef DECOMP_CUT_INCLUDED
#define DECOMP_CUT_INCLUDED

//a cut in terms of the original x variables

//by default it could look for OSI-CGL cuts? for CP procedures
//or the user can provide its own cut generator - force them to return
//a common cut structure?

//the user might have some compact way to store its cut, for example
//with TSP cuts, they just want to store the customer set, and then
//they override a function called expandCutToRow which tells this base
//class how to expand into the LP

#include "UtilHash.h"
#include "DecompPortable.h"

#include "CoinError.hpp"
#include <iostream>
using namespace std;

class DecompCut {
private:
   DecompCut(const DecompCut&);
   DecompCut& operator=(const DecompCut&);

private:
   double           m_lb;        //row lower bound
   double           m_ub;        //row upper bound
   //THINK, or they can stick in as sense/rhs
   double           m_violation; //current violation
   int              m_effCnt;    //effectiveness counter

protected:
   string           m_strHash;
   //TODO - use distance instead of violation? see SAS

public:
   inline double    getLowerBound()   const {
      return m_lb;
   }
   inline double    getUpperBound()   const {
      return m_ub;
   }
   inline double    getViolation()    const {
      return m_violation;
   }
   inline int       getEffCnt()       const {
      return m_effCnt;
   }
   inline string    getStrHash()      const {
      return m_strHash;
   }

public:
   inline void      setLowerBound(const double lb) {
      m_lb = lb;
   }
   inline void      setUpperBound(const double ub) {
      m_ub = ub;
   }
   inline void      setViolation(const double violation) {
      m_violation = violation;
   }

   bool calcViolation(const CoinPackedVector* row,
                      const double*              x);

public:
   //but these should be optional to user! after they show us how to
   //expandCutToRow... then we should be able to handle the hashing and
   //checking isSame, etc... the user can override it, if they can do it
   //faster - but we should not force them

   //now it is essentially a DecompCutOsi
   virtual void     setStringHash(CoinPackedVector* row) {
      //the user can override this if they can do it faster... also
      //should link up with isSame
      char sense;
      double rhs, range;
      UtilBoundToSense(getLowerBound(),
                       getUpperBound(), DecompInf,
                       sense, rhs, range);
      m_strHash = UtilCreateStringHash(row->getNumElements(),
                                       row->getIndices(),
                                       row->getElements(),
                                       sense, rhs);
      //need backup for user
      //throw CoinError("Method was invoked but not overridden.",
      //		    "setStringHash", "DecompCut");1
   }

   virtual void     expandCutToRow(CoinPackedVector* row) {
      throw CoinError("Method was invoked but not overridden.",
                      "expandCutToRow", "DecompCut");
   }

   virtual void     setBounds() {
      throw CoinError("Method was invoked but not overridden.",
                      "setBounds", "DecompCut");
   }

   virtual bool     isSame(const DecompCut* cut) const {
      return false;
   }

   virtual void     print(ostream* os = &cout) const;

public:
   inline void resetEffCnt() {
      m_effCnt = 0;
   }

   /** Increase the effectiveness count by 1 (or to 1 if it was negative).
       Return the new effectiveness count. */
   inline void increaseEffCnt() {
      m_effCnt = m_effCnt <= 0 ? 1 : m_effCnt + 1;
   }

   /** Decrease the effectiveness count by 1 (or to -1 if it was positive).
       Return the new effectiveness count. */
   inline void decreaseEffCnt() {
      m_effCnt = m_effCnt >= 0 ? -1 : m_effCnt - 1;
   }

public:
   DecompCut() :
      m_lb       (0.0),
      m_ub       (0.0),
      m_violation(0.0),
      m_effCnt   (0),
      m_strHash  () {
   };
   virtual ~DecompCut() {};

};

#endif
