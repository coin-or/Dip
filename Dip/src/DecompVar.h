//===========================================================================//
// This file is part of the DIP Solver Framework.                            //
//                                                                           //
// DIP is distributed under the Eclipse Public License as part of the        //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)       //
//                                                                           //
// Conceptual Design: Matthew Galati, SAS Institute Inc.                     //
//                    Ted Ralphs, Lehigh University                          //
//                                                                           //
// Copyright (C) 2002-2013, Lehigh University, Matthew Galati, Ted Ralphs    //
// All Rights Reserved.                                                      //
//===========================================================================//

//this is a lambda_s variable
//assumption throughout is lambda is {0,1}

#ifndef DECOMP_VAR_INCLUDED
#define DECOMP_VAR_INCLUDED

#include "Decomp.h"
#include "UtilHash.h"
#include "UtilMacrosDecomp.h"

class DecompApp;
class DecompAlgoModel;

// --------------------------------------------------------------------- //
class DecompVar {
public:
   //THINK: or user overriden way to store it (like a tree)
   //and a function which says expandVarToCol - just as in cuts
   CoinPackedVector m_s;//this is the var in terms of x-space

private:
   //TODO: lb, ub, "type"?
   double           m_origCost;
   double           m_redCost; //(c - uA'')s - alpha
   int              m_effCnt;  //effectiveness counter
   std::string           m_strHash;
   int              m_blockId;
   int              m_colMasterIndex;
   double           m_norm;

public:
   inline double getOriginalCost()   const {
      return m_origCost;
   }
   inline double getReducedCost()    const {
      return m_redCost;
   }
   inline int    getEffectiveness()  const {
      return m_effCnt;
   }
   inline double getLowerBound()     const {
      return 0.0;         //TODO
   }
   inline double getUpperBound()     const {
      return DecompInf;   //TODO
   }
   inline std::string getStrHash()        const {
      return m_strHash;
   }
   inline int    getBlockId()        const {
      return m_blockId;
   }
   inline int    getColMasterIndex() const {
      return m_colMasterIndex;
   }
   inline double getNorm()           const {
      return m_norm;
   }

   inline void   setColMasterIndex(const int colIndex)  {
      m_colMasterIndex = colIndex;
   }
   inline void   setBlockId (const int blockId)  {
      m_blockId = blockId;
   }
   inline void   setReducedCost (const double redCost)  {
      m_redCost = redCost;
   }
   inline void   setOriginalCost(const double origCost) {
      m_origCost = origCost;
   }

   inline void resetEffectiveness() {
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

   inline double calcNorm() {
      return m_norm = m_s.twoNorm();
   }

   inline void sortVar() {
      m_s.sortIncrIndex();
   }

   bool   isEquivalent(const DecompVar& dvar) {
      return m_s.isEquivalent(dvar.m_s);
   }

   bool isDuplicate(const DecompVarList& vars) {
      DecompVarList::const_iterator vi;

      for (vi = vars.begin(); vi != vars.end(); vi++) {
         if ((*vi)->getStrHash() == this->getStrHash()) {
            return true;
         }
      }

      return false;
   }

   bool doesSatisfyBounds(int                     denseLen,
                          double*                 denseArr,
                          const DecompAlgoModel& model,
                          const double*           lbs,
                          const double*           ubs);

   void fillDenseArr(int      len,
                     double* arr);

public:
   virtual void  print(std::ostream*    os  = &std::cout,
                       DecompApp* app = 0) const;
   virtual void  print(std::ostream*               os,
                       const std::vector<std::string>& colNames,
                       const double*          value    = NULL) const;

public:
   /** @name Copy Constructors */
   DecompVar(const DecompVar& source) :
      m_s       (source.m_s),
      m_origCost(source.m_origCost),
      m_effCnt  (source.m_effCnt),
      m_strHash (source.m_strHash),
      m_blockId (source.m_blockId),
      m_colMasterIndex (source.m_colMasterIndex),
      m_norm    (source.m_norm) {
   }

   DecompVar& operator=(const DecompVar& rhs) {
      if (this != &rhs) {
         m_s        = rhs.m_s;
         m_origCost = rhs.m_origCost;
         m_redCost  = rhs.m_redCost;
         m_effCnt   = rhs.m_effCnt;
         m_strHash  = rhs.m_strHash;
         m_blockId  = rhs.m_blockId;
         m_colMasterIndex = rhs.m_colMasterIndex;
      }

      return *this;
   }

   DecompVar():
      m_s       (),
      m_origCost(0.0),
      m_redCost (0.0),
      m_effCnt  (0),
      m_strHash (),
      m_blockId (0),
      m_colMasterIndex(-1),
      m_norm    (0.0) {
   }

   DecompVar(const std::vector<int>&     ind,
             const double           els,
             const double           redCost,
             const double           origCost) :
      m_s       (),
      m_origCost(origCost),
      m_redCost (redCost),
      m_effCnt  (0),
      m_strHash (),
      m_blockId (0),
      m_colMasterIndex(-1),
      m_norm    (0.0) {
      if (ind.size() > 0) {
         m_s.setConstant(static_cast<int>(ind.size()),
                         &ind[0], els, DECOMP_TEST_DUPINDEX);
         m_strHash = UtilCreateStringHash(static_cast<int>(ind.size()),
                                          &ind[0], els);
         m_norm    = calcNorm();
         sortVar();
      }
   }

   DecompVar(const std::vector<int>&     ind,
             const std::vector<double>& els,
             const double           redCost,
             const double           origCost) :
      m_s       (),
      m_origCost(origCost),
      m_redCost (redCost),
      m_effCnt  (0),
      m_strHash (),
      m_blockId (0),
      m_colMasterIndex(-1),
      m_norm    (0.0) {
      if (ind.size() > 0) {
         m_s.setVector(static_cast<int>(ind.size()),
                       &ind[0], &els[0], DECOMP_TEST_DUPINDEX);
         m_strHash = UtilCreateStringHash(static_cast<int>(ind.size()),
                                          &ind[0], &els[0]);
         m_norm    = calcNorm();
         sortVar();
      }
   }

   DecompVar(const int              len,
             const int*             ind,
             const double*          els,
             const double           origCost) :
      m_s       (),
      m_origCost(origCost),
      m_redCost (0.0),
      m_effCnt  (0),
      m_strHash (),
      m_blockId (0),
      m_colMasterIndex(-1),
      m_norm    (0.0) {
      if (len > 0) {
         m_s.setVector(len, ind, els, DECOMP_TEST_DUPINDEX);
         m_strHash = UtilCreateStringHash(len, ind, els);
         m_norm    = calcNorm();
         sortVar();
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
      m_strHash (),
      m_blockId (0),
      m_colMasterIndex(-1),
      m_norm    (0.0) {
      if (len > 0) {
         m_s.setVector(len, ind, els, DECOMP_TEST_DUPINDEX);
         m_strHash = UtilCreateStringHash(len, ind, els);
         m_norm    = calcNorm();
         sortVar();
      }
   }

   DecompVar(const int              len,
             const int*             ind,
             const double           els,
             const double           origCost) :
      m_s       (),
      m_origCost(origCost),
      m_redCost (0.0),
      m_effCnt  (0),
      m_strHash (),
      m_blockId (0),
      m_colMasterIndex(-1),
      m_norm    (0.0) {
      if (len > 0) {
         m_s.setConstant(len, ind, els, DECOMP_TEST_DUPINDEX);
         m_strHash = UtilCreateStringHash(len, ind, els);
         m_norm    = calcNorm();
         sortVar();
      }
   }

   DecompVar(const int              len,
             const int*             ind,
             const double           els,
             const double           redCost,
             const double           origCost) :
      m_s       (),
      m_origCost(origCost),
      m_redCost (redCost),
      m_effCnt  (0),
      m_strHash (),
      m_blockId (0),
      m_colMasterIndex(-1),
      m_norm    (0.0) {
      if (len > 0) {
         m_s.setConstant(len, ind, els, DECOMP_TEST_DUPINDEX);
         m_strHash = UtilCreateStringHash(len, ind, els);
         m_norm    = calcNorm();
         sortVar();
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
      m_strHash (),
      m_blockId (0),
      m_colMasterIndex(-1),
      m_norm    (0.0) {
      UtilPackedVectorFromDense(denseLen, denseArray, DecompEpsilon, m_s);

      if (m_s.getNumElements() > 0) {
         m_strHash = UtilCreateStringHash(denseLen, denseArray);
         m_norm    = calcNorm();
         sortVar();
      }
   }

   virtual ~DecompVar() {};
};

#endif
