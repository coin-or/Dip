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


#ifndef DECOMP_WAITING_ROW_INCLUDE
#define DECOMP_WAITING_ROW_INCLUDE

//class DecompCut;
#include "DecompCut.h"

//TODO? use to have DecompRow = CoinPackedVector

// ---------------------------------------------------------------------- //
class DecompWaitingRow {
private:
   //THINK
   //DecompWaitingRow & operator=(const DecompWaitingRow &);

private:
   DecompCut*         m_cut;        //the cut
   CoinPackedVector* m_row;         //the row (in terms of x)
   CoinPackedVector* m_rowReform;   //the row (in terms of reformulation)

public:
   inline DecompCut*           getCutPtr() const       {
      return m_cut;
   }
   inline CoinPackedVector* getRowPtr() const       {
      return m_row;
   }
   inline CoinPackedVector* getRowReformPtr() const {
      return m_rowReform;
   }
   inline const double getViolation() const  {
      return m_cut->getViolation();
   }
   inline const double getLowerBound() const {
      return m_cut->getLowerBound();
   }
   inline const double getUpperBound() const {
      return m_cut->getUpperBound();
   }

   inline void   deleteCut()                      {
      //cout << "\ndelete m_cut: " << m_cut << "\n";
      UTIL_DELPTR(m_cut);
   }
   inline void   deleteRow()                      {
      UTIL_DELPTR(m_row);
   }
   inline void   deleteRowReform()                {
      UTIL_DELPTR(m_rowReform);
   }
   inline void   clearCut()                       {
      m_cut = 0;
   }

   inline void   setRow(CoinPackedVector* row) {
      m_row = row;
   }
   inline void   setRowReform(CoinPackedVector* rowReform) {
      m_rowReform = rowReform;
   }

   bool setViolation(const double* x);

public:
   DecompWaitingRow(const DecompWaitingRow& x) {
      m_cut       = x.m_cut;
      m_row       = x.m_row;
      m_rowReform = x.m_rowReform;
   }

   DecompWaitingRow(DecompCut*           cut,
                    CoinPackedVector* row,
                    CoinPackedVector* rowReform = 0) :
      m_cut(cut),
      m_row(row),
      m_rowReform(rowReform) {}

   ~DecompWaitingRow() {}
};

#endif

