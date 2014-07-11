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


#ifndef DECOMP_CUT_POOL_INCLUDE
#define DECOMP_CUT_POOL_INCLUDE

#include "DecompWaitingRow.h"

#include <functional>
using namespace std;

class DecompConstraintSet;

// --------------------------------------------------------------------- //
//TODO: switch to distance
class is_greater_thanD { //member of class instead??
public:
   //TODO: design, waitingcol, rc is member of var, not waiting col,
   //but for waitingrow, distance is member of wr, not of cut - why?
   bool operator()( const DecompWaitingRow& x,
                    const DecompWaitingRow& y) {
      return x.getViolation() > y.getViolation();
   }
};

// --------------------------------------------------------------------- //
class DecompCutPool : public std::vector<DecompWaitingRow> {
private:
   DecompCutPool(const DecompCutPool&);
   DecompCutPool& operator=(const DecompCutPool&);

private:
   static const char* classTag;
   bool m_rowsAreValid;

public:
   const inline bool rowsAreValid() const {
      return m_rowsAreValid;
   }
   inline void setRowsAreValid(bool rowsAreValid) {
      m_rowsAreValid = rowsAreValid;
   }

   void print(ostream* os = &cout) const;  //THINK: virtual??
   void reExpand(const DecompVarList& vars,
                 const int             n_corecols);

   CoinPackedVector* createRowReform(const int                  n_corecols,
                                     const CoinPackedVector* row,
                                     const DecompVarList&       vars);

   //THINK
   //bool isDuplicate(const DecompWaitingRow & wcol);

   bool calcViolations(const double*             x,
                       DecompCutPool::iterator   first,
                       DecompCutPool::iterator   last);
   bool calcViolations(const double*             x) {
      return calcViolations(x, begin(), end());
   }

public:
   DecompCutPool() :
      m_rowsAreValid(true) {}

   ~DecompCutPool() {
      //---
      //--- delete any memory that is left in the waiting rows
      //---
      vector<DecompWaitingRow>::iterator vi;

      for (vi = begin(); vi != end(); vi++) {
         (*vi).deleteCut();
         (*vi).deleteRow();
         (*vi).deleteRowReform();
      }
   }

};

#endif
