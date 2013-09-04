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


#ifndef DECOMP_CUT_POOL_INCLUDE
#define DECOMP_CUT_POOL_INCLUDE

#include "DecompWaitingRow.h"

#include <functional>

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

   void print(std::ostream* os = &std::cout) const;  //THINK: virtual??
   void reExpand(const DecompVarList& vars,
                 const int             n_coreCols,
                 const int             n_artCols);

   CoinPackedVector* createRowReform(const int                n_coreCols,
                                     //const int                n_artCols,
                                     const CoinPackedVector* row,
                                     const DecompVarList&     vars);

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
      std::vector<DecompWaitingRow>::iterator vi;

      for (vi = begin(); vi != end(); vi++) {
         (*vi).deleteCut();
         (*vi).deleteRow();
         (*vi).deleteRowReform();
      }
   }

};

#endif
