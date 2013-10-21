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

#ifndef DECOMP_CONSTRAINTSET_INCLUDED
#define DECOMP_CONSTRAINTSET_INCLUDED

#include "UtilMacros.h"
#include "DecompPortable.h"

// --------------------------------------------------------------------- //
class DecompConstraintSet {
private:
   //need these if vector?
   //DecompConstraintSet(const DecompConstraintSet &);
   //DecompConstraintSet & operator=(const DecompConstraintSet &);

   //THINK - better overall to store as sense!
public:
   //do we want this dependence on CoinPackedMatrix?
   //do we want all this depenence on STL vectors?

   //how does this all work for removing cuts?? if using vector!?

   //must flip to row ordered?
   CoinPackedMatrix* M;
   int                  nBaseRowsOrig;
   int                  nBaseRows;
   vector<string>       rowHash;
   vector<char>         rowSense;
   vector<double>       rowRhs; //ugh... have to carry around
   //ranges?

   vector<double>       rowLB; //vector or double *?
   vector<double>       rowUB;
   vector<double>       colLB;
   vector<double>       colUB;
   vector<int>          integerVars;
   //TODO: why not vector if rest are... nice if consistent
   //double             * objCoeff; //only used for modelCore?

   //TODO: colNames, rowNames


public:
   inline const int getNumRows() const {
      return M->getNumRows();
   }
   inline const int getNumCols() const {
      return M->getNumCols();
   }

public:
   void createRowHash();
   void checkSenseAndBound();
   void sensesToBounds();
   void boundsToSenses();

public:
   DecompConstraintSet() :
      M(0),
      nBaseRowsOrig(0),
      nBaseRows(0),
      rowSense(),
      rowRhs(),
      rowLB(),
      rowUB(),
      colLB(),
      colUB(),
      integerVars()
      //objCoeff(0)
   {};

   ~DecompConstraintSet() {
      UTIL_DELPTR(M);
   };
};

#endif
