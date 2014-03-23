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


#include "UtilHash.h"
#include "DecompConstraintSet.h"

// --------------------------------------------------------------------- //
void DecompConstraintSet::createRowHash()
{
   int    r;
   string strHash;
   const int*     rmat_ind = M->getIndices();
   const double* rmat_els = M->getElements();
   const int*     rmat_beg = M->getVectorStarts();
   const int*     rmat_len = M->getVectorLengths();
   printf("\nNUM COLS: %d", getNumCols());
   printf("\nNUM ROWS: %d", getNumRows());

   for (r = 0; r < getNumRows(); r++) {
      //printf("\n");
      strHash = UtilCreateStringHash(rmat_len[r],
                                     rmat_ind + rmat_beg[r],
                                     rmat_els + rmat_beg[r],
                                     rowSense[r],
                                     rowRhs[r]);
      rowHash.push_back(strHash);
   }
}

// --------------------------------------------------------------------- //
void DecompConstraintSet::checkSenseAndBound()
{
   assert(rowLB.size() + rowRhs.size() > 0);
   assert(rowUB.size() + rowRhs.size() > 0);

   if (rowLB.size() > 0 && rowRhs.size() == 0) {
      boundsToSenses();
   } else if (rowLB.size() == 0 && rowRhs.size() > 0) {
      sensesToBounds();
   }

   assert(rowLB.size() == rowUB.size());
   assert(rowLB.size() == rowRhs.size());
   assert(rowLB.size() == rowSense.size());
}

// --------------------------------------------------------------------- //
void DecompConstraintSet::sensesToBounds()
{
   double rlb, rub;
   int    n_rows = static_cast<int>(rowSense.size());
   rowLB.reserve(n_rows);
   rowUB.reserve(n_rows);

   for (int r = 0; r < n_rows; r++) {
      UtilSenseToBound(rowSense[r], rowRhs[r], 0.0,//TODO
                       DecompInf, rlb, rub);
      rowLB.push_back(rlb);
      rowUB.push_back(rub);
   }
}

// --------------------------------------------------------------------- //
void DecompConstraintSet::boundsToSenses()
{
   char   sense;
   double rhs, range;//not used
   int    n_rows = static_cast<int>(rowLB.size());
   rowRhs.reserve(n_rows);
   rowSense.reserve(n_rows);

   for (int r = 0; r < n_rows; r++) {
      UtilBoundToSense(rowLB[r], rowUB[r], DecompInf,
                       sense, rhs, range);
      rowRhs.push_back(rhs);
      rowSense.push_back(sense);
   }
}

