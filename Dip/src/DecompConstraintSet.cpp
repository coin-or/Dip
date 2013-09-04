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

//===========================================================================//
#include "UtilHash.h"
#include "UtilMacrosDecomp.h"
#include "DecompConstraintSet.h"

using namespace std;

//===========================================================================//
void DecompConstraintSet::prepareModel(bool modelIsCore)
{
   //---
   //--- For each model:
   //---   1.) set row senses and/or bounds
   //---   2.) create row hash
   //---   3.) set nBaseRows
   //---   4.) flip to row ordered, if neccessary (for relaxed too)
   //---   5.) mark integers
   //---   6.) if sparse, set active columns
   //---
   if (!M) {
      return;
   }

   UtilPrintMemUsage(&cout, 2, 2);

   //TODO: needed for relax?
   if (M->isColOrdered()) {
      M->reverseOrdering();
   }

   int numRows     = getNumRows();
   int numCols     = getNumCols();
   int numColsOrig = getNumColsOrig();
   UtilPrintMemUsage(&cout, 2, 2);
   checkSenseAndBound();

   if (modelIsCore) {
      createRowHash();
   }

   nBaseRows = getNumRows();
   //TODO: make this an option
   //---
   //--- if row/col names are not given, make up default ones
   //---
   int i, j;

   if (rowNames.size() == 0) {
      for (i = 0; i < numRows; i++) {
         rowNames.push_back("r(" + UtilIntToStr(i) + ")");
      }
   }

   if (colNames.size() == 0) {
      for (j = 0; j < numCols; j++) {
         colNames.push_back("x(" + UtilIntToStr(j) + ")");
      }
   }

   prepHasRun = true;

   //---
   //--- if active columns were not set (or sparse), set to all columns
   //---   note: this is in terms of the original indices (not sparse)
   //---
   if (isSparse()) {
      //---
      //--- is this case, the user might have set this
      //---   or might not have, so we want to find the
      //---   set based on the mapping, but need to check
      //---   for duplicates, in case the user already
      //---   provided this set
      //---
      set<int> activeColumnsSet(activeColumns.begin(), activeColumns.end());
      map<int, int>::const_iterator mcit;
      activeColumns.reserve(m_sparseToOrig.size());

      for (mcit  = m_sparseToOrig.begin();
            mcit != m_sparseToOrig.end(); mcit++) {
         activeColumnsSet.insert(mcit->second);
      }

      set<int>::iterator sit;
      activeColumns.clear();

      for (sit  = activeColumnsSet.begin();
            sit != activeColumnsSet.end(); sit++) {
         activeColumns.push_back(*sit);
      }
   } else {
      int nActiveColumns = static_cast<int>(activeColumns.size());

      if (!nActiveColumns) {
         UtilIotaN(activeColumns, numColsOrig, 0);
      }
   }

   //---
   //--- if dense format, fix non-active columns
   //---
   if (!isSparse()) {
      fixNonActiveColumns();
   }

   //---
   //--- create set from vector - easier to check overlap, etc
   //---
   vector<int>::iterator vit;

   for (vit = activeColumns.begin(); vit != activeColumns.end(); vit++) {
      activeColumnsS.insert(*vit);
   }

   //---
   //--- set column markers (original number of cols)
   //---
   //UtilFillN(columnMarker, numColsOrig, (int)DecompColNonActive);
   //for(vit = activeColumns.begin(); vit != activeColumns.end(); vit++)
   // columnMarker[*vit] = DecompColActive;
   //for(vit = masterOnlyCols.begin(); vit != masterOnlyCols.end(); vit++)
   // columnMarker[*vit] = DecompColMasterOnly;

   //---
   //--- mark integers (original number of cols)
   //---    only do this for core
   if (modelIsCore) {
      UtilFillN(integerMark, numColsOrig, 'C');

      for (vit = integerVars.begin(); vit != integerVars.end(); vit++) {
         integerMark[*vit] = 'I';
      }
   }
}

//===========================================================================//
void DecompConstraintSet::createRowHash()
{
   int    r;
   string strHash;
   const int*     rmat_ind = M->getIndices();
   const double* rmat_els = M->getElements();
   const int*     rmat_beg = M->getVectorStarts();
   const int*     rmat_len = M->getVectorLengths();

   for (r = 0; r < getNumRows(); r++) {
      strHash = UtilCreateStringHash(rmat_len[r],
                                     rmat_ind + rmat_beg[r],
                                     rmat_els + rmat_beg[r],
                                     rowSense[r],
                                     rowRhs[r]);
      rowHash.push_back(strHash);
   }
}

//===========================================================================//
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

//===========================================================================//
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

//===========================================================================//
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

//===========================================================================//
void DecompConstraintSet::fixNonActiveColumns()
{
   const int numCols     = getNumCols();
   const int nActiveCols = static_cast<int>(activeColumns.size());

   if (nActiveCols == numCols) {
      return;
   }

   int* marker = new int[numCols];

   if (!marker) {
      UtilExceptionMemory("fixNonActiveColumns", "DecompConstraintSet");
   }

   UtilFillN(marker, numCols, 0);
   vector<int>::iterator vi;

   for (vi = activeColumns.begin(); vi != activeColumns.end(); vi++) {
      marker[*vi] = 1;
   }

   int i;

   for (i = 0; i < numCols; i++) {
      if (marker[i]) {
         continue;
      }

      colLB[i] = 0.0;
      colUB[i] = 0.0;
   }

   UTIL_DELARR(marker);
}

//===========================================================================//
CoinPackedMatrix* DecompConstraintSet::sparseToOrigMatrix()
{
   assert(m_isSparse);
   //---
   //--- create a dense row-majored version of the sparse matrix M
   //---
   bool                 colOrdered = M->isColOrdered();
   int                  nCols      = m_numColsOrig;
   int                  nRows      = M->getNumRows();
   CoinPackedMatrix* MRow = NULL;

   if (colOrdered) {
      //---
      //--- first create a row-ordered version
      //---
      MRow = new CoinPackedMatrix();
      CoinAssertHint(MRow, "Error: Out of Memory");
      MRow->reverseOrderedCopyOf(*M);
   } else {
      MRow = new CoinPackedMatrix(*M);
      CoinAssertHint(MRow, "Error: Out of Memory");
   }

   int                  i;
   int                  nElems     = MRow->getNumElements();
   const int*           matInd     = MRow->getIndices();
   const int*           matLen     = MRow->getVectorLengths();
   const double*        matVal     = MRow->getElements();
   const CoinBigIndex* matBeg     = MRow->getVectorStarts();
   int* matIndOrig = new int[nElems];
   CoinAssertHint(matIndOrig, "Error: Out of Memory");

   for (i = 0; i < nElems; i++) {
      matIndOrig[i] = m_sparseToOrig[matInd[i]];
   }

   CoinPackedMatrix* MOrig
      = new CoinPackedMatrix(false,
                             nCols,
                             nRows,
                             nElems,
                             matVal,
                             matIndOrig,
                             matBeg,
                             matLen,
                             0.0, 0.0);
   CoinAssertHint(MOrig, "Error: Out of Memory");
   UTIL_DELPTR(MRow);
   UTIL_DELARR(matIndOrig);
   return MOrig;
}
