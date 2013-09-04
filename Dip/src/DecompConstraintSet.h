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
// Copyright (C) 2002-2013, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

#ifndef DECOMP_CONSTRAINTSET_INCLUDED
#define DECOMP_CONSTRAINTSET_INCLUDED

// --------------------------------------------------------------------- //
#include "Decomp.h"
#include "UtilMacros.h"

// --------------------------------------------------------------------- //
enum ColMarkerType {
   DecompColNonActive = 0,
   DecompColActive    = 1,
   DecompColMasterOnly = 2
};

// --------------------------------------------------------------------- //
class DecompConstraintSet {
public:
   CoinPackedMatrix*    M;
   int                  nBaseRowsOrig;
   int                  nBaseRows;
   std::vector<std::string>       rowHash;
   std::vector<char>         rowSense;
   std::vector<double>       rowRhs;
   std::vector<double>       rowLB;
   std::vector<double>       rowUB;
   std::vector<double>       colLB;
   std::vector<double>       colUB;
   std::vector<int>          integerVars;
   std::vector<char>         integerMark; //'C' = continuous, 'I' = integral
   std::vector<std::string>       colNames;
   std::vector<std::string>       rowNames;
   std::vector<int>          activeColumns; //if block, define the active columns
   std::set<int>             activeColumnsS;//if block, define the active columns
   std::vector<int>          masterOnlyCols;
   bool                 prepHasRun;

   //for storage of several rows of row-majored sparse matrix
   //  to be used with appendRows
   std::vector<CoinBigIndex> m_rowBeg;
   std::vector<int         > m_rowInd;
   std::vector<double      > m_rowVal;

   //for special case of master-only vars
   bool                 m_masterOnly;
   int                  m_masterOnlyIndex;
   double               m_masterOnlyLB;
   double               m_masterOnlyUB;
   bool                 m_masterOnlyIsInt;

   //for special case of sparse representation
   bool          m_isSparse;
   int           m_numColsOrig;
   std::map<int, int> m_origToSparse;
   std::map<int, int> m_sparseToOrig;

public:
   inline void setSparse(const int numColsOrig) {
      m_numColsOrig = numColsOrig;
      m_isSparse    = true;
   }
   inline const bool isSparse() const {
      return m_origToSparse.size() ? true : false;
   };
   inline const CoinPackedMatrix* getMatrix() const {
      return M;
   };
   inline const int getNumRows() const {
      return M ? M->getNumRows() : static_cast<int>(rowLB.size());
   }
   inline const int getNumCols() const {
      return M ? M->getNumCols() : static_cast<int>(colLB.size());
   }
   inline const int getNumColsOrig() const {
      return isSparse() ? m_numColsOrig : getNumCols();
   };
   inline const int getNumInts() const {
      return static_cast<int>(integerVars.size());
   }
   inline const std::vector<int>&     getActiveColumns() const {
      return activeColumns;
   }
   inline const std::vector<std::string>& getRowNames() const {
      return rowNames;
   }
   inline const std::vector<std::string>& getColNames() const {
      return colNames;
   }
   inline std::vector<std::string>& getRowNamesMutable() {
      return rowNames;
   }
   inline std::vector<std::string>& getColNamesMutable() {
      return colNames;
   }
   inline const char*    getIntegerMark() {
      return &integerMark[0];
   }
   inline const int*     getIntegerVars() {
      return &integerVars[0];
   }
   inline const double* getColLB() const {
      return &colLB[0];
   };
   inline const double* getColUB() const {
      return &colUB[0];
   };
   inline const double* getRowLB() const {
      return &rowLB[0];
   };
   inline const double* getRowUB() const {
      return &rowUB[0];
   };
   inline const bool     hasPrepRun() const {
      return prepHasRun;
   };
   inline const bool     isMasterOnly() const {
      return m_masterOnly;
   };
   inline const std::map<int, int>& getMapOrigToSparse() const {
      return m_origToSparse;
   };
   inline const std::map<int, int>& getMapSparseToOrig() const {
      return m_sparseToOrig;
   };
   inline const std::vector<int>& getMasterOnlyCols() const {
      return masterOnlyCols;
   }


public:
   void prepareModel(bool modelIsCore = false);
   void createRowHash();
   void checkSenseAndBound();
   void sensesToBounds();
   void boundsToSenses();
   void fixNonActiveColumns();
   CoinPackedMatrix* sparseToOrigMatrix();

   inline void appendRow(CoinPackedVector& row,
                         double             loBound,
                         double             upBound) {
      M->appendRow(row);
      rowLB.push_back(loBound);
      rowUB.push_back(upBound);
   }
   //inline void appendRow(CoinPackedVector & row,
   //		 double             loBound,
   //		 double             upBound,
   //		 std::string      & rowName){
   //   appendRow(row, loBound, upBound);
   // rowNames.push_back(rowName);
   //}
   inline void appendRow(CoinPackedVector& row,
                         double             loBound,
                         double             upBound,
                         std::string        rowName) {
      appendRow(row, loBound, upBound);
      rowNames.push_back(rowName);
   }

   inline void pushCol(const double loBound,
                       const double upBound,
                       const bool   isInteger   = false,
                       const int    origIndex   = -1) {
      int index = static_cast<int>(colLB.size());
      colLB.push_back(loBound);
      colUB.push_back(upBound);

      if (isInteger) {
         integerVars.push_back(index);
      }

      assert(!(origIndex == -1 && m_isSparse));

      if (origIndex >= 0) {
         m_origToSparse.insert(std::make_pair(origIndex, index));
         m_sparseToOrig.insert(std::make_pair(index, origIndex));
      }
   }

   inline void reserve(const int nCols,
                       const int nRows) {
      M->reserve(nRows, nCols);
      rowLB.reserve(nRows);
      rowUB.reserve(nRows);
      colLB.reserve(nCols);
      colUB.reserve(nCols);
   }
public:
   DecompConstraintSet() :
      M                (0),
      nBaseRowsOrig    (0),
      nBaseRows        (0),
      prepHasRun       (false),
      m_masterOnly     (false),
      m_masterOnlyIndex(0),
      m_masterOnlyLB   (0.0),
      m_masterOnlyUB   (0.0),
      m_masterOnlyIsInt(false),
      m_isSparse       (false),
      m_numColsOrig    (0) {
   };

   ~DecompConstraintSet() {
      UTIL_DELPTR(M);
   };
};

#endif
