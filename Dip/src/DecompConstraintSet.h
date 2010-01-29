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
// Conceptual Design: Matthew Galati, SAS Institute Inc.                     //
//                    Ted Ralphs, Lehigh University                          //
//                                                                           //
// Copyright (C) 2002-2009, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

#ifndef DECOMP_CONSTRAINTSET_INCLUDED
#define DECOMP_CONSTRAINTSET_INCLUDED

// --------------------------------------------------------------------- //
#include "Decomp.h"
#include "UtilMacros.h"

enum ColMarkerType {
   DecompColNonActive = 0,
   DecompColActive    = 1,
   DecompColMasterOnly= 2
};

// --------------------------------------------------------------------- //
class DecompConstraintSet{
public:
   CoinPackedMatrix   * M;
   int                  nBaseRowsOrig;
   int                  nBaseRows;
   vector<string>       rowHash;
   vector<char>         rowSense;
   vector<double>       rowRhs;
   vector<double>       rowLB;
   vector<double>       rowUB;
   vector<double>       colLB;
   vector<double>       colUB;
   vector<int>          integerVars;
   vector<string>       colNames;
   vector<string>       rowNames;
   vector<int>          activeColumns; //if block, define the active columns
   set<int>             activeColumnsS;//if block, define the active columns
   vector<int>          masterOnlyCols;
   vector<int>          columnMarker; //active, non-active, master-only
   bool                 prepHasRun;

public:
   inline const CoinPackedMatrix * getMatrix() const { return M; };
   inline const int getNumRows() const { return M->getNumRows(); }
   inline const int getNumCols() const { return M->getNumCols(); }
   inline const int getNumInts() const {
      return static_cast<int>(integerVars.size());}
   inline const vector<int>    & getActiveColumns() const {
      return activeColumns;}
   inline const vector<string> & getRowNames() const {return rowNames;}
   inline const vector<string> & getColNames() const {return colNames;}   
   inline vector<string> & getRowNamesMutable() {return rowNames;}
   inline vector<string> & getColNamesMutable() {return colNames;}
   inline const int    * getIntegerVars() { return &integerVars[0];}
   inline const double * getColLB() const { return &colLB[0]; };
   inline const double * getColUB() const { return &colUB[0]; };
   inline const double * getRowLB() const { return &rowLB[0]; };
   inline const double * getRowUB() const { return &rowUB[0]; };
   inline const bool     hasPrepRun() const { return prepHasRun; };

public:
   void prepareModel();
   void createRowHash();
   void checkSenseAndBound();
   void sensesToBounds();
   void boundsToSenses();
   void fixNonActiveColumns();
   inline void appendRow(CoinPackedVector & row,
			 double             loBound,
			 double             upBound){
      M->appendRow(row);
      rowLB.push_back(loBound);
      rowUB.push_back(upBound);
   }
   inline void appendRow(CoinPackedVector & row,
			 double             loBound,
			 double             upBound,
			 string           & rowName){
      appendRow(row, loBound, upBound);
      rowNames.push_back(rowName);         
   }
   
   inline void reserve(const int nCols,
		       const int nRows){
      M->reserve(nRows, nCols);
      rowLB.reserve(nRows);
      rowUB.reserve(nRows);
      colLB.reserve(nCols);
      colUB.reserve(nCols);
   }
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
      integerVars(),
      colNames(),
      rowNames(),
      prepHasRun(false)
   {};
  
   ~DecompConstraintSet() {
      UTIL_DELPTR(M);    
   };
};

#endif
