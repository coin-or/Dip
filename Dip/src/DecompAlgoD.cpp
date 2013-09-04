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


#include "DecompApp.h"
#include "DecompVar.h"
#include "DecompAlgoD.h"
#include "DecompCutOsi.h"
#include "DecompSolverResult.h"

using namespace std;

//TODO: generateInitVars should be based on cost = -xhat

// ------------------------------------------------------------------------- //
void DecompAlgoD::phaseUpdate(DecompPhase&   phase,
                              DecompStatus& status)
{
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "phaseUpdate()", m_param.LogDebugLevel, 2);
   DecompAlgo::phaseUpdate(phase, status);

   if (phase == PHASE_DONE && status == STAT_FEASIBLE) {
      //---
      //--- then a decomposition was found, return it
      //---
   }

   //---
   //--- TODO:
   //--- are we stuck?
   //---   11/3/09: I13 PhaseIObj not moving
   //--- TODO: use tailoff function
   //---
   int    changeLen      = m_param.TailoffLength;
   double changePerLimit = m_param.TailoffPercent;

   if (static_cast<int>(m_phaseIObj.size()) > changeLen) {
      vector< double >::reverse_iterator it = m_phaseIObj.rbegin();
      int    len       = 0;
      double prevBound = (*it);
      double diff      = DecompInf;
      double sumDiff   = 0.0;
      double aveDiff   = 0.0;
      double perDiff   = 0.0;

      for ( ; it != m_phaseIObj.rend(); it++) {
         diff = fabs(prevBound - (*it));
         UTIL_DEBUG(m_param.LogDebugLevel, 3,
                    (*m_osLog)
                    << setw(10) << "prevBound="
                    << setw(10) << UtilDblToStr(prevBound, 2)
                    << setw(10) << ", thisBound="
                    << setw(10) << UtilDblToStr((*it)) << endl;
                   );
         sumDiff   += diff;
         prevBound  = (*it);
         len++;

         if (len >= changeLen) {
            break;
         }
      }

      aveDiff = sumDiff / len;

      if (UtilIsZero(prevBound)) {
         perDiff = aveDiff;
      } else {
         perDiff = 100 * aveDiff / fabs(prevBound);
      }

      UTIL_MSG(m_param.LogDebugLevel, 2,
               (*m_osLog)
               << setw(10) << "Percentage difference in obj bound="
               << setw(10) << UtilDblToStr(perDiff, 2) << endl;
              );

      //---
      //--- if the average percentage difference is less than some threshold
      //---    than we are tailing off
      //---
      if (perDiff <= changePerLimit) {
         UTIL_DEBUG(m_param.LogDebugLevel, 3,
                    (*m_osLog) << "DC is tailing off - STOP PROCESS" << endl;
                   );
         phase          = PHASE_DONE;
         m_stopCriteria = DecompStopTailOff;
      }
   }

   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "phaseUpdate()", m_param.LogDebugLevel, 2);
}

//===========================================================================//
void DecompAlgoD::phaseDone()
{
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "phaseDone()", m_param.LogDebugLevel, 1);

   if (m_stopCriteria != DecompStopInfeasible) {
      if (m_param.LogDebugLevel >= 3) {
         printVars(m_osLog);   //use this to warm start DW
      }

      return;
   }

   //---
   //--- decomposition could not be found, this means the
   //---  point we are decomposing is not inside P' and we can
   //---  generate a 'farkas cut'
   //---
   //--- since we use a phase I, our 'proof of infeasibility'
   //---  does not come from the 'dual ray' but rather the objective
   //---  of the oracle
   //---
   //--- by getting here, we have shown that (c=0,A=I)
   //---     (c-uA'')s  - alpha >= 0 for all s in P'
   //--- and
   //---     (c-uA'')s* - alpha <  0
   //---
   //--- in case of many blocks, take the most violated block
   //---
   int            i, b;
   const double* dualSol = m_masterSI->getRowPrice();
   double lhs = 0.0;

   for (i = 0; i < m_numOrigCols; i++) {
      lhs -= dualSol[i] * m_xhatD[i];

      if (m_param.LogDebugLevel >= 3) {
         printf("i:%4d u:%5g x:%5g lhs:%5g\n",
                i, dualSol[i], m_xhatD[i], lhs);
      }
   }

   //---
   //--- pick the alpha that maximizes the violation
   //---
   double alpha = -DecompInf;

   for (b = 0; b < m_numConvexCon; b++) {
      if (dualSol[m_numOrigCols + b] > alpha) {
         alpha = dualSol[m_numOrigCols + b];
      }
   }

   lhs -= alpha;

   if (m_param.LogDebugLevel >= 3) {
      printf("alpha:%5g lhs:%5g\n", alpha, lhs);
   }

   if (lhs < 0) {
      printf(" VIOLATED FARKAS CUT lhs = %g\n", lhs);
      CoinPackedVector cut;
      OsiRowCut        rowCut;

      //---
      //--- Farkas Cut: u*x <= -alpha
      //---
      for (i = 0; i < m_numOrigCols; i++) {
         cut.insert(i, dualSol[i]);
      }

      rowCut.setRow(cut);
      rowCut.setLb(-DecompInf);
      rowCut.setUb(-alpha);
      DecompCutOsi* decompCut = new DecompCutOsi(rowCut);
      decompCut->setStringHash();//constructor should do!
      //decompCut->print(m_osLog);
      (*m_newCuts).push_back(decompCut);
   }

   //---
   //--- comparing to Concorde's LOCALCUT implementation, our
   //---  version correponds to their version (3) in separate.c
   //--- because we have a cost of 1.0 on artificals, the duals (u)
   //---  should be bounded between -1 and 1
   //---
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "phaseDone()", m_param.LogDebugLevel, 1);
}


//===========================================================================//
void DecompAlgoD::masterMatrixAddArtCols(CoinPackedMatrix* masterM,
      double*            colLB,
      double*            colUB,
      double*            objCoeff,
      vector<string>&    colNames,
      int                startRow,
      int                endRow,
      char               origOrBranch)
{
   //---
   //--- min sp + sm
   //---
   //--- ax  = b --> ax + sp - sm  = b, sp >= 0, sm >= 0
   //--- ax <= b --> ax      - sm <= b,          sm >= 0
   //--- ax >= b --> ax + sp      >= b, sp >= 0
   //---
   int              r, colIndex;
   DecompConstraintSet*           modelCore   = m_modelCore.getModel();
   vector<string>& rowNames      = modelCore->colNames;
   bool             hasNames      = rowNames.empty() ? false : true;
   string           colName;
   string           strIndex;
   string           colNameL  = origOrBranch == 'O' ? "sOL(c_" : "sBL(c_";
   string           colNameG  = origOrBranch == 'O' ? "sOG(c_" : "sBG(c_";
   DecompColType    colTypeL  = origOrBranch == 'O' ?
                                DecompCol_ArtForRowL : DecompCol_ArtForBranchL;
   DecompColType    colTypeG  = origOrBranch == 'O' ?
                                DecompCol_ArtForRowG : DecompCol_ArtForBranchG;
   colIndex = masterM->getNumCols();
   vector<CoinBigIndex> colBeg;
   vector<int         > colInd;
   vector<double      > colVal;
   colBeg.push_back(0);

   for (r = startRow; r < endRow; r++) {
      if (hasNames) {
         strIndex = UtilIntToStr(colIndex);
      }

      masterMatrixAddArtCol(colBeg, colInd, colVal,
                            'L', r, colIndex, colTypeL,
                            colLB[colIndex], colUB[colIndex],
                            objCoeff[colIndex]);

      if (hasNames) {
         colName = colNameL + strIndex + "_" + rowNames[r] + ")";
         colNames.push_back(colName);
      }

      colIndex++;
      masterMatrixAddArtCol(colBeg, colInd, colVal,
                            'G', r, colIndex, colTypeG,
                            colLB[colIndex], colUB[colIndex],
                            objCoeff[colIndex]);

      if (hasNames) {
         colName = colNameG + strIndex + "_" + rowNames[r] + ")";
         colNames.push_back(colName);
      }

      colIndex++;
   }

   masterM->appendCols(static_cast<int>(colBeg.size()) - 1,
                       &colBeg[0],
                       &colInd[0],
                       &colVal[0]);
}





//===========================================================================//
void DecompAlgoD::createMasterProblem(DecompVarList& initVars)
{
   //use DecompAlgoPC2?
   //---
   //--- Initialize the solver interface for the master problem.
   //---
   //--- For the master constraint system:
   //---  modelCore  contains [A'', b''], in terms of x.
   //---  m_modelRelax contains [A', b'],   and might contiain multiple blocks.
   //---
   //--- For each block we must add a convexity constraint. Let K be the set
   //--- of blocks.
   //---
   //--- Notation:
   //---  n       = orig number of vars
   //---  m''     = orig number of rows in A'', b''
   //---  |K|     = the number of blocks that defines [A' b']
   //---  s       = a solution (e.p.) to the relaxed problem, size (1xn)
   //---  c       = original cost, size (1xn)
   //---  F'[k]   = the current set of relaxed e.p.'s for block k in K
   //--- a''[i,j] = entry at row i, column j for A'' matrix
   //---  C       = original set of columns (n = |C|)
   //---  R''     = original set of rows in A'' (m''=|R''|)
   //---
   //--- The Dantzig-Wolfe LP:
   //---
   //--- min  sum{k in K, s in F'[k]}
   //---        sum{j in C}(c[j]   * s[j]) * lambda[k][s]
   //--- s.t. sum{k in K, s in F'[k]}
   //---        sum{j in C}(a''[i,j] * s[j])* lambda[k][s] ~ b''[i],  i in R''
   //---      sum{s in F'[k]} lambda[k][s] = 1, k in K
   //---      lambda[k][s]                >= 0, k in K, s in F'[k]
   //---
   //---
   //--- The Dantzig-Wolfe (DECOMP) LP:
   //---
   //--- min  sum{k in K, s in F'[k]}
   //---        sum{j in C}(0   * s[j]) * lambda[k][s]
   //--- s.t. sum{k in K, s in F'[k]}
   //---        sum{j in C}(1 * s[j])* lambda[k][s] = x*[j],  j in C
   //---      sum{s in F'[k]} lambda[k][s] = 1, k in K
   //---      lambda[k][s]                >= 0, k in K, s in F'[k]
   //---
   //---
   //---
   //--- Change for Phase I model.
   //---   Add a slack and/or surplus variable to each master constraint
   //---   including the bounds for branching?? THINK...
   //---
   //--- THINK:
   //--- Do we bother removing these vars once feasible? What about the
   //--- fact that adding cuts could once again cause infeasible....
   //---
   //--- What do we do after a branching? jump back to Phase I?
   //---
   //---
   //--- Phase I:
   //--- min  sum{i in R''} (splus[i] + sminus[i])
   //---
   //--- Phase II:
   //--- min  sum{k in K, s in F'[k]}
   //---        sum{j in C}(c[j]   * s[j]) * lambda[k][s]
   //---
   //--- s.t. sum{k in K, s in F'[k]}
   //---        sum{j in C}(a''[i,j] * s[j])* lambda[k][s] +
   //---          splus[i] - sminus[i] ~ b''[i],  i in R''
   //---      sum{s in F'[k]} lambda[k][s] = 1, k in K
   //---      lambda[k][s]                >= 0, k in K, s in F'[k]
   //---      splus[i]                    >= 0, i in R''
   //---      sminus[i]                   >= 0, i in R''
   //---
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "createMasterProblem()", m_param.LogDebugLevel, 2);
   DecompConstraintSet*           modelCore   = m_modelCore.getModel();
   assert(initVars.size() > 0);
   assert(modelCore);   //TODO: what if core is empty
   int r, c, startRow, endRow;
   int nColsCore = modelCore->getNumCols();
   //int nRowsCore = modelCore->getNumRows();
   //int nIntVars  = static_cast<int>(modelCore->integerVars.size());
   double* dblArrNCoreCols = new double[nColsCore];
   assert(dblArrNCoreCols);
   //---
   //--- set the row counts
   //---
   //m_nRowsOrig   = nRowsCore;
   m_nRowsOrig   = nColsCore;
   m_nRowsBranch = 0;//2 * nIntVars;
   m_nRowsConvex = m_numConvexCon;
   m_nRowsCuts   = 0;

   //---
   //--- set the row types of the rows
   //---   original rows, branch rows, convexity rows
   //---
   for (r = 0; r < m_nRowsOrig; r++) {
      m_masterRowType.push_back(DecompRow_Original);
   }

   //for(r = 0; r < m_nRowsBranch; r++)
   // m_masterRowType.push_back(DecompRow_Branch);
   for (r = 0; r < m_nRowsConvex; r++) {
      m_masterRowType.push_back(DecompRow_Convex);
   }

   //---
   //--- In order to implement simple branching, we are going to
   //--- treat all column bounds as explicit constraints. Then branching
   //--- for DW can be done in the same way it is done for regular CPM.
   //---    NOTE: in D, we don't need to ever branch
   //---
   //coreMatrixAppendColBounds();
   ////////
   //THINK - what need this for?
   //number of original core rows
   modelCore->nBaseRowsOrig = modelCore->nBaseRows;
   //number of original core rows plus branching rows
   modelCore->nBaseRows     = modelCore->getNumRows();
   assert(modelCore->nBaseRowsOrig == modelCore->nBaseRows);
   //---
   //--- create a matrix for the master LP
   //---  make room for original rows and convexity rows
   //---
   int                nRows    = m_nRowsOrig + m_nRowsBranch + m_nRowsConvex;
   int                nColsMax = static_cast<int>(initVars.size())
                                 + 2 * (m_nRowsOrig + m_nRowsBranch);
   double*            colLB    = new double[nColsMax];
   double*            colUB    = new double[nColsMax];
   double*            objCoeff = new double[nColsMax];
   CoinPackedMatrix* masterM  = new CoinPackedMatrix(true, 0, 0);
   vector<string>     colNames;
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 5,
              double* denseCol = new double[nRows];
              CoinAssertHint(colLB && colUB && objCoeff && denseCol && masterM,
                             "Error: Out of Memory");
             );
   //---
   //--- set the number of rows, we will add columns
   //---
   masterM->setDimensions(nRows, 0);
   //---
   //--- create artifical columns in master LP for:
   //---  original rows
   //---
   startRow = 0;
   endRow   = m_nRowsOrig;
   masterMatrixAddArtCols(masterM,
                          colLB,
                          colUB,
                          objCoeff,
                          colNames,
                          startRow, endRow, 'O');
   //startRow = m_nRowsOrig;
   //endRow   = m_nRowsOrig + m_nRowsBranch;
   //masterMatrixAddArtCols(masterM,
   //                     colLB,
   //                     colUB,
   //                     objCoeff,
   //                     colNames,
   //                     startRow, endRow, 'B');
   //TODO: should initVars be in pool and do addVarsFromPool here?
   /*M = new CoinPackedMatrix(true, 0, 0);//col-ordered
   M->setDimensions(nColsCore + m_numConvexCon, 0);

   const int n_colsArt = 2 * nColsCore; //this includes appended... ??
   int n_cols    = static_cast<int>(initVars.size());
   n_cols += n_colsArt;
   double * colLB    = new double[n_cols];
   double * colUB    = new double[n_cols];
   double * obj      = new double[n_cols];
   //double * denseCol = new double[nColsCore + m_numConvexCon];
   CoinAssertHint(colLB && colUB && obj,// && denseCol,
                  "Error: Out of Memory");

   int r, c;
   int row_index     = 0;
   int col_index     = 0;
   DecompVarList::iterator li;


   //TODO: for =, why not just use one art that is free?

   //---
   //--- min sp + sm
   //---
   //--- ax  = b --> ax + sp - sm  = b, sp >= 0, sm >= 0
   //---
   //--- ax <= b --> ax      - sm <= b,          sm >= 0
   //--- ax <= b --> ax + sp - sm <= b, sp  = 0, sm >= 0 (for now)
   //---
   //--- ax >= b --> ax + sp      >= b, sp >= 0
   //--- ax >= b --> ax + sp - sm >= b, sp >= 0, sm  = 0 (for now)
   //---
   //vector<char>   & rowSense = modelCore->rowSense;
   vector<string>   rowNames;
   vector<string>   colNames;
   string colNamePlus, colNameMinus;
   for(c = 0; c < nColsCore; c++){
      CoinPackedVector artColPlus;
      CoinPackedVector artColMinus;
      artColPlus.insert (c,  1.0);
      artColMinus.insert(c, -1.0);

      //---
      //--- append the two artificial columns to the matrix
      //---
      M->appendCol(artColPlus);
      colLB[col_index] = 0.0;
      colUB[col_index] = DecompInf;
      obj[col_index]   = 1.0;
      colNamePlus      = "sP(c_" + UtilIntToStr(col_index)
         + "_" + UtilIntToStr(c) + ")";
      col_index++;

      M->appendCol(artColMinus);
      colLB[col_index] = 0.0;
      colUB[col_index] = DecompInf;
      obj[col_index]   = 1.0;
      colNameMinus = "sM(c_" + UtilIntToStr(col_index)
         + "_" + UtilIntToStr(c) + ")";
      col_index++;

      colNames.push_back(colNamePlus);
      colNames.push_back(colNameMinus);
   }
   */
   int colIndex     = 0;
   int blockIndex   = 0;
   DecompVarList::iterator li;

   //TODO:
   //  this should be calling a function to add var to lp so don't dup code
   for (li = initVars.begin(); li != initVars.end(); li++) {
      //---
      //--- appending these variables (lambda) to end of matrix
      //---   after the artificials
      //---
      colIndex         = masterM->getNumCols();
      m_colIndexUnique = colIndex;
      //---
      //--- store the col index for this var in the master LP
      //---   NOTE: if we remove columns, this will be wrong
      //---
      (*li)->setColMasterIndex(colIndex);
      //---
      //--- we expect the user to define the block id in the var object
      //---
      blockIndex = (*li)->getBlockId();
      //---
      //--- give the column a name
      //---
      string colName = "lam(c_" + UtilIntToStr(m_colIndexUnique)
                       + ",b_" + UtilIntToStr(blockIndex) + ")";
      colNames.push_back(colName);
      UTIL_DEBUG(m_param.LogDebugLevel, 5,
                 (*li)->print(m_osLog, m_app);
                );
      //---
      //--- the column is just the vector s
      //---
      CoinPackedVector* sparseCol = 0;

      if ((*li)->m_s.getNumElements() > 0) {
         sparseCol = new CoinPackedVector((*li)->m_s);
      } else {
         sparseCol = new CoinPackedVector();
      }

      //---
      //--- append the coeff for the approriate convexity constraint
      //---
      sparseCol->insert(nColsCore + blockIndex, 1.0);
      UTIL_DEBUG(m_param.LogDebugLevel, 5,
                 (*m_osLog) << "\nSparse Col: \n";
                 UtilPrintPackedVector(*sparseCol, m_osLog);
                );
      //TODO: check for duplicates (against m_vars)
      //      or force initVars to be sent in with no dups?
      //TODO: do in const blocks
      //---
      //--- append the sparse column to the matrix
      //---
      masterM->appendCol(*sparseCol);
      colLB[colIndex]    = 0.0;
      colUB[colIndex]    = DecompInf;
      objCoeff[colIndex] = 0.0; //for D, we are in PHASEI the whole time
      //---
      //--- set master column type
      //---
      m_masterColType.push_back(DecompCol_Structural);
      //---
      //--- clean-up
      //---
      UTIL_DELPTR(sparseCol); //THINK
   } //END: for(li = initVars.begin(); li != initVars.end(); li++)

   //---
   //--- insert the initial set of variables into the master variable list
   //---
   //THINK: now doing in loop, so can check for dups
   appendVars(initVars);
   //---
   //--- row bounds from core inclding
   //---   original rows (= x*)
   //---
   vector<double> masterRowLB;
   vector<double> masterRowUB;

   for (c = 0; c < nColsCore; c++) {
      masterRowLB.push_back(m_xhatD[c]);
      masterRowUB.push_back(m_xhatD[c]);
   }

   //---
   //--- row bounds for convexity constraints
   //---
   for (r = 0; r < m_numConvexCon; r++) {
      masterRowLB.push_back(1.0);
      masterRowUB.push_back(1.0);
   }

   //---
   //--- load the problem into master's solver interface
   //---
   assert(masterM->getNumRows() == static_cast<int>(masterRowLB.size()));
   assert(masterM->getNumRows() == static_cast<int>(masterRowUB.size()));
   assert(masterM->getNumRows() == static_cast<int>(m_masterRowType.size()));
   assert(masterM->getNumCols() == static_cast<int>(m_masterColType.size()));
   m_masterSI->loadProblem(*masterM,
                           colLB, colUB, objCoeff,
                           &masterRowLB[0],
                           &masterRowUB[0]);

   //---
   //--- load column and row names to OSI
   //---
   if (modelCore->colNames.size() > 0) {
      m_masterSI->setIntParam(OsiNameDiscipline, 2);   //Full-Names
   }

   if (modelCore->colNames.size() > 0) {
      assert(static_cast<int>(modelCore->colNames.size()) ==
             modelCore->getNumCols());
      m_masterSI->setRowNames(modelCore->colNames,
                              0,
                              static_cast<int>(modelCore->colNames.size()),
                              0);
      vector<string> conRowNames;

      for (r = 0; r < m_numConvexCon; r++) {
         string rowName = "conv(b_" + UtilIntToStr(r) + ")";
         conRowNames.push_back(rowName);
      }

      m_masterSI->setRowNames(conRowNames,
                              0,
                              static_cast<int>(conRowNames.size()),
                              static_cast<int>(modelCore->colNames.size()));
   }

   if (colNames.size() > 0)
      m_masterSI->setColNames(colNames,
                              0,
                              static_cast<int>(colNames.size()),
                              0);

   UTIL_DEBUG(m_param.LogDebugLevel, 4,

   for (r = 0; r < m_masterSI->getNumRows(); r++) {
   const string rowN = m_masterSI->getRowName(r);
      printf("Row[%4d] Name: %30s Type: %20s\n",
             r,
             rowN.c_str(),
             DecompRowTypeStr[m_masterRowType[r]].c_str());
   }
   for (c = 0; c < m_masterSI->getNumCols(); c++) {
   const string colN = m_masterSI->getColName(c);
      printf("Col[%4d] Name: %30s Type: %20s\n",
             c,
             colN.c_str(),
             DecompColTypeStr[m_masterColType[c]].c_str());
   }
             );
   //---
   //--- reset unique col index id
   //---
   m_colIndexUnique = masterM->getNumCols();
   //---
   //--- free local memory
   //---
   UTIL_DELPTR(masterM);
   UTIL_DELARR(colLB);
   UTIL_DELARR(colUB);
   UTIL_DELARR(objCoeff);
   UTIL_DELARR(dblArrNCoreCols);
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "createMasterProblem()", m_param.LogDebugLevel, 2);
}
