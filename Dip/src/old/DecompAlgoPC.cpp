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

#include "DecompApp.h"
#include "DecompVar.h"
#include "DecompAlgoPC.h"
#include "DecompPortable.h"

// ------------------------------------------------------------------------- //
const char* DecompAlgoPC::classTag         = "\nD-ALGOPC      : ";

// --------------------------------------------------------------------- //
void DecompAlgoPC::setMasterBounds(const double* lbs,
                                   const double* ubs)
{
   int c;
   const int n_cols = m_modelCore->getNumCols();
   const int beg    = m_modelCore->nBaseRowsOrig;
   //keep in memory
   int*     index  = new int[n_cols];
   double* bounds = new double[2*n_cols];

   for (c = 0; c < n_cols; c++) {
      index[c]      = beg + c;
      bounds[2*c]   = lbs[c];
      bounds[2*c+1] = ubs[c];
   }

   //for(c = 0; c < n_cols; c++){
   //set all at once would be faster
   //m_masterSI->setRowBounds(beg + c, lbs[c], ubs[c]);
   //}
   m_masterSI->setRowSetBounds(index, index + n_cols, bounds);
   UTIL_DELARR(index);
   UTIL_DELARR(bounds);
}

// ------------------------------------------------------------------------- //
void DecompAlgoPC::createMasterProblem(DecompVarList& initVars)
{
   //---
   //--- Initialize the solver interface for the master problem.
   //--- PC: min c(s lam)
   //---     A''s   lam   >= b'',
   //---     sum{s} lam_s  = 1  ,
   //---            lam_s >= 0  , s in F'[0]
   //--- m_modelCore contains [A'', b''], from the original model in
   //--- terms of x. In this function we create the DW-LP in terms of
   //--- lambda, [[A''s, 1], [b'', 1]] and load that into the OSI
   //--- interface m_masterSI.
   //---
   //if 0 is feasible to subproblem, we can relax convexity to <= 1
   int i;
   char sense;
   double rhs, range;
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << classTag << " <---- createMasterProblem() ---- ";
             );
   assert(initVars.size() > 0);
   //NEW idea - THINK
   //change m_modelCore to have bounds as explicit rows?
   //that is NOT true... you change modelCore for cuts anyway - right?
   //do it for now... but realistically... i think app has its copy
   //of core, and algo needs to make a copy of it for self...
   //the working copy... then adding cuts, etc is fine -- UGH
   //but adjusting core goes against the fact that CPM and PC use same core
   //need to treat these special anyway... ? and multiple polytopes idea?
   const int        nColsCore = m_modelCore->getNumCols();
   vector<double> & colLBCore = m_modelCore->colLB; //access methods?
   vector<double> & colUBCore = m_modelCore->colUB;
#ifdef EXPLICIT_BOUNDS
   //need to do in blocks, else slow - if append one at a time (64,000 cols)
   const int numRows = nColsCore;
   //TODO: save in memory
   int* rowStarts = new int[numRows + 1];
   int* rowInd    = new int[numRows];
   double* rowEls    = new double[numRows];
   rowStarts[0] = 0;

   for (i = 0; i < nColsCore; i++) {
      rowStarts[i+1] = rowStarts[i] + 1;
      rowInd[i]      = i;
      rowEls[i]      = 1.0;
   }

   m_modelCore->M->appendRows(numRows, rowStarts, rowInd, rowEls);

   for (i = 0; i < nColsCore; i++) {
      //CoinPackedVector row;
      //row.insert(i, 1.0);
      //m_modelCore->M->appendRow(row);
      //need? shouldn't check these for dups anyway - right?
      //need to treat somewhat differently
      m_modelCore->rowLB.push_back(colLBCore[i]);
      m_modelCore->rowUB.push_back(colUBCore[i]);
      UtilBoundToSense(colLBCore[i], colUBCore[i], DecompInf,
                       sense, rhs, range);
      m_modelCore->rowRhs.push_back(rhs);
      m_modelCore->rowSense.push_back(sense);
      //THINK
      m_modelCore->rowHash.push_back(
         UtilCreateStringHash(1,
                              rowInd + i,
                              rowEls + i,
                              sense, rhs)
      );
   }

   UTIL_DELARR(rowStarts);
   UTIL_DELARR(rowInd);
   UTIL_DELARR(rowEls);
#endif
   int nRowsCore = m_modelCore->getNumRows();
   m_modelCore->nBaseRowsOrig = m_modelCore->nBaseRows;
   m_modelCore->nBaseRows = nRowsCore;
   //THINK? should initVars be in pool and do addVarsFromPool here?
   CoinPackedMatrix* M = new CoinPackedMatrix(true, 0, 0);
   M->setDimensions(nRowsCore + 1, 0);
   //THINK? cutpool contains cut and expanded row - is that correct?
   //there should be a function that "addsVars" but from an argument sent
   //in... this time directly, later from the "pool"??? THINK??
   const int n_cols = static_cast<int>(initVars.size());
   double* colLB = new double[n_cols];
   double* colUB = new double[n_cols];
   double* obj   = new double[n_cols];
   double* denseCol = new double[nRowsCore + 1];
   int col_index     = 0;
   DecompVarList::iterator li;

   for (li = initVars.begin(); li != initVars.end(); li++) {
      UTIL_DEBUG(m_app->m_param.LogDebugLevel, 5,
                 (*li)->print(m_osLog, m_app);
                );
      // ---
      // --- get dense column = A''s, append convexity constraint on end
      // ---
      m_modelCore->M->times((*li)->m_s, denseCol);
      denseCol[nRowsCore] = 1.0;
#if 0

      for (i = 0; i < nColsCore; i++) {
         denseCol[nRowsCore + 1 + i] = 0.0;//THINK
      }

      const int      n_els = (*li)->m_s.getNumElements();

      const int*     ind   = (*li)->m_s.getIndices();

      const double* els   = (*li)->m_s.getElements();

      for (i = 0; i < n_els; i++) {
         denseCol[nRowsCore + 1 + ind[i]] = els[i];
      }

#endif
      // ---
      // --- create a sparse column from the dense column
      // ---
      // THINK: do i need a DecompCol?
      // THINK: does this allocate memory for coinpackedvec twice?
      CoinPackedVector* sparseCol
      = UtilPackedVectorFromDense(nRowsCore + 1,
                                  denseCol, m_app->m_param.TolZero);
      UTIL_DEBUG(m_app->m_param.LogDebugLevel, 5,
                 (*m_osLog) << "\nSparse Col: \n";
                 UtilPrintPackedVector(*sparseCol, m_osLog);
                );
      //
      // --- check for duplicates (against m_vars)
      // ?? or force initVars to be sent in with no dups?
      //
      //TODO: do in const blocks
      // ---
      // --- append the sparse column to the matrix
      // ---
      M->appendCol(*sparseCol);
      colLB[col_index] = 0.0; //THINK: (*li)->getLowerBound();
      colUB[col_index] = DecompInf; //THINK: (*li)->getUpperBound(); //FIX!!
      obj[col_index]   = (*li)->getOriginalCost(); //c.s
      col_index++;
      UTIL_DELPTR(sparseCol); //THINK
   }

   //---
   //--- insert the initial set of variables into the master variable list
   //---
   //THINK: now doing in loop, so can check for dups
   appendVars(initVars);
   //---
   //--- THINK: do we want to adjust m_modelCore directly here?
   //--- adjust row bounds for convexity constraint
   //---
   double* zeroSol    = new double[nColsCore];
   UtilFillN(zeroSol, nColsCore, 0.0);
   bool     isZeroFeas = isIPFeasible(zeroSol);

   if (isZeroFeas) {
      printf("\nZERO SOLUTION IS FEASIBLE, RELAX CONVEX");
   }

   vector<double> masterRowLB(m_modelCore->rowLB);
   vector<double> masterRowUB(m_modelCore->rowUB);

   if (isZeroFeas) {
      masterRowLB.push_back(-DecompInf);
      masterRowUB.push_back(1.0);
   } else {
      masterRowLB.push_back(1.0);
      masterRowUB.push_back(1.0);
   }

   //---
   //--- load the problem into master's solver interface
   //---
   assert(M->getNumRows() == static_cast<int>(masterRowLB.size()));
   assert(M->getNumRows() == static_cast<int>(masterRowUB.size()));
   m_masterSI->loadProblem(*M,  colLB, colUB, obj,
                           &masterRowLB[0],
                           &masterRowUB[0]);
   // ---
   // --- free local memory
   // ---
   UTIL_DELPTR(M);
   UTIL_DELARR(denseCol);
   UTIL_DELARR(colLB);
   UTIL_DELARR(colUB);
   UTIL_DELARR(obj);
   UTIL_DELARR(zeroSol);
   //
   //createMasterStabilization();
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << classTag << "  ---- createMasterProblem() ----> ";
             );
}

/*---------------------------------------------------------------------------*/
void DecompAlgoPC::createMasterStabilization()
{
   //---
   //--- Append the variables which define the stabilizing factor
   //--- for the DW master.
   //---
   //--- PC: min c(s lam) - pi_ w + pi_ z
   //---     A''s   lam   - w + z >= b'',
   //---     sum{s} lam_s  = 1  ,
   //---            lam_s >= 0  , s in F'[0]
   //---           0 <= w <= eps, i in |A''|
   //---           0 <= z <= eps, i in |A''|
   //---
   //the problem here is now the first few vars are e.pe.'s of
   //relaxation, then new few vars are stablization, then after
   //that more columns get appended... this will screw things up
   //if you ever access m_masterSI->M... etc... or the column index
   //which you use for recomposeSolution, etc...
   int i;
   const double eps0 = 0.1;
   //do addCols for speed
   int r;
   int n_rows = m_modelCore->getNumRows();
   int n_cols = m_masterSI->getNumCols();
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << classTag
              << " <---- createMasterStabilization() ---- ";
             );

   for (i = 0; i < n_cols; i++) {
      isStab.push_back(false);
   }

   for (r = 0; r < n_rows; r++) {
      CoinPackedVector colw;
      CoinPackedVector colz;
      colw.insert(r, -1.0);
      colz.insert(r,  1.0);
      //CAREFUL - n_cols will be different now...
      m_masterSI->addCol(colw, 0.0, eps0, -m_piEstimate[r]);
      m_masterSI->addCol(colz, 0.0, eps0,  m_piEstimate[r]);
      isStab.push_back(true);//colw
      isStab.push_back(true);//colz
   }

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << classTag
              << "  ---- createMasterStabilization() ----> ";
             );
}

/*---------------------------------------------------------------------------*/
//PC only? think about name?
void DecompAlgoPC::recomposeSolution(const double* solution,
                                     double*        rsolution)
{
   //STOP - wrong with stabilization...
   //masterSol not corresponding correctly...
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << classTag << " <---- recomposeSolution()   ---- ";
             );
   UtilFillN(rsolution, m_modelCore->getNumCols(), 0.0);
   //m_setD.clear();
   int col_index = 0;
   DecompVarList::const_iterator li;

   //m_vars will contain the columns that got brought in,
   //but we will also have some extras now - which are stab vars
   //and should not be included in the recomposition - and
   //it will be ordered such that we have
   //init vars, stab vars, new vars... so we have to be careful
   //we don't have a mapping from m_vars to col_index, we might be
   //compressing, etc - so we would lose this unless we keep updating
   //the index?
   for (li = m_vars.begin(); li != m_vars.end(); li++) {
      //while(isStab[col_index] == true)
      //  col_index++;
      assert(col_index < m_masterSI->getNumCols());

      //STOP
      if (solution[col_index] > m_app->m_param.TolZero) {
         CoinPackedVector& v = (*li)->m_s;
         const int*     inds  = v.getIndices();
         const double* els   = v.getElements();

         for (int i = 0; i < v.getNumElements(); i++) {
            rsolution[inds[i]] += els[i] * solution[col_index];
            UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
                       printf("\nrsolution[ind[%d]: %d] = %g, els: %g, sol: %g",
                              i, inds[i], rsolution[inds[i]], els[i],
                              solution[col_index]);
                      );
         }

         //m_setD.push_back(*li);
      }

      col_index++;
   }

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << classTag << "  ---- recomposeSolution()   ---->";
             );
}


/*-------------------------------------------------------------------------*/
//because rowReform, this is very specific to PC
void DecompAlgoPC::addCutsToPool(const double*    x,
                                 DecompCutList& newCuts,
                                 int&            n_newCuts)
{
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << classTag << " <---- addCutsToPool()       ---- ";
             );
   unsigned int i;
   int r;
   int  cutIndex = 0;
   bool isDupCore;//also check relax?
   bool isDupPool;
   bool isViolated; //TODO: do something similiar to check for pos-rc vars
   bool addCut;
   DecompCutPool::iterator ci;
   DecompCutList::iterator li = newCuts.begin();

   while (li != newCuts.end()) {
      CoinPackedVector* row       = new CoinPackedVector();
      //---
      //--- create a row (in terms of original formulation, x), from a cut
      //---
      (*li)->expandCutToRow(row);
      //---
      //--- set the hash string (for quick duplicate checks)
      //---
      (*li)->setStringHash(row);
      bool isOptViolated = false;

      for (i = 0; i < m_optPoint.size(); i++) {
         isOptViolated = (*li)->calcViolation(row, &m_optPoint[i][0]);

         if (isOptViolated) {
            (*m_osLog) << "\n\nCUT VIOLATES OPT POINT";
            (*li)->print();
         }

         assert(!isOptViolated);
      }

      //here we will use a hash table - or just STL map
      addCut    = true;
      isDupCore = false;

      for (r = 0; r < m_modelCore->getNumRows(); r++) {
         //override isSame( )
         //in one case you check hash if expanded
         //in user case you check isSame directly
         //this will become hash lookup code
         if (m_modelCore->rowHash[r] == (*li)->getStrHash()) {
            (*m_osLog) << "\n\nCUT IS DUPLICATE with Core";
            (*li)->print();
            isDupCore = true;
            break;
         }
      }

      if (isDupCore) {
         addCut = false;
      } else {
         //---
         //--- is this cut already in pool
         //---
         isDupPool = false;

         for (ci = m_cutpool.begin(); ci != m_cutpool.end(); ci++) {
            if ((*li)->getStrHash() == (*ci).getCutPtr()->getStrHash()) {
               (*m_osLog) << "\n\nCUT IS DUPLICATE with Pool";
               (*li)->print();
               isDupPool = true;
               break;
            }
         }

         if (isDupPool) {
            addCut = false;
         } else {
            isViolated = (*li)->calcViolation(row, x);//also sets it

            if (!isViolated) {
               addCut = false;
            }
         }
      }

      if (addCut) {
         //---
         //--- create a row (in terms of reformulation, lambda), from row
         //---
         CoinPackedVector* rowReform
         = m_cutpool.createRowReform(m_modelCore->getNumCols(),
                                     row, m_vars);

         if (!rowReform) {
            assert(0);
            //vi = erase(vi);//THINK...
         } else {
            DecompWaitingRow waitingRow(*li, row, rowReform);
            //do this in a separate function so addCutsTo is not dependent
            //on passing in osolution for DecompVar
            //waitingRow.setViolation(x);//always on original solution!
            m_cutpool.push_back(waitingRow);
         }

         li++;
      } else {
         //---
         //--- cut is not violated, do not put in cut pool, delete memory
         //---
         UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
                    (*m_osLog) << "\nCUT " << cutIndex
                    << "is not violated, do not put in pool";
                   );
         UTIL_DELPTR(row);
         UTIL_DELPTR(*li); //need to do?
         li = newCuts.erase(li); //does this call cut destructor?
         n_newCuts--;
         //then don't increment li next iter?? think
      }

      cutIndex++;
   }

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << classTag << "  ---- addCutsToPool()       ---->";
             );
   assert(n_newCuts >= 0);
}



/*---------------------------------------------------------------------------*/
int DecompAlgoPC::addCutsFromPool()
{
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << classTag << " <---- addCutsFromPool()     ---- ";
             );
   //TODO partial sort!
   sort(m_cutpool.begin(),
        m_cutpool.end(),
        is_greater_thanD());
#if 0

   if (m_app->m_param.logLevel > 10) {
      cout << "\nCUT POOL BEFORE:" << endl;
      m_cutpool.print();
      //write a func for this, member of DecompModel
      cout << "\nMODEL CUTS BEFORE: " << endl;
      DecompCutList::iterator it;
      int row_index = 0;

      for (it = m_model.cuts.begin(); it != m_model.cuts.end(); it++) {
         cout << row_index << " : ";
         (*it)->print();
         cout << endl;
         row_index++;
      }

      cout << endl;
   }

#endif
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << "\nCUT POOL BEFORE:\n";
              m_cutpool.print(m_osLog);
              (*m_osLog) << "\nCUTS BEFORE:\n";
              printCuts(m_osLog);
             );
   const int maxcuts_toadd = 100;//m_app->m_param.cut_maxcuts_periter;
   int n_newrows = CoinMin(static_cast<int>(m_cutpool.size()), maxcuts_toadd);
   //since we use a list - find_first won't help as it returns an
   //iterator not an index in the list... UGH
   int index = 0;
   DecompCutPool::iterator li;

   for (li = m_cutpool.begin(); li != m_cutpool.end(); li++) {
      if ((*li).getViolation() < 0.00001) { //PARM
         break;
      }

      index++;
   }

   n_newrows = std::min<int>(n_newrows, index);//never add anything not violated

   if (n_newrows > 0) {
      m_varpool.setColsAreValid(false);
   }

   //TODO: look into coin build...
   double* rlb = new double[n_newrows];
   double* rub = new double[n_newrows];
   const CoinPackedVectorBase** rowReformBlock =
      new const CoinPackedVectorBase*[n_newrows];
   const CoinPackedVectorBase** rowBlock   =
      new const CoinPackedVectorBase*[n_newrows];
   index = 0;

   for (li = m_cutpool.begin(); li != m_cutpool.end(); li++) {
      if (index >= n_newrows) {
         break;
      }

      const CoinPackedVector* rowReform = (*li).getRowReformPtr();

      const CoinPackedVector* row       = (*li).getRowPtr();

      rowReformBlock[index] = rowReform;

      rowBlock[index]       = row;

      rlb[index] = (*li).getLowerBound();

      rub[index] = (*li).getUpperBound();

      m_cuts.push_back((*li).getCutPtr());

      m_modelCore->rowHash.push_back((*li).getCutPtr()->getStrHash());

      index++;
   }

   m_masterSI->addRows(n_newrows, rowReformBlock, rlb, rub);
   m_modelCore->M->appendRows(n_newrows, rowBlock);

   for (index = 0; index < n_newrows; index++) {
      printf("\nadding to rowLB: %g", rlb[index]);
      printf("\nadding to rowUB: %g", rub[index]);
      m_modelCore->rowLB.push_back(rlb[index]);
      m_modelCore->rowUB.push_back(rub[index]);
      double range, rhs;
      char sense;
      UtilBoundToSense(rlb[index], rub[index], DecompInf,
                       sense, rhs, range);
      m_modelCore->rowSense.push_back(sense);
      m_modelCore->rowRhs.push_back(rhs);
   }

   //clean up
   index = 0;

   for (li = m_cutpool.begin(); li != m_cutpool.end(); li++) {
      if (index >= n_newrows) {
         break;
      }

      (*li).deleteRowReform();
      (*li).deleteRow();
      (*li).clearCut();//need to do this?
      index++;
   }

   m_cutpool.erase(m_cutpool.begin(), li);
   UTIL_DELARR(rowReformBlock);
   UTIL_DELARR(rowBlock);
   UTIL_DELARR(rlb);
   UTIL_DELARR(rub);
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << "\nCUT POOL AFTER:\n";
              m_cutpool.print(m_osLog);
              (*m_osLog) << "\nCUTS AFTER:\n";
              printCuts(m_osLog);
             );
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << classTag << "  ---- addCutsFromPool()     ---->";
             );
   return n_newrows;
}


#if 0
// ------------------------------------------------------------------------- //
decompPhase DecompAlgoPC::phaseUpdate(const decompPhase   phase,
                                      const decompStat    stat,
                                      int&                n_newCuts,
                                      int&                n_newVars,
                                      int&                n_cutCalls,
                                      int&                n_priceCalls)
{
   //THINK: rather than call it price,cut - call it inner outer?
   //TODO: need a frequency parameter for how you do things...
   //a) which to start with, b) how many consecutive price, cut
   //or allow it to price out, then cut till none... many strategies to offer
   //need tailing off detection
   //options:
   // just price = (a), x = inf, y = 0
   // just cut   = (a), x = 0,   y = inf
   // (a) price out up to x iters, then cut up to y iters, repeat
   //   for now just look at (a)
   // (b) cut out up to x iters, then price up to y iters
   decompPhase nextPhase = PHASE_UNKNOWN;
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << classTag << " <---- phaseUpdate() ---- ";
             );
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << "\nn_newCuts: "   << n_newCuts;
              (*m_osLog) << " n_newVars: "    << n_newVars;
              (*m_osLog) << " n_cutCalls: "   << n_cutCalls;
              (*m_osLog) << " n_priceCalls: " << n_priceCalls;
              (*m_osLog) << "\nPHASEIN : "    << decompPhaseStr[phase] << "\t";
              (*m_osLog) << "STAT IN: "       << decompStatStr[stat];
             );

   switch (phase) {
   case PHASE_INIT:

      if (n_priceCalls < m_app->m_param.PriceIters) {
         nextPhase = PHASE_PRICE;
      } else if (n_cutCalls < m_app->m_param.CutIters) {
         nextPhase = PHASE_CUT;
      } else {
         nextPhase = PHASE_DONE;
      }

      break;
   case PHASE_PRICE:

      //---
      //--- we are pricing
      //---
      if (stat == STAT_FEASIBLE) {
         if ((n_priceCalls > m_app->m_param.PriceIters) || (n_newVars == 0)) {
            //---
            //--- we priced out
            //---
            if ((n_newCuts == 0) && (n_cutCalls > 0)) {
               //---
               //--- we didn't find any cuts last iteration either, DONE
               //---
               nextPhase = PHASE_DONE;
            } else if (n_cutCalls < m_app->m_param.CutIters) {
               //---
               //--- switch to cut phase, reset the counters
               //---
               nextPhase  = PHASE_CUT;
               n_newCuts  = 0;
               n_cutCalls = 0;
            } else {
               //---
               //--- we priced out, and called too many cut steps in a row
               //---
               nextPhase = PHASE_DONE;
            }
         } else {
            //---
            //--- we found a var last time, try again
            //---
            nextPhase = PHASE_PRICE;
         }
      } else {
         //---
         //--- master is infeasible, we have to price to break infeasiblity
         //---
         nextPhase = PHASE_PRICE;
      }

      //TODO: what about problems that are really infeasible - THINK
      break;
   case PHASE_CUT:

      //---
      //--- we are cutting
      //---
      if (stat == STAT_FEASIBLE) {
         if ((n_cutCalls >= m_app->m_param.CutIters) || (n_newCuts == 0)) {
            //---
            //---  too many cut steps in a row OR we found no violations
            //---
            if ((n_newVars == 0) && (n_priceCalls > 0)) {
               //---
               //--- we didn't find any vars last iteration either, DONE
               //---
               nextPhase = PHASE_DONE;
            } else if (n_priceCalls < m_app->m_param.PriceIters) {
               //---
               //--- switch to price phase, reset the counters
               //---
               nextPhase    = PHASE_PRICE;
               n_newVars    = 0;
               n_priceCalls = 0;
            } else {
               //---
               //--- we found no violations, and called too many price
               //--- steps in a row, DONE
               //---
               nextPhase = PHASE_DONE;
            }
         } else {
            //---
            //--- we found a cut last time, try again
            //---
            nextPhase = PHASE_CUT;
         }
      } else {
         //---
         //--- master is infeasible, we have to price to break infeasiblity
         //---
         nextPhase = PHASE_PRICE;
      }

      //TODO: what about problems that are really infeasible - THINK
      break;
   case PHASE_DONE:
      break;
   case PHASE_UNKNOWN:
   default:
      assert(0);//TODO: real error codes?
   }

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << "\nPHASEOUT: " << decompPhaseStr[nextPhase] << "\t";
             );
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << classTag << " <---- phaseUpdate() ---- ";
             );
   return nextPhase;
}
#endif
