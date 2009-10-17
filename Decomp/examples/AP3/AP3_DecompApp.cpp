//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, Lehigh University                                 //
//                                                                           //
// Copyright (C) 2002-2007, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

#include "DecompCutOsi.h"
#include "AP3_DecompApp.h"

//TODO: brute force solve

// --------------------------------------------------------------------- //
void AP3_DecompApp::initializeApp(UtilParameters & utilParam) 
   throw(CoinError)
{
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "initializeApp()", m_param.LogDebugLevel, 2);
   
   //---
   //--- get application parameters
   //---   
   m_appParam.getSettings(utilParam);   
      
   string fileName = m_appParam.DataDir 
      + UtilDirSlash() + m_appParam.Instance + ".txt";
   m_ap3data.readInstance(fileName.c_str());
   m_ap3data.m_instance = m_appParam.Instance;

   string fileNameSol = m_appParam.DataDir 
      + UtilDirSlash() + m_appParam.Instance + ".sol";
   m_ap3data.readOptimalBound(fileNameSol.c_str());
   
   int d;
   int dimension = m_ap3data.m_dimension;
   m_assigncostMin  = new double*[dimension];
   m_assignindexMin = new int*   [dimension];
   CoinAssertHint(m_assigncostMin && m_assignindexMin, "Error: Out of Memory");
   for(d = 0; d < dimension; d++){
      m_assigncostMin[d]  = new double[dimension];
      m_assignindexMin[d] = new int   [dimension];
      CoinAssertHint(m_assigncostMin[d] && m_assignindexMin[d], 
                 "Error: Out of Memory");
   }

   //TODO: don't bother building if not using LP solver for this
   {
      m_siAP = new OsiLpSolverInterface();
      CoinAssertHint(m_siAP, "Error: Out of Memory");

      //---
      //--- Two-Indexed Assignment Problem Relaxation (e.g., MODEL_I):
      //---
      //--- min {jk in J x K} c'_jk x_jk
      //---  sum {k in K}           x_jk = 1, forall j in J
      //---  sum {j in J}           x_jk = 1, forall k in K
      //---  x_jk in {0,1}, for all jk in J x K
      //--- where c_'jk = min{i in I} c_ijk
      //---
      //--- The model stays the same, only objective changes
      //--- so, construct only once.
      //---
      int      dimensionSq = dimension * dimension;
      int      n_rows      = 2 * dimension;
      int      n_nonzeros  = 2 * dimensionSq;
      int    * rowInd = new int   [n_nonzeros];
      int    * rowBeg = new int   [n_rows + 1];
      int    * rowLen = new int   [n_rows    ];
      double * rowEls = new double[n_nonzeros];
      double * colLB  = new double[dimensionSq];
      double * colUB  = new double[dimensionSq];
      double * rowB   = new double[n_rows];
      CoinAssertHint(rowInd && rowBeg && rowEls && colLB && colUB && rowB,
                     "Error: Out of Memory");
      
      CoinFillN(rowEls, n_nonzeros,  1.0);
      CoinFillN(rowB,   n_rows,      1.0);
      CoinFillN(colLB,  dimensionSq, 0.0);
      CoinFillN(colUB,  dimensionSq, 1.0);
      
      int rowIndex, ind1, ind2, nz_index;      

      rowIndex         = 0;
      rowBeg[rowIndex] = 0;
      nz_index         = 0;

      //---
      //---  sum {k in K}           x_jk = 1, forall j in J
      //---      
      for(ind1 = 0; ind1 < dimension; ind1++){        
         for(ind2 = 0; ind2 < dimension; ind2++){
            rowInd[nz_index++] = index2(ind1,ind2);
         }
         rowBeg[rowIndex+1] = rowBeg[rowIndex] + dimension;
         rowLen[rowIndex]   = dimension;
         rowIndex++;
      }
      
      //---
      //---  sum {j in J}           x_jk = 1, forall k in K
      //---
      for(ind2 = 0; ind2 < dimension; ind2++){
         for(ind1 = 0; ind1 < dimension; ind1++){
            rowInd[nz_index++] = index2(ind1,ind2);
         }
         rowBeg[rowIndex+1] = rowBeg[rowIndex] + dimension;
         rowLen[rowIndex]   = dimension;
         rowIndex++;
      }
      
      CoinPackedMatrix M(false, dimensionSq, n_rows, n_nonzeros,
                         rowEls, rowInd, rowBeg, rowLen);
      m_siAP->loadProblem(M, colLB, colUB, colLB, rowB, rowB);
      
      UTIL_DELARR(rowInd);
      UTIL_DELARR(rowBeg);
      UTIL_DELARR(rowLen);
      UTIL_DELARR(rowEls);
      UTIL_DELARR(colLB);
      UTIL_DELARR(colUB);
      UTIL_DELARR(rowB);      
   }
   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "initializeApp()", m_param.LogDebugLevel, 2);
}

// --------------------------------------------------------------------- //
void
AP3_DecompApp::APPcreateModel(double                        *& objCoeff,
                              map<int, DecompConstraintSet*> & modelCore,
                              map<int, vector<DecompConstraintSet* > > & modelRelax)
{
   
   //---
   //--- createModel is a pure virtual method of DecompApp and must 
   //--- be derived by the application class to define the partitioning
   //--- of constraints into [A,b] = [A',b'] union [A'', b'']
   //---
   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "APPcreateModel()", m_param.LogDebugLevel, 2);

   //---
   //--- Three-Indexed Assignment Problem:
   //---
   //--- min {ijk in I x J x K} c_ijk x_ijk
   //---  sum {jk in J x K}            x_ijk = 1, forall i in I
   //---  sum {ik in I x K}            x_ijk = 1, forall j in J
   //---  sum {ij in I x J}            x_ijk = 1, forall k in K
   //---  x_ijk in {0,1}, for all ijk in I x J x K
   //---
  
   //---
   //--- Partition the problem into I/J/K.
   //---   The modelCore  will be one set of constraints (I, J, or K).
   //---   The modelRelax will be the rest (an instance of 2AP).
   //---
   int dimension   = m_ap3data.m_dimension;
   int dimensionSq = dimension * dimension;
   int n_cols      = m_ap3data.m_ncolsFull;

   //---
   //--- open temporary memory for use in storing rows
   //---
   int    * rowInd = new int   [dimensionSq];
   double * rowEls = new double[dimensionSq];
   //if(!(rowInd && rowEls))
   //  CoinAssert("Error: Out of Memory");
  
   //---
   //--- open memory for the objective coefficients of modelCore
   //---  who is responsible to open this memory? AP3_DecompApp
   //---  who is responsible to free this memory? DecompAlgo
   //---
   objCoeff    = new double[n_cols];
   memcpy(objCoeff, m_ap3data.m_assigncost, n_cols * sizeof(double));
   
   //---
   //--- create Model I
   //---
   vector< DecompConstraintSet* > modelRelaxVI;   
   DecompConstraintSet * modelCoreI  = new DecompConstraintSet();
   DecompConstraintSet * modelRelaxI = new DecompConstraintSet();  
   createModelPart(MODEL_I, rowInd, rowEls, modelCoreI, modelRelaxI);
   modelRelaxVI.push_back(modelRelaxI);
   modelCore.insert(make_pair(MODEL_I, modelCoreI));
   modelRelax.insert(make_pair(MODEL_I, modelRelaxVI));
  
   //---
   //--- create Model J
   //---
   vector< DecompConstraintSet* > modelRelaxVJ;   
   DecompConstraintSet * modelCoreJ  = new DecompConstraintSet();
   DecompConstraintSet * modelRelaxJ = new DecompConstraintSet();
   createModelPart(MODEL_J, rowInd, rowEls, modelCoreJ, modelRelaxJ);
   modelRelaxVJ.push_back(modelRelaxJ);
   modelCore.insert(make_pair(MODEL_J, modelCoreJ));
   modelRelax.insert(make_pair(MODEL_J, modelRelaxVJ));

   //---
   //--- create Model K
   //---
   vector< DecompConstraintSet* > modelRelaxVK;   
   DecompConstraintSet * modelCoreK  = new DecompConstraintSet();
   DecompConstraintSet * modelRelaxK = new DecompConstraintSet();
   createModelPart(MODEL_K, rowInd, rowEls, modelCoreK, modelRelaxK);
   modelRelaxVK.push_back(modelRelaxK);
   modelCore.insert(make_pair(MODEL_K, modelCoreK));
   modelRelax.insert(make_pair(MODEL_K, modelRelaxVK));
  
   //---
   //--- TODO: weird if objCoeff is allocated here, but not deleted
   //--- free local memory
   //--- 
   UTIL_DELARR(rowInd);
   UTIL_DELARR(rowEls);

   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "APPcreateModel()", m_param.LogDebugLevel, 2);
}

//--------------------------------------------------------------------- //
void AP3_DecompApp::createModelPart(const int             modelType,
                                    int                 * rowInd,
                                    double              * rowEls,
                                    DecompConstraintSet * modelCore,
                                    DecompConstraintSet * modelRelax)
   throw(CoinError) 
{

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createModelPart()", m_param.LogDebugLevel, 2);
   
   const int dimension   = m_ap3data.m_dimension;
   const int dimensionSq = dimension * dimension;
   const int n_cols      = m_ap3data.m_ncolsFull;
   const int n_rows      = m_ap3data.m_nrowsFull;
   const int n_rowsThird = static_cast<int>(n_rows / 3);   

   //---
   //--- for multi-polytope, we must choose a fixed A''
   //---     A = [A'', A'[k]]
   //---
   //--- it always must be true that
   //---     A'' inter A[k] contains A
   //---
   //--- so, if A[k] is a partition (rather than nested), as
   //--- it is in AP3, then we are forced to use A'' = A
   //---
      
   //---
   //--- set the constraint matrix of modelCore and modelRelax
   //---  who is responsible to open this memory? AP3_DecompApp
   //---  who is responsible to free this memory? DecompAlgo 
   //---
   modelCore->M = new CoinPackedMatrix(false, 0.0, 0.0);
   CoinAssertHint(modelCore->M, "Error: Out of Memory");
   modelCore->M->setDimensions(0, n_cols);

   if(m_param.PriceMultiPoly)
      modelCore->M->reserve(n_rows, n_rows * dimensionSq);
   else
      modelCore->M->reserve(n_rowsThird, n_rowsThird * dimensionSq);

   modelRelax->M = new CoinPackedMatrix(false, 0.0, 0.0);
   CoinAssertHint(modelRelax->M, "Error: Out of Memory");
   modelRelax->M->setDimensions(0, n_cols);
   modelRelax->M->reserve(2 * n_rowsThird, 2 * n_rowsThird * dimensionSq);

   int (AP3_DecompApp::*indexFunc)(const int, const int, const int) const = 0;
   switch(modelType){
   case MODEL_I:      
      //---
      //---  sum {jk in J x K}            x_ijk = 1, forall i in I [Core ]
      //---  sum {ik in I x K}            x_ijk = 1, forall j in J [Relax]
      //---  sum {ij in I x J}            x_ijk = 1, forall k in K [Relax]
      //---
      indexFunc = &AP3_DecompApp::indexIJK;
      break;
   case MODEL_J:      
      //---
      //---  sum {jk in J x K}            x_ijk = 1, forall i in I [Relax]
      //---  sum {ik in I x K}            x_ijk = 1, forall j in J [Core ]
      //---  sum {ij in I x J}            x_ijk = 1, forall k in K [Relax]
      //---
      indexFunc = &AP3_DecompApp::indexJIK;
      break;
   case MODEL_K:      
      //---
      //---  sum {jk in J x K}            x_ijk = 1, forall i in I [Relax]
      //---  sum {ik in I x K}            x_ijk = 1, forall j in J [Relax]
      //---  sum {ij in I x J}            x_ijk = 1, forall k in K [Core ]
      //---
      indexFunc = &AP3_DecompApp::indexKIJ;
      break;
   default:
      CoinAssertHint(0, "Error: Bad Argument for modelType");            
   }
   
   
   int ind1, ind2, ind3, len;
   CoinFillN(rowEls, dimensionSq, 1.0);
   for(ind1 = 0; ind1 < dimension; ind1++){      
      len = 0;
      for(ind2 = 0; ind2 < dimension; ind2++){
         for(ind3 = 0; ind3 < dimension; ind3++){
            rowInd[len++] = (this->*indexFunc)(ind1,ind2,ind3);
         }
      }
      CoinAssertHint(len == dimensionSq, "Error in construction len != n^2");
      modelCore->M->appendRow(len, rowInd, rowEls);
      
      len = 0;
      for(ind2 = 0; ind2 < dimension; ind2++){
         for(ind3 = 0; ind3 < dimension; ind3++){
            rowInd[len++] = (this->*indexFunc)(ind2,ind1,ind3);
         }
      }
      CoinAssertHint(len == dimensionSq, "Error in construction len != n^2");
      if(m_param.PriceMultiPoly)
         modelCore->M->appendRow(len, rowInd, rowEls);
      modelRelax->M->appendRow(len, rowInd, rowEls);
      
      len = 0;
      for(ind2 = 0; ind2 < dimension; ind2++){
         for(ind3 = 0; ind3 < dimension; ind3++){
            rowInd[len++] = (this->*indexFunc)(ind2,ind3,ind1);
         }
      }
      CoinAssertHint(len == dimensionSq, "Error in construction len != n^2");
      if(m_param.PriceMultiPoly)
         modelCore->M->appendRow(len, rowInd, rowEls);
      modelRelax->M->appendRow(len, rowInd, rowEls);   
   }
   if(m_param.PriceMultiPoly){
      CoinAssert(modelCore->M->getNumRows()  ==   3*dimension);
   }
   else{
      CoinAssert(modelCore->M->getNumRows()  ==   dimension);
   }
   CoinAssert(modelRelax->M->getNumRows() == 2*dimension);
         
   //---
   //--- set the row upper and lower bounds
   //--- set the col upper and lower bounds
   //---
   int n_CoreRows  = modelCore->M->getNumRows();
   int n_RelaxRows = modelRelax->M->getNumRows();
   UtilFillN(modelCore->rowLB,  n_CoreRows,  1.0);
   UtilFillN(modelCore->rowUB,  n_CoreRows,  1.0);
   UtilFillN(modelRelax->rowLB, n_RelaxRows, 1.0);
   UtilFillN(modelRelax->rowUB, n_RelaxRows, 1.0);

   //THINK: is colLB/UB for Core vs Relax ever different?
   UtilFillN(modelCore->colLB,  n_cols, 0.0);
   UtilFillN(modelCore->colUB,  n_cols, 1.0);
   UtilFillN(modelRelax->colLB, n_cols, 0.0);
   UtilFillN(modelRelax->colUB, n_cols, 1.0);

   //#define DEBUG_AP3_10_3
#ifdef DEBUG_AP3_10_3
   modelCore->colUB[372]  = 1;
   modelRelax->colUB[372] = 1;
   modelCore->colLB[889]  = 1;
   modelRelax->colLB[889] = 1;
#endif
 
   //---
   //--- set the indices of the integer variables of modelRelax
   //---
   UtilIotaN(modelRelax->integerVars, n_cols, 0);

   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "createModelPart()", m_param.LogDebugLevel, 2);
}


//TODO:
// rethink even this design. user doesn't need to know DecompVar or any
// of this - user can just provide a vector, a solution to oracle, and
// framework can handle checking rc, etc and pushing that into vars
// the user JUST provides a solver for oracle given a cost vector

// advanced user can do things like check for negRC... the advanced user
// can also derive their own DecompVar (if they know how to store it more
// compactly (see TSP, as an example)... but, in end, the framework must
// expand it anyway... 

//--------------------------------------------------------------------- //
DecompStatus AP3_DecompApp::APPsolveRelaxed(const int             whichModel,
                                            const double        * redCostX,
                                            const double        * origCost,
                                            const double          alpha,
                                            const int             n_origCols,
                                            const bool            checkRC,
                                            const bool            checkDup,
                                            OsiSolverInterface  * m_subprobSI,
                                            list<DecompVar*>    & vars){
   
   //TODO: use lp solver
   //TODO: use ap solver

   //---
   //--- Three-Indexed Assignment Problem:
   //---
   //--- min {ijk in I x J x K} c_ijk x_ijk
   //---  sum {jk in J x K}            x_ijk = 1, forall i in I
   //---  sum {ik in I x K}            x_ijk = 1, forall j in J
   //---  sum {ij in I x J}            x_ijk = 1, forall k in K
   //---  x_ijk in {0,1}, for all ijk in I x J x K
   //---
   //--- Two-Indexed Assignment Problem Relaxation (e.g., MODEL_I):
   //---
   //--- min {jk in J x K} c'_jk x_jk
   //---  sum {k in K}           x_jk = 1, forall j in J
   //---  sum {j in J}           x_jk = 1, forall k in K
   //---  x_jk in {0,1}, for all jk in J x K
   //--- where c_'jk = min{i in I} c_ijk
   //---

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "APPsolveRelaxed()", m_param.LogDebugLevel, 2);


#if 0
   {
      //TEMP
      int i;
      for(i = 0; i < n_origCols; i++){
         cout << "\nREDCOST: ";
         printOriginalColumn(i);
         cout << " " << redCostX[i] 
              << " orig " << origCost[i];
      }
   }
#endif


   //---
   //--- calculate c_'jk = min{i in I} c_ijk, or
   //--- calculate c_'ik = min{i in J} c_ijk, or
   //--- calculate c_'ij = min{i in K} c_ijk
   //---
   int (AP3_DecompApp::*indexFunc)(const int, const int, const int) const = 0;
   switch(whichModel){
   case MODEL_I:
      indexFunc = &AP3_DecompApp::indexIJK;
      break;
   case MODEL_J:
      indexFunc = &AP3_DecompApp::indexJIK;
      break;
   case MODEL_K:      
      indexFunc = &AP3_DecompApp::indexKIJ;
      break;
   default:
      CoinAssertHint(0, "Error: Bad Argument for modelType");      
   }
   
   double min_cost;
   int    index, ind1, ind2, ind3, min_index;
   int    dimension = m_ap3data.m_dimension;
   for(ind2 = 0; ind2 < dimension; ind2++){
      for(ind3 = 0; ind3 < dimension; ind3++){
         min_cost  = DecompInf;
         min_index = 0;
         for(ind1 = 0; ind1 < dimension; ind1++){      
            index = (this->*indexFunc)(ind1,ind2,ind3);
            if(redCostX[index] < min_cost){
               min_cost   = redCostX[index];
               min_index  = index;
            }            
         }
         m_assigncostMin [ind2][ind3] = min_cost;
         m_assignindexMin[ind2][ind3] = min_index;
         CoinAssert(min_cost < (DecompInf/2.0));
      }
   }
   
   //TODO: option to use ap solver vs lp solver

   //---
   //--- assigncostMin double array format is good for ap solver
   //---   but lp solver needs a single array version
   //--- 
   int col_index = 0;
   for(ind2 = 0; ind2 < dimension; ind2++){
      for(ind3 = 0; ind3 < dimension; ind3++){
         m_siAP->setObjCoeff(col_index, m_assigncostMin[ind2][ind3]);
         col_index++;
      }
   }   


   //---
   //--- solve the LP relaxation of the 2AP (integral polytope)
   //---
   m_siAP->messageHandler()->setLogLevel(m_param.LogLpLevel);
   m_siAP->initialSolve();
   CoinAssert(m_siAP->isProvenOptimal());

   //deal with status issues


   //---
   //--- store the solution as a DecompVar and push into list
   //---
   int            i;
   pair<int,int>  p;
   vector<int>    apInd;
   vector<double> apEls(dimension, 1.0);
   apInd.reserve(dimension);
   
   double varRedCost    = 0.0;
   double varOrigCost   = 0.0;
   const double * lpSol = m_siAP->getColSolution();
   for(i = 0; i < m_siAP->getNumCols(); i++){
      CoinAssertDebug( UtilIsZero(lpSol[i],       1.0e-4) || 
                       UtilIsZero(1.0 - lpSol[i], 1.0e-4));
      if(lpSol[i] > 0.5){
         //---
         //--- convert back from 2D to 3D case (stored in m_assignindexMin)
         //---
         p     = index2Inv(i);
         index = m_assignindexMin[p.first][p.second];
         CoinAssertDebug(index < n_origCols);
         
         varRedCost  += redCostX[index];
         varOrigCost += origCost[index];

         apInd.push_back(index);
      }
   }
   varRedCost += alpha; //RC = c-uA''s - alpha


   DecompVar * var = new DecompVar(apInd, apEls, varRedCost, varOrigCost);
   
   //TODO: framework should do all this for the user!
   bool doPush = true;
   if(checkRC && varRedCost > -1.e-10) //THINK: dualTol?
      doPush = false;
   else if(checkDup){
      DecompVarList::iterator it;
      for(it = vars.begin(); it != vars.end(); it++){
         if((*it)->isEquivalent(*var)){
            UTIL_DEBUG(m_param.LogDebugLevel, 3,
                       (*m_osLog) << "\nDuplicate variable, not adding.";
                       );
            doPush = false;
            break;
         }
      }
   }
   
   //just use as sanity check doPush = 0
   //doPush = 0;
   if(doPush){
      UTIL_DEBUG(m_param.LogDebugLevel, 5,
                 var->print();
                 );
      vars.push_back(var);
   }
   else
      UTIL_DELPTR(var);
   
   //vars.push_back(new DecompVar(apInd, apEls, varRedCost, varOrigCost));
   
   
   //TODO: why can't framework do the work of calculating redCost and obj?
   //overload var constructor so user can supply or not?, and why does user
   //need to know alpha? because they need to know what is acceptable... 
   //to push? they might want to filter out nonnegative rc

   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "APPsolveRelaxed()", m_param.LogDebugLevel, 2);
   
   return STAT_FEASIBLE; //think
}


#if 0
//--------------------------------------------------------------------- //
//STOP - see old/AAP/LP for both version, x* and DECOMP on s
int AP3_DecompApp::generateCuts(const double              * x, 
                                const DecompConstraintSet & modelCore,
                                const DecompConstraintSet & modelRelax,
                                DecompCutList             & newCuts){
  

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "generateCuts()", m_param.LogDebugLevel, 2);
   
   int                n_cuts           = 0;

   







   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "generateCuts()", m_param.LogDebugLevel, 2);

   return n_cuts;
  
}
#endif

//--------------------------------------------------------------------- //
void AP3_DecompApp::printOriginalColumn(const int   index,
                                        ostream   * os) const {
   int i, j, k;
   index3Inv(index, i, j, k);
   (*os) << "x[ " << i << " , " << j << " , " << k << "]";
}
