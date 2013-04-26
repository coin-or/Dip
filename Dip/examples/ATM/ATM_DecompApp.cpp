//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, Lehigh University                                 //
//                                                                           //
// Copyright (C) 2002-2013, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

//===========================================================================//
#include "ATM_DecompApp.h"

//===========================================================================//
void ATM_DecompApp::initializeApp(UtilParameters & utilParam) {
   

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "initializeApp()", m_appParam.LogLevel, 2);
   
   //---
   //--- get application parameters
   //---
   m_appParam.getSettings(utilParam);
   if(m_appParam.LogLevel >= 1)
      m_appParam.dumpSettings(m_osLog); 
   
   //---
   //--- read instance
   //
   string fileNameA = m_appParam.DataDir
      + UtilDirSlash() + m_appParam.DataAtm;
   string fileNameD = m_appParam.DataDir
      + UtilDirSlash() + m_appParam.DataDate;
   string fileNameAD = m_appParam.DataDir
      + UtilDirSlash() + m_appParam.DataAtmDate;
   m_instance.readInstance(fileNameA,
                           fileNameD,
                           fileNameAD);

   //---
   //--- create models
   //---
   createModels();

   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "initializeApp()", m_appParam.LogLevel, 2);
}

//===========================================================================//
void ATM_DecompApp::addColumnNamesA(DecompConstraintSet * model,
				    const string          prefix,
				    const int             offset){
   
   int       a, colIndex;
   string    colName;
   const int nAtms  = m_instance.getNAtms();
   colIndex = offset;
   for(a = 0; a < nAtms; a++){      
      colName = prefix + "("
	 + UtilIntToStr(colIndex) + "_" 
         + m_instance.getAtmName(a) + ")";
      model->colNames.push_back(colName);
      colIndex++;
   }
}

//===========================================================================//
void ATM_DecompApp::addColumnNamesAT(DecompConstraintSet * model,
				     const string          prefix,
				     const int             offset){
   
   int       a, t, colIndex;
   string    colName;
   const int nAtms  = m_instance.getNAtms();
   const int nSteps = m_appParam.NumSteps;
   colIndex = offset;
   for(a = 0; a < nAtms; a++){
      if(m_appParam.UseTightModel){
	 for(t = 0; t <= nSteps; t++){
	    colName = prefix + "("
	       + UtilIntToStr(colIndex)    + "_"
	       + m_instance.getAtmName(a) + "," 
	       + UtilIntToStr(t) + ")";
	    model->colNames.push_back(colName);
	    colIndex++;
	 }
      }
      else{
	 for(t = 0; t < nSteps; t++){
	    colName = prefix + "("
	       + UtilIntToStr(colIndex)    + "_"
	       + m_instance.getAtmName(a) + "," 
	       + UtilIntToStr(t+1) + ")";
	    model->colNames.push_back(colName);
	    colIndex++;
	 }
      }
   }
}

//===========================================================================//
void ATM_DecompApp::addColumnNamesAD(DecompConstraintSet * model,
				     const string          prefix,
				     const int             offset){

   int       a, d, colIndex;
   string    colName;
   pair<int,int>                 adP;
   vector<int>::const_iterator   vi;
   const vector<int>           & pairsAD = m_instance.getPairsAD();
   colIndex = offset;
   for(vi = pairsAD.begin(); vi != pairsAD.end(); vi++){
      adP = m_instance.getIndexADInv(*vi);
      a   = adP.first;
      d   = adP.second;            
      colName  = prefix + "("
	 + UtilIntToStr(colIndex)     + "_"
	 + m_instance.getAtmName(a)  + "," 
         + m_instance.getDateName(d) + ")";
      model->colNames.push_back(colName);
      colIndex++;
   }
}

//===========================================================================//
void ATM_DecompApp::createModelColumns(DecompConstraintSet * model,
				       const int             atmIndex,
                                       const int             dateIndex){

   //---
   //--- Column bounds:
   //---    for a in A, d in D:
   //---       f+[a,d], f-[a,d] >= 0
   //---                f-[a,d] <= w[a,d]
   //---       v[a,d] in {0,1}
   //---    for a in A:
   //---       x2[a] in [0,1], x3[a] >= 0
   //--- NOTE: x3 no UB causing unbounded subproblems? try <= 1000
   //---    for a in A, t in T
   //---       x1[a,t] in {0,1}, z[a,t] in [0,1]
   //---
   //---  If this is for block atmIndex(>=0), 
   //---      fix all a!=atmIndex in A columns to 0.
   //---
   //---  If this is for block dateIndex(>=0), 
   //---      fix all d!=dateIndex in D columns to 0.
   //---
   string    colName;
   const int nPairs     = m_instance.getNPairs();
   const int nAtms      = m_instance.getNAtms();
   const int nSteps     = m_appParam.NumSteps;
   const int nAtmsSteps = getNAtmsSteps();
   const int nCols      = numCoreCols();

   const double      * w_ad    = m_instance.get_w_ad(); //dense storage
   const vector<int> & pairsAD = m_instance.getPairsAD();   
   vector<int>::const_iterator vi;

   //---   
   //--- set the col upper and lower bounds
   //---
   UtilFillN(model->colLB, nCols,  0.0);
   UtilFillN(model->colUB, nCols,  1.0);

   int index;
   index = getColOffset_fm();
   for(vi = pairsAD.begin(); vi != pairsAD.end(); vi++){
      model->colUB[index] = w_ad[*vi];
      index++;
   }
   
   UtilFillN(&model->colUB[0] + getColOffset_fp(), nPairs, DecompInf);
   
   //Another case of needing extreme rays??
   //UtilFillN(&model->colUB[0] + getColOffset_x3(), nAtms,  DecompInf);
   UtilFillN(&model->colUB[0] + getColOffset_x3(), nAtms,  1000.0);

   //---
   //--- if this is for block a, fix all others to 0.
   //---
   if(atmIndex >= 0){
      int a, t;
      int index_x1 = getColOffset_x1();
      int index_z  = getColOffset_z();
      int index_x2 = getColOffset_x2();
      int index_x3 = getColOffset_x3();
      int end      = nSteps;
      if(m_appParam.UseTightModel)
	 end++;
      for(a = 0; a < nAtms; a++){
         if(a != atmIndex){
            model->colUB[index_x2] = 0.0;
            model->colUB[index_x3] = 0.0;
            for(t = 0; t < end; t++){
               model->colUB[index_x1++] = 0.0;
               model->colUB[index_z ++] = 0.0;
            }
         }
         else{
            for(t = 0; t < end; t++){
	       model->activeColumns.push_back(index_x1++);
	       model->activeColumns.push_back(index_z++);
	    }
	    model->activeColumns.push_back(index_x2);
	    model->activeColumns.push_back(index_x3);	    
         }
         index_x2++;
         index_x3++;
      }
      int index_fp = getColOffset_fp();
      int index_fm = getColOffset_fm();      
      int index_v  = getColOffset_v();
      for(vi = pairsAD.begin(); vi != pairsAD.end(); vi++){
	 a = m_instance.getIndexADInv(*vi).first; 
	 if(a != atmIndex){
            model->colUB[index_fp] = 0.0;
            model->colUB[index_fm] = 0.0;
            model->colUB[index_v]  = 0.0;
         }
	 else{
	    model->activeColumns.push_back(index_fp);	    
	    model->activeColumns.push_back(index_fm); 
            model->activeColumns.push_back(index_v);	    
	 }
	 index_fp++;
	 index_fm++;	 
         index_v++;
      }
   }

   if(dateIndex >= 0){
      int d;
      int index_fp = getColOffset_fm();
      int index_fm = getColOffset_fm();            
      int index_v  = getColOffset_v();
      for(vi = pairsAD.begin(); vi != pairsAD.end(); vi++){
	 d = m_instance.getIndexADInv(*vi).second; 
	 if(d != dateIndex){
            model->colUB[index_fp] = 0.0;
            model->colUB[index_fm] = 0.0;            
            model->colUB[index_v]  = 0.0;
         }
	 else{
	    model->activeColumns.push_back(index_fp);
	    model->activeColumns.push_back(index_fm); 
            model->activeColumns.push_back(index_v); 
	 } 
	 index_fp++;
	 index_fm++;	 
         index_v++;
      }
   }

   //---
   //--- set the indices of the integer variables of model: x1, v
   //---
   UtilIotaN(model->integerVars, nAtmsSteps, getColOffset_x1());
   UtilIotaN(model->integerVars, nPairs,     getColOffset_v());
   
   //---
   //--- set column names for debugging
   //---
   addColumnNamesAT(model, "x1", getColOffset_x1());
   addColumnNamesAT(model, "z",  getColOffset_z());
   addColumnNamesAD(model, "fp", getColOffset_fp());
   addColumnNamesAD(model, "fm", getColOffset_fm());
   addColumnNamesA (model, "x2", getColOffset_x2());
   addColumnNamesA (model, "x3", getColOffset_x3());
   addColumnNamesAD(model, "v",  getColOffset_v());

#if 0
   {
      int i;
      for(i = 0; i < static_cast<int>(model->colNames.size()); i++){
	 printf("COL[%4d] = %20s LB:%g UB:%g\n", 
		i, 
		model->colNames[i].c_str(),
		model->colLB[i],
		model->colUB[i]);
      }
   }
#endif
}

//===========================================================================//
int ATM_DecompApp::createConZtoX(DecompConstraintSet * model,
				 const int             atmIndex){

   //---
   //---  for a in A, t in T:
   //---     z[a,t] = x1[a,t] * x2[a] 
   //---         <==>
   //--- OLD:
   //---  for a in A:
   //---     sum{t in T}  x1[a,t] <= 1 (where, T = 1..n)
   //---  for a in A, t in T:
   //---     z[a,t] >= 0 <= 1,                  
   //---     z[a,t] <= x1[a,t],            
   //---     z[a,t] <= x2[a],              
   //---     z[a,t] >= x1[a,t] + x2[a] - 1.
   //--- NEW: 
   //---  for a in A:
   //---     sum{t in T}  x1[a,t] = 1 (where, T = 0..n)
   //---     sum{t in T}   z[a,t] = x2[a]
   //---  for a in A, t in T:
   //---     z[a,t] >= 0 <= 1,                  
   //---     z[a,t] <= x1[a,t].
   //---
   int       t;
   int       nRows = 0;
   const int nSteps = m_appParam.NumSteps;
   if(m_appParam.UseTightModel){
      for(t = 0; t <= nSteps; t++){
	 CoinPackedVector row;
	 string strAT = "(a_" + m_instance.getAtmName(atmIndex) 
	    + ",t_" + UtilIntToStr(t) + ")";
	 string rowName = "ztox1" + strAT;	 
	 row.insert(colIndex_z (atmIndex,t),  1.0);
	 row.insert(colIndex_x1(atmIndex,t), -1.0);
	 model->appendRow(row, -DecompInf, 0.0, rowName);
	 nRows++;
      }
      CoinPackedVector row2;
      string rowName2 = "ztox2(a_" 
	 + m_instance.getAtmName(atmIndex) + ")"; 
      row2.insert(colIndex_x2(atmIndex), -1.0);
      for(t = 0; t <= nSteps; t++){
	 row2.insert(colIndex_z (atmIndex,t),  1.0);
      }
      model->appendRow(row2, 0.0, 0.0, rowName2);	 
      nRows++;         
   }
   else{
      for(t = 0; t < nSteps; t++){
	 CoinPackedVector row1, row2, row3;
	 string strAT = "(a_" + m_instance.getAtmName(atmIndex) 
	    + ",t_" + UtilIntToStr(t+1) + ")";
	 string rowName1 = "ztox1" + strAT;
	 string rowName2 = "ztox2" + strAT;
	 string rowName3 = "ztox3" + strAT;
	 
	 row1.insert(colIndex_z (atmIndex,t),  1.0);
	 row1.insert(colIndex_x1(atmIndex,t), -1.0);
	 model->appendRow(row1, -DecompInf, 0.0, rowName1);
	 
	 row2.insert(colIndex_z (atmIndex,t),  1.0);
	 row2.insert(colIndex_x2(atmIndex)  , -1.0);
	 model->appendRow(row2, -DecompInf, 0.0, rowName2);
	 
	 row3.insert(colIndex_z (atmIndex,t),  1.0);
	 row3.insert(colIndex_x1(atmIndex,t), -1.0);
	 row3.insert(colIndex_x2(atmIndex)  , -1.0);
	 model->appendRow(row3, -1.0, DecompInf, rowName3);
	 
	 nRows+=3;         
      }
   }
   return nRows;
}


//===========================================================================//
int ATM_DecompApp::createConPickOne(DecompConstraintSet * model,
				    const int             atmIndex){
   //---
   //---       sum{t in T} x1[a,t] <= 1
   //---      
   int              t;
   CoinPackedVector row;
   string rowName = "pickone_x1(a_" + m_instance.getAtmName(atmIndex) + ")";
   if(m_appParam.UseTightModel){
      for(t = 0; t <= m_appParam.NumSteps; t++)
	 row.insert(colIndex_x1(atmIndex,t), 1.0);
      model->appendRow(row, 1.0, 1.0, rowName);
   }
   else{
      for(t = 0; t < m_appParam.NumSteps; t++)
	 row.insert(colIndex_x1(atmIndex,t), 1.0);
      model->appendRow(row, -DecompInf, 1.0, rowName);
   }
   return 1;
}

//===========================================================================//
int ATM_DecompApp::createConCount(DecompConstraintSet * model,
				  const int             atmIndex){
  
   //---
   //---       sum{d in D} v[a,d]  <= K[a] (for count)      
   //---
   CoinPackedVector row;
   string rowName = "count(a_" + m_instance.getAtmName(atmIndex) + ")";
   
   pair<int,int>       adP;
   int                 pairIndex = 0;
   const vector<int> & pairsAD   = m_instance.getPairsAD();
   const double      * K_a       = m_instance.get_K_a();
   vector<int>::const_iterator vi;
   //TODO: this can be faster by storing the incident d's for each a
   for(vi = pairsAD.begin(); vi != pairsAD.end(); vi++){
      adP = m_instance.getIndexADInv(*vi);
      if(atmIndex != adP.first){
	 pairIndex++;
	 continue;      
      }      
      row.insert(colIndex_v(pairIndex), 1.0);
      pairIndex++;
   }
   model->appendRow(row, -DecompInf, K_a[atmIndex], rowName);
   return 1;
}

//===========================================================================//
DecompConstraintSet * ATM_DecompApp::createModelCore1(bool includeCount){

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createModelCore1()", m_appParam.LogLevel, 2);

   int   a, d, pairIndex;

   pair<int,int>               adP;
   vector<int>::const_iterator vi;

   const int           nAtms   = m_instance.getNAtms();
   const int           nDates  = m_instance.getNDates();
   const vector<int> & pairsAD = m_instance.getPairsAD();   
   const double      * B_d     = m_instance.get_B_d();
   const int           nCols   = numCoreCols();
   int                 nRows   = nDates; 

   if(includeCount)
      nRows += nAtms;

   DecompConstraintSet * model = new DecompConstraintSet();
   CoinAssertHint(model, "Error: Out of Memory");

   model->M = new CoinPackedMatrix(false, 0.0, 0.0);
   CoinAssertHint(model->M, "Error: Out of Memory");
   model->M->setDimensions(0, nCols);
   model->reserve(nRows, nCols);

   //---
   //---    for d in D:
   //---      sum{a in A} (f+[a,d] - f-[a,d]) <= B[d]
   //---
   CoinPackedVector * rowsD = new CoinPackedVector[nDates];   
   pairIndex = 0;
   for(vi = pairsAD.begin(); vi != pairsAD.end(); vi++){
      adP = m_instance.getIndexADInv(*vi);
      a   = adP.first;
      d   = adP.second;
      rowsD[d].insert(getColOffset_fp() + pairIndex,  1.0);
      rowsD[d].insert(getColOffset_fm() + pairIndex, -1.0);
      pairIndex++;
   }

   for(d = 0; d < nDates; d++){
      string rowName = "budget(d_" 
	 + m_instance.getDateName(d) + ")";
      model->appendRow(rowsD[d], -DecompInf, B_d[d], rowName);
   }
   UTIL_DELARR(rowsD);

   if(includeCount){
      //---  for a in A:
      //---    sum{d in D} v[a,d] <= K[a]
      //---   
      for(a = 0; a < nAtms; a++){
	 createConCount(model, a);
      }
   }
   
   //---
   //--- create model columns
   //---
   createModelColumns(model);
   
   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "createModelCore1()", m_appParam.LogLevel, 2);

   return model;
}

//===========================================================================//
DecompConstraintSet * ATM_DecompApp::createModelCore2(){

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createModelCore2()", m_appParam.LogLevel, 2);

   int                 a;
   const int           nAtms   = m_instance.getNAtms();
   const int           nCols   = numCoreCols();
   const int           nRows   = (2*nAtms);
   int                 rowCnt  = 0;

   //---
   //--- A'' (core):
   //---    for a in A:
   //---       sum{t in T} x1[a,t] <= 1
   //---       sum{d in D} v[a,d]  <= K[a]
   //---    for a in A, t in T:
   //---       z[a,t] <= x1[a,t],            
   //---       z[a,t] <= x2[a],              
   //---       z[a,t] >= x1[a,t] + x2[a] - 1.
   //---
   
   DecompConstraintSet * model = new DecompConstraintSet();
   CoinAssertHint(model, "Error: Out of Memory");

   //---
   //--- create new matrix
   //---
   model->M = new CoinPackedMatrix(false, 0.0, 0.0);
   CoinAssertHint(model->M, "Error: Out of Memory");
   model->M->setDimensions(0, nCols);
   model->reserve(nRows, nCols);

   //---
   //--- create model constraints
   //---
   for(a = 0; a < nAtms; a++){
      rowCnt += createConPickOne(model, a);
      rowCnt += createConCount  (model, a);
   }
   printf("nRows = %d, rowCnt = %d\n", nRows, rowCnt);
   assert(rowCnt == nRows);
   
   //---
   //--- create model variables/columns
   //---
   createModelColumns(model);

   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "createModelCore2()", m_appParam.LogLevel, 2);
   
   return model;
}

//===========================================================================//
DecompConstraintSet * ATM_DecompApp::createModelCoreCount(){

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createModelCoreCount()", m_appParam.LogLevel, 2);

   int                 a, nRows;
   const int           nAtms      = m_instance.getNAtms();
   const int           nAtmsSteps = getNAtmsSteps();
   const int           nCols      = numCoreCols();
   if(m_appParam.UseTightModel)
      nRows = 3*nAtms +   nAtmsSteps;
   else
      nRows = 2*nAtms + 3*nAtmsSteps;
   int                 rowCnt     = 0;

   //---
   //--- A'' (core):
   //---    for a in A:
   //---       sum{d in D} v[a,d]  <= K[a]
   //---
   
   DecompConstraintSet * model = new DecompConstraintSet();
   CoinAssertHint(model, "Error: Out of Memory");

   //---
   //--- create new matrix
   //---
   model->M = new CoinPackedMatrix(false, 0.0, 0.0);
   CoinAssertHint(model->M, "Error: Out of Memory");
   model->M->setDimensions(0, nCols);
   model->reserve(nRows, nCols);

   //---
   //--- create model constraints
   //---
   for(a = 0; a < nAtms; a++){
      rowCnt += createConCount  (model, a);
      rowCnt += createConZtoX   (model, a);//OPT
      rowCnt += createConPickOne(model, a);//OPT
   }
   //THINK - add conZtoX to core and remove from relax?
   printf("nRows = %d, rowCnt = %d\n", nRows, rowCnt);
   assert(rowCnt == nRows);
   
   //---
   //--- create model variables/columns
   //---
   createModelColumns(model);

   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "createModelCoreCount()", m_appParam.LogLevel, 2);
   
   return model;
}

//===========================================================================//
DecompConstraintSet * 
ATM_DecompApp::createModelRelax1(const int a,
				 bool      includeCount){
   
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createModelRelax1()", m_appParam.LogLevel, 2);

   int                t, d, nRows, pairIndex;
   double             rhs, coefA, coefC;
   pair<int,int>      adP;
   vector<int>::const_iterator vi;
   
   const int           nSteps    = m_appParam.NumSteps;
   const int           nDates    = m_instance.getNDates();
   const vector<int> & pairsAD   = m_instance.getPairsAD();   

   const double      * a_ad      = m_instance.get_a_ad(); //dense storage
   const double      * b_ad      = m_instance.get_b_ad(); //dense storage
   const double      * c_ad      = m_instance.get_c_ad(); //dense storage
   const double      * d_ad      = m_instance.get_d_ad(); //dense storage
   const double      * e_ad      = m_instance.get_e_ad(); //dense storage
   const double      * w_ad      = m_instance.get_w_ad(); //dense storage

   const int           nCols     = numCoreCols();
   int                 nRowsMax  = nDates + 1 + (3*nSteps);
   if(includeCount)
      nRowsMax += nDates + 1;

   //---
   //---    for a in A, d in D:
   //---       f+[a,d] - f-[a,d] =  
   //---         a[a,d] sum{t in T} (t/n) x1[a,t] +
   //---         b[a,d] x2[a]                     +     
   //---         c[a,d] sum{t in T} (t/n) z[a,t]  + 
   //---         d[a,d] x3[a]                     +
   //---         e[a,d]
   //---       f-[a,d] <= w[a,d] * v[a,d]  (for count)
   //---    for a in A:
   //---       sum{t in T} x1[a,t] <= 1
   //---       sum{d in D} v[a,d]  <= K[a] (for count) [OPT]

   //Probably not... 
   //THINK: can't we just post-process this? z=x1*x2
   //---    for a in A, t in T:
   //---       z[a,t] <= x1[a,t],            
   //---       z[a,t] <= x2[a],              
   //---       z[a,t] >= x1[a,t] + x2[a] - 1.
   //---
   
   DecompConstraintSet * model = new DecompConstraintSet();
   CoinAssertHint(model, "Error: Out of Memory");

   model->M = new CoinPackedMatrix(false, 0.0, 0.0);
   CoinAssertHint(model->M, "Error: Out of Memory");
   model->M->setDimensions(0, nCols);
   model->M->reserve(nRowsMax, nCols);
   model->rowLB.reserve(nRowsMax);
   model->rowUB.reserve(nRowsMax);

   //---
   //---    for a in A, d in D:
   //---       f+[a,d] - f-[a,d] =  
   //---         a[a,d] sum{t in T} (t/n) x1[a,t] +
   //---         b[a,d] x2[a]                     +     
   //---         c[a,d] sum{t in T} (t/n) z[a,t]  + 
   //---         d[a,d] x3[a]                     +
   //---         e[a,d]
   //---       f-[a,d] <= w[a,d] * v[a,d]  (for count)
   //---
   nRows     = 0;
   pairIndex = 0;
   for(vi = pairsAD.begin(); vi != pairsAD.end(); vi++){
      adP = m_instance.getIndexADInv(*vi);
      if(a != adP.first){
	 pairIndex++;
	 continue;      
      }
      d = adP.second;

      CoinPackedVector row;
      string rowName = "demand_def(a_" 
	 + m_instance.getAtmName(a) + ",d_" 
	 + m_instance.getDateName(d) + ")";

      row.insert(colIndex_fp(pairIndex),  1.0);
      row.insert(colIndex_fm(pairIndex), -1.0);
      coefA = -a_ad[*vi] / nSteps;
      coefC = -c_ad[*vi] / nSteps;
      if(m_appParam.UseTightModel){
	 for(t = 0; t <= nSteps; t++){ 
	    row.insert(colIndex_x1(a,t), t * coefA);
	    row.insert(colIndex_z (a,t), t * coefC);
	 }
      }
      else{
	 for(t = 0; t < nSteps; t++){ 
	    row.insert(colIndex_x1(a,t), (t+1) * coefA);
	    row.insert(colIndex_z (a,t), (t+1) * coefC);
	 }
      }
      row.insert(colIndex_x2(a), -b_ad[*vi]);
      row.insert(colIndex_x3(a), -d_ad[*vi]);

      rhs = e_ad[*vi];
      model->M->appendRow(row);
      model->rowLB.push_back(rhs);
      model->rowUB.push_back(rhs);
      model->rowNames.push_back(rowName);      
      nRows++;
      
      CoinPackedVector rowLink;
      string rowNameLink = "linkv(a_" 
         + m_instance.getAtmName(a) + ",d_" 
         + m_instance.getDateName(d) + ")";	 
      rowLink.insert(colIndex_fm(pairIndex), 1.0);
      rowLink.insert(colIndex_v (pairIndex), -w_ad[*vi]);
      model->M->appendRow(rowLink);
      model->rowLB.push_back(-DecompInf);
      model->rowUB.push_back(0.0);
      model->rowNames.push_back(rowNameLink);	 
      nRows++;	 
   
      pairIndex++;
   }

   //---
   //---    for a in A:
   //---       sum{t in T} x1[a,t] <= 1
   //---       sum{d in D} v[a,d]  <= K[a] (for count)
   //---
   nRows += createConPickOne(model, a);
   if(includeCount)
      nRows += createConCount(model, a);
   
   //---
   //---    for a in A, t in T:
   //---       z[a,t] <= x1[a,t],            
   //---       z[a,t] <= x2[a],              
   //---       z[a,t] >= x1[a,t] + x2[a] - 1.
   //---
   nRows += createConZtoX(model, a);

   //---
   //--- create model columns
   //---
   createModelColumns(model, a);
   
   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "createModelRelax1()", m_appParam.LogLevel, 2);
   
   return model;
}

//===========================================================================//
DecompConstraintSet * ATM_DecompApp::createModelRelax2(const int d){
   
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createModelRelax2()", m_appParam.LogLevel, 2);
   
   int                a, t, nRows, pairIndex, nRowsMax;
   double             rhs, coefA, coefC;
   pair<int,int>      adP;
   vector<int>::const_iterator vi;
   
   const int           nSteps    = m_appParam.NumSteps;
   const int           nAtms     = m_instance.getNAtms();
   const vector<int> & pairsAD   = m_instance.getPairsAD();   

   const double      * a_ad      = m_instance.get_a_ad(); //dense storage
   const double      * b_ad      = m_instance.get_b_ad(); //dense storage
   const double      * c_ad      = m_instance.get_c_ad(); //dense storage
   const double      * d_ad      = m_instance.get_d_ad(); //dense storage
   const double      * e_ad      = m_instance.get_e_ad(); //dense storage
   const double      * w_ad      = m_instance.get_w_ad(); //dense storage
   const double      * B_d       = m_instance.get_B_d();

   const int           nCols     = numCoreCols();
   if(m_appParam.UseTightModel)
      nRowsMax  = 3*nAtms + 1 + getNAtmsSteps();
   else
      nRowsMax  = 2*nAtms + 1 + (3*getNAtmsSteps());

   //---
   //--- A'[d] for d in D (independent blocks)
   //---    for a in A
   //---       f+[a,d] - f-[a,d] =  
   //---         a[a,d] sum{t in T} (t/n) x1[a,t] +
   //---         b[a,d] x2[a]                     +     
   //---         c[a,d] sum{t in T} (t/n) z[a,t]  + 
   //---         d[a,d] x3[a]                     +
   //---         e[a,d]
   //---       f-[a,d] <= w[a,d] * v[a,d]
   //---    sum{a in A} (f+[a,d] - f-[a,d]) <= B[d]
   //---    for a in A, t in T:
   //---       z[a,t] <= x1[a,t],            
   //---       z[a,t] <= x2[a],              
   //---       z[a,t] >= x1[a,t] + x2[a] - 1.
   //---
   
   DecompConstraintSet * model = new DecompConstraintSet();
   CoinAssertHint(model, "Error: Out of Memory");

   model->M = new CoinPackedMatrix(false, 0.0, 0.0);
   CoinAssertHint(model->M, "Error: Out of Memory");
   model->M->setDimensions(0, nCols);
   model->reserve(nRowsMax, nCols);

   CoinPackedVector rowBudget;
   nRows     = 0;
   pairIndex = 0;
   for(vi = pairsAD.begin(); vi != pairsAD.end(); vi++){
      adP = m_instance.getIndexADInv(*vi);
      if(d != adP.second){
	 pairIndex++;
	 continue;      
      }
      a = adP.first;

      //---
      //--- sum{a in A} (f+[a,d] - f-[a,d]) <= B[d]
      //---
      rowBudget.insert(getColOffset_fp() + pairIndex,  1.0);
      rowBudget.insert(getColOffset_fm() + pairIndex, -1.0);


      //---
      //---       f+[a,d] - f-[a,d] =  
      //---         a[a,d] sum{t in T} (t/n) x1[a,t] +
      //---         b[a,d] x2[a]                     +     
      //---         c[a,d] sum{t in T} (t/n) z[a,t]  + 
      //---         d[a,d] x3[a]                     +
      //---         e[a,d]
      //---
      CoinPackedVector row;
      string rowName = "demand_def(a_" 
	 + m_instance.getAtmName(a) + ",d_" 
	 + m_instance.getDateName(d) + ")";
      row.insert(colIndex_fp(pairIndex),  1.0);
      row.insert(colIndex_fm(pairIndex), -1.0);
      coefA = -a_ad[*vi] / nSteps;
      coefC = -c_ad[*vi] / nSteps;
      if(m_appParam.UseTightModel){
	 for(t = 0; t <= nSteps; t++){ 
	    row.insert(colIndex_x1(a,t), t * coefA);
	    row.insert(colIndex_z (a,t), t * coefC);
	 }
      }
      else{
	 for(t = 0; t < nSteps; t++){ 
	    row.insert(colIndex_x1(a,t), (t+1) * coefA);
	    row.insert(colIndex_z (a,t), (t+1) * coefC);
	 }
      }
      row.insert(colIndex_x2(a), -b_ad[*vi]);
      row.insert(colIndex_x3(a), -d_ad[*vi]);

      rhs = e_ad[*vi];
      model->appendRow(row, rhs, rhs, rowName);
      nRows++;

      //---
      //---       f-[a,d] <= w[a,d] * v[a,d]      
      //---      
      CoinPackedVector rowLink;
      string rowNameLink = "linkv(a_" 
         + m_instance.getAtmName(a) + ",d_" 
         + m_instance.getDateName(d) + ")";	 
      rowLink.insert(colIndex_fm(pairIndex), 1.0);
      rowLink.insert(colIndex_v (pairIndex), -w_ad[*vi]);
      model->appendRow(rowLink, -DecompInf, 0.0, rowNameLink);
      nRows++;	 
   
      pairIndex++;
   }
   
   string rowNameBudget = "budget(d_" + m_instance.getDateName(d) + ")";
   model->appendRow(rowBudget, -DecompInf, B_d[d], rowNameBudget);
   nRows++;

  for(a = 0; a < nAtms; a++)
     nRows += createConZtoX(model, a);
   assert(nRows <= nRowsMax);
   
   //---
   //--- create model columns
   //---
   createModelColumns(model, -1, d);
   
   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "createModelRelax2()", m_appParam.LogLevel, 2);

   return model;
}


//===========================================================================//
DecompConstraintSet * ATM_DecompApp::createModelRelaxCount(){
   
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createModelRelaxCount()", m_appParam.LogLevel, 2);
   
   int                a, d, t, nRows, pairIndex, nRowsMax;
   double             rhs, coefA, coefC;
   pair<int,int>      adP;
   vector<int>::const_iterator vi;
   
   const int           nSteps     = m_appParam.NumSteps;
   const int           nAtms      = m_instance.getNAtms();
   const int           nDates     = m_instance.getNDates();
   const int           nPairs     = m_instance.getNPairs();
   const int           nAtmsSteps = getNAtmsSteps();
   const vector<int> & pairsAD    = m_instance.getPairsAD();   

   const double      * a_ad      = m_instance.get_a_ad(); //dense storage
   const double      * b_ad      = m_instance.get_b_ad(); //dense storage
   const double      * c_ad      = m_instance.get_c_ad(); //dense storage
   const double      * d_ad      = m_instance.get_d_ad(); //dense storage
   const double      * e_ad      = m_instance.get_e_ad(); //dense storage
   const double      * w_ad      = m_instance.get_w_ad(); //dense storage
   const double      * B_d       = m_instance.get_B_d();

   const int           nCols     = numCoreCols();
   if(m_appParam.UseTightModel)
      nRowsMax  = nDates + 2*nAtms +   nAtmsSteps + 2*nPairs;
   else
      nRowsMax  = nDates +   nAtms + 3*nAtmsSteps + 2*nPairs;

   //TODO:
   // try the nested idea - so master has ConZtoX and we price without that
   //   but we solve with gap the harder oracle with those constraints

   //---
   //--- A' (relax):
   //---    for d in D:
   //---       sum{a in A} (f+[a,d] - f-[a,d]) <= B[d]  <---- makes it hard
   //OPT
   //---    for a in A:
   //---       sum{t in T}  x1[a,t] <= 1
   //OPT
   //---    for a in A, t in T:
   //---       z[a,t] <= x1[a,t],            
   //---       z[a,t] <= x2[a],              
   //---       z[a,t] >= x1[a,t] + x2[a] - 1.      

   //---    for a in A, d in D:
   //---       f+[a,d] - f-[a,d] =  
   //---         a[a,d] sum{t in T} (t/n) x1[a,t] +
   //---         b[a,d] x2[a]                     +     
   //---         c[a,d] sum{t in T} (t/n) z[a,t]  + 
   //---         d[a,d] x3[a]                     +
   //---         e[a,d]
   //---       f-[a,d] <= w[a,d] * v[a,d]
   //---
   
   DecompConstraintSet * model = new DecompConstraintSet();
   CoinAssertHint(model, "Error: Out of Memory");

   model->M = new CoinPackedMatrix(false, 0.0, 0.0);
   CoinAssertHint(model->M, "Error: Out of Memory");
   model->M->setDimensions(0, nCols);
   model->reserve(nRowsMax, nCols);

   //---
   //--- create nDates empty rowBudget packed vectors
   //---
   vector<CoinPackedVector> rowBudget;
   for(d = 0; d < nDates; d++){
      CoinPackedVector row;
      rowBudget.push_back(row);
   }

   nRows     = 0;
   pairIndex = 0;
   for(vi = pairsAD.begin(); vi != pairsAD.end(); vi++){
      adP = m_instance.getIndexADInv(*vi);
      a = adP.first;
      d = adP.second;

      //---
      //--- sum{a in A} (f+[a,d] - f-[a,d]) <= B[d]
      //---
      rowBudget[d].insert(getColOffset_fp() + pairIndex,  1.0);
      rowBudget[d].insert(getColOffset_fm() + pairIndex, -1.0);


      //---
      //---       f+[a,d] - f-[a,d] =  
      //---         a[a,d] sum{t in T} (t/n) x1[a,t] +
      //---         b[a,d] x2[a]                     +     
      //---         c[a,d] sum{t in T} (t/n) z[a,t]  + 
      //---         d[a,d] x3[a]                     +
      //---         e[a,d]
      //---
      CoinPackedVector row;
      string rowName = "demand_def(a_" 
	 + m_instance.getAtmName(a) + ",d_" 
	 + m_instance.getDateName(d) + ")";
      row.insert(colIndex_fp(pairIndex),  1.0);
      row.insert(colIndex_fm(pairIndex), -1.0);
      coefA = -a_ad[*vi] / nSteps;
      coefC = -c_ad[*vi] / nSteps;
      if(m_appParam.UseTightModel){
	 for(t = 0; t <= nSteps; t++){ 
	    row.insert(colIndex_x1(a,t), t * coefA);
	    row.insert(colIndex_z (a,t), t * coefC);
	 }
      }
      else{
	 for(t = 0; t < nSteps; t++){ 
	    row.insert(colIndex_x1(a,t), (t+1) * coefA);
	    row.insert(colIndex_z (a,t), (t+1) * coefC);
	 }
      }
      row.insert(colIndex_x2(a), -b_ad[*vi]);
      row.insert(colIndex_x3(a), -d_ad[*vi]);
      
      rhs = e_ad[*vi];
      model->appendRow(row, rhs, rhs, rowName);
      nRows++;

      //---
      //---       f-[a,d] <= w[a,d] * v[a,d]      
      //---
      CoinPackedVector rowLink;
      string rowNameLink = "linkv(a_" 
         + m_instance.getAtmName(a) + ",d_" 
         + m_instance.getDateName(d) + ")";	 
      rowLink.insert(colIndex_fm(pairIndex), 1.0);
      rowLink.insert(colIndex_v (pairIndex), -w_ad[*vi]);
      model->appendRow(rowLink, -DecompInf, 0.0, rowNameLink);
      nRows++;	 

      pairIndex++;
   }

   for(d = 0; d < nDates; d++){   
      string rowNameBudget = "budget(d_" + m_instance.getDateName(d) + ")";
      model->appendRow(rowBudget[d], -DecompInf, B_d[d], rowNameBudget);
      nRows++;
   }

   //---
   //---    for a in A:
   //---       sum{t in T}  x1[a,t] <= 1
   //---    for a in A, t in T:
   //---       z[a,t] <= x1[a,t],            
   //---       z[a,t] <= x2[a],              
   //---       z[a,t] >= x1[a,t] + x2[a] - 1.      
   //---   
   //THINK: if you use this, shouldn't we postprocess z=x1*x2?
   //  putting those constaints in master are unneccessary... 
   //  in fact, why not just do that in block version too... 
   //for(a = 0; a < nAtms; a++){
      //nRows += createConZtoX(model, a);
      //nRows += createConPickOne(model, a);
   //}
   assert(nRows <= nRowsMax);
   
   //---
   //--- create model columns
   //---
   createModelColumns(model);
   
   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "createModelRelaxCount()", m_appParam.LogLevel, 2);

   return model;
}

//===========================================================================//
void ATM_DecompApp::createModels(){

   //---
   //--- This function does the work to create the different models 
   //---  that will be used. This memory is owned by the user. It will
   //---  be passed to the application interface and used by the algorithms.
   //---
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "APPcreateModel()", m_appParam.LogLevel, 2);
   
   //---
   //--- The original problem is in the form of a MINLP:
   //---    min  sum{a in A, d in D} | f[a,d] | 
   //---    s.t.
   //---    for a in A, d in D:
   //---       f[a,d] =  
   //---         a[a,d] x1[a]       +
   //---         b[a,d] x2[a]       +     
   //---         c[a,d] x1[a] x2[a] + 
   //---         d[a,d] x3[a]       +
   //---         e[a,d]
   //---    for d in D:
   //---       sum{a in A} f[a,d]      <= B[d]
   //---    for a in A:
   //---      |{d in D : f[a,d] < 0} | <= K[a] (notice strict ineq)
   //---    for a in A, d in D:
   //---      f[a,d] free
   //---    for a in A:
   //---      x1[a], x2[a] in [0,1], x3[a] >= 0
   //---

   //---
   //--- Discretization of continuous variable.
   //---      n     = number of steps
   //---      T     = {1, ..., n}
   //---      sum{t in T} (t/n) * x1[a,t]  = x1[a], for a in A
   //---      sum{t in T}         x1[a,t] <= 1    , for a in A
   //---   so, if n=10, x1[a] in {0,0.1,0.2,...,1.0}.
   //---

   //---
   //--- Linearization of product of a binary and a continuous.
   //---  For a in A, t in T:
   //---     z[a,t] = x1[a,t] * x2[a] 
   //---         <==>
   //---     z[a,t] >= 0 <= 1,                  
   //---     z[a,t] <= x1[a,t],            
   //---     z[a,t] <= x2[a],              
   //---     z[a,t] >= x1[a,t] + x2[a] - 1.
   //---

   //---
   //--- Linearization of the absolute value.
   //---  For a in A, d in D: 
   //---     f+[a,d], f-[a,d] >= 0
   //---     f[a,d] = f+[a,d] - f-[a,d]
   //---    |f[a,d]|= f+[a,d] + f-[a,d]
   //--- 

   //---
   //--- To model the count constraints.
   //---  For a in A:
   //---    |{d in D : f[a,d]            < 0}| <= K[a]
   //---      <==>
   //---    |{d in D : f+[a,d] - f-[a,d] < 0}| <= K[a]
   //---
   //---  At optimality, if f[a,d] != 0, then either
   //---    f+[a,d] > 0 and f-[a,d] = 0, or
   //---    f+[a,d] = 0 and f-[a,d] > 0.
   //---  So, to count the cases when f[a,d] < 0, we can restrict
   //---  attention to the cases where f-[a,d] > 0.
   //---
   //---  With some application specific stuff, 
   //---      we know that f-[a,d] <= w[a,d].
   //---
   //---  So, to count the number of cases we can use the following:
   //---      f-[a,d] >  0 ==> v[a,d] = 1
   //---      f-[a,d] <= 0 <== v[a,d] = 0
   //---        <==>
   //---      f-[a,d] <= w[a,d] * v[a,d]
   //---
   //---  and, then
   //---    |{d in D : f+[a,d] - f-[a,d] <  0}| <= K[a]
   //---        <==>
   //---    sum{d in D} v[a,d]                  <= K[a]
   //---

   //---
   //--- (Approximate) MILP Reformulation
   //---
   //---    min  sum{a in A, d in D} f+[a,d] + f-[a,d] 
   //---    s.t.
   //---    for a in A, d in D:
   //---       f+[a,d] - f-[a,d] =  
   //---         a[a,d] sum{t in T} (t/n) x1[a,t] +
   //---         b[a,d] x2[a]                     +     
   //---         c[a,d] sum{t in T} (t/n) z[a,t]  + 
   //---         d[a,d] x3[a]                     +
   //---         e[a,d]
   //---       f-[a,d] <= w[a,d] * v[a,d]
   //---    for a in A:
   //---       sum{t in T}  x1[a,t] <= 1
   //---    for a in A, t in T:
   //---       z[a,t] <= x1[a,t],            
   //---       z[a,t] <= x2[a],              
   //---       z[a,t] >= x1[a,t] + x2[a] - 1.
   //---    for d in D:
   //---       sum{a in A} (f+[a,d] - f-[a,d]) <= B[d]
   //---    for a in A:
   //---       sum{d in D} v[a,d]              <= K[a]
   //---    for a in A, d in D:
   //---       f+[a,d], f-[a,d] >= 0
   //---                f-[a,d] <= w[a,d]
   //---       v[a,d] in {0,1}
   //---    for a in A:
   //---       x2[a], x3[a] in [0,1]
   //---    for a in A, t in T
   //---       x1[a,t] in {0,1}, z[a,t] in [0,1]
   //---

   //---
   //--- UPDATE (April 2010). 
   //---
   //--- A tighter formulation of the 
   //---   linearization of product of a binary and a continuous 
   //---   is possible due to the constraint: sum{t in T}  x1[a,t] <= 1
   //---
   //---  For a in A, t in T:
   //---     z[a,t] = x1[a,t] * x2[a] 
   //---         <==>
   //--- OLD:
   //---  for a in A:
   //---     sum{t in T}  x1[a,t] <= 1 (where, T = 1..n)
   //---  for a in A, t in T:
   //---     z[a,t] >= 0 <= 1,                  
   //---     z[a,t] <= x1[a,t],            
   //---     z[a,t] <= x2[a],              
   //---     z[a,t] >= x1[a,t] + x2[a] - 1.
   //--- NEW: 
   //---  for a in A:
   //---     sum{t in T}  x1[a,t] = 1 (where, T = 0..n)
   //---     sum{t in T}   z[a,t] = x2[a]
   //---  for a in A, t in T:
   //---     z[a,t] >= 0 <= 1,                  
   //---     z[a,t] <= x1[a,t].
   //---


   //---
   //--- Columns
   //---   x1[a,t] (binary)
   //---   z [a,t]
   //---   f+[a,d] 
   //---   f-[a,d]
   //---   x2[a]
   //---   x3[a]
   //----  v[a,d] (binary)
   //---
   int i, a;
   int nAtms   = m_instance.getNAtms();
   int numCols = numCoreCols();

   //---
   //--- Construct the objective function.
   //---    Coefficients of f+ and f- are 1.0, the rest are 0.
   //---
   m_objective = new double[numCols];
   if(!m_objective)
      throw UtilExceptionMemory("createModels", "MMKP_DecompApp");
   UtilFillN(m_objective, numCols, 0.0);
   for(i = getColOffset_fp(); i < getColOffset_x2(); i++)
      m_objective[i] = 1.0;
   setModelObjective(m_objective);
   
   
   //---
   //--- A'' (core):
   //---    for d in D:
   //---       sum{a in A} (f+[a,d] - f-[a,d]) <= B[d]
   //---    for a in A: <--- option ( core or relax )
   //---       sum{d in D} v[a,d] <= K[a]
   //---
   //--- A'[a] for a in A (independent blocks)
   //---    for d in D:
   //---       f+[a,d] - f-[a,d] =  
   //---         a[a,d] sum{t in T} (t/n) x1[a,t] +
   //---         b[a,d] x2[a]                     +     
   //---         c[a,d] sum{t in T} (t/n) z[a,t]  + 
   //---         d[a,d] x3[a]                     +
   //---         e[a,d]
   //---       f-[a,d] <= w[a,d] * v[a,d]
   //---    sum{t in T}  x1[a,t] <= 1
   //---    for t in T:
   //---       z[a,t] <= x1[a,t],            
   //---       z[a,t] <= x2[a],              
   //---       z[a,t] >= x1[a,t] + x2[a] - 1.
   //---    sum{d in D} v[a,d] <= K[a] <--- option ( core or relax )
   //---
   if(m_appParam.ModelNameCore == "BUDGET"){
      DecompConstraintSet * modelCore = createModelCore1(false);
      m_models.push_back(modelCore);
      setModelCore(modelCore, "BUDGET");
   }
   if(m_appParam.ModelNameCore == "BUDGET_COUNT"){
      DecompConstraintSet * modelCore = createModelCore1();
      m_models.push_back(modelCore);
      setModelCore(modelCore, "BUDGET_COUNT");
   }   
   if(m_appParam.ModelNameRelax == "CASH"){
      for(a = 0; a < nAtms; a++){
	 DecompConstraintSet * modelRelax = createModelRelax1(a, false);
         m_models.push_back(modelRelax);
         setModelRelax(modelRelax, "CASH" + UtilIntToStr(a), a);
      }
   }
   if(m_appParam.ModelNameRelax == "CASH_COUNT"){
      for(a = 0; a < nAtms; a++){
	 DecompConstraintSet * modelRelax = createModelRelax1(a);
         m_models.push_back(modelRelax);
         setModelRelax(modelRelax, "CASH_COUNT" + UtilIntToStr(a), a);
      }
   }
   if(m_appParam.ModelNameRelaxNest == "CASH_COUNT"){
      for(a = 0; a < nAtms; a++){
	 DecompConstraintSet * modelRelax = createModelRelax1(a);
         m_models.push_back(modelRelax);
         setModelRelaxNest(modelRelax, "CASH_COUNT" + UtilIntToStr(a), a);
      }
   }


   //TODO: solve this A' with gap as relaxNest - but
   //  we are doing blocks - so find a column then break it out
   //  tell framework to do that? or do as user? show how we can 
   //  get stuff back from framework to setup and solve... 
   /*{
      //---
      //--- Version MODEL_MASTER_COUNT
      //---    is relaxation too hard?
      //---
      //--- A'' (core):
      //---    for a in A:
      //---       sum{d in D} v[a,d] <= K[a]
      //---
      //--- A' (relax):
      //---    for d in D:
      //---       sum{a in A} (f+[a,d] - f-[a,d]) <= B[d]
      //---    for a in A:
      //---       sum{t in T}  x1[a,t] <= 1
      //---    for a in A, t in T:
      //---       z[a,t] <= x1[a,t],            
      //---       z[a,t] <= x2[a],              
      //---       z[a,t] >= x1[a,t] + x2[a] - 1.      
      //---    for a in A, d in D:
      //---       f+[a,d] - f-[a,d] =  
      //---         a[a,d] sum{t in T} (t/n) x1[a,t] +
      //---         b[a,d] x2[a]                     +     
      //---         c[a,d] sum{t in T} (t/n) z[a,t]  + 
      //---         d[a,d] x3[a]                     +
      //---         e[a,d]
      //---       f-[a,d] <= w[a,d] * v[a,d]
      //---
     vector< DecompConstraintSet* > modelRelaxV;   
     DecompConstraintSet * modelCoreCount  = createModelCoreCount();
     DecompConstraintSet * modelRelaxCount = createModelRelaxCount();
     modelRelaxV.push_back(modelRelaxCount);   
     modelCore.insert (make_pair(MODEL_COUNT, modelCoreCount));
     modelRelax.insert(make_pair(MODEL_COUNT, modelRelaxV));
     }*/
   
   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "APPcreateModel()", m_appParam.LogLevel, 2);
}

//===========================================================================//
bool ATM_DecompApp::APPisUserFeasible(const double * x,
				      const int      nCols,
				      const double   tolZero){
   //---
   //--- sanity check - is this solution feasible?
   //---   since we provide a full matrix, DECOMP will also check 
   //---   that the algebra gives a feasible point
   //---

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "APPisUserFeasible()", m_appParam.LogLevel, 2);
   

   double lhs, rhs, coeff;
   int    a, d, t;
   int    pairIndex = 0;
   bool   isFeas    = true;

   const int           nSteps    = m_appParam.NumSteps;
   const int           nAtms     = m_instance.getNAtms();
   const vector<int> & pairsAD   = m_instance.getPairsAD();   
   const double      * K_a       = m_instance.get_K_a();
   const double      * a_ad      = m_instance.get_a_ad(); //dense storage
   const double      * b_ad      = m_instance.get_b_ad(); //dense storage
   const double      * c_ad      = m_instance.get_c_ad(); //dense storage
   const double      * d_ad      = m_instance.get_d_ad(); //dense storage
   const double      * e_ad      = m_instance.get_e_ad(); //dense storage
   double            * count     = new double[nAtms];
    
   //---
   //--- is the flow variable matching x,z variables
   //---    for a in A, d in D:
   //---       f+[a,d] - f-[a,d] =  
   //---         a[a,d] sum{t in T} (t/n) x1[a,t] +
   //---         b[a,d] x2[a]                     +     
   //---         c[a,d] sum{t in T} (t/n) z[a,t]  + 
   //---         d[a,d] x3[a]                     +
   //---         e[a,d]
   //---
   pair<int,int>                 adP;
   vector<int>::const_iterator   vi;
   double                        actViol = 0.0;
   double                        relViol = 0.0;
   for(vi = pairsAD.begin(); vi != pairsAD.end(); vi++){
      adP = m_instance.getIndexADInv(*vi);
      a   = adP.first;

      lhs  = x[getColOffset_fp() + pairIndex];
      lhs -= x[getColOffset_fm() + pairIndex];
      
      rhs  = e_ad[*vi]
	 + b_ad[*vi] * x[colIndex_x2(a)]
	 + d_ad[*vi] * x[colIndex_x3(a)];
      if(m_appParam.UseTightModel){
	 for(t = 0; t <= nSteps; t++){
	    coeff = (double)(t) / (double)nSteps;
	    rhs += a_ad[*vi] * coeff * x[colIndex_x1(a,t)];
	    rhs += c_ad[*vi] * coeff * x[colIndex_z (a,t)];
	 }
      }
      else{
	 for(t = 0; t < nSteps; t++){
	    coeff = (double)(t+1) / (double)nSteps;
	    rhs += a_ad[*vi] * coeff * x[colIndex_x1(a,t)];
	    rhs += c_ad[*vi] * coeff * x[colIndex_z (a,t)];
	 }
      }
      actViol = fabs(lhs-rhs);
      if(UtilIsZero(lhs,1.0e-3))
	relViol = actViol;
      else
	relViol = actViol / std::fabs(lhs);
      if(relViol > 0.05){
	printf("NOT FEASIBLE lhs=%12.10f, rhs=%12.10f\n", lhs, rhs);
	 isFeas = false;
      }      
      pairIndex++;
   }


   //---
   //--- did the indicators work correctly?
   //---    f+[a,d] - f-[a,d] < 0 ==> v[a,d] = 1
   //---
   pairIndex = 0;
   for(vi = pairsAD.begin(); vi != pairsAD.end(); vi++){
      adP  = m_instance.getIndexADInv(*vi);
      a    = adP.first;
      lhs  = x[getColOffset_fp() + pairIndex];
      lhs -= x[getColOffset_fm() + pairIndex];
      //if(lhs < -DecompEpsilon){
      if(lhs < -0.01){
	 if(fabs(x[getColOffset_v() + pairIndex] - 1.0) > tolZero){
	    printf("NOT FEASIBLE fp=%10.5f fm=%10.5f f=%10.5f < 0, but v=%g not 1\n",
		   x[getColOffset_fp() + pairIndex],
		   x[getColOffset_fm() + pairIndex],
		   lhs, x[getColOffset_v() + pairIndex]);	 
	    isFeas = false;
	 }
      }
      pairIndex++;
   }
   
   //---
   //--- did the surrogate z work correctly?
   //---   z[a,t] = sum{t in T} x1[a,t] * x2[a]
   //---

   //---
   //--- is budget satisfied?
   //---    for d in D:
   //---       sum{a in A} (f+[a,d] - f-[a,d]) <= B[d]
   //---

   //---
   //--- is the count constraint satisifed
   //---    for a in A:
   //---       sum{d in D} v[a,d]              <= K[a]
   //---
   pairIndex = 0;
   UtilFillN(count, nAtms, 0.0);
   for(vi = pairsAD.begin(); vi != pairsAD.end(); vi++){
      adP  = m_instance.getIndexADInv(*vi);
      a    = adP.first;
      d    = adP.second;
      count[a] += x[getColOffset_v() + pairIndex];
      pairIndex++;
   }
   for(a = 0; a < nAtms; a++){
#if 0
      printf("COUNT[a=%3d->%10s]: %5g <= K=%5g\n",
	     a, 
	     m_instance.getAtmName(a).c_str(),
	     count[a], 
	     K_a[a]);
#endif	     
      if(count[a] > (K_a[a] + 0.01)){
#if 0
	 printf("NOT FEASIBLE a:%d count=%g K=%g\n", a, count[a], K_a[a]);
#endif
	 isFeas = false;
      }
   }

#if 0
   printf("IsUserFeas = %d\n", isFeas);
#endif

   //---
   //--- free local memory
   //---
   UTIL_DELARR(count);

   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "APPisUserFeasible()", m_appParam.LogLevel, 2);
   
   return isFeas;
}
