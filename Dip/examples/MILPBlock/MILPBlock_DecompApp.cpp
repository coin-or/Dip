//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, Lehigh University                                 //
//                                                                           //
// Copyright (C) 2002-2009, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

//===========================================================================//
#include "DecompAlgo.h"
#include "MILPBlock_DecompApp.h"

//===========================================================================//
void MILPBlock_DecompApp::initializeApp(UtilParameters & utilParam)  {
   
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "initializeApp()", m_appParam.LogLevel, 2);
   
   //---
   //--- get application parameters
   //---   
   m_appParam.getSettings(utilParam);   
   if(m_appParam.LogLevel >= 1)
      m_appParam.dumpSettings();

   //---
   //--- read MILPBlock instance (mps format)
   //---
   string fileName = m_appParam.DataDir 
      + UtilDirSlash() + m_appParam.Instance;   

   m_mpsIO.messageHandler()->setLogLevel(m_param.LogLpLevel);

   int rstatus = m_mpsIO.readMps(fileName.c_str());   
   if(rstatus < 0){
      cerr << "Error: Filename = " << fileName << " failed to open." << endl;
      throw UtilException("I/O Error.", "initalizeApp", "MILPBlock_DecompApp");
   }
   if(m_appParam.LogLevel >= 2)
      (*m_osLog) << "Objective Offset = " 
                 << UtilDblToStr(m_mpsIO.objectiveOffset()) << endl;

   //---
   //--- set best known lb/ub
   //---
   double offset = m_mpsIO.objectiveOffset();
   setBestKnownLB(m_appParam.BestKnownLB + offset);
   setBestKnownUB(m_appParam.BestKnownUB + offset);

   //---
   //--- read block file
   //---
   readBlockFile();

   //---
   //--- create models
   //---
   createModels();

   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "initializeApp()", m_appParam.LogLevel, 2);
}

//===========================================================================//
void MILPBlock_DecompApp::readBlockFile(){

   ifstream is;
   string   fileName = m_appParam.DataDir 
      + UtilDirSlash() + m_appParam.BlockFile;
   
   //---
   //--- open file streams
   //---
   UtilOpenFile(is, fileName.c_str());

   int i, rowId, numRowsInBlock, blockId;
   if(m_appParam.LogLevel >= 1)
      (*m_osLog) << "Reading " << fileName << endl;
   if(m_appParam.BlockFileFormat == "List" ||
      m_appParam.BlockFileFormat == "LIST"){

      //---
      //--- The block file defines those rows in each block.
      //---   <block id>  <num rows in block>
      //---     <row ids...>
      //---   <block id>  <num rows in block>
      //---     <row ids...>
      //---
      
      //---
      //--- TODO: also allow row names (instead of row ids)
      //---      
      while(!is.eof()){
	 is >> blockId;
	 is >> numRowsInBlock;
	 if(is.eof()) break;
	 vector<int> rowsInBlock;
	 for(i = 0; i < numRowsInBlock; i++){
	    is >> rowId;
	    rowsInBlock.push_back(rowId);
	 }
	 m_blocks.insert(make_pair(blockId, rowsInBlock));
	 if(is.eof()) break;      
      }
   }
   else if(m_appParam.BlockFileFormat == "Pair" ||
	   m_appParam.BlockFileFormat == "PAIR"){
      //---
      //--- <block id> <row id> 
      //---  ...
      //---
      int blockIdCurr = 0;
      is >> blockId;
      while(!is.eof()){
	 vector<int> rowsInBlock;
	 while(blockId == blockIdCurr && !is.eof()){
	    is >> rowId;
	    rowsInBlock.push_back(rowId);
	    is >> blockId;
	 }
	 m_blocks.insert(make_pair(blockIdCurr, rowsInBlock));	 
	 blockIdCurr = blockId;
	 if(is.eof()) break;
      }      
   } else{
      assert(0);
   }

   if(m_appParam.LogLevel >= 3){
      map<int, vector<int> >::iterator mit;
      vector<int>           ::iterator vit;
      for(mit = m_blocks.begin(); mit != m_blocks.end(); mit++){
         (*m_osLog) << "Block " << (*mit).first << " : ";
         for(vit = (*mit).second.begin(); vit != (*mit).second.end(); vit++)
            (*m_osLog) << (*vit) << " ";
         (*m_osLog) << endl;
      }
   }
   
   is.close();
}

/*
//===========================================================================//
void MILPBlock_DecompApp::readInitSolutionFile(){

   ifstream is;
   string   fileName = m_appParam.DataDir 
      + UtilDirSlash() + m_appParam.InitSolutionFile;
   ////////////STOP
   
   //---
   //--- open file streams
   //---
   UtilOpenFile(is, fileName.c_str());

   int i, rowId, numRowsInBlock, blockId;
   if(m_appParam.LogLevel >= 1)
      (*m_osLog) << "Reading " << fileName << endl;
   if(m_appParam.BlockFileFormat == "List" ||
      m_appParam.BlockFileFormat == "LIST"){

      //---
      //--- The block file defines those rows in each block.
      //---   <block id>  <num rows in block>
      //---     <row ids...>
      //---   <block id>  <num rows in block>
      //---     <row ids...>
      //---
      
      //---
      //--- TODO: also allow row names (instead of row ids)
      //---      
      while(!is.eof()){
	 is >> blockId;
	 is >> numRowsInBlock;
	 if(is.eof()) break;
	 vector<int> rowsInBlock;
	 for(i = 0; i < numRowsInBlock; i++){
	    is >> rowId;
	    rowsInBlock.push_back(rowId);
	 }
	 m_blocks.insert(make_pair(blockId, rowsInBlock));
	 if(is.eof()) break;      
      }
   }
   else if(m_appParam.BlockFileFormat == "Pair" ||
	   m_appParam.BlockFileFormat == "PAIR"){
      //---
      //--- <block id> <row id> 
      //---  ...
      //---
      int blockIdCurr = 0;
      is >> blockId;
      while(!is.eof()){
	 vector<int> rowsInBlock;
	 while(blockId == blockIdCurr && !is.eof()){
	    is >> rowId;
	    rowsInBlock.push_back(rowId);
	    is >> blockId;
	 }
	 m_blocks.insert(make_pair(blockIdCurr, rowsInBlock));	 
	 blockIdCurr = blockId;
	 if(is.eof()) break;
      }      
   } else{
      assert(0);
   }

   if(m_appParam.LogLevel >= 3){
      map<int, vector<int> >::iterator mit;
      vector<int>           ::iterator vit;
      for(mit = m_blocks.begin(); mit != m_blocks.end(); mit++){
         (*m_osLog) << "Block " << (*mit).first << " : ";
         for(vit = (*mit).second.begin(); vit != (*mit).second.end(); vit++)
            (*m_osLog) << (*vit) << " ";
         (*m_osLog) << endl;
      }
   }
   
   is.close();
   }*/

//===========================================================================//
void 
MILPBlock_DecompApp::findActiveColumns(const vector<int> & rowsPart,
                                       set<int>          & activeColsSet){

   const CoinPackedMatrix * M    = m_mpsIO.getMatrixByRow();
   const int              * ind  = M->getIndices();
   const int              * beg  = M->getVectorStarts();
   const int              * len  = M->getVectorLengths();
   const int              * indR = NULL;
   
   //---
   //--- which columns are present in this part's rows
   //---
   int k, r;
   vector<int>::const_iterator it;
   for(it = rowsPart.begin(); it != rowsPart.end(); it++){
      r    = *it;
      indR = ind + beg[r];
      for(k = 0; k < len[r]; k++){
         activeColsSet.insert(indR[k]);
      }
   }
}

//===========================================================================//
void
MILPBlock_DecompApp::createModelMasterOnlys(vector<int> & masterOnlyCols){

   int            nBlocks     = static_cast<int>(m_blocks.size());
   const int      nCols       = m_mpsIO.getNumCols();
   const double * colLB       = m_mpsIO.getColLower();
   const double * colUB       = m_mpsIO.getColUpper();
   const char   * integerVars = m_mpsIO.integerColumns();
   int            nMasterOnlyCols =
      static_cast<int>(masterOnlyCols.size());

   if(m_appParam.LogLevel >= 1){
      (*m_osLog) << "nCols           = " << nCols << endl;
      (*m_osLog) << "nMasterOnlyCols = " << nMasterOnlyCols << endl;
   }

   if(nMasterOnlyCols == 0)
      return;


   int i, j;
   vector<int>::iterator vit;
   for(vit = masterOnlyCols.begin(); vit != masterOnlyCols.end(); vit++){
      i = *vit;
            
      DecompConstraintSet * model = new DecompConstraintSet();
      model->M = new CoinPackedMatrix(false, 0.0, 0.0);
      if(!model->M)
         throw UtilExceptionMemory("createModels", "MILPBlock_DecompApp");
      //TODO: Ouch - memory wise - this creates a block per
      //  master-only var of size nCols! that is dense space??
      //similar issue for SSR... when size of nCols is big
      //need to treat these special and not create the explicit
      //block or constraint system at all
      model->M->setDimensions(0, nCols);
      model->reserve(1, nCols);
      
      //---
      //--- fix all non-active columns to 0, this is everything
      //---  except the master-only columns
      //---
      UtilFillN(model->colLB, nCols, 0.0);
      UtilFillN(model->colUB, nCols, 0.0);
      
      //---
      //--- set the master-only vars but watch for unbounded
      //---
      model->colLB[i] = colLB[*vit];
      model->colUB[i] = colUB[*vit];
      if(m_appParam.ColumnUB <  1.0e15)
	 if(colUB[i] >  1.0e15)
	    model->colUB[i] = m_appParam.ColumnUB;
      if(m_appParam.ColumnLB > -1.0e15)
	 if(colLB[i] < -1.0e15)
	    model->colLB[i] = m_appParam.ColumnLB;

      //---
      //--- the master-only columns are the only active columns
      //---
      model->activeColumns.push_back(i);
      
      //---
      //--- set the indices of the integer variables of modelRelax
      //---  also set the column names, if they exist
      //---
      for(j = 0; j < nCols; j++){
         const char * colName = m_mpsIO.columnName(j);
         if(colName)
            model->colNames.push_back(colName);
         if(integerVars && integerVars[j]){
            model->integerVars.push_back(j);            
         }
      }   

      //---
      //--- to avoid any issues with an empty constraint matrix
      //---   just add one column bound as an explicity row
      //---
      CoinPackedVector row;
      string           rowName  = "fake_row";
      int              colIndex = i;
      if(m_appParam.LogLevel >= 2){
         (*m_osLog) << "Masteronly colindex = " << colIndex << endl;
         (*m_osLog) << "  colUB now rowUB   = " << model->colUB[colIndex] 
                    << endl;
      }
      row.insert(colIndex, 1.0);
      model->appendRow(row, -DecompInf, model->colUB[colIndex], rowName);

      if(m_appParam.LogLevel >= 2){
         (*m_osLog) << "model numcols= " << model->getNumCols() << endl;
         (*m_osLog) << "model numrows= " << model->getNumRows() << endl;
      }

      m_modelR.insert(make_pair(nBlocks, model));
      setModelRelax(model, 
                    "master_only" + UtilIntToStr(i), nBlocks);
      nBlocks++;
   }

   return;   
}

//===========================================================================//
void
MILPBlock_DecompApp::createModelMasterOnlys2(vector<int> & masterOnlyCols){

   int            nBlocks     = static_cast<int>(m_blocks.size());
   const int      nCols       = m_mpsIO.getNumCols();
   const double * colLB       = m_mpsIO.getColLower();
   const double * colUB       = m_mpsIO.getColUpper();
   const char   * integerVars = m_mpsIO.integerColumns();
   int            nMasterOnlyCols =
      static_cast<int>(masterOnlyCols.size());

   if(m_appParam.LogLevel >= 1){
      (*m_osLog) << "nCols           = " << nCols << endl;
      (*m_osLog) << "nMasterOnlyCols = " << nMasterOnlyCols << endl;
   }

   if(nMasterOnlyCols == 0)
      return;


   int i, j;
   vector<int>::iterator vit;
   for(vit = masterOnlyCols.begin(); vit != masterOnlyCols.end(); vit++){
      i = *vit;

      //THINK:
      //  what-if master-only var is integer and bound is not at integer
            
      DecompConstraintSet * model = new DecompConstraintSet();
      model->m_masterOnly      = true;
      model->m_masterOnlyIndex = i;
      model->m_masterOnlyLB    = colLB[i];
      model->m_masterOnlyUB    = colUB[i];
      if(m_appParam.ColumnUB <  1.0e15)
	 if(colUB[i] >  1.0e15)
	    model->m_masterOnlyUB = m_appParam.ColumnUB;
      if(m_appParam.ColumnLB > -1.0e15)
	 if(colLB[i] < -1.0e15)
	    model->m_masterOnlyLB = m_appParam.ColumnLB;

      m_modelR.insert(make_pair(nBlocks, model));
      setModelRelax(model, 
                    "master_only" + UtilIntToStr(i), nBlocks);
      nBlocks++;
   }

   return;   
}

//===========================================================================//
DecompConstraintSet * 
MILPBlock_DecompApp::createModelMasterOnly(vector<int> & masterOnlyCols){
					  
   const int      nCols       = m_mpsIO.getNumCols();
   const double * colLB       = m_mpsIO.getColLower();
   const double * colUB       = m_mpsIO.getColUpper();
   const char   * integerVars = m_mpsIO.integerColumns();
   int            nMasterOnlyCols =
      static_cast<int>(masterOnlyCols.size());

   if(m_appParam.LogLevel >= 1){
      (*m_osLog) << "nCols           = " << nCols << endl;
      (*m_osLog) << "nMasterOnlyCols = " << nMasterOnlyCols << endl;
   }

   if(nMasterOnlyCols == 0)
      return NULL;
   
   DecompConstraintSet * model = new DecompConstraintSet();
   model->M = new CoinPackedMatrix(false, 0.0, 0.0);
   if(!model->M)
      throw UtilExceptionMemory("createModels", "MILPBlock_DecompApp");
   model->M->setDimensions(0, nCols);
   model->reserve(1, nCols);

   //---
   //--- fix all non-active columns to 0, this is everything
   //---  except the master-only columns
   //---
   UtilFillN(model->colLB, nCols, 0.0);
   UtilFillN(model->colUB, nCols, 0.0);
   
   //---
   //--- set the master-only vars but watch for unbounded
   //---
   int i;
   vector<int>::iterator vit;
   for(vit = masterOnlyCols.begin(); vit != masterOnlyCols.end(); vit++){
      i = *vit;
      model->colLB[i] = colLB[*vit];
      model->colUB[i] = colUB[*vit];
      if(m_appParam.ColumnUB <  1.0e15)
	 if(colUB[i] >  1.0e15)
	    model->colUB[i] = m_appParam.ColumnUB;
      if(m_appParam.ColumnLB > -1.0e15)
	 if(colLB[i] < -1.0e15)
	    model->colLB[i] = m_appParam.ColumnLB;
   }

   //---
   //--- the master-only columns are the only active columns
   //---
   model->activeColumns.insert(model->activeColumns.end(),
			       masterOnlyCols.begin(),
			       masterOnlyCols.end());

   
   //---
   //--- set the indices of the integer variables of modelRelax
   //---  also set the column names, if they exist
   //---
   for(i = 0; i < nCols; i++){
      const char * colName = m_mpsIO.columnName(i);
      if(colName)
         model->colNames.push_back(colName);
      if(integerVars && integerVars[i]){
         model->integerVars.push_back(i);            
      }
   }   

   //---
   //--- to avoid any issues with an empty constraint matrix
   //---   just add one column bound as an explicity row
   //---
   CoinPackedVector row;
   string           rowName  = "fake_row";
   int              colIndex = masterOnlyCols[0];
   if(m_appParam.LogLevel >= 2){
      (*m_osLog) << "Masteronly colindex = " << colIndex << endl;
      (*m_osLog) << "  colUB now rowUB   = " << model->colUB[colIndex] << endl;
   }
   row.insert(colIndex, 1.0);
   model->appendRow(row, -DecompInf, model->colUB[colIndex], rowName);

   if(m_appParam.LogLevel >= 2){
      (*m_osLog) << "model numcols= " << model->getNumCols() << endl;
      (*m_osLog) << "model numrows= " << model->getNumRows() << endl;
   }

   return model;   
}

//===========================================================================//
DecompConstraintSet * 
MILPBlock_DecompApp::createModelPart(const int   nRowsPart,
                                     const int * rowsPart){

   const int      nCols       = m_mpsIO.getNumCols();
   const double * rowLB       = m_mpsIO.getRowLower();
   const double * rowUB       = m_mpsIO.getRowUpper();
   const double * colLB       = m_mpsIO.getColLower();
   const double * colUB       = m_mpsIO.getColUpper();
   const char   * integerVars = m_mpsIO.integerColumns();
   
   DecompConstraintSet * model = new DecompConstraintSet();
   model->M = new CoinPackedMatrix(false, 0.0, 0.0);
   if(!model->M)
      throw UtilExceptionMemory("createModels", "MILPBlock_DecompApp");
   model->reserve(nRowsPart, nCols);
   model->M->submatrixOf(*m_mpsIO.getMatrixByRow(), nRowsPart, rowsPart);   
   
   //---
   //--- set the row upper and lower bounds
   //--- set the col upper and lower bounds
   //---
   int i, r;
   for(i = 0; i < nRowsPart; i++){
      r = rowsPart[i];
      const char * rowName = m_mpsIO.rowName(r);
      if(rowName)
         model->rowNames.push_back(rowName);
      model->rowLB.push_back(rowLB[r]);
      model->rowUB.push_back(rowUB[r]);
   }
   copy(colLB, colLB + nCols, back_inserter( model->colLB)  );
   copy(colUB, colUB + nCols, back_inserter( model->colUB)  );

   //---
   //--- big fat hack... we don't deal with dual rays yet,
   //---  so, we assume subproblems are bounded
   //---
   //--- NOTE: might also need to tighten LBs
   //---
   //--- Too small - ATM infeasible!
   //--- Too big   - round off issues with big coeffs in 
   //---             master-only vars
   //---
   //--- TODO: need extreme rays or bounded subproblems from user
   //---
   if(m_appParam.ColumnUB < 1.0e15){
      for(i = 0; i < nCols; i++){
	 if(colUB[i] > 1.0e15){
	    model->colUB[i] = m_appParam.ColumnUB;
	 }
      }
   }
   if(m_appParam.ColumnLB > -1.0e15){
      for(i = 0; i < nCols; i++){
	 if(colLB[i] < -1.0e15){
	    model->colLB[i] = m_appParam.ColumnLB;
	 }
      }
   }

   //---
   //--- set the indices of the integer variables of modelRelax
   //---  also set the column names, if they exist
   //---
   for(i = 0; i < nCols; i++){
      const char * colName = m_mpsIO.columnName(i);
      if(colName)
         model->colNames.push_back(colName);
      if(integerVars && integerVars[i]){
         model->integerVars.push_back(i);            
      }
   }   
   return model;   
}


//===========================================================================//
void MILPBlock_DecompApp::createModels(){

   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "createModels()", m_appParam.LogLevel, 2);
   
   //---
   //--- how many rows to put into relaxation
   //---
   int            i, nRowsRelax, nRowsCore;
   const int      nRows       = m_mpsIO.getNumRows();
   const int      nCols       = m_mpsIO.getNumCols();
   int            nBlocks     = static_cast<int>(m_blocks.size());

   map<int, vector<int> >::iterator mit;
   nRowsRelax = 0;
   for(mit = m_blocks.begin(); mit != m_blocks.end(); mit++)
      nRowsRelax += static_cast<int>((*mit).second.size());   
   nRowsCore = nRows - nRowsRelax; 

   UTIL_MSG(m_appParam.LogLevel, 2,
            (*m_osLog) << "Instance    = " << m_appParam.Instance << endl;
            (*m_osLog) << " nRows      = " << nRows     << endl;
            (*m_osLog) << " nCols      = " << nCols     << endl;
            (*m_osLog) << " nBlocks    = " << nBlocks   << endl;
            (*m_osLog) << " nRowsCore  = " << nRowsCore << endl;
            (*m_osLog) << " nRowsRelax = " << nRowsRelax 
            << " [ " << 100*nRowsRelax/nRows << " % ]" << endl;
            );

   //---
   //--- setup markers for core and relax rows
   //---   
   int * rowsMarker = new int[nRows];
   int * rowsCore   = new int[nRowsCore];   
   UtilFillN(rowsMarker, nRows, -1);//-1 will mark core rows

   for(mit = m_blocks.begin(); mit != m_blocks.end(); mit++){
      vector<int> & rowsRelax = (*mit).second;
      vector<int>::iterator vit;
      for(vit = rowsRelax.begin(); vit != rowsRelax.end(); vit++)
         rowsMarker[*vit] = (*mit).first;
   }
   
   int nRowsCoreTmp  = 0;
   for(i = 0; i < nRows; i++){
      if(rowsMarker[i] == -1)
         rowsCore[nRowsCoreTmp++]   = i;
   }
   assert(nRowsCoreTmp == nRowsCore);
   
   UTIL_MSG(m_appParam.LogLevel, 3,
            (*m_osLog) << "Core  Rows:";
            for(i = 0; i < nRowsCore; i++)
               (*m_osLog) << rowsCore[i] << " ";
            (*m_osLog) << "\n";
            );
      
   //---
   //--- Construct the objective function.
   //---
   m_objective = new double[nCols];
   if(!m_objective)
      throw UtilExceptionMemory("createModels", "MILPBlock_DecompApp");
   memcpy(m_objective, 
          m_mpsIO.getObjCoefficients(), nCols * sizeof(double));
   setModelObjective(m_objective);

   //---
   //--- Construct the core matrix.
   //---
   DecompConstraintSet * modelCore = createModelPart(nRowsCore, rowsCore);
      
   //---
   //--- save a pointer so we can delete it later
   //---
   m_modelC = modelCore;

   //---
   //--- Construct the relaxation matrices.
   //---
   for(mit = m_blocks.begin(); mit != m_blocks.end(); mit++){
      vector<int> & rowsRelax  = (*mit).second;
      int           nRowsRelax = static_cast<int>(rowsRelax.size());
      if(m_appParam.LogLevel >= 1)
         (*m_osLog) << "Create model part nRowsRelax = " 
                    << nRowsRelax << " (Block=" << (*mit).first << ")" << endl;
      DecompConstraintSet * modelRelax 
         = createModelPart(nRowsRelax, &rowsRelax[0]);

      //---
      //--- find and set active columns
      //---
      set<int>::iterator sit;
      set<int> activeColsSet;
      findActiveColumns(rowsRelax, activeColsSet);
      for(sit = activeColsSet.begin(); sit != activeColsSet.end(); sit++)
         modelRelax->activeColumns.push_back(*sit);
      if(m_appParam.LogLevel >= 3){
         (*m_osLog) << "Active Columns:" << endl;
         UtilPrintVector(modelRelax->activeColumns, m_osLog);
         UtilPrintVector(modelRelax->activeColumns, 
                         modelCore->getColNames(), m_osLog);
      }

      //---
      //--- save a pointer so we can delete it later
      //---
      m_modelR.insert(make_pair((*mit).first, modelRelax));
   }

   //---
   //--- figure out which columns are not active in any subprobs
   //---  we refer to these as "master-only" variables
   //---
   int * colMarker = new int[nCols];
   if(!colMarker)
      throw UtilExceptionMemory("createModels", "MILPBlock_DecompApp");
   UtilFillN(colMarker, nCols, 0);
   
   vector<int>                      ::iterator vi;
   map   <int, DecompConstraintSet*>::iterator mdi;
   for(mdi = m_modelR.begin(); mdi != m_modelR.end(); mdi++){
      vector<int> & activeColumns = (*mdi).second->activeColumns;
      for(vi = activeColumns.begin(); vi != activeColumns.end(); vi++){
	 colMarker[*vi] = 1;
      }
   }

   for(i = 0; i < nCols; i++){
      if(!colMarker[i]){
         if(m_appParam.LogLevel >= 3){
            (*m_osLog) << "Column " << setw(5) << i << " -> "
                       << setw(25) << modelCore->colNames[i]
                       << " is not in union of blocks." << endl;
         }
	 modelCore->masterOnlyCols.push_back(i);
      }
   }
   if(m_appParam.LogLevel >= 3){
      (*m_osLog) << "Master only columns:" << endl;
      UtilPrintVector(modelCore->masterOnlyCols, m_osLog);
      UtilPrintVector(modelCore->masterOnlyCols, 
                      modelCore->colNames, m_osLog);
   }

   //---
   //--- set core and system in framework
   //---
   setModelCore(modelCore, "core");
   
   for(mdi = m_modelR.begin(); mdi != m_modelR.end(); mdi++){
      //---
      //--- append master-only cols to each block
      //---
      //vector<int> & activeColumns = (*mdi).second->activeColumns;
      //activeColumns.insert(activeColumns.end(),
      //		   modelCore->masterOnlyCols.begin(),
      //		   modelCore->masterOnlyCols.end());
      
      //---
      //--- fix column bounds on non-active columns
      //---
      (*mdi).second->fixNonActiveColumns();
      
      //---
      //--- set system in framework
      //---
      setModelRelax((*mdi).second,
                    "relax" + UtilIntToStr((*mdi).first),
                    (*mdi).first);
   }

   //---
   //--- create an extra "empty" block for the master-only vars
   //---   since I don't know what OSI will do with empty problem
   //---   we will make column bounds explicity rows
   //---
   int nMasterOnlyCols = static_cast<int>(modelCore->masterOnlyCols.size());
   if(nMasterOnlyCols){
      if(m_appParam.LogLevel >= 1)
         (*m_osLog) << "Create model part Master-Only." << endl;

      if(m_appParam.MasterOnlyOneBlock){
         DecompConstraintSet * modelMasterOnly
            = createModelMasterOnly(modelCore->masterOnlyCols);
         int nBlocks = static_cast<int>(m_blocks.size());
         m_modelR.insert(make_pair(nBlocks, modelMasterOnly));

         //---
         //--- set system in framework
         //---
         setModelRelax(modelMasterOnly, "master_only", nBlocks);
         
         if(m_appParam.LogLevel >= 3){
            (*m_osLog) << "Active Columns:" << endl;
            UtilPrintVector(modelMasterOnly->activeColumns, m_osLog);
            UtilPrintVector(modelMasterOnly->activeColumns, 
                            modelCore->getColNames(), m_osLog);
         }
      }
      else{
         //createModelMasterOnlys(modelCore->masterOnlyCols);
         createModelMasterOnlys2(modelCore->masterOnlyCols);
      }
   }
      
   //---
   //--- free up local memory
   //---
   UTIL_DELARR(rowsMarker);
   UTIL_DELARR(rowsCore);
   UTIL_DELARR(colMarker);
   
   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "createModels()", m_appParam.LogLevel, 2);   
}

/*
#if 0
//===========================================================================//
DecompSolverStatus 
MILPBlock_DecompApp::solveRelaxedNest(const int          whichBlock,
				      const double     * redCostX,
				      const double       convexDual,
				      DecompVarList    & varList){
   
   //---
   //--- solve full model heuristically  as IP
   //---   if get incumbent, break them out into approriate blocks
   //---   and return those partial columns
   //---


   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "solveRelaxedNest()", m_appParam.LogLevel, 2);
   
   /////STOP
   //---
   //--- this allows user direct access to access methods in 
   //---   algorithm interface (in case they want to use any 
   //---   of its data)
   //---
   //--- for example, if the user wants to enforce the branching
   //---   decisions in the oracle
   //--- TODO: can this be done using mcknap solver?
   //---
   //const DecompAlgo * decompAlgo = getDecompAlgo();
   //const double     * colLBNode  = decompAlgo->getColLBNode();
   //const double     * colUBNode  = decompAlgo->getColUBNode();

   //---
   //--- in the case where oracle=MCKP0, we have a specialized solver
   //--- in the case where oracle=MDKP,  we do not have a specialized solver
   //---   so, we just return with no solution and the framework will 
   //---   attempt to use the built-in MILP solver
   //---
   DecompSolverStatus solverStatus = DecompSolStatNoSolution;
   if(m_appParam.ModelNameRelax == "MCKP0"){
      vector<int>           solInd;
      vector<double>        solEls;
      double                varRedCost    = 0.0;
      double                varOrigCost   = 0.0;
      MMKP_MCKnap         * mcknapK       = m_mcknap[whichBlock];
      //TODO: check status return codes here
      mcknapK->solveMCKnap(redCostX, m_objective,
                           solInd, solEls, varRedCost, varOrigCost);
      assert(static_cast<int>(solInd.size()) == m_instance.getNGroupRows());
      
      UTIL_DEBUG(m_param.LogDebugLevel, 4,
                 printf("PUSH var with k = %d RC = %g origCost = %g\n", 
                        whichBlock, varRedCost - convexDual, varOrigCost);
                 );
      
      //the user should not have to calculate orig cost too
      //  the framework can do this... in fact the framework shoudl
      //  calculate red-cost too... but the user might want to check this stuff
      
      //this is way too confusing for user to remember they need -alpha!
      //  let framework do that - also setting the block id - framework!
      DecompVar * var = new DecompVar(solInd, solEls, 
                                      varRedCost - convexDual, varOrigCost);
      var->setBlockId(whichBlock);
      varList.push_back(var);      
      solverStatus = DecompSolStatOptimal;
   }
           
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "solveRelaxedNest()", m_appParam.LogLevel, 2);
   return solverStatus;
}
#endif
*/
