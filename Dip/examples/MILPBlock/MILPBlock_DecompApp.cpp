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
   printf("Reading %s\n", fileName.c_str());
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
	    //printf("block=%d rowId=%d\n", blockId, rowId);
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

   map<int, vector<int> >::iterator mit;
   vector<int>           ::iterator vit;
   for(mit = m_blocks.begin(); mit != m_blocks.end(); mit++){
      (*m_osLog) << "Block " << (*mit).first << " : ";
      for(vit = (*mit).second.begin(); vit != (*mit).second.end(); vit++)
	 (*m_osLog) << (*vit) << " ";
      (*m_osLog) << endl;
   }
   
   is.close();
}

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
DecompConstraintSet * 
MILPBlock_DecompApp::createModelMasterOnly(vector<int> & masterOnlyCols){
					  
   const int      nCols       = m_mpsIO.getNumCols();
   const double * colLB       = m_mpsIO.getColLower();
   const double * colUB       = m_mpsIO.getColUpper();
   const char   * integerVars = m_mpsIO.integerColumns();
   int            nMasterOnlyCols =
      static_cast<int>(masterOnlyCols.size());

   printf("nCols           = %d\n", nCols);
   printf("nMasterOnlyCols = %d\n", nMasterOnlyCols);

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
      if(m_appParam.ColumnUB < 1.0e15)
	 if(colUB[i] > 1.0e15)
	    model->colUB[i] = m_appParam.ColumnUB;
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
   printf("masteronly colindex=%d\n", colIndex);
   printf("colUB now rowUB    =%g\n", model->colUB[colIndex]);
   row.insert(colIndex, 1.0);
   model->appendRow(row, -DecompInf, model->colUB[colIndex], rowName);

   printf("model numcols= %d\n", model->getNumCols());
   printf("model numrows= %d\n", model->getNumRows());

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
	 //printf("i: %5d lb: %8.5f ub:%8.5f\n", i, colLB[i], colUB[i]);
	 if(colUB[i] > 1.0e15){
	    //printf("colUB[%d]: %g\n", i, colUB[i]);
	    model->colUB[i] = m_appParam.ColumnUB;
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
      throw UtilExceptionMemory("createModels", "MMKP_DecompApp");
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
      printf("Create model part nRowsRelax = %d (Block=%d)\n", 
	     nRowsRelax, (*mit).first); 
      fflush(stdout);
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
      printf("Active Columns:\n");
      UtilPrintVector(modelRelax->activeColumns);
      UtilPrintVector(modelRelax->activeColumns, modelCore->getColNames());

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
      throw UtilExceptionMemory("createModels", "MMKP_DecompApp");
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
	 (*m_osLog) << "Column " << setw(5) << i << " -> "
		    << setw(25) << modelCore->colNames[i]
		    << " is not in union of blocks." << endl;
	 modelCore->masterOnlyCols.push_back(i);
      }
   }
   (*m_osLog) << "Master only columns:" << endl;
   UtilPrintVector(modelCore->masterOnlyCols);
   UtilPrintVector(modelCore->masterOnlyCols, modelCore->colNames);

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
      printf("Create model part Master-Only.\n");
      DecompConstraintSet * modelMasterOnly
	 = createModelMasterOnly(modelCore->masterOnlyCols);
      UtilPrintVector(modelMasterOnly->activeColumns);
      UtilPrintVector(modelMasterOnly->activeColumns, 
		      modelCore->getColNames());
      int nBlocks = static_cast<int>(m_blocks.size());
      m_modelR.insert(make_pair(nBlocks, modelMasterOnly));
      //---
      //--- set system in framework
      //---
      setModelRelax(modelMasterOnly, "master_only", nBlocks);
      
      printf("Active Columns:\n");
      UtilPrintVector(modelMasterOnly->activeColumns);
      UtilPrintVector(modelMasterOnly->activeColumns, 
		      modelCore->getColNames());
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



