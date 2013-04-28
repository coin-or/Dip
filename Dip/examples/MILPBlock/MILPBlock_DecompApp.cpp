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
   //--- is there a permutation file?
   //---  this file just remaps the row ids
   //--- (for use in submission of atm to MIPLIB2010 and debugging)
   //---
   map<int,int>           permute;
   map<int,int>::iterator mit;
   string       fileNameP = m_appParam.DataDir 
      + UtilDirSlash() + m_appParam.PermuteFile;
   
   if(m_appParam.PermuteFile.size() > 0){
      ifstream isP;
      int      rowIdOld, rowIdNew;
      //---
      //--- open file streams
      //---
      UtilOpenFile(isP, fileName.c_str());
      while(!isP.eof()){
	 if(isP.eof()) break;
	 isP >> rowIdOld >> rowIdNew;
	 permute.insert(make_pair(rowIdOld, rowIdNew));
      }
      isP.close();
   }

   
   //---
   //--- open file streams
   //---
   UtilOpenFile(is, fileName.c_str());

   int i, rowId, rowIdP, numRowsInBlock, blockId;
   if(m_appParam.LogLevel >= 1)
      (*m_osLog) << "Reading " << fileName << endl;

   map<int, vector<int> > blocks;
   map<int, vector<int> >::iterator blocksIt;
   if(m_appParam.BlockFileFormat == "List" ||
      m_appParam.BlockFileFormat == "LIST"){

      //---
      //--- The block file defines those rows in each block.
      //---   <block id>  <num rows in block>
      //---     <row ids...>
      //---   <block id>  <num rows in block>
      //---     <row ids...>
      //---      
      while(!is.eof()){
	 is >> blockId;
	 is >> numRowsInBlock;
	 if(is.eof()) break;
	 vector<int> rowsInBlock;
	 for(i = 0; i < numRowsInBlock; i++){
	    is >> rowId;
	    mit = permute.find(rowId);
	    if(mit != permute.end())
	       rowsInBlock.push_back(mit->second);
	    else
	       rowsInBlock.push_back(rowId);
	 }
	 blocks.insert(make_pair(blockId, rowsInBlock));
	 if(is.eof()) break;      
      }
   }
   else if(m_appParam.BlockFileFormat == "Pair" ||
	   m_appParam.BlockFileFormat == "PAIR"){
      //---
      //--- <block id> <row id> 
      //---  ...
      //---
      is >> blockId;
      while(!is.eof()){
	 is >> rowId;
	 mit = permute.find(rowId);
	 if(mit != permute.end())
	    rowIdP = mit->second;
	 else
	    rowIdP = rowId;
	 blocksIt = blocks.find(blockId);
	 if(blocksIt != blocks.end())
	    blocksIt->second.push_back(rowIdP);
	 else{
	    vector<int> rowsInBlocks;
	    rowsInBlocks.push_back(rowIdP);
	    blocks.insert(make_pair(blockId, rowsInBlocks));	 
	 }
	 is >> blockId;
	 if(is.eof()) break;

      }      
   } else if(m_appParam.BlockFileFormat == "PairName" ||
	     m_appParam.BlockFileFormat == "PAIRNAME"){
      //---
      //--- <block id> <row name> 
      //---  ...
      //---

      //---
      //--- first create a map from row name to row id from mps
      //---   CHECK: mps to OSI guaranteed to keep order of rows?
      //---
      map<string, int>           rowNameToId;
      map<string, int>::iterator rowNameToIdIt;
      for(i = 0; i < m_mpsIO.getNumRows(); i++){
	 rowNameToId.insert(make_pair(m_mpsIO.rowName(i), i));
      }
      
      string rowName     = "";
      is >> blockId;
      while(!is.eof()){	 
	 is >> rowName;
	 if(is.eof()) 
	    break;
	 rowNameToIdIt = rowNameToId.find(rowName);
	 if(rowNameToIdIt != rowNameToId.end()){
	    rowId = rowNameToIdIt->second;
	    //printf("rowName=%s rowId=%d\n", rowName.c_str(), rowId);
	 }
	 else{
	    //---
	    //--- NOTE: this can happen if we use a presolved mps file
	    //---  with an original blocks file
	    //---
	    if(m_appParam.LogLevel >= 3){
	       (*m_osLog) << "Warning: Row name ("
			  << rowName << " in block file " 
			  << "is not found in mps file" << endl;
	    }
	    //throw UtilException("Invalid Input.", 
	    //		"readBlockFile", "MILPBlock_DecompApp");
	    rowId = -1;
	 }	    
	 if(rowId != -1){
	    mit = permute.find(rowId);
	    if(mit != permute.end())
	       rowIdP = mit->second;
	    else
	       rowIdP = rowId;	 
	    blocksIt = blocks.find(blockId);
	    if(blocksIt != blocks.end())
	       blocksIt->second.push_back(rowIdP);
	    else{
	       vector<int> rowsInBlocks;
	       rowsInBlocks.push_back(rowIdP);
	       blocks.insert(make_pair(blockId, rowsInBlocks));	 
	    }
	 }
	 is >> blockId;
	 if(is.eof()) 
	    break;
      }      
   } else{
      cerr << "Error: BlockFileFormat = " 
	   << m_appParam.BlockFileFormat 
	   << " is an invalid type. Valid types = (List,Pair,PairName)." 
	   << endl;
      throw UtilException("Invalid Parameter.", 
			  "readBlockFile", "MILPBlock_DecompApp");
   }

   //---
   //--- after presolve, some blocks might have been completely
   //---  removed - renumber the block ids - it is arbitrary anyway
   //--- and copy into class object m_blocks
   //---
   blockId = 0;
   for(blocksIt = blocks.begin(); blocksIt != blocks.end(); blocksIt++){
      m_blocks.insert(make_pair(blockId, blocksIt->second));
      blockId++;
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
   //exit(1);
   is.close();
}

//===========================================================================//
void MILPBlock_DecompApp::readInitSolutionFile(DecompVarList & initVars){

   //TODO: is this ok for sparse?

   ifstream is;
   string   fileName = m_appParam.DataDir 
      + UtilDirSlash() + m_appParam.InitSolutionFile;
   if(m_appParam.InitSolutionFile == "")
      return;

   //---
   //--- create map from col name to col index
   //---
   int                    i;
   map<string,int>        colNameToIndex;
   const vector<string> & colNames = m_modelC->getColNames();
   for(i = 0; i < m_modelC->getNumCols(); i++)
      colNameToIndex.insert(make_pair(colNames[i], i));
   
   //---
   //--- create a map from col index to block index
   //---
   map<int,int> colIndexToBlockIndex;
   map<int, DecompConstraintSet*>::iterator mit;
   const double * colLB = m_modelC->getColLB();
   const double * colUB = m_modelC->getColUB();
   for(mit = m_modelR.begin(); mit != m_modelR.end(); mit++){
      int                   blockIndex = mit->first;
      DecompConstraintSet * model      = mit->second;
      if(model->m_masterOnly){         
         colIndexToBlockIndex.insert(make_pair(model->m_masterOnlyIndex,
                                               blockIndex));
      }
      else{         
         const vector<int> & activeColumns = model->getActiveColumns();
         vector<int>::const_iterator vit;
         for(vit = activeColumns.begin(); vit != activeColumns.end(); vit++){
            colIndexToBlockIndex.insert(make_pair(*vit, blockIndex));
         }
      }
   }

   //---
   //--- open file streams
   //---
   UtilOpenFile(is, fileName.c_str());
   if(m_appParam.LogLevel >= 1)
      (*m_osLog) << "Reading " << fileName << endl;

   //---
   //--- create variables for each block of each solution
   //---
   int    solutionIndex, colIndex, blockIndex;
   string colName;
   double colValue;
   char   line[1000];
   map< pair<int,int>, pair< vector<int>,vector<double> > > varTemp;
   map< pair<int,int>, pair< vector<int>,vector<double> > >::iterator it;
   is.getline(line, 1000);


   //TODO? master-only
   // 1. if user gives lb, then add lb only
   //    if 0, add 0-col? or just let it take care of from PI?   
   // 2. if user gives ub, then add ub only
   // 3. if user gives betwen bounds, then add lb and ub
   //    unless it is general integer
   while(!is.eof()){
      is >> solutionIndex >> colName >> colValue;
      if(is.eof()) break;
      colIndex        = colNameToIndex[colName];
      blockIndex      = colIndexToBlockIndex[colIndex];
      DecompConstraintSet * model = m_modelR[blockIndex];
      if(model->m_masterOnly){
         printf("MasterOnly col=%s value=%g lb=%g ub=%g",
                colName.c_str(), colValue, colLB[colIndex], colUB[colIndex]);
         if(colValue < (colUB[colIndex]-1.0e-5) &&
            colValue > (colLB[colIndex]+1.0e-5)){
            printf(" --> in between bounds");
            //TODO: if so, should add both lb and ub
         }
         printf("\n");
      }
      pair<int,int> p = make_pair(solutionIndex, blockIndex);
      it = varTemp.find(p);
      if(it == varTemp.end()){         
         vector<int>    ind;
         vector<double> els;
         ind.push_back(colIndex);
         els.push_back(colValue);
         varTemp.insert(make_pair(p, make_pair(ind, els)));
      }
      else{
         vector<int>    & ind = it->second.first;
         vector<double> & els = it->second.second;
         ind.push_back(colIndex);
         els.push_back(colValue);         
      }
   }

   //---
   //--- create DecompVar's from varTemp
   //---
   for(it = varTemp.begin(); it != varTemp.end(); it++){
      const pair<int,int>                 & indexPair  = it->first;
      pair< vector<int>, vector<double> > & columnPair = it->second;
      double      origCost = 0.0;
      for(i = 0; i < static_cast<int>(columnPair.first.size()); i++){
         origCost += columnPair.second[i] *
            m_objective[columnPair.first[i]];            
      }
      DecompVar * var = new DecompVar(columnPair.first,
                                      columnPair.second,
                                      -1.0,
                                      origCost);
      var->setBlockId(indexPair.second);

      var->print(m_osLog, colNames);

      initVars.push_back(var);
      printf("Adding initial variable with origCost = %g\n", origCost);
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


   int i;
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
      //0=cont, 1=integer
      model->m_masterOnlyIsInt = 
         (integerVars && integerVars[i]) ? true : false;
      if(colUB[i]            >  1.0e15 &&
	 m_appParam.ColumnUB >= 1.0e15)
	 (*m_osLog) << "WARNING: Master-only column " << i 
		    << " has unbounded upper bound. DIP does not"
		    << " yet support extreme rays. Please bound all"
		    << " variables or use the ColumnUB parameter." << endl;
      if(colLB[i]            <  -1.0e15 && 
	 m_appParam.ColumnLB <= -1.0e15)
	 (*m_osLog) << "WARNING: Master-only column " << i 
		    << " has unbounded lower bound. DIP does not"
		    << " yet support extreme rays. Please bound all"
		    << " variables or use the ColumnLB parameter." << endl;
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
void
MILPBlock_DecompApp::createModelPart(DecompConstraintSet * model,
				     const int             nRowsPart,
                                     const int           * rowsPart){
   
   const int      nCols       = m_mpsIO.getNumCols();
   const double * rowLB       = m_mpsIO.getRowLower();
   const double * rowUB       = m_mpsIO.getRowUpper();
   const double * colLB       = m_mpsIO.getColLower();
   const double * colUB       = m_mpsIO.getColUpper();
   const char   * integerVars = m_mpsIO.integerColumns();
   
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
      if(m_appParam.UseNames){
	 const char * rowName = m_mpsIO.rowName(r);
	 if(rowName)
	    model->rowNames.push_back(rowName);
      }
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
      if(m_appParam.UseNames){
	 const char * colName = m_mpsIO.columnName(i);
	 if(colName)
	    model->colNames.push_back(colName);
      }
      if(integerVars && integerVars[i]){
         model->integerVars.push_back(i);            
      }
   }   
}

//===========================================================================//
void
MILPBlock_DecompApp::createModelPartSparse(DecompConstraintSet * model,
					   const int             nRowsPart,
					   const int           * rowsPart){

   const int      nColsOrig   = m_mpsIO.getNumCols();
   const double * rowLB       = m_mpsIO.getRowLower();
   const double * rowUB       = m_mpsIO.getRowUpper();
   const double * colLB       = m_mpsIO.getColLower();
   const double * colUB       = m_mpsIO.getColUpper();
   const char   * integerVars = m_mpsIO.integerColumns();

   //---
   //--- set model as sparse
   //---
   model->setSparse(nColsOrig);

   bool                  isInteger;
   int                   nCols, origIndex, newIndex;
   vector<int>::iterator vit;
   newIndex = 0;
   for(vit  = model->activeColumns.begin();
       vit != model->activeColumns.end(); vit++){
      origIndex = *vit;
      if(integerVars && integerVars[origIndex])
         isInteger  = true;
      else
         isInteger  = false;      
      model->pushCol(colLB[origIndex], 
		     colUB[origIndex],
		     isInteger,
		     origIndex);	     

      //---
      //--- big fat hack... we don't deal with dual rays yet,
      //---  so, we assume subproblems are bounded
      //---
      if(m_appParam.ColumnUB < 1.0e15){
	 if(colUB[origIndex] > 1.0e15){
	    model->colUB[newIndex] = m_appParam.ColumnUB;
	 }
      }
      if(m_appParam.ColumnLB > -1.0e15){
	 if(colLB[origIndex] < -1.0e15){
	    model->colLB[newIndex] = m_appParam.ColumnLB;
	 }
      }

      if(m_appParam.UseNames){
	 const char * colName = m_mpsIO.columnName(origIndex);
	 if(colName)
	    model->colNames.push_back(colName);
      }
      newIndex++;
   }

   nCols    = static_cast<int>(model->activeColumns.size());
   assert(static_cast<int>(model->colLB.size()) == nCols);
   assert(static_cast<int>(model->colUB.size()) == nCols);
   
   model->M = new CoinPackedMatrix(false, 0.0, 0.0);
   if(!model->M)
      throw UtilExceptionMemory("createModels", "MILPBlock_DecompApp");
   model->M->setDimensions(0, nCols);
   model->reserve(nRowsPart, nCols);

   //---
   //--- for each row in rowsPart, create the row using sparse mapping
   //---
   int                      i, k, r, begInd;
   const map<int,int>     & origToSparse   = model->getMapOrigToSparse();   
   const CoinPackedMatrix * M              = m_mpsIO.getMatrixByRow(); 
   const int              * matInd         = M->getIndices();
   const CoinBigIndex     * matBeg         = M->getVectorStarts();
   const int              * matLen         = M->getVectorLengths();
   const double           * matVal         = M->getElements();
   const int              * matIndI        = NULL;
   const double           * matValI        = NULL;

   vector<CoinBigIndex>   & rowBeg         = model->m_rowBeg;//used as temp
   vector<int         >   & rowInd         = model->m_rowInd;//used as temp
   vector<double      >   & rowVal         = model->m_rowVal;//used as temp
   map<int,int>::const_iterator mit;

   begInd = 0;
   rowBeg.push_back(0);
   for(i = 0; i < nRowsPart; i++){
      r = rowsPart[i];
      if(m_appParam.UseNames){
	 const char * rowName = m_mpsIO.rowName(r);
	 if(rowName)
	    model->rowNames.push_back(rowName);
      }
      model->rowLB.push_back(rowLB[r]);
      model->rowUB.push_back(rowUB[r]);

      matIndI = matInd + matBeg[r];
      matValI = matVal + matBeg[r];
      for(k = 0; k < matLen[r]; k++){
	 origIndex = matIndI[k];
	 mit       = origToSparse.find(origIndex);
	 assert(mit != origToSparse.end());
	 rowInd.push_back(mit->second);
	 rowVal.push_back(matValI[k]);			  
      }
      begInd += matLen[r];
      rowBeg.push_back(begInd);
   }
   model->M->appendRows(nRowsPart,
			&rowBeg[0],
			&rowInd[0],
			&rowVal[0]);
   rowBeg.clear();
   rowInd.clear();
   rowVal.clear();
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
   if(m_appParam.ObjectiveSense == -1){
      for(i = 0; i < nCols; i++)
	 m_objective[i] *= -1;
   }
   setModelObjective(m_objective);

   //---
   //--- Construct the core matrix.
   //---
   DecompConstraintSet * modelCore = new DecompConstraintSet();
   createModelPart(modelCore, nRowsCore, rowsCore);
      
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

      DecompConstraintSet * modelRelax = new DecompConstraintSet();
      CoinAssertHint(modelRelax, "Error: Out of Memory");

      //---
      //--- find and set active columns
      //---
      set<int>::iterator sit;
      set<int> activeColsSet;
      findActiveColumns(rowsRelax, activeColsSet);
      for(sit = activeColsSet.begin(); sit != activeColsSet.end(); sit++)
         modelRelax->activeColumns.push_back(*sit);      	 
      
      if(m_appParam.UseSparse){
         //---
         //--- create model part (using sparse API)
         //---         
	 createModelPartSparse(modelRelax, nRowsRelax, &rowsRelax[0]);
      }
      else{
         //---
         //--- create model part (using dense API)
         //---         
	 createModelPart(modelRelax, nRowsRelax, &rowsRelax[0]);
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
            if(modelCore->getColNames().size() > 0)
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
      if(modelCore->getColNames().size() > 0)
	 UtilPrintVector(modelCore->masterOnlyCols, 
			 modelCore->getColNames(), m_osLog);
   }

   //---
   //--- set core and system in framework
   //---
   setModelCore(modelCore, "core");
   
   for(mdi = m_modelR.begin(); mdi != m_modelR.end(); mdi++){
      DecompConstraintSet * modelRelax = (*mdi).second;
      //---
      //--- set system in framework
      //---
      setModelRelax((*mdi).second,
                    "relax" + UtilIntToStr((*mdi).first),
                    (*mdi).first);

      if(m_appParam.LogLevel >= 3){
	 (*m_osLog) << "Active Columns:" << endl;
	 UtilPrintVector(modelRelax->activeColumns, m_osLog);
	 if(modelCore->getColNames().size() > 0)
	    UtilPrintVector(modelRelax->activeColumns, 
			    modelCore->getColNames(), m_osLog);
      }
   }

   //---
   //--- create an extra "empty" block for the master-only vars
   //---   since I don't know what OSI will do with empty problem
   //---   we will make column bounds explicity rows
   //---
   ///////////STOP - don't need anymore if DECOMP_MASTERONLY_DIRECT
#if 1
   int nMasterOnlyCols = static_cast<int>(modelCore->masterOnlyCols.size());
   if(nMasterOnlyCols){
      if(m_appParam.LogLevel >= 1)
         (*m_osLog) << "Create model part Master-Only." << endl;
      createModelMasterOnlys(modelCore->masterOnlyCols);
   }
#endif
      
   //---
   //--- free up local memory
   //---
   UTIL_DELARR(rowsMarker);
   UTIL_DELARR(rowsCore);
   UTIL_DELARR(colMarker);
   
   UtilPrintFuncEnd(m_osLog, m_classTag,
		    "createModels()", m_appParam.LogLevel, 2);   
   //exit(1);
}


//===========================================================================//
int MILPBlock_DecompApp::generateInitVars(DecompVarList & initVars){	
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "generateInitVars()", m_appParam.LogLevel, 2);

   readInitSolutionFile(initVars);
   
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "generateInitVars()", m_appParam.LogLevel, 2);  
   return static_cast<int>(initVars.size());
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
