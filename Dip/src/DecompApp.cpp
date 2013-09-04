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
#include "DecompConfig.h"
#include <vector>
#include <set>
#include <fstream>
#include <string>
#include "iterator"
//#if defined(autoDecomp) && defined(PaToH)

#include <vector>
#include <set>
#include <fstream>
#include <string>
#include "iterator"
#if defined(_OPENMP)
#include "omp.h"
#endif
//#if defined(autoDecomp) && defined(PaToH)
#if  defined(PaToH)

#include "patoh.h"

//#elif defined(autoDecomp)
#else

extern "C" {
#if defined (COIN_HAS_HMETIS)
#include "hmetis.h"
#endif
}

#endif

using namespace std;

// --------------------------------------------------------------------- //
void DecompApp::startupLog()
{
   if (m_param.LogLevel >= 0) {
      (*m_osLog)
            << "\n========================================================"
            << "\n========================================================"
            <<   "\nWelcome to the DIP Decomposition Framework"
            <<   "\nCopyright 2002-2013 Lehigh University and others"
            <<   "\nAll Rights Reserved"
            <<   "\nDistributed under the Eclipse Public License 1.0"
            <<   "\nVersion: " << DIP_VERSION
            <<   "\nBuild Date: " << __DATE__
#ifdef DIP_SVN_REV
            <<   "\nRevision Number: " << DIP_SVN_REV
#endif
            << "\n========================================================"
            << "\n========================================================"
            << "\n";
   }

   if (m_param.LogLevel > 1) {
      //m_param.dumpSettings(m_osLog);
   }
}

// --------------------------------------------------------------------- //
int DecompApp::generateInitVars(DecompVarList& initVars)
{
   // ---
   // --- this function does nothing by default
   // ---
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "generateInitVars()", m_param.LogDebugLevel, 2);
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "generateInitVars()", m_param.LogDebugLevel, 2);
   return 0;
}

// --------------------------------------------------------------------- //
int DecompApp::generateCuts(const double*   x,
                            DecompCutList& newCuts)
{
   // ---
   // --- this function does nothing by default
   // ---
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "generateCuts()", m_param.LogDebugLevel, 2);
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "generateCuts()", m_param.LogDebugLevel, 2);
   return 0;
}

/*-------------------------------------------------------------------------*/
void DecompApp::printOriginalColumn(const int index,
                                    ostream* os) const
{
   (*os) << index << " ";
}

/*-------------------------------------------------------------------------*/
void DecompApp::printOriginalSolution(const int              n_cols,
                                      const vector<string>& colNames,
                                      const double*          solution,
                                      ostream*               os) const
{
   int i;
   bool hasNames = false;

   //---
   //--- do we have column names?
   //---
   if (colNames.size() > 0) {
      hasNames = true;
   }

   (*os) << setiosflags(ios::fixed | ios::showpoint);

   for (i = 0; i < n_cols; i++) {
      if (!UtilIsZero(solution[i])) {
         printOriginalColumn(i, os);

         if (hasNames)
            (*os) << "\t" << colNames[i]
                  << "\t" << solution[i] << endl;
         else {
            (*os) << "\t" << solution[i] << endl;
         }
      }
   }

   (*os) << resetiosflags(ios::fixed | ios::showpoint | ios::scientific);
}

/*
 *  The following methods are from MILPBlock_DecompApp
 */

void DecompApp::initializeApp(UtilParameters& utilParam)
{
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "initializeApp()", m_param.LogLevel, 2);

   //---
   //--- get application parameters
   //---

   if (m_param.LogLevel >= 1) {
      m_param.dumpSettings();
   }

   if (!m_param.Concurrent) {
      //---
      //--- read block file
      //---
      readBlockFile();
   } else
      // automatic structure detection
   {
      #pragma omp critical
      singlyBorderStructureDetection();
   }

   /*
    * After identifying the strucuture either through files or
    * automatic structure detection, call the method below to
    * create models
    *
    */
   createModels();
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "initializeApp()", m_param.LogLevel, 2);
}


const CoinPackedMatrix* DecompApp::readProblem(UtilParameters& utilParam)
{
   //---
   //--- get application parameters
   //---
   m_param.getSettings(utilParam);

   if (m_param.LogLevel >= 1) {
      m_param.dumpSettings();
   }

   //---
   //--- read MILP instance (mps format)
   //---
   string fileName;

   if (m_param.DataDir != "") {
      fileName = m_param.DataDir + UtilDirSlash() + m_param.Instance;
   } else {
      fileName = m_param.Instance;
   }

   if (m_param.Instance.empty()) {
      cerr << "================================================" << std::endl
           << "Usage:"
           << "./dip  --MILP:BlockFileFormat List" << std::endl
           << "       --MILP:Instance /FilePath/ABC.mps" << std::endl
           << "       --MILP:BlockFile /FilePath/ABC.block" << std::endl
           << "================================================" << std::endl
           << std::endl;
      exit(0);
   }

   int rstatus = 0;
   bool foundFormat = false;

   if (m_param.InstanceFormat == "") {
      string::size_type idx = fileName.rfind('.');

      if (idx != string::npos) {
         string extension = fileName.substr(idx + 1);

         if (extension == "MPS" || extension == "mps") {
            m_param.InstanceFormat = "MPS";
         } else if (extension == "LP" || extension == "lp") {
            m_param.InstanceFormat = "LP";
         }
      } else {
         cerr << "File format not specified and no file extension" << endl;
         throw UtilException("I/O Error.", "initializeApp", "DecompApp");
      }
   }

   if (m_param.InstanceFormat == "MPS") {
      m_mpsIO.messageHandler()->setLogLevel(m_param.LogLpLevel);
   } else if (m_param.InstanceFormat == "LP") {
      m_lpIO.messageHandler()->setLogLevel(m_param.LogLpLevel);
   }

   if (m_param.InstanceFormat == "MPS") {
      rstatus = m_mpsIO.readMps(fileName.c_str());
      foundFormat = true;
   } else if (m_param.InstanceFormat == "LP") {
      m_lpIO.readLp(fileName.c_str());
      foundFormat = true;
   }

   if (!foundFormat) {
      cerr << "Error: Format = " << m_param.InstanceFormat << " unknown."
           << endl;
      throw UtilException("I/O Error.", "initalizeApp", "DecompApp");
   }

   if (rstatus < 0) {
      cerr << "Error: Filename = " << fileName << " failed to open." << endl;
      throw UtilException("I/O Error.", "initalizeApp", "DecompApp");
   }

   if (m_param.LogLevel >= 2)
      if (m_param.InstanceFormat == "MPS") {
         (*m_osLog) << "Objective Offset = "
                    << UtilDblToStr(m_mpsIO.objectiveOffset()) << endl;
      } else if (m_param.InstanceFormat == "LP") {
         (*m_osLog) << "Objective Offset = "
                    << UtilDblToStr(m_lpIO.objectiveOffset()) << endl;
      }

   //---
   //--- set best known lb/ub
   //---
   double offset = 0;

   if (m_param.InstanceFormat == "MPS") {
      offset = m_mpsIO.objectiveOffset();
   } else if (m_param.InstanceFormat == "LP") {
      offset = m_lpIO.objectiveOffset();
   }

   setBestKnownLB(m_param.BestKnownLB + offset);
   setBestKnownUB(m_param.BestKnownUB + offset);
   preprocess();

   if (m_param.InstanceFormat == "MPS") {
      return m_mpsIO.getMatrixByRow();
   } else if (m_param.InstanceFormat == "LP") {
      return m_lpIO.getMatrixByRow();
   }
}



void DecompApp::preprocess() {}

void DecompApp::readBlockFile()
{
   ifstream is;
   string fileName;

   if (m_param.DataDir != "") {
      fileName = m_param.DataDir + UtilDirSlash() + m_param.Instance;
   } else {
      fileName = m_param.BlockFile;
   }

   //---
   //--- is there a permutation file?
   //---  this file just remaps the row ids
   //--- (for use in submission of atm to MIPLIB2010 and debugging)
   //---
   map<int, int>               permute;
   map<int, int>::iterator     mit;

   if (m_param.PermuteFile.size() > 0) {
      if (m_param.DataDir != "") {
         fileName = m_param.DataDir + UtilDirSlash() + m_param.PermuteFile;
      } else {
         fileName = m_param.PermuteFile;
      }

      ifstream isP;
      int      rowIdOld, rowIdNew;
      //---
      //--- open file streams
      //---
      UtilOpenFile(isP, fileName.c_str());

      while (!isP.eof()) {
         if (isP.eof()) {
            break;
         }

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
   //---
   //--- first create a map from row name to row id from mps
   //---   CHECK: mps to OSI guaranteed to keep order of rows?
   //---
   map<string, int>           rowNameToId;
   map<string, int>::iterator rowNameToIdIt;
   int numRows = 0;

   if (m_param.InstanceFormat == "MPS") {
      numRows = m_mpsIO.getNumRows();

      for (i = 0; i < numRows; i++) {
         rowNameToId.insert(make_pair(m_mpsIO.rowName(i), i));
      }
   } else if (m_param.InstanceFormat == "LP") {
      numRows = m_lpIO.getNumRows();

      for (i = 0; i < numRows; i++) {
         rowNameToId.insert(make_pair(m_lpIO.rowName(i), i));
      }
   }

   if (m_param.LogLevel >= 1) {
      (*m_osLog) << "Reading " << fileName << endl;
   }

   map<int, vector<int> > blocks;
   map<int, vector<int> >::iterator blocksIt;

   if (m_param.BlockFileFormat == "") {
      string::size_type idx = fileName.rfind('.');

      if (idx != string::npos) {
         string extension = fileName.substr(idx + 1);

         if (extension == "DEC" || extension == "dec") {
            m_param.BlockFileFormat = "ZIBList";
         } else if (extension == "block" || extension == "blk") {
            m_param.BlockFileFormat = "List";
         }
      } else {
         cerr << "File format not specified and no file extension" << endl;
         throw UtilException("I/O Error.", "initializeApp", "DecompApp");
      }
   }

   if (m_param.BlockFileFormat == "List" ||
         m_param.BlockFileFormat == "LIST") {
      //---
      //--- The block file defines those rows in each block.
      //---   <block id>  <num rows in block>
      //---     <row ids...>
      //---   <block id>  <num rows in block>
      //---     <row ids...>
      //---
      string rowName;

      while (!is.eof()) {
         is >> blockId;
         is >> numRowsInBlock;

         if (is.eof()) {
            break;
         }

         vector<int> rowsInBlock;

         for (i = 0; i < numRowsInBlock; i++) {
            is >> rowId;
            mit = permute.find(rowId);

            if (mit != permute.end()) {
               rowsInBlock.push_back(mit->second);
            } else {
               rowsInBlock.push_back(rowId);
            }
         }

         blocks.insert(make_pair(blockId, rowsInBlock));

         if (is.eof()) {
            break;
         }
      }
   } else if (m_param.BlockFileFormat == "ZIBList" ||
              m_param.BlockFileFormat == "ZIBLIST") {
      //-- The block file defines those rows in each block.
      //--   NBLOCKS
      //--   <numBlocks>
      //--   BLOCK <block id>
      //--   <row names...>
      //--   BLOCK  <block id>
      //--   <row names...>
      int numBlocks = 0;
      string tmp, rowName;

      while (!numBlocks) {
         is >> tmp;

         if (tmp == "NBLOCKS") {
            is >> numBlocks;
         }
      }

      while (tmp != "BLOCK") {
         is >> tmp;
      }

      while (!is.eof() && rowName != "MASTERCONSS") {
         is >> blockId;
         vector<int> rowsInBlock;

         while (true) {
            is >> rowName;

            if (is.eof()) {
               break;
            }

            if (rowName == "BLOCK" || rowName == "MASTERCONSS") {
               break;
            }

            rowNameToIdIt = rowNameToId.find(rowName);

            if (rowNameToIdIt != rowNameToId.end()) {
               rowId = rowNameToIdIt->second;
            } else {
               std::cout << "Warning: Unrecognized row name" << rowName;
               std::cout << "in block file" << std::endl;
            }

            rowsInBlock.push_back(rowId);
         }

         blocks.insert(make_pair(blockId, rowsInBlock));

         if (is.eof()) {
            break;
         }
      }
   } else if (m_param.BlockFileFormat == "Pair" ||
              m_param.BlockFileFormat == "PAIR") {
      //---
      //--- <block id> <row id>
      //---  ...
      //---
      is >> blockId;

      while (!is.eof()) {
         is >> rowId;
         mit = permute.find(rowId);

         if (mit != permute.end()) {
            rowIdP = mit->second;
         } else {
            rowIdP = rowId;
         }

         blocksIt = blocks.find(blockId);

         if (blocksIt != blocks.end()) {
            blocksIt->second.push_back(rowIdP);
         } else {
            vector<int> rowsInBlocks;
            rowsInBlocks.push_back(rowIdP);
            blocks.insert(make_pair(blockId, rowsInBlocks));
         }

         is >> blockId;

         if (is.eof()) {
            break;
         }
      }
   } else if (m_param.BlockFileFormat == "PairName" ||
              m_param.BlockFileFormat == "PAIRNAME") {
      //---
      //--- <block id> <row name>
      //---  ...
      //---
      string rowName     = "";
      is >> blockId;

      while (!is.eof()) {
         is >> rowName;

         if (is.eof()) {
            break;
         }

         rowNameToIdIt = rowNameToId.find(rowName);

         if (rowNameToIdIt != rowNameToId.end()) {
            rowId = rowNameToIdIt->second;
            //printf("rowName=%s rowId=%d\n", rowName.c_str(), rowId);
         } else {
            //---
            //--- NOTE: this can happen if we use a presolved mps file
            //---  with an original blocks file
            //---
            if (m_param.LogLevel >= 3) {
               (*m_osLog) << "Warning: Row name ("
                          << rowName << " in block file "
                          << "is not found in instance file" << endl;
            }

            //throw UtilException("Invalid Input.",
            //		"readBlockFile", "DecompApp");
            rowId = -1;
         }

         if (rowId != -1) {
            mit = permute.find(rowId);

            if (mit != permute.end()) {
               rowIdP = mit->second;
            } else {
               rowIdP = rowId;
            }

            blocksIt = blocks.find(blockId);

            if (blocksIt != blocks.end()) {
               blocksIt->second.push_back(rowIdP);
            } else {
               vector<int> rowsInBlocks;
               rowsInBlocks.push_back(rowIdP);
               blocks.insert(make_pair(blockId, rowsInBlocks));
            }
         }

         is >> blockId;

         if (is.eof()) {
            break;
         }
      }
   } else {
      cerr << "Error: BlockFileFormat = "
           << m_param.BlockFileFormat
           << " is an invalid type. Valid types = (List,ZIBlist,Pair,PairName)."
           << endl;
      throw UtilException("Invalid Parameter.",
                          "readBlockFile", "DecompApp");
   }

   //---
   //--- after presolve, some blocks might have been completely
   //---  removed - renumber the block ids - it is arbitrary anyway
   //--- and copy into class object m_blocks
   //---
   blockId = 0;

   for (blocksIt = blocks.begin(); blocksIt != blocks.end(); blocksIt++) {
      m_blocks.insert(make_pair(blockId, blocksIt->second));
      blockId++;
   }

   if (m_param.LogLevel >= 3) {
      map<int, vector<int> >::iterator mit;
      vector<int>           ::iterator vit;

      for (mit = m_blocks.begin(); mit != m_blocks.end(); mit++) {
         (*m_osLog) << "Block " << (*mit).first << " : ";

         for (vit = (*mit).second.begin(); vit != (*mit).second.end(); vit++) {
            (*m_osLog) << (*vit) << " ";
         }

         (*m_osLog) << endl;
      }
   }

   //exit(1);
   is.close();
}


void DecompApp::readInitSolutionFile(DecompVarList& initVars)
{
   //TODO: is this ok for sparse?
   ifstream is;
   string   fileName = m_param.DataDir
                       + UtilDirSlash() + m_param.InitSolutionFile;

   if (m_param.InitSolutionFile == "") {
      return;
   }

   //---
   //--- create map from col name to col index
   //---
   int                    i;
   map<string, int>        colNameToIndex;
   const vector<string>& colNames = m_modelC->getColNames();

   for (i = 0; i < m_modelC->getNumCols(); i++) {
      colNameToIndex.insert(make_pair(colNames[i], i));
   }

   //---
   //--- create a map from col index to block index
   //---
   map<int, int> colIndexToBlockIndex;
   map<int, DecompConstraintSet*>::iterator mit;
   const double* colLB = m_modelC->getColLB();
   const double* colUB = m_modelC->getColUB();

   for (mit = m_modelR.begin(); mit != m_modelR.end(); mit++) {
      int                   blockIndex = mit->first;
      DecompConstraintSet* model      = mit->second;

      if (model->m_masterOnly) {
         colIndexToBlockIndex.insert(make_pair(model->m_masterOnlyIndex,
                                               blockIndex));
      } else {
         const vector<int>& activeColumns = model->getActiveColumns();
         vector<int>::const_iterator vit;

         for (vit = activeColumns.begin(); vit != activeColumns.end(); vit++) {
            colIndexToBlockIndex.insert(make_pair(*vit, blockIndex));
         }
      }
   }

   //---
   //--- open file streams
   //---
   UtilOpenFile(is, fileName.c_str());

   if (m_param.LogLevel >= 1) {
      (*m_osLog) << "Reading " << fileName << endl;
   }

   //---
   //--- create variables for each block of each solution
   //---
   int    solutionIndex, colIndex, blockIndex;
   string colName;
   double colValue;
   char   line[1000];
   map< pair<int, int>, pair< vector<int>, vector<double> > > varTemp;
   map< pair<int, int>, pair< vector<int>, vector<double> > >::iterator it;
   is.getline(line, 1000);

   //TODO? master-only
   // 1. if user gives lb, then add lb only
   //    if 0, add 0-col? or just let it take care of from PI?
   // 2. if user gives ub, then add ub only
   // 3. if user gives betwen bounds, then add lb and ub
   //    unless it is general integer
   while (!is.eof()) {
      is >> solutionIndex >> colName >> colValue;

      if (is.eof()) {
         break;
      }

      colIndex        = colNameToIndex[colName];
      blockIndex      = colIndexToBlockIndex[colIndex];
      DecompConstraintSet* model = m_modelR[blockIndex];

      if (model->m_masterOnly) {
         printf("MasterOnly col=%s value=%g lb=%g ub=%g",
                colName.c_str(), colValue, colLB[colIndex], colUB[colIndex]);

         if (colValue < (colUB[colIndex] - 1.0e-5) &&
               colValue > (colLB[colIndex] + 1.0e-5)) {
            printf(" --> in between bounds");
            //TODO: if so, should add both lb and ub
         }

         printf("\n");
      }

      pair<int, int> p = make_pair(solutionIndex, blockIndex);
      it = varTemp.find(p);

      if (it == varTemp.end()) {
         vector<int>    ind;
         vector<double> els;
         ind.push_back(colIndex);
         els.push_back(colValue);
         varTemp.insert(make_pair(p, make_pair(ind, els)));
      } else {
         vector<int>&     ind = it->second.first;
         vector<double>& els = it->second.second;
         ind.push_back(colIndex);
         els.push_back(colValue);
      }
   }

   //---
   //--- create DecompVar's from varTemp
   //---
   for (it = varTemp.begin(); it != varTemp.end(); it++) {
      const pair<int, int>&                  indexPair  = it->first;
      pair< vector<int>, vector<double> >& columnPair = it->second;
      double      origCost = 0.0;

      for (i = 0; i < static_cast<int>(columnPair.first.size()); i++) {
         origCost += columnPair.second[i] *
                     m_objective[columnPair.first[i]];
      }

      DecompVar* var = new DecompVar(columnPair.first,
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


void DecompApp::findActiveColumns(const vector<int>& rowsPart,
                                  set<int>&           activeColsSet)
{
   const CoinPackedMatrix* M = NULL;

   if (m_param.InstanceFormat == "MPS") {
      M    = m_mpsIO.getMatrixByRow();
   } else if (m_param.InstanceFormat == "LP") {
      M    = m_lpIO.getMatrixByRow();
   }

   const int*               ind  = M->getIndices();
   const int*               beg  = M->getVectorStarts();
   const int*               len  = M->getVectorLengths();
   const int*               indR = NULL;
   //---
   //--- which columns are present in this part's rows
   //---
   int k, r;
   vector<int>::const_iterator it;

   for (it = rowsPart.begin(); it != rowsPart.end(); it++) {
      r    = *it;
      indR = ind + beg[r];

      for (k = 0; k < len[r]; k++) {
         activeColsSet.insert(indR[k]);
      }
   }
}

void DecompApp::createModelMasterOnlys(vector<int>& masterOnlyCols)
{
   int            nBlocks     = static_cast<int>(m_blocks.size());
   int      nCols       = 0;
   double* colLB       = NULL;
   double* colUB       = NULL;
   char*    integerVars = NULL;

   if (m_param.InstanceFormat == "MPS") {
      nCols       = m_mpsIO.getNumCols();
      colLB       = const_cast <double*>(m_mpsIO.getColLower());
      colUB       = const_cast <double*>(m_mpsIO.getColUpper());
      integerVars = const_cast <char*>  (m_mpsIO.integerColumns());
   } else if (m_param.InstanceFormat == "LP") {
      nCols       = m_lpIO.getNumCols();
      colLB       = const_cast <double*>(m_lpIO.getColLower());
      colUB       = const_cast <double*>(m_lpIO.getColUpper());
      integerVars = const_cast <char*>  (m_lpIO.integerColumns());
   }

   int            nMasterOnlyCols =
      static_cast<int>(masterOnlyCols.size());

   if (m_param.LogLevel >= 1) {
      (*m_osLog) << "nCols           = " << nCols << endl;
      (*m_osLog) << "nMasterOnlyCols = " << nMasterOnlyCols << endl;
   }

   if (nMasterOnlyCols == 0) {
      return;
   }

   int i;
   vector<int>::iterator vit;

   for (vit = masterOnlyCols.begin(); vit != masterOnlyCols.end(); vit++) {
      i = *vit;
      //THINK:
      //  what-if master-only var is integer and bound is not at integer
      DecompConstraintSet* model = new DecompConstraintSet();
      model->m_masterOnly      = true;
      model->m_masterOnlyIndex = i;
      model->m_masterOnlyLB    = colLB[i];
      model->m_masterOnlyUB    = colUB[i];
      //0=cont, 1=integer
      model->m_masterOnlyIsInt =
         (integerVars && integerVars[i]) ? true : false;

      if (colUB[i]            >  1.0e15 &&
            m_param.ColumnUB >= 1.0e15)
         (*m_osLog) << "WARNING: Master-only column " << i
                    << " has unbounded upper bound. DIP does not"
                    << " yet support extreme rays. Please bound all"
                    << " variables or use the ColumnUB parameter." << endl;

      if (colLB[i]            <  -1.0e15 &&
            m_param.ColumnLB <= -1.0e15)
         (*m_osLog) << "WARNING: Master-only column " << i
                    << " has unbounded lower bound. DIP does not"
                    << " yet support extreme rays. Please bound all"
                    << " variables or use the ColumnLB parameter." << endl;

      if (m_param.ColumnUB <  1.0e15)
         if (colUB[i] >  1.0e15) {
            model->m_masterOnlyUB = m_param.ColumnUB;
         }

      if (m_param.ColumnLB > -1.0e15)
         if (colLB[i] < -1.0e15) {
            model->m_masterOnlyLB = m_param.ColumnLB;
         }

      m_modelR.insert(make_pair(nBlocks, model));
      setModelRelax(model,
                    "master_only" + UtilIntToStr(i), nBlocks);
      nBlocks++;
   }

   return;
}

void DecompApp::createModelPart(DecompConstraintSet* model,
                                const int             nRowsPart,
                                const int*            rowsPart)
{
   int      nCols       = 0;
   double* rowLB       = NULL;
   double* rowUB       = NULL;
   double* colLB       = NULL;
   double* colUB       = NULL;
   char*    integerVars = NULL;

   if (m_param.InstanceFormat == "MPS") {
      nCols = m_mpsIO.getNumCols();
      rowLB       = const_cast<double*>(m_mpsIO.getRowLower());
      rowUB       = const_cast<double*>(m_mpsIO.getRowUpper());
      colLB       = const_cast<double*>(m_mpsIO.getColLower());
      colUB       = const_cast<double*>(m_mpsIO.getColUpper());
      integerVars = const_cast<char*>  (m_mpsIO.integerColumns());
   } else if (m_param.InstanceFormat == "LP") {
      nCols = m_lpIO.getNumCols();
      rowLB       = const_cast<double*>(m_lpIO.getRowLower());
      rowUB       = const_cast<double*>(m_lpIO.getRowUpper());
      colLB       = const_cast<double*>(m_lpIO.getColLower());
      colUB       = const_cast<double*>(m_lpIO.getColUpper());
      integerVars = const_cast<char*>  (m_lpIO.integerColumns());
   }

   model->M = new CoinPackedMatrix(false, 0.0, 0.0);

   if (!model->M) {
      throw UtilExceptionMemory("createModels", "DecompApp");
   }

   model->reserve(nRowsPart, nCols);

   if (m_param.InstanceFormat == "MPS") {
      model->M->submatrixOf(*m_mpsIO.getMatrixByRow(), nRowsPart, rowsPart);
   } else if (m_param.InstanceFormat == "LP") {
      model->M->submatrixOf(*m_lpIO.getMatrixByRow(), nRowsPart, rowsPart);
   }

   //---
   //--- set the row upper and lower bounds
   //--- set the col upper and lower bounds
   //---
   int i, r;

   for (i = 0; i < nRowsPart; i++) {
      r = rowsPart[i];

      if (m_param.UseNames) {
         const char* rowName = NULL;

         if (m_param.InstanceFormat == "MPS") {
            rowName = m_mpsIO.rowName(r);
         } else if (m_param.InstanceFormat == "LP") {
            rowName = m_lpIO.rowName(r);
         }

         if (rowName) {
            model->rowNames.push_back(rowName);
         }
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
   if (m_param.ColumnUB < 1.0e15) {
      for (i = 0; i < nCols; i++) {
         if (colUB[i] > 1.0e15) {
            model->colUB[i] = m_param.ColumnUB;
         }
      }
   }

   if (m_param.ColumnLB > -1.0e15) {
      for (i = 0; i < nCols; i++) {
         if (colLB[i] < -1.0e15) {
            model->colLB[i] = m_param.ColumnLB;
         }
      }
   }

   //---
   //--- set the indices of the integer variables of modelRelax
   //---  also set the column names, if they exist
   //---
   for (i = 0; i < nCols; i++) {
      if (m_param.UseNames) {
         const char* colName = NULL;

         if (m_param.InstanceFormat == "MPS") {
            colName = m_mpsIO.columnName(i);
         } else if (m_param.InstanceFormat == "LP") {
            colName = m_lpIO.columnName(i);
         }

         if (colName) {
            model->colNames.push_back(colName);
         }
      }

      if (integerVars && integerVars[i]) {
         model->integerVars.push_back(i);
      }
   }
}



void DecompApp::createModelPartSparse(DecompConstraintSet* model,
                                      const int             nRowsPart,
                                      const int*            rowsPart)
{
   int      nColsOrig   = 0;
   double* rowLB       = NULL;
   double* rowUB       = NULL;
   double* colLB       = NULL;
   double* colUB       = NULL;
   char*    integerVars = NULL;

   if (m_param.InstanceFormat == "MPS") {
      nColsOrig   = m_mpsIO.getNumCols();
      rowLB       = const_cast<double*>(m_mpsIO.getRowLower());
      rowUB       = const_cast<double*>(m_mpsIO.getRowUpper());
      colLB       = const_cast<double*>(m_mpsIO.getColLower());
      colUB       = const_cast<double*>(m_mpsIO.getColUpper());
      integerVars = const_cast<char*>  (m_mpsIO.integerColumns());
   } else if (m_param.InstanceFormat == "LP") {
      nColsOrig   = m_lpIO.getNumCols();
      rowLB       = const_cast<double*>(m_lpIO.getRowLower());
      rowUB       = const_cast<double*>(m_lpIO.getRowUpper());
      colLB       = const_cast<double*>(m_lpIO.getColLower());
      colUB       = const_cast<double*>(m_lpIO.getColUpper());
      integerVars = const_cast<char*>  (m_lpIO.integerColumns());
   }

   //---
   //--- set model as sparse
   //---
   model->setSparse(nColsOrig);
   bool                  isInteger;
   int                   nCols, origIndex, newIndex;
   vector<int>::iterator vit;
   newIndex = 0;

   for (vit  = model->activeColumns.begin();
         vit != model->activeColumns.end(); vit++) {
      origIndex = *vit;

      if (integerVars && integerVars[origIndex]) {
         isInteger  = true;
      } else {
         isInteger  = false;
      }

      model->pushCol(colLB[origIndex],
                     colUB[origIndex],
                     isInteger,
                     origIndex);

      //---
      //--- big fat hack... we don't deal with dual rays yet,
      //---  so, we assume subproblems are bounded
      //---
      if (m_param.ColumnUB < 1.0e15) {
         if (colUB[origIndex] > 1.0e15) {
            model->colUB[newIndex] = m_param.ColumnUB;
         }
      }

      if (m_param.ColumnLB > -1.0e15) {
         if (colLB[origIndex] < -1.0e15) {
            model->colLB[newIndex] = m_param.ColumnLB;
         }
      }

      if (m_param.UseNames) {
         const char* colName = NULL;

         if (m_param.InstanceFormat == "MPS") {
            colName = m_mpsIO.columnName(origIndex);
         } else if (m_param.InstanceFormat == "LP") {
            colName = m_lpIO.columnName(origIndex);
         }

         if (colName) {
            model->colNames.push_back(colName);
         }
      }

      newIndex++;
   }

   nCols    = static_cast<int>(model->activeColumns.size());
   assert(static_cast<int>(model->colLB.size()) == nCols);
   assert(static_cast<int>(model->colUB.size()) == nCols);
   model->M = new CoinPackedMatrix(false, 0.0, 0.0);

   if (!model->M) {
      throw UtilExceptionMemory("createModels", "DecompApp");
   }

   model->M->setDimensions(0, nCols);
   model->reserve(nRowsPart, nCols);
   //---
   //--- for each row in rowsPart, create the row using sparse mapping
   //---
   int                      i, k, r, begInd;
   const map<int, int>&      origToSparse   = model->getMapOrigToSparse();
   const CoinPackedMatrix* M              = NULL;

   if (m_param.InstanceFormat == "MPS") {
      M              = m_mpsIO.getMatrixByRow();
   } else if (m_param.InstanceFormat == "LP") {
      M              = m_lpIO.getMatrixByRow();
   }

   const int*               matInd         = M->getIndices();
   const CoinBigIndex*      matBeg         = M->getVectorStarts();
   const int*               matLen         = M->getVectorLengths();
   const double*            matVal         = M->getElements();
   const int*               matIndI        = NULL;
   const double*            matValI        = NULL;
   vector<CoinBigIndex>&    rowBeg         = model->m_rowBeg;//used as temp
   vector<int         >&    rowInd         = model->m_rowInd;//used as temp
   vector<double      >&    rowVal         = model->m_rowVal;//used as temp
   map<int, int>::const_iterator mit;
   begInd = 0;
   rowBeg.push_back(0);

   for (i = 0; i < nRowsPart; i++) {
      r = rowsPart[i];

      if (m_param.UseNames) {
         const char* rowName = NULL;

         if (m_param.InstanceFormat == "MPS") {
            rowName = m_mpsIO.rowName(r);
         } else if (m_param.InstanceFormat == "LP") {
            rowName = m_lpIO.rowName(r);
         }

         if (rowName) {
            model->rowNames.push_back(rowName);
         }
      }

      model->rowLB.push_back(rowLB[r]);
      model->rowUB.push_back(rowUB[r]);
      matIndI = matInd + matBeg[r];
      matValI = matVal + matBeg[r];

      for (k = 0; k < matLen[r]; k++) {
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


void DecompApp::createModels()
{
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "createModels()", m_param.LogLevel, 2);
   //---
   //--- how many rows to put into relaxation
   //---
   int            i, nRowsRelax, nRowsCore;
   int      nRows       = 0;
   int      nCols       = 0;

   if (m_param.InstanceFormat == "MPS") {
      nRows       = m_mpsIO.getNumRows();
      nCols       = m_mpsIO.getNumCols();
   } else if (m_param.InstanceFormat == "LP") {
      nRows       = m_lpIO.getNumRows();
      nCols       = m_lpIO.getNumCols();
   }

   int            nBlocks     = static_cast<int>(m_blocks.size());
   map<int, vector<int> >::iterator mit;
   nRowsRelax = 0;

   for (mit = m_blocks.begin(); mit != m_blocks.end(); mit++) {
      nRowsRelax += static_cast<int>((*mit).second.size());
   }

   nRowsCore = nRows - nRowsRelax;
   UTIL_MSG(m_param.LogLevel, 2,
            (*m_osLog) << "Instance    = " << m_param.Instance << endl;
            (*m_osLog) << " nRows      = " << nRows     << endl;
            (*m_osLog) << " nCols      = " << nCols     << endl;
            (*m_osLog) << " nBlocks    = " << nBlocks   << endl;
            (*m_osLog) << " nRowsCore  = " << nRowsCore << endl;
            (*m_osLog) << " nRowsRelax = " << nRowsRelax
            << " [ " << 100 * nRowsRelax / nRows << " % ]" << endl;
           );
   //---
   //--- setup markers for core and relax rows
   //---
   int* rowsMarker = new int[nRows];
   int* rowsCore   = new int[nRowsCore];
   UtilFillN(rowsMarker, nRows, -1);//-1 will mark core rows

   for (mit = m_blocks.begin(); mit != m_blocks.end(); mit++) {
      vector<int>& rowsRelax = (*mit).second;
      vector<int>::iterator vit;

      for (vit = rowsRelax.begin(); vit != rowsRelax.end(); vit++) {
         rowsMarker[*vit] = (*mit).first;
      }
   }

   int nRowsCoreTmp  = 0;

   for (i = 0; i < nRows; i++) {
      if (rowsMarker[i] == -1) {
         rowsCore[nRowsCoreTmp++]   = i;
      }
   }

   assert(nRowsCoreTmp == nRowsCore);
   UTIL_MSG(m_param.LogLevel, 3,
            (*m_osLog) << "Core  Rows:";

            for (i = 0; i < nRowsCore; i++)
            (*m_osLog) << rowsCore[i] << " ";
            (*m_osLog) << "\n";
           );

   //---
   //--- Construct the objective function.
   //---
   double* objective = new double[nCols];

   if (!objective) {
      throw UtilExceptionMemory("createModels", "DecompApp");
   }

   if (m_param.InstanceFormat == "MPS") {
      memcpy(objective,
             m_mpsIO.getObjCoefficients(), nCols * sizeof(double));
   } else if (m_param.InstanceFormat == "LP") {
      memcpy(objective,
             m_lpIO.getObjCoefficients(), nCols * sizeof(double));
   }

   if (m_param.ObjectiveSense == -1) {
      for (i = 0; i < nCols; i++) {
         objective[i] *= -1;
      }
   }

   setModelObjective(objective);
   //---
   //--- Construct the core matrix.
   //---
   DecompConstraintSet* modelCore = new DecompConstraintSet();
   createModelPart(modelCore, nRowsCore, rowsCore);
   //---
   //--- save a pointer so we can delete it later
   //---
   m_modelC = modelCore;

   //---
   //--- Construct the relaxation matrices.
   //---
   for (mit = m_blocks.begin(); mit != m_blocks.end(); mit++) {
      vector<int>& rowsRelax  = (*mit).second;
      int           nRowsRelax = static_cast<int>(rowsRelax.size());

      if (m_param.LogLevel >= 1)
         (*m_osLog) << "Create model part nRowsRelax = "
                    << nRowsRelax << " (Block=" << (*mit).first << ")" << endl;

      DecompConstraintSet* modelRelax = new DecompConstraintSet();
      CoinAssertHint(modelRelax, "Error: Out of Memory");
      //---
      //--- find and set active columns
      //---
      set<int>::iterator sit;
      set<int> activeColsSet;
      findActiveColumns(rowsRelax, activeColsSet);

      for (sit = activeColsSet.begin(); sit != activeColsSet.end(); sit++) {
         modelRelax->activeColumns.push_back(*sit);
      }

      if (m_param.UseSparse) {
         //---
         //--- create model part (using sparse API)
         //---
         createModelPartSparse(modelRelax, nRowsRelax, &rowsRelax[0]);
      } else {
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
   int* colMarker = new int[nCols];

   if (!colMarker) {
      throw UtilExceptionMemory("createModels", "DecompApp");
   }

   UtilFillN(colMarker, nCols, 0);
   vector<int>                      ::iterator vi;
   map   <int, DecompConstraintSet*>::iterator mdi;

   for (mdi = m_modelR.begin(); mdi != m_modelR.end(); mdi++) {
      vector<int>& activeColumns = (*mdi).second->activeColumns;

      for (vi = activeColumns.begin(); vi != activeColumns.end(); vi++) {
         colMarker[*vi] = 1;
      }
   }

   for (i = 0; i < nCols; i++) {
      if (!colMarker[i]) {
         if (m_param.LogLevel >= 3) {
            if (modelCore->getColNames().size() > 0)
               (*m_osLog) << "Column " << setw(5) << i << " -> "
                          << setw(25) << modelCore->colNames[i]
                          << " is not in union of blocks." << endl;
         }

         modelCore->masterOnlyCols.push_back(i);
      }
   }

   if (m_param.LogLevel >= 3) {
      (*m_osLog) << "Master only columns:" << endl;
      UtilPrintVector(modelCore->masterOnlyCols, m_osLog);

      if (modelCore->getColNames().size() > 0)
         UtilPrintVector(modelCore->masterOnlyCols,
                         modelCore->getColNames(), m_osLog);
   }

   //---
   //--- set core and system in framework
   //---
   setModelCore(modelCore, "core");

   for (mdi = m_modelR.begin(); mdi != m_modelR.end(); mdi++) {
      DecompConstraintSet* modelRelax = (*mdi).second;
      //---
      //--- set system in framework
      //---
      setModelRelax((*mdi).second,
                    "relax" + UtilIntToStr((*mdi).first),
                    (*mdi).first);

      if (m_param.LogLevel >= 3) {
         (*m_osLog) << "Active Columns:" << endl;
         UtilPrintVector(modelRelax->activeColumns, m_osLog);

         if (modelCore->getColNames().size() > 0)
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
   int nMasterOnlyCols = static_cast<int>(modelCore->masterOnlyCols.size());

   if (nMasterOnlyCols) {
      if (m_param.LogLevel >= 1) {
         (*m_osLog) << "Create model part Master-Only." << endl;
      }

      if (m_param.LogLevel >= 2) {
         (*m_osLog) << "The number of master Only Cols is: "
                    << nMasterOnlyCols << std::endl;
         (*m_osLog) << "(kappa) The percentage of Master Only Cols over total columns is "
                    << double(nMasterOnlyCols) / nCols << std::endl;
      }

      createModelMasterOnlys(modelCore->masterOnlyCols);
   }

   //---
   //--- free up local memory
   //---
   UTIL_DELARR(rowsMarker);
   UTIL_DELARR(rowsCore);
   UTIL_DELARR(colMarker);
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "createModels()", m_param.LogLevel, 2);
   //exit(1);
}

/*
int DecompApp::generateInitVars(DecompVarList & initVars){
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "generateInitVars()", m_param.LogLevel, 2);

   readInitSolutionFile(initVars);

   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "generateInitVars()", m_param.LogLevel, 2);
   return static_cast<int>(initVars.size());
}

*/

void DecompApp::singlyBorderStructureDetection()
{
   //======================================================================
   // Using Row-net hypergraph model for automatic matrix decomposition
   //======================================================================
   int numRows = 0;
   int numCols = 0;
   int numElements = 0;

   if (m_param.InstanceFormat == "MPS") {
      numRows = m_mpsIO.getNumRows();
      numCols = m_mpsIO.getNumCols();
      numElements = m_mpsIO.getNumElements();
      m_matrix = m_mpsIO.getMatrixByRow();
   } else if (m_param.InstanceFormat == "LP") {
      numRows = m_mpsIO.getNumRows();
      numCols = m_mpsIO.getNumCols();
      numElements = m_mpsIO.getNumElements();
      m_matrix = m_mpsIO.getMatrixByRow();
   }

   // get the column/row index for by-row matrix
   const int* minorIndex = m_matrix->getIndices();
   const int* majorIndex = m_matrix->getMajorIndices();
   UTIL_MSG(m_param.LogLevel, 2,
            (*m_osLog)
            << "The number of rows is " << numRows << "\n"
            << "The number of columns is " << numCols << "\n"
            << "The number of elements is " << numElements
            << "\n";
           );
   // Here assume the matrix to be partitioned into singly-bordered
   // diagonal matrix. ToDo: complete the Doubly-bordered diagonal
   // matrix mapping
   // the number of vertices
   int numVertices ;
   // The number of hyperedges
   int numHyperedges;
#ifdef doublyBordered
   numVertices = numElements;
   numHyperedges = numRows + numCols;
#else
   // default singlyBordered
   numVertices = numElements ;
   numHyperedges = numRows;
#endif
   /*
     throw UtilException("Please provide a valid partitioning model" ,
   		 "initializeApp", "DecompApp");
   */
   //======================================================================
   // The code below  prepares the input parameters of the hypergraph
   // partitioning. for details , please refer to
   //
   //  http://glaros.dtc.umn.edu/gkhome/fetch/sw/hmetis/manual.pdf
   //
   //======================================================================
   // pointer to hyperedge
   int* eptr = new int[numHyperedges + 1];
   // the length vector for each row (number of nonzero elements
   // in each row)
   const int* lengthRows = m_matrix->getVectorLengths();
   /*
     assigning the pointer to hyperedges, indicating the number of
     vertices in each hyperedge
   */
   eptr[0] = 0;

   for (int k = 1 ; k < numHyperedges + 1; k ++ ) {
      eptr[k] = eptr[k - 1] + lengthRows[k - 1];
   }

   assert(eptr[numHyperedges] == numVertices);
   // declaring the hyperedges, which correspond to the rows in the matrix
   int* eind = new int[numVertices];

   for (int i = 0; i < numVertices; i ++) {
      eind[i] = minorIndex[i];
   }

   // declaring the number of partitions
   int nparts = m_param.NumBlocks;
   // weights of vertices and hyperedges
   int* vwgts = new int[numVertices] ;
   int* hewgts = new int[numHyperedges] ;
   /*
    * declare boolean variables indicating whether the nonzero elements
    * are integer or not (corresponding to the vertices), the row in
    * the matrix has integer columns (corresponding to the hyperedges)
    * or not
    */
   bool* intVertices = new bool[numVertices];
   bool* intHyperedges = new bool[numHyperedges];

   // initialization
   for (int i = 0 ; i < numHyperedges; i ++) {
      intHyperedges[i] = false ;
   }

   // the index of the original consstraint matrix
   int index_base = 0 ;
   int index = 0;
   // counter of integer variables
   int intCounter = 0 ;
   // vector containing the number of integer elements in each row
   int* intLengthRows = new int[numRows];

   for (int i = 0 ; i < numRows ; i ++) {
      intLengthRows[i] = 0 ;
   }

   bool isInteger = false;

   for (int i = 0; i < numRows ; i ++) {
      intCounter = 0;

      for ( int j = 0 ; j < lengthRows[i] ; j ++ ) {
         index = index_base + j ;
         // determine whether the corresponding column is
         // integer or not

         if (m_param.InstanceFormat == "MPS") {
            isInteger = m_mpsIO.isInteger(minorIndex[index]);
         } else if (m_param.InstanceFormat == "LP") {
            isInteger = m_lpIO.isInteger(minorIndex[index]);
         }

         if (isInteger) {
            intVertices[index] = true;
            intHyperedges[majorIndex[index]] = true;
            intCounter ++;
         } else {
            intVertices[index] = false;
         }
      }

      index_base = index_base + lengthRows[i];
      intLengthRows[i] = intCounter ;
   }

   /*
    *  define the weight parameter in the hypergraph
    */

   // assign the weights on vertices

   for (int i = 0 ; i < numVertices ; i ++) {
#ifdef VARIABLE_WEIGHT

      if (intVertices[i]) {
         vwgts[i] = 2 ;
      } else {
         vwgts[i] = 1;
      }

#else
      vwgts[i] = 1;
#endif
   }

   // assign the weights on hyperedges

   for (int i = 0 ; i < numHyperedges; i ++) {
#ifdef VARIABLE_WEIGHT

      if (intHyperedges[i]) {
         hewgts[i] = 2 * lengthRows[i];
      } else {
         hewgts[i] = 1;
      }

#else
      hewgts[i] = 1;
#endif
   }

   // part is an array of size nvtxs that returns the computed partition
   int* part = new int[numElements];
   // edgecut is the number of hyperedge cut
   int* edgecut = new int[1];
   edgecut[0] = 0 ;
   int* options = new int[1];
   // 0 indicates the default paraemter value, 1 otherwise;
   options[0] = 0 ;
   int* partweights = new int[nparts];

   // initialization
   for (int i = 0; i < nparts ; i ++) {
      partweights[i] = 1;
   }

   // calling HMETIS_PartKway API to perform the hypergraph partitioning
#ifdef PaToH
   PaToH_Parameters args;
   args._k = nparts;
   //  PaToH_Initialize_Parameters(&args, PATOH_CUTPART, PATOH_SUGPARAM_QUALITY);
   PaToH_Initialize_Parameters(&args, PATOH_CONPART,
                               PATOH_SUGPARAM_DEFAULT);
   // the number of constraint in the multilevel algorithm
   int nconst = 1; // single constraint
   PaToH_Alloc(&args, numVertices, numHyperedges, nconst,
               vwgts, hewgts, eptr, eind);
   int cut = 0;
   PaToH_Part(&args, numVertices, numHyperedges, nconst , 0 , vwgts,
              hewgts, eptr, eind, NULL, part, partweights, &cut);
   edgecut[0] = cut ;
   int computedCut = PaToH_Compute_Cut(nparts, PATOH_CONPART, numVertices,
                                       numHyperedges, hewgts, eptr, eind,
                                       part);
   UTIL_MSG(m_param.LogDebugLevel, 2,
            (*m_osLog)
            << "The computedCut is "
            << computedCut << "\n";
           );
#else
   clock_t begin = clock();
#if defined(COIN_HAS_HMETIS)
   // maximum load imbalance (%)
   int ubfactor = 5;
   HMETIS_PartRecursive(numVertices, numHyperedges, vwgts, eptr,
                        eind, hewgts, nparts, ubfactor, options, part, edgecut);
#endif
   clock_t end = clock();
   UTIL_MSG(m_param.LogLevel, 2,
            (*m_osLog)
            << "********************************************" << "\n"
            << "The time elapse for hypergraph partitioning is "
            << (end - begin) / CLOCKS_PER_SEC << " seconds" << "\n"
            << "********************************************"
            << "\n";
           );
#endif
   /*
    *The following codes try to find the hyperedges in the
    * vertex separator set by traversing the hyperedges.
    *If a hyperedge has vertices in more than one part,
    *then it is a cut hyperedge then it is in separator.
    */
   /*
    * define a set to store the coupling rows (hyperedges in
    * the vertex separator set)
    */
   std::set<int> netSet;
   std::set<int> :: iterator netIter;
   // initilizations for global index
   index = 0 ;
   index_base = 0 ;
   int tempBase = 0 ;

   /*
    * Identify the coupling row in the matrix by storing
    * them in a net set
    */

   for ( int i = 0 ; i < numRows ; i ++) {
      for ( int j = 0 ; j < lengthRows[i] ; j ++ ) {
         index = index_base + j ;

         if ( j == 0) {
            tempBase = part[minorIndex[index_base]];
         } else {
            if (tempBase != part[minorIndex[index]]) {
               netSet.insert(i);
               j = lengthRows[i];
            }
         }
      }

      //update index_base
      index_base = index_base + lengthRows[i];
   }

   UTIL_MSG(m_param.LogLevel, 2,
            (*m_osLog)
            << "The size of the net set after finding the coupling"
            << "row is "
            << static_cast<int>(netSet.size())
            << "\n";
           );
   // Eliminate the coupling rows from each partition set
   std::set<int> numRowIndex;
   std::set<int> :: iterator rowIter;
   std::vector<int> rowsBlock;
   /*
    * truePartNum indicates the true partition number after eliminating the partitioned
    * blocks where there is no element in the block
    */
   int truePartNum = 0 ;

   for (int part_index = 0 ; part_index < nparts; part_index ++) {
      // first, store the rows in different nets
      for (int j = 0; j < numElements; j ++) {
         if (part[minorIndex[j]] == part_index ) {
            numRowIndex.insert(majorIndex[j]);
         }
      }

      //second, temp stores the row index that is duplicating in the
      // coupling net set, and then removes it
      std::vector<int> temp;

      for (rowIter = numRowIndex.begin(); rowIter != numRowIndex.end(); rowIter ++) {
         for (netIter = netSet.begin(); netIter != netSet.end(); netIter ++) {
            if ((*rowIter) == (*netIter)) {
               temp.push_back(*rowIter);
            }
         }
      }

      for (int s = 0 ; s < static_cast<int>(temp.size()); s ++) {
         numRowIndex.erase(temp.at(s));
      }

      if (numRowIndex.size() != 0) {
         for (rowIter = numRowIndex.begin(); rowIter != numRowIndex.end();
               rowIter++) {
            rowsBlock.push_back(*rowIter);
         }

         m_blocks.insert(make_pair(truePartNum, rowsBlock));
         truePartNum ++;
      }

      numRowIndex.clear();
      rowsBlock.clear();
      temp.clear();
   }

   UTIL_DELARR(eptr);
   UTIL_DELARR(eind);
   //UTIL_DELARR(minorIndex);
   //UTIL_DELARR(lengthRows);
   UTIL_DELARR(majorIndex);
   UTIL_DELARR(intLengthRows);
   UTIL_DELARR(edgecut);
   UTIL_DELARR(part);
   UTIL_DELARR(vwgts);
   UTIL_DELARR(hewgts);
   UTIL_DELARR(options);
   UTIL_DELARR(intVertices);
   UTIL_DELARR(intHyperedges);
   UTIL_DELARR(partweights);
#ifdef PaToH
   PaToH_Free();
#endif
   std::cout << "The number of blocks is " << nparts << std::endl;
}




#if 0
// --------------------------------------------------------------------- //
void DecompApp::setOptimalPoint(vector< vector<double> >& optPoint)
{
   //
   // ---
   // --- this function does nothing by default
   // ---
   //
   UTIL_DEBUG(m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << " <---- setOptimalPoint()  ----\n";
             );
   UTIL_DEBUG(m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << "  ---- setOptimalPoint()  ---->\n";
             );
}
#endif
