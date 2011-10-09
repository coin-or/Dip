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

//===========================================================================//
#include "UtilMacros.h"

//===========================================================================//
#include "BcpsDecompModel.h"
#include "BcpsDecompTreeNode.h"
#include "BcpsDecompNodeDesc.h"

//===========================================================================//
using namespace std;

//===========================================================================//
const char * BcpsDecompModel::m_classTag         = "\nBCPSD         : ";

#define LOGLEVEL 10

//===========================================================================//
//=== pure virtual methods from BcpsModel or AlpsModel ======================//
//===========================================================================//

//===========================================================================//
void BcpsDecompModel::readInstance(const char * dataFile){
  //TODO: why is readInstance virtual?
  //TODO: log levels
  UTIL_DEBUG(LOGLEVEL, 3,
             cout << m_classTag << " <---- readInstance() ---- ";
             );
  
  UTIL_DEBUG(LOGLEVEL, 3,
             cout << m_classTag << "  ---- readInstance() ----> ";
             );
}

//===========================================================================//
AlpsTreeNode * BcpsDecompModel::createRoot(){

  //---
  //--- NOTE: Root will be deleted by ALPS. Root is an explicit node.
  //---
  BcpsDecompTreeNode * root = new BcpsDecompTreeNode();
  BcpsDecompNodeDesc * desc 
    = new BcpsDecompNodeDesc(this,
                             &decompAlgo_->m_modelCore->colLB[0],
                             &decompAlgo_->m_modelCore->colUB[0]);
  root->setDesc(desc);

  //??

  //STOP
#if 0
  BlisTreeNode* root = new BlisTreeNode;
  BlisNodeDesc* desc = new BlisNodeDesc(this);
  root->setDesc(desc);

  //-------------------------------------------------------------
  // NOTE: Although original data are stored in model when reading. 
  //   Root desc still store a full copy of col and row bounds when creating.
  //   It will store soft differences after finding a branching object. 
  //   The soft difference are due to reduced cost fixing and probing.
  //   Also the added cols and rows will be stored.
  //-------------------------------------------------------------
  int k;

  BcpsVariable ** vars = getCoreVariables();
  BcpsConstraint ** cons = getCoreConstraints();

  int *varIndices1 = new int [numCoreVariables_];
  int *varIndices2 = new int [numCoreVariables_];
  int *varIndices3 = NULL; //new int [numCoreVariables_];
  int *varIndices4 = NULL; //new int [numCoreVariables_];
  double *vlhe = new double [numCoreVariables_];
  double *vuhe = new double [numCoreVariables_];
  double *vlse = NULL; //new double [numCoreVariables_];
  double *vuse = NULL; //new double [numCoreVariables_];

  int *conIndices1 = new int [numCoreConstraints_];
  int *conIndices2 = new int [numCoreConstraints_];
  int *conIndices3 = NULL; //new int [numCoreConstraints_];
  int *conIndices4 = NULL; //new int [numCoreConstraints_];
  double *clhe = new double [numCoreConstraints_];
  double *cuhe = new double [numCoreConstraints_];
  double *clse = NULL; //new double [numCoreConstraints_];
  double *cuse = NULL; //new double [numCoreConstraints_];

  //-------------------------------------------------------------
  // Get var bounds and indices.
  //-------------------------------------------------------------

  for (k = 0; k < numCoreVariables_; ++k) {
    vlhe[k] = vars[k]->getLbHard();
    vuhe[k] = vars[k]->getUbHard();
    //vlse[k] = vars[k]->getLbSoft();
    //vuse[k] = vars[k]->getUbSoft();
    varIndices1[k] = k;
    varIndices2[k] = k;
    
    //varIndices3[k] = k;
    //varIndices4[k] = k;

#ifdef BLIS_DEBUG_MORE
    std::cout << "BLIS: createRoot(): var "<< k << ": hard: lb=" << vlhe[k]
	      << ", ub=" << vuhe[k] << std::endl;
#endif  
    
  }

  //-------------------------------------------------------------  
  // Get con bounds and indices.
  //-------------------------------------------------------------

  for (k = 0; k < numCoreConstraints_; ++k) {
    clhe[k] = cons[k]->getLbHard();
    cuhe[k] = cons[k]->getUbHard();
    //clse[k] = cons[k]->getLbSoft();
    //cuse[k] = cons[k]->getUbSoft();
    conIndices1[k] = k;
    conIndices2[k] = k;
    //conIndices3[k] = k;
    //conIndices4[k] = k;
  }

  int *tempInd = NULL;
  BcpsObject **tempObj = NULL;
  
  desc->assignVars(0 /*numRem*/, tempInd,
		   0 /*numAdd*/, tempObj,
		   false, numCoreVariables_, varIndices1, vlhe, /*Var hard lb*/
		   false, numCoreVariables_, varIndices2, vuhe, /*Var hard ub*/
		   false, 0, varIndices3, vlse, /*Var soft lb*/
		   false, 0, varIndices4, vuse);/*Var soft ub*/
  desc->assignCons(0 /*numRem*/, tempInd,
		   0 /*numAdd*/, tempObj,
		   false, numCoreConstraints_,conIndices1,clhe, /*Con hard lb*/
		   false, numCoreConstraints_,conIndices2,cuhe, /*Con hard ub*/
		   false, 0,conIndices3,clse, /*Con soft lb*/
		   false, 0,conIndices4,cuse);/*Con soft ub*/

  //-------------------------------------------------------------  
  // Mark it as an explicit node.
  //-------------------------------------------------------------


#endif  
  root->setExplicit(1);
  return root;
}

//===========================================================================//
//=== setup functions =======================================================//
//===========================================================================//
void BcpsDecompModel::init(){
  UTIL_DEBUG(LOGLEVEL, 3,
             cout << m_classTag << " <---- init()         ---- ";
             );
  
  UTIL_DEBUG(LOGLEVEL, 3,
             cout << m_classTag << "  ---- init()         ----> ";
             );
}
