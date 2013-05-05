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
// Copyright (C) 2002-2011, Lehigh University, Matthew Galati, Ted Ralphs    //
// All Rights Reserved.                                                      //
//===========================================================================//

//===========================================================================//
#include "AlpsDecompModel.h"
#include "AlpsDecompNodeDesc.h"
#include "AlpsDecompTreeNode.h"
#ifdef COIN_HAS_MPI
#include "AlpsKnowledgeBrokerMPI.h"
#else
#include "AlpsKnowledgeBrokerSerial.h"
#endif

using namespace std;

//===========================================================================//
void AlpsDecompModel::setAlpsSettings(){

   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "setAlpsSettings()", m_param.msgLevel, 3);

   AlpsPar()->setEntry(AlpsParams::logFileLevel,    m_param.logFileLevel);
   AlpsPar()->setEntry(AlpsParams::printSolution,   m_param.printSolution);
   AlpsPar()->setEntry(AlpsParams::checkMemory,     m_param.checkMemory);
   AlpsPar()->setEntry(AlpsParams::msgLevel,        m_param.msgLevel);
   AlpsPar()->setEntry(AlpsParams::nodeLimit,       m_param.nodeLimit);
   AlpsPar()->setEntry(AlpsParams::nodeLogInterval, m_param.nodeLogInterval);
   AlpsPar()->setEntry(AlpsParams::searchStrategy, m_param.searchStrategy); 

   double timeLimit = m_decompAlgo->getParam().LimitTime;
   AlpsPar()->setEntry(AlpsParams::timeLimit,       timeLimit);

   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "setAlpsSettings()", m_param.msgLevel, 3);
}


//===========================================================================//
AlpsTreeNode * AlpsDecompModel::createRoot(){
      
   //---
   //--- Create the root node description and set explicit (no diff'ing)
   //---    NOTE: Alps will delete this memory;
   //---
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "createRoot()", m_param.msgLevel, 3);

   AlpsDecompTreeNode * root = new AlpsDecompTreeNode();
   assert(root);


   CoinAssert(m_decompAlgo);
   const DecompAlgoModel & modelCore = m_decompAlgo->getModelCore();
   CoinAssert(modelCore.getModel()->getColLB());
   CoinAssert(modelCore.getModel()->getColUB());
   
   AlpsDecompNodeDesc * desc 
      = new AlpsDecompNodeDesc(this,
			       modelCore.getModel()->getColLB(),
                               modelCore.getModel()->getColUB());
   assert(desc);
   
   root->setDesc(desc);
   //root->setExplicit(1);  
   
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "setAlpsSettings()", m_param.msgLevel, 3);

   return root;
}

//===========================================================================//
bool AlpsDecompModel::fathomAllNodes() {
   double feasBound   = ALPS_OBJ_MAX;
   double relBound    = ALPS_OBJ_MAX;
   double gapVal      = ALPS_OBJ_MAX;   
   double currAbsGap_ = ALPS_OBJ_MAX;
   double currRelGap_ = ALPS_OBJ_MAX;
   AlpsTreeNode * bestNode = NULL;
    
   // Compute gap
   feasBound = broker_->getIncumbentValue();
   bestNode  = broker_->getBestNode();    

   //printf("feasBound= %12.10f\n", feasBound);

   if (bestNode){
      relBound = bestNode->getQuality();
      m_bestLB = relBound;
      //printf("bestNode m_bestLB= %12.10f\n", m_bestLB);
   }
   else{
      m_bestLB = getKnowledgeBroker()->getBestQuality();
      //printf("no bestNode m_bestLB= %12.10f\n", m_bestLB);
   }

   if (relBound > ALPS_OBJ_MAX_LESS) {
      currAbsGap_ = currRelGap_ = 0.0;
   }
   else if (feasBound < ALPS_OBJ_MAX_LESS) {
      gapVal      = ALPS_MAX(0, feasBound - relBound);
      currAbsGap_ = ALPS_MAX(0, gapVal);
      currRelGap_ = 100.0 * UtilCalculateGap(relBound, feasBound);
   }
   //printf("+++ Process %d: currAbsGap_ %g, currRelGap_%g\n",
   //     broker_->getProcRank(), currAbsGap_,  currRelGap_);
   //TODO: make option
   double optimalAbsGap_ = 1.0e-6;
   double optimalRelGap_ = 0.01;//0.01%
   //TODO: cutoffIncrement (currentUB-cutoffIncrement)
   if ( (currAbsGap_ <= optimalAbsGap_ + ALPS_ZERO) ||
        (currRelGap_ <= optimalRelGap_ + ALPS_ZERO) ) {
      m_bestLB = feasBound;
      return true;
   }
   else {
      return false;
   }
}




//===========================================================================//
AlpsExitStatus AlpsDecompModel::solve(){
   /** \todo Parallel version. */

#ifdef UTIL_USE_TIMERS
   globalTimer.reset();
#endif

   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "solve()", m_param.msgLevel, 3);

   //---
   //--- Since the setup phase for DECOMP includes generating initial
   //---   columns and creating the master problem it could be 
   //---   a significant amount of time. So, we need to adjust the 
   //---   time limit for the residual.
   //---
   DecompAlgo  * decompAlgo  = getDecompAlgo();
   DecompStats & decompStats = decompAlgo->getStats();
   DecompParam & decompParam = decompAlgo->getMutableParam();

   double timeLimit = decompParam.LimitTime;
   double timeLeft  = timeLimit - decompStats.timerOverall.getRealTime();   
   AlpsPar()->setEntry(AlpsParams::timeLimit, timeLeft);

   //---
   //--- copy relevant parameters to DecompParam from AlpsParam
   //---
   decompParam.LimitNodes = m_param.nodeLimit;
   


#ifdef COIN_HAS_MPI
   //---
   // declare an AlpsKnowledgeBroker for parallel application
   //---
   AlpsKnowledgeBrokerMPI alpsBroker(0, NULL, *this);
#else

   //---
   // declare an AlpsKnowledgeBroker for serial application
   //---
   AlpsKnowledgeBrokerSerial alpsBroker(0, NULL, *this);
#endif

   
   //---
   //--- search for the best solution 
   //---
   alpsBroker.search(this);
   if(m_param.msgLevel > 0)
      m_decompAlgo->getDecompStats().printOverallStats();

   //---
   //--- store best LB/UB objective found
   //---
   m_bestUB         = alpsBroker.getBestQuality();
   m_nodesProcessed = alpsBroker.getNumNodesProcessed();
   if(alpsBroker.getSolStatus() != AlpsExitStatusOptimal){
      AlpsTreeNode * bestNode = NULL;
      //if stops on time, have the nodes been free'd?
      bestNode = alpsBroker.getBestNode();
      if(bestNode)
	 m_bestLB = bestNode->getQuality();
      else
	 m_bestLB =-ALPS_OBJ_MAX;
   }
   m_alpsStatus = alpsBroker.getSolStatus();

   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "solve()", m_param.msgLevel, 3);

   return alpsBroker.getSolStatus();
}

void
AlpsDecompModel::encodeDecompAlgo(AlpsEncoded*& encoded) const{


} 

void
AlpsDecompModel::decodeDecompAlgo(AlpsEncoded& encoded){

  


}


void 
AlpsDecompModel::encodeAlpsDecompParam(AlpsEncoded*& encoded){


  encoded->writeRep(getParam().logFileLevel); 
  
  encoded->writeRep(getParam().printSolution); 

  encoded->writeRep(getParam().checkMemory); 

  encoded->writeRep(getParam().msgLevel); 

  encoded->writeRep(getParam().nodeLimit); 

  encoded->writeRep(getParam().nodeLogInterval); 

  encoded->writeRep(getParam().searchStrategy); 
  
}

AlpsDecompParam 
AlpsDecompModel::decodeAlpsDecompParam(AlpsEncoded& encoded){

  int logFileLevel; 

  bool printSolution; 

  bool checkMemory; 

  int msgLevel; 

  int nodeLimit; 

  int nodeLogInterval; 

  int searchStrategy; 

  
  encoded.readRep(logFileLevel); 
  
  encoded.writeRep(printSolution); 

  encoded.writeRep(checkMemory); 

  encoded.writeRep(msgLevel); 

  encoded.writeRep(nodeLimit); 

  encoded.writeRep(nodeLogInterval); 

  encoded.writeRep(searchStrategy); 

  AlpsDecompParam param(logFileLevel, printSolution,
			checkMemory, msgLevel,
			nodeLimit, nodeLogInterval,
			searchStrategy); 

  return param; 

}




//#########################################################################
// Send model and root so that initial solve
AlpsEncoded* AlpsDecompModel::encode()
{

  AlpsReturnStatus status = AlpsReturnStatusOk; 

  AlpsEncoded* encoded = new AlpsEncoded(AlpsKnowledgeTypeModel); 

  //*****************************************
  // Encode Alps part
  //*****************************************

  status = encodeAlps(encoded); 


  //*****************************************
  // Encode DecompModel part
  //*****************************************
  

  //write the model data into representation

  encoded->writeRep(m_bestLB);
  
  encoded->writeRep(m_bestUB); 

  encoded->writeRep(m_nodesProcessed);

  encoded->writeRep(m_alpsStatus);

  
  const CoinPackedMatrix* matrixByCol = m_decompAlgo->getMasterOSI()->getMatrixByCol();

  int numRows = matrixByCol->getNumRows();

  encoded->writeRep(numRows);

  int numCols = matrixByCol->getNumCols();

  encoded->writeRep(numCols); 


  const double* collb = m_decompAlgo->getMasterOSI()->getColLower();

  encoded->writeRep(collb, numCols);

  const double* colub = m_decompAlgo->getMasterOSI()->getColUpper();

  encoded->writeRep(colub, numCols);

  const double* obj = m_decompAlgo->getMasterOSI()->getObjCoefficients();

  encoded->writeRep(obj, numCols);

  const double objSense = m_decompAlgo->getMasterOSI()->getObjSense();

  encoded->writeRep(objSense);

  const double* rowlb = m_decompAlgo->getMasterOSI()->getRowLower();

  encoded->writeRep(rowlb, numRows);

  const double* rowub = m_decompAlgo->getMasterOSI()->getRowUpper();

  encoded->writeRep(rowub, numRows);

  int numElements = m_decompAlgo->getMasterOSI()->getNumElements();

  encoded->writeRep(numElements);

  const double* elementValue = matrixByCol->getElements();

  encoded->writeRep(elementValue, numElements);

  const CoinBigIndex* colStart = matrixByCol->getVectorStarts();

  int numStart = numCols + 1;

  encoded->writeRep(colStart, numStart);

  const int* index = matrixByCol->getIndices();

  encoded->writeRep(index, numElements);
  
  int numberIntegers = 0;


  int* integerVariable = new int[numElements];
  for (int iColumn=0; iColumn<numCols; iColumn++) {
    if( m_decompAlgo->getMasterOSI()->isInteger(iColumn))
      integerVariable[numberIntegers++]=iColumn;
    
  }
  if (!numberIntegers){
    integerVariable = NULL;
  }

  encoded->writeRep(numberIntegers);
  encoded->writeRep(integerVariable, numberIntegers);
  

  delete [] integerVariable; 


  encodeDecompAlgo(encoded);

  encodeAlpsDecompParam(encoded);
  
  return encoded;

}


//#########################################################################
//

void
AlpsDecompModel::decodeToSelf(AlpsEncoded& encoded)
{

  AlpsReturnStatus status = AlpsReturnStatusOk;

  double m_bLB; 
  double m_bUB; 
  int m_nodesP; 
  int m_alpsStat; 

  encoded.readRep(m_bLB); 
  encoded.readRep(m_bUB); 
  encoded.readRep(m_nodesP);
  encoded.readRep(m_alpsStat);

  decodeDecompAlgo(encoded); 

  decodeAlpsDecompParam(encoded); 
  

}
