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

//===========================================================================//
#include "AlpsDecompModel.h"
#include "AlpsDecompNodeDesc.h"
#include "AlpsDecompTreeNode.h"
#include "AlpsDecompSolution.h"

#ifdef COIN_HAS_MPI
#include "AlpsKnowledgeBrokerMPI.h"
#else
#include "AlpsKnowledgeBrokerSerial.h"
#endif

using namespace std;

//===========================================================================//
void AlpsDecompModel::setAlpsSettings()
{
   //TODO: use stream not cout
   UtilPrintFuncBegin(&cout, m_classTag,
                      "setAlpsSettings()", m_param.msgLevel, 3);
   AlpsPar()->setEntry(AlpsParams::logFileLevel,    m_param.logFileLevel);
   AlpsPar()->setEntry(AlpsParams::printSolution,   m_param.printSolution);
   AlpsPar()->setEntry(AlpsParams::checkMemory,     m_param.checkMemory);
   AlpsPar()->setEntry(AlpsParams::msgLevel,        m_param.msgLevel);
   AlpsPar()->setEntry(AlpsParams::nodeLimit,       m_param.nodeLimit);
   AlpsPar()->setEntry(AlpsParams::nodeLogInterval, m_param.nodeLogInterval);
   double timeLimit = m_decompAlgo->getParam().LimitTime;
   AlpsPar()->setEntry(AlpsParams::timeLimit,       timeLimit);
   UtilPrintFuncEnd(&cout, m_classTag,
                    "setAlpsSettings()", m_param.msgLevel, 3);
}


//===========================================================================//
AlpsTreeNode* AlpsDecompModel::createRoot()
{
   //---
   //--- Create the root node description and set explicit (no diff'ing)
   //---    NOTE: Alps will delete this memory;
   //---
   UtilPrintFuncBegin(&cout, m_classTag,
                      "createRoot()", m_param.msgLevel, 3);
   AlpsDecompTreeNode* root = new AlpsDecompTreeNode();
   assert(root);
   CoinAssert(m_decompAlgo);
   const DecompAlgoModel& modelCore = m_decompAlgo->getModelCore();
   CoinAssert(modelCore.getModel()->getColLB());
   CoinAssert(modelCore.getModel()->getColUB());
   AlpsDecompNodeDesc* desc
   = new AlpsDecompNodeDesc(this,
                            modelCore.getModel()->getColLB(),
                            modelCore.getModel()->getColUB());
   assert(desc);
   root->setDesc(desc);
   //root->setExplicit(1);
   UtilPrintFuncEnd(&cout, m_classTag,
                    "setAlpsSettings()", m_param.msgLevel, 3);
   return root;
}

//===========================================================================//
bool AlpsDecompModel::fathomAllNodes()
{
   double feasBound   = ALPS_OBJ_MAX;
   double relBound    = ALPS_OBJ_MAX;
   double gapVal      = ALPS_OBJ_MAX;
   double currAbsGap_ = ALPS_OBJ_MAX;
   double currRelGap_ = ALPS_OBJ_MAX;
   AlpsTreeNode* bestNode = NULL;
   // Compute gap
   feasBound = broker_->getIncumbentValue();
   bestNode  = broker_->getBestNode();

   //printf("feasBound= %12.10f\n", feasBound);
   if (bestNode) {
      relBound = bestNode->getQuality();
      m_bestLB = relBound;
      //printf("bestNode m_bestLB= %12.10f\n", m_bestLB);
   } else {
      m_bestLB = getKnowledgeBroker()->getBestQuality();
      //printf("no bestNode m_bestLB= %12.10f\n", m_bestLB);
   }

   if (relBound > ALPS_OBJ_MAX_LESS) {
      currAbsGap_ = currRelGap_ = 0.0;
   } else if (feasBound < ALPS_OBJ_MAX_LESS) {
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
   } else {
      return false;
   }
}




//===========================================================================//
AlpsExitStatus AlpsDecompModel::solve(int argc, char* argv[])
{
   /** \todo Parallel version. */
#ifdef UTIL_USE_TIMERS
   globalTimer.reset();
#endif
   UtilPrintFuncBegin(&cout, m_classTag,
                      "solve()", m_param.msgLevel, 3);
   //---
   //--- Since the setup phase for DECOMP includes generating initial
   //---   columns and creating the master problem it could be
   //---   a significant amount of time. So, we need to adjust the
   //---   time limit for the residual.
   //---
   DecompAlgo*   decompAlgo  = getDecompAlgo();
   DecompStats& decompStats = decompAlgo->getStats();
   DecompParam& decompParam = decompAlgo->getMutableParam();
   double timeLimit = decompParam.LimitTime;
   double timeLeft  = timeLimit - decompStats.timerOverall.getRealTime();
   AlpsPar()->setEntry(AlpsParams::timeLimit, timeLeft);
   //---
   //--- copy relevant parameters to DecompParam from AlpsParam
   //---
   decompParam.LimitNodes = m_param.nodeLimit;
#ifdef COIN_HAS_MPI
   //---
   //--- declare an AlpsKnowledgeBroker for paralell application
   //---
   AlpsKnowledgeBrokerMPI alpsBroker(argc, argv, *this);
#else
   //---
   //--- declare an AlpsKnowledgeBroker for serial application
   //---
   AlpsKnowledgeBrokerSerial alpsBroker(0, NULL, *this);
#endif
   //---
   //--- search for the best solution
   //---
   alpsBroker.search(this);

   if (m_param.msgLevel > 0) {
      m_decompAlgo->getDecompStats().printOverallStats();
   }

   //---
   //--- store best LB/UB objective found
   //---
   m_bestUB         = alpsBroker.getBestQuality();
   m_nodesProcessed = alpsBroker.getNumNodesProcessed();

   if (alpsBroker.getSolStatus() != AlpsExitStatusOptimal) {
      AlpsTreeNode* bestNode = NULL;
      //if stops on time, have the nodes been free'd?
      bestNode = alpsBroker.getBestNode();

      if (bestNode) {
         m_bestLB = bestNode->getQuality();
      } else {
         m_bestLB = -ALPS_OBJ_MAX;
      }
   }

   m_alpsStatus = alpsBroker.getSolStatus();
   UtilPrintFuncEnd(&cout, m_classTag,
                    "solve()", m_param.msgLevel, 3);
   return alpsBroker.getSolStatus();
}

AlpsEncoded*
AlpsDecompModel::encode() const
{
   AlpsReturnStatus status = AlpsReturnStatusOk;
   // NOTE: "AlpsKnowledgeTypeModel" is the type name.
   AlpsEncoded* encoded = new AlpsEncoded(AlpsKnowledgeTypeModel);
   status = encodeAlps(encoded);
   status = encodeAlpsDecomp(encoded);
   return encoded;
}

void
AlpsDecompModel::decodeToSelf(AlpsEncoded& encoded)
{
   AlpsReturnStatus status = AlpsReturnStatusOk;
   status = decodeAlps(encoded);
   status = decodeAlpsDecomp(encoded);
}

AlpsReturnStatus
AlpsDecompModel::encodeAlpsDecomp(AlpsEncoded* encoded) const
{
   AlpsReturnStatus status = AlpsReturnStatusOk;
   // write the model data into representation
   encoded->writeRep(m_bestLB);
   encoded->writeRep(m_bestUB);
   encoded->writeRep(m_nodesProcessed);
   encoded->writeRep(m_alpsStatus);
   //   m_param.pack(encoded);
   return status;
}

AlpsReturnStatus
AlpsDecompModel::decodeAlpsDecomp(AlpsEncoded& encoded)
{
   AlpsReturnStatus status = AlpsReturnStatusOk;
   double d_m_bestLB, d_m_bestUB;
   int d_m_nodesProcessed, d_m_alpsStatus;
   encoded.readRep(d_m_bestLB);
   encoded.readRep(d_m_bestUB);
   encoded.readRep(d_m_nodesProcessed);
   encoded.readRep(d_m_alpsStatus);
   //   m_param.unpack(encoded);
   m_bestLB = d_m_bestLB;
   m_bestUB = d_m_bestUB;
   m_nodesProcessed = d_m_nodesProcessed;
   m_alpsStatus = d_m_alpsStatus;
   return status;
}


AlpsEncoded*
AlpsDecompModel::packSharedKnowlege()
{
   return NULL;
}

void
AlpsDecompModel::unpackSharedKnowledge(AlpsEncoded&)
{
}


/** Register knowledge. */
void
AlpsDecompModel::registerKnowledge()
{
   // Register model, solution, and tree node
   assert(broker_);
   broker_->registerClass(AlpsKnowledgeTypeModel, new AlpsDecompModel);

   if (broker_->getMsgLevel() > 100) {
      std::cout << "AlpsDecomp: Register Alps model." << std::endl;
   }

   broker_->registerClass(AlpsKnowledgeTypeNode, new AlpsDecompTreeNode(this));

   if (broker_->getMsgLevel() > 100) {
      std::cout << "AlpsDecomp: Register Alps node." << std::endl;
   }

   broker_->registerClass(AlpsKnowledgeTypeSolution, new AlpsDecompSolution());

   if (broker_->getMsgLevel() > 100) {
      std::cout << "AlpsDecomp: Register Alps solution." << std::endl;
   }
}
