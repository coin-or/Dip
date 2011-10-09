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

/**
 * \todo cleanup decomp status codes - scope them like Alps
 */

//===========================================================================//
#include "DecompAlgo.h"
#include "DecompApp.h"

//===========================================================================//
#include "AlpsKnowledgeBroker.h"
#include "AlpsDecompTreeNode.h"
#include "AlpsDecompNodeDesc.h"
#include "AlpsDecompSolution.h"
#include "AlpsDecompModel.h"

//===========================================================================//
#include "CoinUtility.hpp"

using namespace std;

//===========================================================================//
AlpsTreeNode * 
AlpsDecompTreeNode::createNewTreeNode(AlpsNodeDesc *& desc) const {

   //---
   //--- Create a new tree node, set node description.
   //---    NOTE: we are not using differencing, constructs node from scratch
   //---
   AlpsDecompModel    * model 
      = dynamic_cast<AlpsDecompModel*>(desc->getModel());
   AlpsDecompParam    & param = model->getParam();
   
   UtilPrintFuncBegin(&cout, m_classTag,
                      "createNewTreeNode()", param.msgLevel, 3);   
   AlpsDecompTreeNode * node = new AlpsDecompTreeNode();
   node->desc_ = desc;
   UtilPrintFuncEnd(&cout, m_classTag,
                    "createNewTreeNode()", param.msgLevel, 3);
   return node;  
}

//===========================================================================//
void AlpsDecompTreeNode::checkIncumbent(AlpsDecompModel      * model,
                                        const DecompSolution * decompSol){
   
   AlpsDecompParam    & param      = model->getParam();
   DecompAlgo         * decompAlgo = model->getDecompAlgo();
   //---
   //--- decompAlgo found an IP (and user) feasible point
   //---	
   double currentUB   = getKnowledgeBroker()->getIncumbentValue(); 
   double candidateUB = decompSol->getQuality();
   UTIL_DEBUG(param.msgLevel, 3,
              cout
              << "DecompAlgo found IP incum = " 
              << UtilDblToStr(candidateUB)
              << " currentUB " << UtilDblToStr(currentUB) << endl;
              );
   if(candidateUB < currentUB){
      //---
      //--- create a new solution and add to alps knowledge 
      //---      
      AlpsDecompSolution * alpsDecompSol = 
         new AlpsDecompSolution(decompSol->getSize(),
                                decompSol->getValues(),
                                decompSol->getQuality(),
                                decompAlgo->getDecompApp());	    
      getKnowledgeBroker()->addKnowledge(AlpsKnowledgeTypeSolution, 
                                         alpsDecompSol, 
                                         candidateUB);

      //---
      //--- print the new solution (if debugging)
      //---
      const DecompApp            * app         
         = decompAlgo->getDecompApp();
      const DecompConstraintSet  * modelCore 
         = decompAlgo->getModelCore().getModel();
      UTIL_DEBUG(param.msgLevel, 3,
                 app->printOriginalSolution(decompSol->getSize(),
                                            modelCore->getColNames(),
                                            decompSol->getValues());
                 );
   }
}

//===========================================================================//
int AlpsDecompTreeNode::process(bool isRoot,
				bool rampUp){

   //---
   //--- get pointer / reference to model, node description, decomp algo(s)
   //---
   AlpsDecompNodeDesc * desc 
      = dynamic_cast<AlpsDecompNodeDesc*>(desc_);
   AlpsDecompModel    * model 
      = dynamic_cast<AlpsDecompModel*>(desc->getModel());
   AlpsDecompParam    & param      = model->getParam();
   DecompAlgo         * decompAlgo = model->getDecompAlgo();
   CoinAssertDebug(desc && model);

   UtilPrintFuncBegin(&cout, m_classTag,
                      "process()", param.msgLevel, 3);
   UTIL_DEBUG(param.msgLevel, 3,
	      cout 
	      << "Start process of node: " << getIndex()
	      << " (parent = " << getParentIndex() << ")" << endl;
	      );
   
   int            status       = AlpsReturnStatusOk;
   bool           doFathom     = false;
   DecompStatus   decompStatus = STAT_FEASIBLE;
   double         relTolerance = 0.0001; //0.01% means optimal (make param)
   double         gap;

   //---
   //--- check if this can be fathomed based on parent by objective cutoff
   //---          
   double currentUB       = getKnowledgeBroker()->getIncumbentValue();
   double parentObjValue  = getQuality();   
   double primalTolerance = 1.0e-6;
   double globalLB        = -DecompInf;
   double globalUB        =  DecompInf;
   double thisQuality;

   AlpsTreeNode        * bestNode  = NULL;
   const double        * lbs       = desc->lowerBounds_;
   const double        * ubs       = desc->upperBounds_;      
   const DecompApp     * app       = decompAlgo->getDecompApp();
   DecompConstraintSet * modelCore = decompAlgo->getModelCore().getModel();   
   const int             n_cols    = modelCore->getNumCols();


   /** \todo get primalTolerance from parameter */
   if ((parentObjValue - primalTolerance) > currentUB) {
      doFathom = true;
      UTIL_DEBUG(param.msgLevel, 3,
                 cout << "Fathom since parentObjValue=" 
                 << setw(10) << UtilDblToStr(parentObjValue)
                 << " currentUB = " << setw(10) << UtilDblToStr(currentUB) << endl;
                 );
      goto TERM_PROCESS;
   }
   
   //---
   //--- the destructor initializes quality_ = infinity
   //---   we really want -infinity
   //---
   if(isRoot)
     quality_ = -ALPS_OBJ_MAX;
   
   //---
   //--- reset user-currentUB (if none given, this will have no effect)
   //---
   decompAlgo->setObjBoundIP(decompAlgo->getCutoffUB());

   if(!isRoot){
      //---
      //--- set the master column bounds (for this node in tree)
      //---
      
      //---
      //--- for debugging, print column bounds that differ from original
      //---
      UTIL_MSG(param.msgLevel, 3,
	       int              c;
	       double           diffLB;
	       double           diffUB;
	       vector<double> & colLBCore = modelCore->colLB; 
	       vector<double> & colUBCore = modelCore->colUB;
	       for(c = 0; c < n_cols; c++){
		  diffLB = lbs[c] - colLBCore[c];
		  diffUB = ubs[c] - colUBCore[c]; 
		  if(!UtilIsZero(diffLB) || !UtilIsZero(diffUB)){    
		     cout << "bound-diffs c: " << c << " -> ";
		     app->printOriginalColumn(c, &cout);
		     cout << "\t(lb,ub): (" << colLBCore[c] << ","
			  << colUBCore[c] << ")\t->\t(" << lbs[c]
			  << "," << ubs[c] << ")" << endl;
		  }
		  //else{
		  // cout << "\t(lb,ub): (" << colLBCore[c] << ","
		  //  << colUBCore[c] << ")\t->\t(" << lbs[c]
		  //  << "," << ubs[c] << ")" << endl;
		  //}
	       }
	       );   
      decompAlgo->setMasterBounds(lbs, ubs);
      decompAlgo->setSubProbBounds(lbs,ubs);
   }
   else{
      //---
      //--- check to see if we got lucky in generating init vars
      //---
      if(decompAlgo->getXhatIPBest())
         checkIncumbent(model, decompAlgo->getXhatIPBest());
      
      //---
      //--- This is a first attempt at a redesign of branching rows.
      //---  We still have all of them explicitly defined, but we
      //---  relax them and only explicitly enforce things as we branch.
      //---
      //--- A more advanced attempt would treat branching rows as cuts
      //---  and add them dynamically.
      //---
      //--- In root node, set all branching rows to "free" by relaxing
      //---   lb and ub.
      //---
      //--- NOTE: this should also be done for all nodes except for the
      //---   rows that represent bounds that we have branched on.
      //---      
      if(decompAlgo->getAlgo() == PRICE_AND_CUT){
	 int      c;
	 double * lbsInf = new double[n_cols];
	 double * ubsInf = new double[n_cols];
	 for(c = 0; c < n_cols; c++){
	    lbsInf[c] = -DecompInf;
	    ubsInf[c] =  DecompInf;
	    //printf("root c:%d lb=%g ub=%g\n",
	    //   c, lbs[c], ubs[c]);
	 } 
         decompAlgo->setMasterBounds(lbsInf, ubsInf);
	 UTIL_DELARR(lbsInf);
	 UTIL_DELARR(ubsInf);
	 
	 //actually, don't need to do this - these should already be set
	 decompAlgo->setSubProbBounds(lbs,ubs);
      }
   }
   
      
   //---
   //--- update the currentUB value for decomp algo
   //---
   currentUB = getKnowledgeBroker()->getIncumbentValue();
   decompAlgo->setObjBoundIP(currentUB);//??

   gap      = DecompInf;
   globalUB = getKnowledgeBroker()->getIncumbentValue();
   if(!isRoot){
      bestNode = getKnowledgeBroker()->getBestNode();
      globalLB = bestNode->getQuality();

      //---
      //--- if the overall gap is tight enough, fathom whatever is left
      //---
      gap = UtilCalculateGap(globalLB, globalUB);
      if(gap <= relTolerance){	 
	 doFathom = true;
	 UTIL_MSG(param.msgLevel, 3,
		    cout << "Fathom Node " << getIndex() << " since globalLB= "
		    << setw(10) << UtilDblToStr(globalLB) 		    
		    << " globalUB = "  << setw(10) << UtilDblToStr(globalUB)
		    << " gap = "     << setw(10) << UtilDblToStr(gap) << endl;
		    );
	 goto TERM_PROCESS;
      }      
   }
   







   
   //---
   //--- solve the bounding problem (DecompAlgo)
   //---
   decompStatus = decompAlgo->processNode(getIndex(), globalLB, globalUB);
      
   //---
   //--- during processNode, did we find any IP feasible points?
   //--- 
   if(decompAlgo->getXhatIPBest()){
      checkIncumbent(model, decompAlgo->getXhatIPBest());
      //---
      //--- update the local currentUB value and the decomp global UB
      //---
      currentUB = getKnowledgeBroker()->getIncumbentValue();
      decompAlgo->setObjBoundIP(currentUB);
   }
   
   switch(decompStatus){
   case STAT_FEASIBLE:
      //---
      //--- the relaxation is feasible
      //---   if the new bound is > current currentUB, fathom 
      //---   else                                , branch
      //---
      thisQuality = decompAlgo->getObjBestBoundLB();           //LB (min)
      currentUB      = getKnowledgeBroker()->getIncumbentValue(); //UB (min)
      if(thisQuality > quality_)
         quality_ = thisQuality;

      //watch tolerance here... if quality is close enough, fathom it
      gap = UtilCalculateGap(thisQuality, currentUB);
      //if(gap <= relTolerance){
      if (quality_ >= currentUB){
         doFathom = true;
         UTIL_DEBUG(param.msgLevel, 3,
                    cout << "Fathom since thisQuality= "
                    << setw(10) << UtilDblToStr(thisQuality) 
                    << " quality_= " << setw(10) << UtilDblToStr(quality_)
                    << " currentUB = "  << setw(10) << UtilDblToStr(currentUB)
                    << " gap = "     << setw(10) << UtilDblToStr(gap) << endl;
                    );
      }      
      UTIL_MSG(param.msgLevel, 3,
               cout << "Node " << getIndex()
               << " quality " << UtilDblToStr(quality_)
               << " currentUB "  << UtilDblToStr(currentUB)
               << " doFathom " << doFathom << endl;
               );      
      break;
   case STAT_INFEASIBLE:
      //---
      //--- the relaxation is infeasible, fathom 
      //---
      thisQuality = -ALPS_OBJ_MAX;
      doFathom    = true;	 
      UTIL_MSG(param.msgLevel, 3,
               cout << "Fathom since node infeasible\n";
               );
      break;
   default:
      assert(0);
   }

   //TODO: control by decomp log level?
   UTIL_MSG(param.msgLevel, 3,	
	    cout << "Node " << getIndex()
	    << " bestQuality "  << UtilDblToStr(quality_)
	    << " bestFeasible " << UtilDblToStr(currentUB) << endl;
	    );
   
 TERM_PROCESS:
   //STOP: if do fathom when node limit hit, then it gives wrong LB
   //  what is the proper status setting if node limit is hit to stop
   //  but not fathom so as to lose the proper bound
   //if(param.nodeLimit == 0)
   //   status = AlpsExitStatusNodeLimit;

   //---
   //--- for nodeLimit == 0, we do not want it to look for 
   //---   branching candidates since in some cases we stop due to 
   //---   gap without a branching candidate and do not want to have to 
   //---   return (since we are not evaluating any more nodes anyway)
   //--- so, we fake it by acting like a branching candidate was found
   //---
   if(param.nodeLimit == 0)
      setStatus(AlpsNodeStatusPregnant);
   else if(doFathom) // || param.nodeLimit == 0)
      setStatus(AlpsNodeStatusFathomed);
   else
      status = chooseBranchingObject(model);
   
   UtilPrintFuncEnd(&cout, m_classTag,
                    "process()", param.msgLevel, 3);
   return status;
}

//===========================================================================//
int AlpsDecompTreeNode::chooseBranchingObject(AlpsModel * model) { 

   

   AlpsDecompNodeDesc * desc = 
     dynamic_cast<AlpsDecompNodeDesc*>(desc_);
   AlpsDecompModel * m     = dynamic_cast<AlpsDecompModel*>(desc->getModel());
   AlpsDecompParam & param = m->getParam();
   
   UtilPrintFuncBegin(&cout, m_classTag, "chooseBranchingObject()", 
                      param.msgLevel, 3);


   bool gotBranch = m->getDecompAlgo()->chooseBranchSet(downBranchLB_, 
                                                        downBranchUB_, 
                                                        upBranchLB_, 
                                                        upBranchUB_);
   if(!gotBranch){
      setStatus(AlpsNodeStatusEvaluated);
      //---
      //--- but if we can't branch on this and it DID finish pricing out
      //---   that means DW_LB=DW_UB for that node, then we are done 
      //---   processing it and we should fathom(?)
      //--- all we have to check is that LB=UB, since LB is updated 
      //--- despite the tailoff - so their should be a gap... 
      //--- the UB for this node, not the global UB...
      //---
      //printf("BestLB at this Node = %g\n", decompAlgo->getObjBestBoundLB());
      //printf("BestLB at this Node = %g\n", decompAlgo->getObjBestBoundUB();)
   }
   else{
     //---
     //--- we can go ahead and branch on this variable
     //---   meaning we will produce children (hence, the name pregnant)
     //---
      setStatus(AlpsNodeStatusPregnant);
   }
   UtilPrintFuncEnd(&cout, m_classTag, "chooseBranchingObject()", 
                    param.msgLevel, 3);

   return AlpsReturnStatusOk;
}

//===========================================================================//
std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > 
AlpsDecompTreeNode::branch()  { 

   AlpsDecompNodeDesc * desc  
      = dynamic_cast<AlpsDecompNodeDesc*>(desc_);
   AlpsDecompModel    * m     
      = dynamic_cast<AlpsDecompModel*>(desc->getModel());
   AlpsDecompParam    & param = m->getParam();

   UtilPrintFuncBegin(&cout, m_classTag, "branch()", param.msgLevel, 3);

   std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > newNodes;

   double *  oldLbs  = desc->lowerBounds_;
   double *  oldUbs  = desc->upperBounds_;
   const int numCols = desc->numberCols_;  
   CoinAssert(oldLbs && oldUbs && numCols);

   if((downBranchLB_.size() + downBranchUB_.size() == 0) || 
      (upBranchLB_.size()   + upBranchUB_.size()   == 0)) {
      std::cout << "AlpsDecompError: " 
                << "downBranch_.size() = "
                << downBranchLB_.size() + downBranchUB_.size() 
                << "; upBranch_.size() = "
                << upBranchLB_.size() + upBranchUB_.size() 
                << "; index_ = " << index_ << std::endl;
      throw CoinError("empty branch variable set(s)",
                      "branch", "AlpsDecompTreeNode");
   }
   
   /*if((branchedOn_ < 0) || (branchedOn_ >= numCols)) {
      std::cout << "AlpsDecompError: branchedOn_ = " 
                << branchedOn_ << "; numCols = " 
                << numCols << "; index_ = " << index_ << std::endl;
      throw CoinError("branch index is out of range", 
                      "branch", "AlpsDecompTreeNode");
                      }*/
  

   double * newLbs = new double[numCols];
   double * newUbs = new double[numCols];
   std::copy(oldLbs, oldLbs + numCols, newLbs);
   std::copy(oldUbs, oldUbs + numCols, newUbs);
   
   //printf("Start with old bounds:\n");
   //for(i = 0; i < numCols; i++){
   // printf("ind:%d -> lb:%g ub:%g\n", i, newLbs[i], newUbs[i]);
   //}
   
   
   double objVal(getQuality());

   
   //newLbs[branchedOn_] = oldLbs[branchedOn_];
   //newUbs[branchedOn_] = floor(branchedOnVal_);

   AlpsDecompNodeDesc* child = 0;

   //---
   //--- Branch down
   //---
   for (unsigned i = 0; i < downBranchLB_.size(); i++) {
      if((downBranchLB_[i].first < 0) || 
         (downBranchLB_[i].first >= numCols)) {
         std::cout << "AlpsDecompError: downBranchLB_[" << i << "] variable = "
                   << downBranchLB_[i].first << "; numCols = "
                   << numCols << "; index_ = " << index_ << std::endl;
         throw CoinError("branch index is out of range",
                         "branch", "AlpsDecompTreeNode");
      }
      newLbs[downBranchLB_[i].first] = downBranchLB_[i].second;
   }
   
   for (unsigned i = 0; i < downBranchUB_.size(); i++) {
      if((downBranchUB_[i].first < 0) || 
         (downBranchUB_[i].first >= numCols)) {
         std::cout << "AlpsDecompError: downBranchUB_[" << i << "] variable = "
                   << downBranchUB_[i].first << "; numCols = "
                   << numCols << "; index_ = " << index_ << std::endl;
         throw CoinError("branch index is out of range",
                         "branch", "AlpsDecompTreeNode");
      }
      newUbs[downBranchUB_[i].first] = downBranchUB_[i].second;
   }
   assert(downBranchLB_.size() +downBranchUB_.size() > 0);
   //assert(branchedOn_ >= 0);
   child = new AlpsDecompNodeDesc(m, newLbs, newUbs);
   //child->setBranchedOn(branchedOn_);
   //child->setBranchedVal(branchedOnVal_);
   child->setBranchedDir(-1);
   newNodes.push_back(CoinMakeTriple(static_cast<AlpsNodeDesc *>(child),
                                     AlpsNodeStatusCandidate,
                                     objVal));
       

   //---
   //--- Branch up
   //---
   //TODO: this can be done more cheaply than a full copy
   std::copy(oldLbs, oldLbs + numCols, newLbs);
   std::copy(oldUbs, oldUbs + numCols, newUbs);
   //newUbs[branchedOn_] = oldUbs[branchedOn_];
   //newLbs[branchedOn_] = ceil(branchedOnVal_);
   for (unsigned i = 0; i < upBranchLB_.size(); i++) {
      if((upBranchLB_[i].first < 0) || 
         (upBranchLB_[i].first >= numCols)) {
         std::cout << "AlpsDecompError: upBranchLB_[" << i << "] variable = "
                   << upBranchLB_[i].first << "; numCols = "
                   << numCols << "; index_ = " << index_ << std::endl;
         throw CoinError("branch index is out of range",
                         "branch", "AlpsDecompTreeNode");
      }
      newLbs[upBranchLB_[i].first] = upBranchLB_[i].second;
   }

   for (unsigned i = 0; i < upBranchUB_.size(); i++) {
      if((upBranchUB_[i].first < 0) || 
         (upBranchUB_[i].first >= numCols)) {
         std::cout << "AlpsDecompError: upBranchUB_[" << i << "] variable = "
                   << upBranchUB_[i].first << "; numCols = "
                   << numCols << "; index_ = " << index_ << std::endl;
         throw CoinError("branch index is out of range",
                         "branch", "AlpsDecompTreeNode");
      }      
      newUbs[upBranchUB_[i].first] = upBranchUB_[i].second;
   }
   
   child = 0;
   assert(upBranchLB_.size() + upBranchUB_.size() > 0);
   child = new AlpsDecompNodeDesc(m, newLbs, newUbs);
   //child->setBranchedOn(branchedOn_);
   //child->setBranchedVal(branchedOnVal_);
   child->setBranchedDir(1);
   newNodes.push_back(CoinMakeTriple(static_cast<AlpsNodeDesc *>(child),
                                     AlpsNodeStatusCandidate,
                                     objVal));
   //---
   //--- clean-up
   //---
   if (newLbs != 0) {
      delete [] newLbs;
      newLbs = 0;
   }
   if (newUbs != 0) {
      delete [] newUbs;
      newUbs = 0;
   }

   //---
   //--- change this node's status to branched
   //---    
   setStatus(AlpsNodeStatusBranched);

   UtilPrintFuncEnd(&cout, m_classTag, "branch()", param.msgLevel, 3);  
   return newNodes;   
}

