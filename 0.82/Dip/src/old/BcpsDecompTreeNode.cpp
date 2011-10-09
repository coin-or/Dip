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
#include "DecompAlgo.h"
#include "DecompApp.h"

//really starting to look like an Alps derivation, not Bcps

//===========================================================================//
#include "AlpsKnowledgeBroker.h"
#include "BcpsDecompTreeNode.h"
#include "BcpsDecompNodeDesc.h"
#include "BcpsDecompSolution.h"
#include "BcpsDecompModel.h"
#include "CoinUtility.hpp"
 

//class BcpsDecompModel;
//===========================================================================//
AlpsTreeNode * 
BcpsDecompTreeNode::createNewTreeNode(AlpsNodeDesc *& desc) const {

   //---
   //--- create a new tree node, set node description
   //---
   BcpsDecompTreeNode * node = new BcpsDecompTreeNode();
   node->desc_ = desc;
   return node;  
}

//===========================================================================//
int BcpsDecompTreeNode::installSubProblem(BcpsModel * model) { 

   //is this pure virtual? yes

   /** Extract node information (bounds, constraints, variables) from 
       this node and load the information into the relaxation solver,
       such as linear programming solver.*/
   //what about cuts that have been added - how do we add them back in here?
  
   printf("\nBcpsDecompTreeNode::installSubProblem");
   return AlpsReturnStatusOk; 
}
  
//===========================================================================//
int BcpsDecompTreeNode::bound(BcpsModel * model) { 
   //TOOD: do nothing, I am rewriting process, but bound is pure virtual
   //should I really be deriving from ALPs? should DECOMP really be the
   //Bcps layer? we'll see when we get to differencing?
   printf("\nInside BcpsDecompTreeNode::bound() - do nothing");
   return AlpsReturnStatusOk; 
}

//===========================================================================//
int BcpsDecompTreeNode::process(bool isRoot, bool rampUp){
  
   //TODO: driven by DECOMP parameters?
   int    status = AlpsReturnStatusOk;
   bool   needBranch;
   double primalTolerance = 1.0e-6;
   decompStat stat = STAT_FEASIBLE;

   printf("\n\n\nSTART to process NODE: %d, at DEPTH: %d",
          getIndex(), getDepth());
  
  
   //---
   //--- TODO: who controls this? will i rampup? rampdown?
   //--- possible AlpsPhase: 
   //---   ALPS_PHASE_RAMPUP
   //---   ALPS_PHASE_SEARCH
   //---   ALPS_PHASE_RAMPDOWN
   //---
   AlpsPhase phase = getKnowledgeBroker()->getPhase();

   //---
   //--- get pointer to model, node description, decomp algo
   //---
   BcpsDecompNodeDesc * desc 
      = dynamic_cast<BcpsDecompNodeDesc*>(desc_);
   BcpsDecompModel    * model 
      = dynamic_cast<BcpsDecompModel*>(desc->getModel());
   DecompAlgo         * decompAlgo = model->getDecompAlgo();

   const double * lbs = desc->lowerBounds_;
   const double * ubs = desc->upperBounds_;
   if(!isRoot){
      //for PC easier if add bounds as cuts??
      //what do i need to update in decomp?
      //set decompAlgo's bounds

      //C
      //m_modelCore (need?) and m_masterSI?
      decompAlgo->setMasterBounds(lbs, ubs);
   }
    

   //---
   //--- check if this can be fathomed by objective cutoff
   //--- if not, set node active and add 1 to number of nodes
   //---
   double cutoff         = getKnowledgeBroker()->getIncumbentValue();
   double parentObjValue = getQuality();
   if ((parentObjValue - primalTolerance) > cutoff) {


    
      setStatus(AlpsNodeStatusFathomed);    
      goto TERM_PROCESS;
   }
   else {
      model->setActiveNode(this);
      model->addNumNodes();
   }

   //---
   //--- Restore, load and solve the subproblem.
   //---  (1) LP infeasible
   //---      a. set status to be fathom.
   //---  (2) LP feasible
   //---      a. MILP feasible. Check whether need update incumbent.
   //---      b. LP feasible but not MIP feasible. Check whether can be 
   //---         fathomed, if not, choose a branch variable.
   //---
  
   //---
   //--- extract info from this node and load subproblem into DECOMP solver
   //---
   //STOP - here is where we have to enforce the branching rules
   //inside decomp - resetting col bounds??
   //installSubProblem(model);

   //-------------------------------------------------------------------------
   //what about cuts?
   //const double *lbs = desc->lowerBounds_;
   //const double *ubs = desc->upperBounds_;

   //for(i = 0; i < numCols; ++i) {
   //  m->solver()->setColBounds(i, lbs[i], ubs[i]);
   //}


   //---
   //--- let DECOMP process the node (bounding step)
   //---
   //TODO: need to pass into decomp the UB, to see if can stop and fathom
   if(getIndex() > 0)
      printf("\nIn Tree");

   //getKnowledgeBroker()->printQualityActiveNodes();

   
   stat = decompAlgo->processNode(getIndex());
  

   //did we get any incumbents - if so, update with the best
   if(decompAlgo->getXhatIPBest()){
      printf("\nIP Feasible found in Decomp\n");
      // IP feasible       
      quality_ = decompAlgo->getXhatIPBest()->getQuality();
      if (quality_ < cutoff) {  
         //getKnowledgeBroker()->getNodeSelection()->setWeight(0.0);

         assert(desc->numberCols_ == decompAlgo->getXhatIPBest()->getSize());
         BcpsDecompSolution* ksol = 
            new BcpsDecompSolution(desc->numberCols_, 
                                   decompAlgo->getXhatIPBest()->getValues(),
                                   decompAlgo->getXhatIPBest()->getQuality()
                                   );
         printf("\nadd knowledge quality_ = %g", quality_);
         getKnowledgeBroker()->addKnowledge(AlpsKnowledgeTypeSolution, 
                                            ksol, 
                                            quality_); 
         // Update cutoff
         cutoff = getKnowledgeBroker()->getIncumbentValue();
         decompAlgo->setTrueUpperBound(cutoff);
      }
   }
   
      
      
      
   //THINK!
   //quality_ = decompAlgo->getTrueLowerBound();
   
   //decomp needs to returns status info
   
   //TODO: what to do next?
   //if LP feasible
   // if IP feasible and improving, create a BlisSolution
   // else quality_ (the LB) was also set in DECOMP... 
   
   
   //decomp needs a status here
   //quality_ = decompAlgo->getTrueLowerBound();//THINK! what if inf
   //quality_ = decompAlgo->getMasterSolverInterface()->getObjValue();
   //THINK
   
   
   
   needBranch = true;
   if(stat == STAT_FEASIBLE){
      //lpfeasible -> this can be false - need to have decomp return code
      //STOP add ip feas check
      
      bool isIPFeas = decompAlgo->isIPFeasible(decompAlgo->getX(),
					       1.0e-4,
                                               1.0e-4);
      if(isIPFeas){//ipfeasible -> should also check user feas check
         printf("\nIP Feasible - but is it user feasible?\n");
         // IP feasible
	 decompAlgo->getApp()->APPisUserFeasible(decompAlgo->getX(),
						 desc->numberCols_,
						 1.0e-5);
	 quality_ = decompAlgo->getMasterSolverInterface()->getObjValue();
         if (quality_ < cutoff) {  

            getKnowledgeBroker()->getNodeSelection()->setWeight(0.0);
            //who deletes this memory?
            BcpsDecompSolution* ksol = 
               new BcpsDecompSolution(desc->numberCols_, 
                                      decompAlgo->getX(),
                                      quality_);
            printf("\nadd knowledge quality_ = %g", quality_);
            getKnowledgeBroker()->addKnowledge(AlpsKnowledgeTypeSolution, 
                                               ksol, 
                                               quality_); 
            // Update cutoff
            cutoff = getKnowledgeBroker()->getIncumbentValue();
            decompAlgo->setTrueUpperBound(cutoff);
         }
         setStatus(AlpsNodeStatusFathomed);
         goto TERM_PROCESS;
      }
      else{
	 quality_ = decompAlgo->getTrueLowerBound();//THINK! what if inf
         //tailoff detection is done by DECOMP
         cutoff   = getKnowledgeBroker()->getIncumbentValue();
         if (quality_ > cutoff) {
            setStatus(AlpsNodeStatusFathomed);
            goto TERM_PROCESS;
         }
         needBranch = true;
         //reducedCostFix(model);//TODO? 
         //check for tailoff
         //remove non-core slack constraints??
      }
   }
   else if(stat == STAT_INFEASIBLE){
      //what are other possible cases?
      setStatus(AlpsNodeStatusFathomed);
      quality_ = -ALPS_OBJ_MAX;
      needBranch = false;
   }
  
   //TODO: of delete slack rows - in decomp doing this
   //TODO: primal huristics?
  
   //if(fathomed || !keepOn) {
   //	    // Infeasible, fathomed, tailing off or some user's criteriall.
   //    break;
   //}
  

  
   //------------------------------------------------------
   // Select branching object
   //------------------------------------------------------
   //STOP
   if (needBranch) { 
      status = chooseBranchingObject(model);
   }
  
 TERM_PROCESS:
  
   return status;
}


//===========================================================================//
int BcpsDecompTreeNode::chooseBranchingObject(BcpsModel * model) { 
   //int    branchedOnIndex;
   //double branchedOnValue;

   BcpsDecompNodeDesc* desc = 
      dynamic_cast<BcpsDecompNodeDesc*>(desc_);
   BcpsDecompModel* m = dynamic_cast<BcpsDecompModel*>(desc->getModel());

   m->getDecompAlgo()->chooseBranchVar(branchedOn_,
                                       branchedOnVal_);

   if(branchedOn_ == -1)
      setStatus(AlpsNodeStatusFathomed);
   else{
      setStatus(AlpsNodeStatusPregnant);
      //and set desc? the current lower bounds?    
   }
         

   //who deletes? 
   //what if already exists?
   //later on, use BcpsBranchStrategy?
  
   //member of BcpsTreeNode
   //??
   /* pure virtual
      if(branchObject_) {
      delete branchObject_;
      branchObject_ = NULL;
      }
      branchObject_ = BcpsBranchObject(model,
      branchedOnIndex,
      -1, //direction?
      branchedOnValue);
   */
   return AlpsReturnStatusOk;//??
}

//===========================================================================//
std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > 
BcpsDecompTreeNode::branch()  { 
   std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > newNodes;

   //what does it mean to branch up or down - this will depend on which
   //decomp algo you are using - let's just make it work for CPM first

   //int branchedOn_    = branchObject_->getObjectIndex();
   //int branchedOnVal_ = branchObject_->getValue();

   BcpsDecompNodeDesc* desc = 
      dynamic_cast<BcpsDecompNodeDesc*>(desc_);
   BcpsDecompModel* m = dynamic_cast<BcpsDecompModel*>(desc->getModel());
  
   //no differencing - that ok??? - for decomp, sure
   //not going to brnach that much anyway - we hope
   double* oldLbs = desc->lowerBounds_;
   double* oldUbs = desc->upperBounds_;
   const int numCols = desc->numberCols_;
  
   //int numInts(m->getNumInterdictVars());
   const double * sol = m->getDecompAlgo()->getX();
   assert(oldLbs && oldUbs && numCols);
  
   if((branchedOn_ < 0) || (branchedOn_ >= numCols)) {
      std::cout << "BcpsDecompError: branchedOn_ = " 
                << branchedOn_ << "; numCols = " 
                << numCols << "; index_ = " << index_ << std::endl;
      throw CoinError("branch index is out of range", 
                      "branch", "BcpsDecompTreeNode");
   }
  

   double* newLbs = new double[numCols];
   double* newUbs = new double[numCols];
   std::copy(oldLbs, oldLbs + numCols, newLbs);
   std::copy(oldUbs, oldUbs + numCols, newUbs);


   double objVal(getQuality());

   // Branch down
   newLbs[branchedOn_] = oldLbs[branchedOn_];
   newUbs[branchedOn_] = floor(branchedOnVal_);//floor(branchedOnVal_+1.0e-5);
  
   BcpsDecompNodeDesc* child;
   assert(branchedOn_ >= 0);
   child = new BcpsDecompNodeDesc(m, newLbs, newUbs);
   child->setBranchedOn(branchedOn_);
   child->setBranchedVal(branchedOnVal_);
   child->setBranchedDir(-1);
   newNodes.push_back(CoinMakeTriple(static_cast<AlpsNodeDesc *>(child),
                                     AlpsNodeStatusCandidate,
                                     objVal));
       
   // Branch up
   newUbs[branchedOn_] = oldUbs[branchedOn_];
   newLbs[branchedOn_] = ceil(branchedOnVal_);//ceil(branchedOnVal_ - 1.0e-5);
   child = 0;
   child = new BcpsDecompNodeDesc(m, newLbs, newUbs);
   child->setBranchedOn(branchedOn_);
   child->setBranchedVal(branchedOnVal_);
   child->setBranchedDir(1);
   newNodes.push_back(CoinMakeTriple(static_cast<AlpsNodeDesc *>(child),
                                     AlpsNodeStatusCandidate,
                                     objVal));

   if (newLbs != 0) {
      delete [] newLbs;
      newLbs = 0;
   }
   if (newUbs != 0) {
      delete [] newUbs;
      newUbs = 0;
   }
  
   return newNodes;  
}

