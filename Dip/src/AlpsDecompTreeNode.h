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

#ifndef AlpsDecompTreeNode_h_
#define AlpsDecompTreeNode_h_

//===========================================================================//
#include "Alps.h"
#include "AlpsTreeNode.h"
#include "AlpsDecompNodeDesc.h"
#include "AlpsNodeDesc.h"

//class AlpsNodeDesc;
class AlpsDecompModel;

//===========================================================================//
class AlpsDecompTreeNode : public AlpsTreeNode {
private:
   /** class tag for debugging */
   std::string m_classTag;

   std::vector< std::pair<int, double> > downBranchLB_;
   std::vector< std::pair<int, double> > downBranchUB_;
   std::vector< std::pair<int, double> > upBranchLB_;
   std::vector< std::pair<int, double> > upBranchUB_;

public:
   /** Default constructor. */
   AlpsDecompTreeNode() :
      AlpsTreeNode(),
      m_classTag  ("ALPSTN") {
      //quality_ = -ALPS_OBJ_MAX;//MVG?
   }

   AlpsDecompTreeNode(AlpsNodeDesc *&desc){
     desc_ = desc; 
   }


   AlpsDecompTreeNode(AlpsDecompModel* m ) {
      AlpsDecompTreeNode();
   }


   /** Destructor */
   virtual ~AlpsDecompTreeNode() {
   }

   bool checkIncumbent(AlpsDecompModel*       model,
                       const DecompSolution* decompSol);

   //---
   //--- pure virtual functions from AlpsTreeNode or AlpsTreeNode
   //---

   /** Create a new node based on given desc. */
   AlpsTreeNode* createNewTreeNode(AlpsNodeDesc*& desc) const;

   /** To be defined.?? */
   int chooseBranchingObject(AlpsModel* model);



   /** Select a branching object based on give branching strategy. */
   int selectBranchObject(AlpsDecompModel* model,
                          bool& foundSol,
                          int numPassesLeft);

   //point direct to DECOMP?
   /** Performing the bounding operation. */
   int process(bool isRoot = false, bool rampUp = false);

   /** Takes the explicit description of the current active node and
       creates the children's descriptions, which contain information
       about how the branching is to be done. The stati of the children
       are AlpsNodeStatusCandidate. */
   std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > branch();

   //For now, we assume there is only one variable being branched on.
   int getBranchedVar() {
      if (!downBranchLB_.empty()) {
         return downBranchLB_[0].first;
      } else {
         return -1;
      }
   }


   inline std::vector< std::pair<int, double> >  getDownBranchLB() const {
      return downBranchLB_;
   }
   inline std::vector< std::pair<int, double> >  getDownBranchUB() const {
      return downBranchUB_;
   }
   inline std::vector< std::pair<int, double> >  getUpBranchLB() const {
      return upBranchLB_;
   }
   inline std::vector< std::pair<int, double> >  getUpBranchUB() const {
      return upBranchUB_;
   }

   inline void setBranchBound_upUB(int& BranchSize, int*& BranchIndices,
                                   double*& BranchValues) {
      if (BranchSize) {
         upBranchUB_.clear();

         for (int i = 0 ; i < BranchSize; ++i) {
            upBranchUB_.push_back(std::make_pair(BranchIndices[i], BranchValues[i]));
         }
      } else {
         upBranchUB_.clear();
      }
   }

   inline void setBranchBound_upLB(int& BranchSize, int*& BranchIndices,
                                   double*& BranchValues) {
      if (BranchSize) {
         upBranchLB_.clear();

         for (int i = 0 ; i < BranchSize; ++i) {
            upBranchLB_.push_back(std::make_pair(BranchIndices[i], BranchValues[i]));
         }
      } else {
         upBranchLB_.clear();
      }
   }

   inline void setBranchBound_downLB(int& BranchSize, int*& BranchIndices,
                                     double*& BranchValues) {
      if (BranchSize) {
         downBranchLB_.clear();

         for (int i = 0 ; i < BranchSize; ++i) {
            downBranchLB_.push_back(std::make_pair(BranchIndices[i], BranchValues[i]));
         }
      } else {
         downBranchLB_.clear();
      }
   }

   inline void setBranchBound_downUB(int& BranchSize, int*& BranchIndices,
                                     double*& BranchValues) {
      if (BranchSize) {
         downBranchUB_.clear();

         for (int i = 0 ; i < BranchSize; ++i) {
            downBranchUB_.push_back(std::make_pair(BranchIndices[i], BranchValues[i]));
         }
      } else {
         downBranchUB_.clear();
      }
   }


   /*
   inline void setBranchBound(int& BranchSize, int*& BranchIndices,
                              double*& BranchValues,
                              std::vector< std::pair<int, double> > & branchInfo) {
      if (BranchSize) {
         for (int i = 0 ; i < BranchSize; ++i) {
      branchInfo.push_back(std::make_pair(BranchIndices[i], BranchValues[i]));
    }
      } else {
         branchInfo.clear();
      }
   }
   */

   /** Encode this node for message passing. *\/ */
   virtual AlpsEncoded* encode() const;

   /* /\** Decode a node from an encoded object. *\/ */
   virtual AlpsKnowledge* decode(AlpsEncoded&) const;

};

#endif
