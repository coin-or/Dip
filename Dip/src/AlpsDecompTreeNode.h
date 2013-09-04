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

//===========================================================================//
class AlpsNodeDesc;
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
      m_classTag  ("ALPSTN")

   {
      //quality_ = -ALPS_OBJ_MAX;//MVG?
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
};

#endif
