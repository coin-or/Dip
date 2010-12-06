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
// Copyright (C) 2002-2009, Lehigh University, Matthew Galati, Ted Ralphs    //
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
//---
//--- BcpsDecompTreeNode is derived from BcpsTreeNode
//---   BcpsTreeNode virtual methods that must be derived:
//---      (1) chooseBranchingObject (pure virtual)
//---      (2) installSubProblem     (pure virtual)
//---      (3) bound                 (pure virtual)
//---      (4) branch                (pure virtual)
//---   BcpsTreeNode virtual methods that should be derived:
//---      (1) generateConstraints   (does nothing)
//---      (2) generateVariables     (does nothing)
//---      (3) handleBoundingStatus  (does nothing)
//---      (4) process               (the meat of the algo?)
//--- BcpsTreeNode       is derived from AlpsTreeNode
//---   AlpsTreeNode virtual methods that must be derived:
//---      (1) createNewTreeNode     (pure virtual)
//---      (2) process               (pure virtual)
//---      (3) branch                (pure virtual)
//---   AlpsTreeNode virtual methods that should be derived:
//---      (1) convertToExplicit()   (does nothing)
//---      (2) convertToRelative()   (does nothing)
//---      

//===========================================================================//
class AlpsDecompTreeNode : public AlpsTreeNode {
private:
   /** class tag for debugging */
   string m_classTag;

   std::vector< std::pair<int, double> > downBranchLB_;
   std::vector< std::pair<int, double> > downBranchUB_;
   std::vector< std::pair<int, double> > upBranchLB_;
   std::vector< std::pair<int, double> > upBranchUB_;
  
   //int branchedOn_;
   //double branchedOnVal_;
  
public:
   /** Default constructor. */
   AlpsDecompTreeNode() :
      AlpsTreeNode(),
      m_classTag  ("ALPSTN")
     
   {
      //quality_ = -ALPS_OBJ_MAX;//MVG?
   }
    
   /** Destructor */
   virtual ~AlpsDecompTreeNode(){
   }

   void checkIncumbent(AlpsDecompModel      * model,
                       const DecompSolution * decompSol);
    
   //---
   //--- pure virtual functions from AlpsTreeNode or AlpsTreeNode
   //---  
    
   /** Create a new node based on given desc. */
   AlpsTreeNode * createNewTreeNode(AlpsNodeDesc *& desc) const;

   /** To be defined.?? */
   int chooseBranchingObject(AlpsModel * model);


   //point direct to DECOMP?
   /** Performing the bounding operation. */
   int process(bool isRoot = false, bool rampUp = false);
  
   /** Takes the explicit description of the current active node and 
       creates the children's descriptions, which contain information 
       about how the branching is to be done. The stati of the children
       are AlpsNodeStatusCandidate. */
   std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > branch(); 
};

#endif
