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

#ifndef BcpsDecompTreeNode_h_
#define BcpsDecompTreeNode_h_

//THINK: gardner discussion about interface class to hide the details
//only show use the virtual methods 

//===========================================================================//
#include "Alps.h"
#include "BcpsTreeNode.h"

//===========================================================================//
class AlpsTreeNode;
class AlpsNodeDesc;

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
class BcpsDecompTreeNode : public BcpsTreeNode {
 private:
  /** class tag for debugging */
  static const char * m_classTag;
  
  int branchedOn_;
  double branchedOnVal_;
  
 private:
  //---
  //--- disable default copy constructor
  //---
  BcpsDecompTreeNode(const BcpsDecompTreeNode&);
  BcpsDecompTreeNode& operator=(const BcpsDecompTreeNode&);
  
 public:
  /** Default constructor. */
  BcpsDecompTreeNode()
    : BcpsTreeNode() 
    {
    }
    
    /** Destructor */
    virtual ~BcpsDecompTreeNode(){
    }
    
    //---
    //--- pure virtual functions from BcpsTreeNode or AlpsTreeNode
    //---  
    
    /** Create a new node based on given desc. */
    AlpsTreeNode * createNewTreeNode(AlpsNodeDesc *& desc) const;

  /** To be defined.?? */
  int chooseBranchingObject(BcpsModel * model);

  /** intall subproblem */
  int installSubProblem(BcpsModel * model);
  
  /** Bounding procedure */
  int bound(BcpsModel * model);

  //point direct to DECOMP?
  /** Performing the bounding operation. */
  int process(bool isRoot, bool rampUp);
  
  /** Takes the explicit description of the current active node and 
      creates the children's descriptions, which contain information 
      about how the branching is to be done. The stati of the children
      are AlpsNodeStatusCandidate. */
  std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > branch(); 
};

#endif
