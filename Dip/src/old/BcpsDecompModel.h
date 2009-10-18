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

//copyright

//===========================================================================//
#ifndef BcpsDecompModel_h_
#define BcpsDecompModel_h_

//===========================================================================//
#include "DecompAlgo.h"
#include "BcpsModel.h"


//===========================================================================//
//---
//--- BcpsDecompModel is derived from BcpsModel
//---   BcpsModel has no virtual functions
//--- BcpsModel       is derived from AlpsModel
//---   AlpsModel virtual methods that must be derived:
//---      (1) readInstance   ?? why ??
//---      (2) createRoot
//---      

//===========================================================================//
class AlpsTreeNode;
//class BcpsModel;
//class DecompAlgo;

//===========================================================================//
class BcpsDecompModel : public BcpsModel {
 
 private:
  /** class tag for debugging   */
  static const char * m_classTag;

  /** ptr to decomp algo        */
  DecompAlgo        * decompAlgo_;

  /** ptr to active node        */
  AlpsTreeNode      * activeNode_;

  /** number of processed nodes */
  int                 numNodes_;

 private:
  //---
  //--- disable default copy constructor
  //---
  BcpsDecompModel(const BcpsDecompModel&);
  BcpsDecompModel& operator=(const BcpsDecompModel&);

 public:
  /** Default constructor. */
  BcpsDecompModel()
    : BcpsModel() 
    {
      init();
    }

  /** Default constructor. */
  BcpsDecompModel(DecompAlgo * decompAlgo) : 
    BcpsModel(), 
    decompAlgo_(decompAlgo)
    {
      init();
    }
  
  /** Destructor */
  virtual ~BcpsDecompModel(){
  }
  
  //---
  //--- pure virtual functions from BcpsModel or AlpsModel
  //---  

  //TODO: this will be decomp's read-in?
  /** Read in the instance data */
  void readInstance(const char* dataFile);

  /** create the root node */
  AlpsTreeNode * createRoot();
  
  //---
  //--- functions for setup or access
  //---
  /** initialize the model data */
  void init();

  /** get a ptr to the decomp algo */
  inline DecompAlgo * getDecompAlgo() const { return decompAlgo_; }

  //TODO: should this be at Bcps layer?
  /** set active node */
  inline void setActiveNode(AlpsTreeNode * node) { activeNode_ = node; }
  
  /** increment node count */
  inline void addNumNodes(int newNodes = 1) { numNodes_ += newNodes; }

  inline int getNumRows(){
    return decompAlgo_->m_modelCore->getNumRows();
  }

  inline int getNumCols(){
    return decompAlgo_->m_modelCore->getNumCols();
  }

};

#endif
