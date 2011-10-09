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
#ifndef BcpsDecompNodeDesc_h_
#define BcpsDecompNodeDesc_h_

//===========================================================================//
#include "AlpsEncoded.h"
#include "BcpsNodeDesc.h"
#include "BcpsDecompModel.h"
#include "UtilMacrosAlps.h"

//===========================================================================//
class CoinWarmStartBasis;
//class BcpsDecompModel;

//===========================================================================//
//---
//--- BcpsDecompNodeDesc is derived from BcpsNodeDesc
//---   BcpsNodeDesc has no virtual functions
//--- BcpsNodeDesc       is derived from AlpsNodeDesc
//---   AlpsNodeDesc virtual methods that must be derived:
//---      (1) encode ?? needed for serial ??
//---      (2) decode
//---      

//===========================================================================//
class BcpsDecompNodeDesc : public BcpsNodeDesc {
  
 private:
    
 public:
  /* Here, we need to fill in what the node description will look
     like. For now, we will not use differencing -- just explicitly
     represent it. Probably this means that we will just store the
       original problem data and a list of the variables that have been
       fixed. */
  
  /** */
  double* lowerBounds_;
  /** */
  double* upperBounds_;
  
  /** Number of rows in problem (before these cuts).  This
      means that for top of chain it must be rows at continuous */
  int numberRows_;
  ///
  int numberCols_;


    /** Branched direction to create it. */
    int branchedDir_;

    /** Branched object index to create it. */
    int branchedInd_;

    /** Branched value to create it. */
    double branchedVal_;

    //THINK: different derivations for different algos?
    /** Warm start. */
    CoinWarmStartBasis *basis_;
    
 public:

    /** Default constructor. */
    BcpsDecompNodeDesc() :
        BcpsNodeDesc(),
        branchedDir_(0),
        branchedInd_(-1),
        branchedVal_(0.0),
	basis_(NULL)
        {}

    /** Useful constructor. */
    //TODO: BlisModel as arg? or should be BcpsModel?
    //BcpsDecompNodeDesc(BcpsDecompModel * m) 
    BcpsDecompNodeDesc(BcpsModel * m) 
	:
	BcpsNodeDesc(m),
        branchedDir_(0),
        branchedInd_(-1),
        branchedVal_(0.0),
	basis_(NULL)
	{}

    BcpsDecompNodeDesc(BcpsDecompModel * m, 
                       const double    * lb, 
                       const double    * ub) 
	:
      BcpsNodeDesc(m),
      branchedDir_(0),
      branchedInd_(-1),
      branchedVal_(0.0),
      basis_(NULL)
      {
        numberRows_ = m->getNumRows();
        numberCols_ = m->getNumCols();
        assert(numberRows_ && numberCols_);
        lowerBounds_ = new double [numberCols_];
        upperBounds_ = new double [numberCols_];
        memcpy(lowerBounds_, lb, sizeof(double)*numberCols_);
        memcpy(upperBounds_, ub, sizeof(double)*numberCols_);
      }

    /** Destructor. */
    virtual ~BcpsDecompNodeDesc() { 
      if (lowerBounds_ != 0) {
        delete [] lowerBounds_;
        lowerBounds_ = 0;
      }
      if (upperBounds_ != 0) {
        delete [] upperBounds_;
        upperBounds_ = 0;
      }
      delete basis_; 
    }

    /** Set basis. */ 
    void setBasis(CoinWarmStartBasis *&ws) { 
        if (basis_) { delete basis_; }
        basis_= ws;
        ws = NULL; 
    }

    /** Get warm start basis. */
    CoinWarmStartBasis * getBasis() const { return basis_; }

    void setBranchedOn(int b) { branchedInd_ = b; }
    int getBranchedOn() const { return branchedInd_; }

    /** Set branching direction. */
    void setBranchedDir(int d) { branchedDir_ = d; }

    /** Get branching direction. */
    int getBranchedDir() const { return branchedDir_; }

    /** Set branching object index. */
    void setBranchedInd(int d) { branchedInd_ = d; }

    /** Get branching object index. */
    int getBranchedInd() const { return branchedInd_; }

    /** Set branching value. */
    void setBranchedVal(double d) { branchedVal_ = d; }

    /** Get branching direction. */
    double getBranchedVal() const { return branchedVal_; }

 protected:

    //---
    //--- helper functions for encode/decode
    //---

    /** Pack blis portion of node description into an encoded. */
    AlpsReturnStatus encodeBcpsDecomp(AlpsEncoded *encoded) const {
	AlpsReturnStatus status = AlpsReturnStatusOk;

	encoded->writeRep(branchedDir_);
	encoded->writeRep(branchedInd_);
	encoded->writeRep(branchedVal_);

	// Basis
	int ava = 0;
	if (basis_) {
	    ava = 1;
	    encoded->writeRep(ava);
            //should this be a util func or blis func?
            //seems pretty standard, alps/coin util type stuff
	    UtilAlpsEncodeWarmStart(encoded, basis_);
	}
	else {
	    encoded->writeRep(ava);
	}
	
	return status;
    }

    /** Unpack blis portion of node description from an encoded. */
    AlpsReturnStatus decodeBcpsDecomp(AlpsEncoded &encoded) {
	AlpsReturnStatus status = AlpsReturnStatusOk;
	
	encoded.readRep(branchedDir_);
	encoded.readRep(branchedInd_);
	encoded.readRep(branchedVal_);
	
	// Basis
	int ava;
	encoded.readRep(ava);
	if (ava == 1) {
	    basis_ = UtilAlpsDecodeWarmStart(encoded, &status);
	}
	else {
	    basis_ = NULL;
	}
	
	return status;
    }

 public:

    //---
    //--- pure virtual functions from BcpsNodeDesc or AlpsNodeDesc
    //---  
    
    /** Pack node description into an encoded. */
    virtual AlpsReturnStatus encode(AlpsEncoded *encoded) const {
    	AlpsReturnStatus status = AlpsReturnStatusOk;
	
	status = encodeBcps(encoded);
	status = encodeBcpsDecomp(encoded);
	
	return status;
    }

    /** Unpack a node description from an encoded. Fill member data. */
    virtual AlpsReturnStatus decode(AlpsEncoded &encoded) {
	
    	AlpsReturnStatus status = AlpsReturnStatusOk;
	
	status = decodeBcps(encoded);
	status = decodeBcpsDecomp(encoded);

	return status;
    }
    
};
#endif
