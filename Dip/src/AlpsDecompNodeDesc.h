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
#ifndef AlpsDecompNodeDesc_h_
#define AlpsDecompNodeDesc_h_

//===========================================================================//
#include "AlpsEncoded.h"
#include "AlpsNodeDesc.h"
#include "AlpsDecompModel.h"
#include "UtilMacrosAlps.h"

//===========================================================================//
class CoinWarmStartBasis;

//===========================================================================//
/**
 * \class AlpsDecompNodeDesc
 * \brief
 * Derivation of AlpsNodeDesc for DECOMP.
 *
 * An object derived from AlpsNodeDesc. This stores the description
 * of a search tree node. For DECOMP, we are not using differencing,
 * so, we only need to store the bounds set during branching.
 *
 * AlpsDecompNodeDesc is derived from AlpsNodeDesc
 *    AlpsModel has no pure virtual functions
 *
 * Virtual methods that should are derived here:
 *    encode
 *    decode
 *
 * \see
 * AlpsNodeDesc
 *
 * \todo
 * Invent a way to lose weight on a donut diet.
 * Use differencing scheme.
 */
//===========================================================================//

//===========================================================================//
class AlpsDecompNodeDesc : public AlpsNodeDesc {

 private:
   
   //----------------------------------------------------------------------//
   /**      
    * @name Data.
    * @{      
    */
   //----------------------------------------------------------------------//
   
   /**
    * Store the name of the class (for logging/debugging) - "who am I?"
    */
   string m_classTag;

   //STOP
   //what here do i need?
 public:  
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
   AlpsDecompNodeDesc() :
      AlpsNodeDesc(),
      branchedDir_(0),
      branchedInd_(-1),
      branchedVal_(0.0),
      basis_(NULL)
      {}

   /** Useful constructor. */
   AlpsDecompNodeDesc(AlpsModel * m) 
      :
      AlpsNodeDesc(m),
      branchedDir_(0),
      branchedInd_(-1),
      branchedVal_(0.0),
      basis_(NULL)
      {}

   AlpsDecompNodeDesc(AlpsDecompModel * m, 
		      const double    * lb, 
		      const double    * ub) 
      :
      AlpsNodeDesc(m),
      branchedDir_(0),
      branchedInd_(-1),
      branchedVal_(0.0),
      basis_(NULL)
      {
	 numberRows_ = m->getNumCoreRows();
	 numberCols_ = m->getNumCoreCols();
	 assert(numberRows_ && numberCols_);
	 lowerBounds_ = new double [numberCols_];
	 upperBounds_ = new double [numberCols_];
	 memcpy(lowerBounds_, lb, sizeof(double)*numberCols_);
	 memcpy(upperBounds_, ub, sizeof(double)*numberCols_);
      }

   /** Destructor. */
   virtual ~AlpsDecompNodeDesc() { 
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
   AlpsReturnStatus encodeAlpsDecomp(AlpsEncoded *encoded) const {
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
   AlpsReturnStatus decodeAlpsDecomp(AlpsEncoded &encoded) {
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
   //--- pure virtual functions from AlpsNodeDesc or AlpsNodeDesc
   //---  
    
   /** Pack node description into an encoded. */
   virtual AlpsReturnStatus encode(AlpsEncoded *encoded) const {
      AlpsReturnStatus status = AlpsReturnStatusOk;
	
      //status = encodeAlps(encoded);
      status = encodeAlpsDecomp(encoded);
	
      return status;
   }

   /** Unpack a node description from an encoded. Fill member data. */
   virtual AlpsReturnStatus decode(AlpsEncoded &encoded) {
	
      AlpsReturnStatus status = AlpsReturnStatusOk;
	
      //status = decodeAlps(encoded);
      status = decodeAlpsDecomp(encoded);

      return status;
   }
    
};
#endif
