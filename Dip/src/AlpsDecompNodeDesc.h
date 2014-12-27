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
#ifndef AlpsDecompNodeDesc_h_
#define AlpsDecompNodeDesc_h_

//===========================================================================//
#include "AlpsEncoded.h"
#include "AlpsNodeDesc.h"
#include "AlpsDecompModel.h"
#include "UtilMacrosAlps.h"

//===========================================================================//
//class AlpsDecompModel;
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
   std::string m_classTag;

public:

   /** number of columns in original space */
   int numberCols_;
   /** lower bounds in original space */
   double* lowerBounds_;
   /** upper bounds in original space */
   double* upperBounds_;
   /** Branched direction to create it. */
   int branchedDir_;
   /** Branched set of indices/values to create it. */
   std::vector< std::pair<int, double> > branched_;

   //THINK: different derivations for different algos? need?
   /** Warm start. */
   // not needed I think
   //   CoinWarmStartBasis* basis_;

public:

   /** Default constructor. */
   AlpsDecompNodeDesc() :
      AlpsNodeDesc(),
      lowerBounds_(NULL),
      upperBounds_(NULL),
      branchedDir_(0) {
   }

   /** Useful constructor. */
   AlpsDecompNodeDesc(AlpsDecompModel* m)
      :
      AlpsNodeDesc(m),
      branchedDir_(0)
      //      basis_(NULL)
   {
      numberCols_ = m->getNumCoreCols();
      assert(numberCols_);
      lowerBounds_ = new double[numberCols_];
      upperBounds_ = new double[numberCols_];

      for (int i = 0 ; i < numberCols_ ; i ++) {
         lowerBounds_[i] = -DecompInf;
         upperBounds_[i] = DecompInf;
      }
   }

   AlpsDecompNodeDesc(AlpsDecompModel* m,
                      const double*     lb,
                      const double*     ub)
      :
      AlpsNodeDesc(m)
      //    basis_(NULL)
   {
      numberCols_ = m->getNumCoreCols();
      assert(numberCols_);
      lowerBounds_ = new double [numberCols_];
      upperBounds_ = new double [numberCols_];
      memcpy(lowerBounds_, lb, sizeof(double)*numberCols_);
      memcpy(upperBounds_, ub, sizeof(double)*numberCols_);
      branchedDir_ = 0;
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

      //      delete basis_;
   }

   /** Set basis. */
   /*
      void setBasis(CoinWarmStartBasis*& ws) {

         if (basis_) {
            delete basis_;
         }

         basis_ = ws;
         ws = NULL;
      }
   */
   /** Get warm start basis. */
   /*
      CoinWarmStartBasis* getBasis() const {
         return basis_;
      }
   */
   /** Set branching direction. */
   void setBranchedDir(int d) {
      branchedDir_ = d;
   }

   /** Get branching direction. */
   int getBranchedDir() const {
      return branchedDir_;
   }

   /** Set branching set. */
   void setBranched(std::vector< std::pair<int, double> > b) {
      branched_ = b;
   }

   /** Get branching set. */
   std::vector< std::pair<int, double> > getBranched() const {
      return branched_;
   }

public:
   //---
   //--- helper functions for encode/decode
   //---

   /** Pack AlpsDecomp portion of node description into an encoded. */
   AlpsReturnStatus encodeAlpsDecomp(AlpsEncoded* encoded) const {
      AlpsReturnStatus status = AlpsReturnStatusOk;
      double* lb = new double[numberCols_];
      double* ub = new double[numberCols_];
      memcpy( lb, lowerBounds_, numberCols_ * sizeof(double) );
      memcpy( ub, upperBounds_, numberCols_ * sizeof(double) );
      encoded->writeRep(numberCols_);
      encoded->writeRep(lb, numberCols_);
      encoded->writeRep(ub, numberCols_);
      encoded->writeRep(branchedDir_);
      delete [] lb;
      delete [] ub;
      return status;
   }

   /** Unpack AlpsDecompNodeDesc portion of node description from an encoded. */
   AlpsReturnStatus decodeAlpsDecomp(AlpsEncoded& encoded) {
      AlpsReturnStatus status = AlpsReturnStatusOk;
      encoded.readRep(numberCols_);
      encoded.readRep(lowerBounds_, numberCols_);
      encoded.readRep(upperBounds_, numberCols_);
      encoded.readRep(branchedDir_);
      return status;
   }

   //---
   //--- pure virtual functions from AlpsNodeDesc or AlpsNodeDesc
   //---

   /** Pack node description into an encoded. */
   virtual AlpsReturnStatus encode(AlpsEncoded* encoded) const {
      AlpsReturnStatus status = AlpsReturnStatusOk;
      status = encodeAlpsDecomp(encoded);
      return status;
   }

   /** Unpack a node description from an encoded. Fill member data. */
   virtual AlpsReturnStatus decode(AlpsEncoded& encoded) {
      AlpsReturnStatus status = AlpsReturnStatusOk;
      status = decodeAlpsDecomp(encoded);
      return status;
   }

};
#endif