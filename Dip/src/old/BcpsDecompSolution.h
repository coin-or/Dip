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

/*===========================================================================*
 * This file is part of the Abstract Library for Parallel Search (ALPS).     *
 *                                                                           *
 * ALPS is distributed under the Common Public License as part of the        *
 * COIN-OR repository (http://www.coin-or.org).                              *
 *                                                                           *
 * Authors: Yan Xu, SAS Institute Inc.                                       *
 *          Ted Ralphs, Lehigh University                                    *
 *          Laszlo Ladanyi, IBM T.J. Watson Research Center                  *
 *          Matthew Saltzman, Clemson University                             *
 *                                                                           * 
 *                                                                           *
 * Copyright (C) 2001-2006, Lehigh University, Yan Xu, and Ted Ralphs.       *
 *===========================================================================*/

#ifndef BcpsDecompSolution_h
#define BcpsDecompSolution_h

#include "AlpsSolution.h"
#include "BcpsDecompModel.h"

/** This class holds a MIP feasible primal solution. */
class BcpsDecompSolution : public AlpsSolution {
 private:
    int size_;
    double* value_;
    double objective_;
    
 public:
    BcpsDecompSolution() 
	: 
	size_(0), 
	value_(0), 
	objective_() 
	{}
    BcpsDecompSolution(const int s, const double* val, const double obj) 
	: 
	size_(s)
	{ 
	    if (size_ >= 0) {
		value_ = new double [size_];
		memcpy(value_, val, sizeof(double) * size_);
	    }
	}

    ~BcpsDecompSolution() { 
	if (value_ != 0) {
	    delete [] value_; 
	    value_ = 0;
	}
    }
  
    /** Get the objective value value */
    double getObjValue() const { return objective_; }

    virtual double getQuality() const { return getObjValue(); }
  
    /** Get the size of the solution */
    int getSize() const { return size_; }
    
    /** Get the column solution */
    const double* getColSolution() const 
	{ return value_; }
    
    /** Get item i in the solution vector */
    double getColSolution(int i) const { return value_[i]; }
    
    /** Print out the solution.*/
    virtual void print(std::ostream& os) const;
  
#if 0
    /** The method that encodes the solution into a encoded object. */
    virtual AlpsEncoded* encode() const;
  
    /** The method that decodes the solution from a encoded object. */
    virtual AlpsKnowledge* decode(AlpsEncoded&) const;
#endif
};

#endif
