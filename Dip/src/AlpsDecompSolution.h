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

#ifndef AlpsDecompSolution_h
#define AlpsDecompSolution_h

#include "AlpsSolution.h"
#include "AlpsDecompModel.h"

class AlpsDecompSolution : public AlpsSolution {
 private:
    int size_;
    double* value_;
    double objective_;
    
 public:
    AlpsDecompSolution() 
	: 
	size_(0), 
	value_(0), 
	objective_() 
	{}
    AlpsDecompSolution(const int s, const double* val, const double obj) 
	: 
	size_(s)
	{ 
	    if (size_ >= 0) {
		value_ = new double [size_];
		memcpy(value_, val, sizeof(double) * size_);
	    }
	}

    ~AlpsDecompSolution() { 
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
    //virtual void print(std::ostream& os) const;
  

};

#endif
