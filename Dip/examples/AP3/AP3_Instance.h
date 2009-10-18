//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, Lehigh University                                 //
//                                                                           //
// Copyright (C) 2002-2007, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

#ifndef AP3_INSTANCE_INCLUDED
#define AP3_INSTANCE_INCLUDED

// --------------------------------------------------------------------- //
#include "UtilMacros.h"
#include "CoinError.hpp"
using namespace std;
// --------------------------------------------------------------------- //

/*!
 * \class AP3_Instance
 * Stores an instance of 3-Indexed Assignment Problem (AP3)
 *
 * \todo Instance reader is not very robust.
 */

// --------------------------------------------------------------------- //
class AP3_Instance {
 public:
   /** The dimension of the instance. */
   int      m_dimension;

   /** The number of cols in full formulation.  */
   int      m_ncolsFull;

   /** The number of rows in full formulation.  */
   int      m_nrowsFull;
          
   /** Cost of assignment (size = dimension^3). */
   double * m_assigncost;

   /** Name of AP3 instance */
   string   m_instance;

   /** Known optimal bound (for debugging)      */
   double   m_optBound;
   
 public:
   /** Read the optimal bound from file */
   void readOptimalBound(const char * fileName) {
      string   buffer;
      ifstream is;      

      UtilOpenFile(is, fileName);
      
      is >> buffer 
         >> m_optBound;

      cout << "\nOPTIMAL BOUND = " << m_optBound << "\n";
      
      is.close();
   }



   //should this be in instance or decompapp?
   inline void index3Inv(const int   index,
                         int       & ind1, 
                         int       & ind2,
                         int       & ind3) const { 
      int n   = m_dimension;
      int nsq = n * n;
      int indexMod;
      ind1     = index / nsq;
      indexMod = index % nsq; 
      ind2     = indexMod / n;
      ind3     = indexMod % n;
   }


    /** Read in an instance from data file */
   void readInstance(const char * fileName) {
      int      c, n_cols, dimension;
      string   buffer;
      ifstream is;

      UtilOpenFile(is, fileName);

      //---
      //--- Data Format (Example) - Grundel
      //---   overall dimension             3
      //---   dimension of each vector      5 5 5
      //---   cost of assignment (i,j,k)    18 15 14 16 16
      //---                                 19 ...
      //---
      //--- NOTE: assume balanced (each vector has same dimension)
      //---
      is >> dimension;
      CoinAssert(dimension == 3);

      is >> m_dimension;
      getline(is, buffer);
      
      n_cols       = m_dimension * m_dimension * m_dimension;
      m_assigncost = new double[n_cols];
      CoinAssertHint(m_assigncost, "Error: Out of Memory");
      for(c = 0; c < n_cols; c++){
         is >> m_assigncost[c];
         
         //int ind1,ind2,ind3;
         //index3Inv(c, ind1, ind2, ind3);         
         //printf("(%d,%d,%d): %g\n", ind1+1, ind2+1, ind3+1, m_assigncost[c]);
      }      
      is.close();      

      //---
      //--- initialize other data members
      //---
      m_ncolsFull = m_dimension * m_dimension * m_dimension;
      m_nrowsFull = m_dimension * m_dimension * 3;
   }
    
 public:
   /** @name Constructor and Destructor */

   /** Default constructor. */
   AP3_Instance() :
      m_dimension(0),
      m_ncolsFull(0),
      m_nrowsFull(0),
      m_assigncost(0),
      m_optBound(DecompInf)
      {}
   
   /** Constructor. Construct an instance from data file */
   AP3_Instance(const char * fileName) :
      m_dimension(0),
      m_ncolsFull(0),
      m_nrowsFull(0),
      m_assigncost(0),
      m_optBound(DecompInf)
     {
       readInstance(fileName);
     }
   
   virtual ~AP3_Instance() {
     UTIL_DELARR(m_assigncost);
   };
};

#endif
