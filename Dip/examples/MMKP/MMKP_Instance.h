//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, Lehigh University                                 //
//                                                                           //
// Copyright (C) 2002-2013, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

#ifndef MMKP_INSTANCE_INCLUDED
#define MMKP_INSTANCE_INCLUDED

//===========================================================================//
#include "UtilMacros.h"
using namespace std;
//===========================================================================//
class MMKP_Param;
//===========================================================================//

//===========================================================================//
/*!
 * \class MMKP_Instance
 * A class to store an instance of the 
 *     Multi-Dimensional Mulit-Choice Knapsack Problem (MMKP).
 *
 *     max  sum{i in 1..n, j in 1..l[i]} v[i,j]   x[i,j]
 *     s.t. sum{i in 1..n, j in 1..l[i]} r[k,i,j] x[i,j] <= b[k], k in 1..m
 *          sum{j in 1..l[i]}                     x[i,j]  = 1   , i in 1..n
 *          x[i,j] in {0,1}, i in 1..n, j in 1..l[i]
 *
 */

//===========================================================================//
class MMKP_Instance {
private:
   /** MMKP_Instance problem instance data */   
   int       m_nKnapRows;        //m
   int       m_nGroupRows;       //n
   int       m_nGroupCols;       //l
   double *  m_capacity;         //b[k = 1..m]
   double *  m_value;            //v[i,j]
   double ** m_weight;           //r[k,i,j]

   /** MMKP_Instance best known LB/UB */
   bool      m_isProvenOptimal;
   double    m_bestKnownLB;
   double    m_bestKnownUB;
   

public:
   /** @name Access methods. */
   inline const int       getNKnapRows ()  const {return m_nKnapRows; }
   inline const int       getNGroupRows()  const {return m_nGroupRows;}
   inline const int       getNGroupCols()  const {return m_nGroupCols;}
   inline const double *  getCapacity  ()  const {return m_capacity;  }
   inline const double *  getValue     ()  const {return m_value;     }
   inline const double * const* getWeight()const {return m_weight;    }

public:
   /** @name Helper Methods. */
   void readInstance     (string & fileName,
                          string & dataFormat);
   void readInstanceSimon(string & fileName);
   void readBestKnown(string & fileName,
                      string & instanceName);

   inline void initMembers(){
      m_nKnapRows  = 0;
      m_nGroupRows = 0;
      m_nGroupCols = 0;
      m_capacity   = NULL;
      m_value      = NULL;
      m_weight     = NULL;      
      m_isProvenOptimal =  false;
      m_bestKnownLB     = -1.e20;
      m_bestKnownUB     =  1.e20;
   }

   inline const int getIndexIJ(const int i,
                               const int j) const{
      return (i * m_nGroupCols) + j;
   }
   
   inline pair<int,int> getIndexInv(const int index) const {      
      return make_pair(index / m_nGroupCols, index % m_nGroupCols);
   }

   inline const double getBestKnownLB() const {return m_bestKnownLB;}
   inline const double getBestKnownUB() const {return m_bestKnownUB;}
   
public:
   /** @name Constructor and Destructor */

   /** Default constructor. */
   MMKP_Instance(){
      initMembers();
   };
   
   /** Default constructor. Takes an instance of UtilParameters */
   MMKP_Instance(string & fileName) {
      string dataFormat = "hifi";
      initMembers();
      readInstance(fileName, dataFormat);
   }
   
   /** Default destructor. */
   ~MMKP_Instance() {
      int k;
      UTIL_DELARR(m_capacity);
      UTIL_DELARR(m_value);
      for(k = 0; k < m_nKnapRows; k++)
         UTIL_DELARR(m_weight[k]);
      UTIL_DELARR(m_weight);
   };
};

#endif
