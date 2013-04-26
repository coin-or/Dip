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

#ifndef ATM_INSTANCE_INCLUDED
#define ATM_INSTANCE_INCLUDED

// --------------------------------------------------------------------- //
#include "UtilMacros.h"

using namespace std;

// --------------------------------------------------------------------- //

/*!
 * \class ATM_Instance
 * A class to store an instance of the 
 *     ATM Cash Management Problem (ATM).
 * 
 * The original problem is in the form of a MINLP:
 *    min  sum{a in A, d in D} | f[a,d] | 
 *    s.t.
 *    for a in A, d in D:
 *       f[a,d] =  
 *         a[a,d] x1[a]       +
 *         b[a,d] x2[a]       +     
 *         c[a,d] x1[a] x2[a] + 
 *         d[a,d] x3[a]       +
 *         e[a,d]
 *    for d in D:
 *       sum{a in A} f[a,d]     <= B[d]
 *    for a in A:
 *      |{d in D : f[a,d] <= 0} <= K[a]
 *    for a in A, d in D:
 *      f[a,d] free
 *    for a in A:
 *      x1[a], x2[a] in [0,1] x3[a] >= 0
 *
 */

// --------------------------------------------------------------------- //
class ATM_Instance {

private:
   /** ATM_Instance problem instance data */   
   int                 m_nAtms;         //number of atms  (20  of these)
   int                 m_nDates;        //number of dates (272 of these)
   map<string, int>    m_strToIntAtms;  //map from name to index (atms)
   map<string, int>    m_strToIntDates; //map from date to index (dates)
   vector<string>      m_intToStrAtms;  //map from index to name (atms)
   vector<string>      m_intToStrDates; //map from index to date (dates)
   vector<int>         m_pairsAD;       //pairs of indices with observations

   double            * m_a_ad;          //a[a,d]
   double            * m_b_ad;          //b[a,d]
   double            * m_c_ad;          //c[a,d]
   double            * m_d_ad;          //d[a,d]
   double            * m_e_ad;          //e[a,d]
   double            * m_w_ad;          //w[a,d]
   double            * m_B_d;           //B[d]
   double            * m_K_a;           //K[a]

public:
   /** @name Access methods. */
   inline const int    getNAtms    ()  const {return m_nAtms; }
   inline const int    getNDates   ()  const {return m_nDates;}
   inline const int    getNPairs   ()  const {
      return static_cast<int>(m_pairsAD.size());
   }
   inline const vector<int> & getPairsAD() const {
      return m_pairsAD;
   }
   inline const int getIndexAD(int a, int d) const {
      return a * m_nDates + d;
   }
   inline const pair<int,int> getIndexADInv(int ad) const {
      return make_pair(ad / m_nDates, ad % m_nDates);
   }
   inline const double * get_a_ad() const {return m_a_ad;}
   inline const double * get_b_ad() const {return m_b_ad;}
   inline const double * get_c_ad() const {return m_c_ad;}
   inline const double * get_d_ad() const {return m_d_ad;}
   inline const double * get_e_ad() const {return m_e_ad;}
   inline const double * get_w_ad() const {return m_w_ad;}
   inline const double * get_B_d()  const {return m_B_d;}
   inline const double * get_K_a()  const {return m_K_a;}
   inline const string & getAtmName(const int a) const{
      return m_intToStrAtms[a];
   }
   inline const string & getDateName(const int d) const{
      return m_intToStrDates[d];
   }

   
public:
   /** @name Helper Methods. */
   void readInstance(string & fileNameA,
		     string & fileNameD,
		     string & fileNameAD);
   void generateRandom(const int nAtms,
		       const int nDates,
		       const int seed);
   void initMembers(){
      m_nAtms = 0;
      m_nDates= 0;
      m_a_ad  = NULL;
      m_b_ad  = NULL;
      m_c_ad  = NULL;
      m_d_ad  = NULL;
      m_e_ad  = NULL;
      m_w_ad  = NULL;
      m_B_d   = NULL;
      m_K_a   = NULL;
   }

public:
   /** @name Constructor and Destructor */
   
   /** Default constructor. Takes an instance of UtilParameters */
   ATM_Instance(){initMembers();};
 
 ATM_Instance(string & fileNameA,
              string & fileNameD,
              string & fileNameAD)
   {
      initMembers();
      readInstance(fileNameA,
		   fileNameD,
		   fileNameAD);
   }
   
   ~ATM_Instance() {
      UTIL_DELARR(m_a_ad);
      UTIL_DELARR(m_b_ad);
      UTIL_DELARR(m_c_ad);
      UTIL_DELARR(m_d_ad);
      UTIL_DELARR(m_e_ad);
      UTIL_DELARR(m_w_ad);
      UTIL_DELARR(m_B_d);
      UTIL_DELARR(m_K_a);
   };
};

#endif
