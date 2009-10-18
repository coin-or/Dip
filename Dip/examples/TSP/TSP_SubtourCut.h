// $Id: TSP_SubtourCut.hpp,v 1.9 2004/08/10 03:43:15 magh Exp $

/*-------------------------------------------------------------------------
  Author: Matthew Galati (magh@lehigh.edu)

  (c) Copyright 2004 Lehigh University. All Rights Reserved.

  This software is licensed under the Common Public License. Please see 
  accompanying file for terms.
  ---------------------------------------------------------------------------*/

#ifndef TSP_SUBTOUR_CUT_INCLUDED
#define TSP_SUBTOUR_CUT_INCLUDED

/*----------------------------------------------------------------------
  TSP_SubtourCut: TSP subtour elimination constraint
  
  (i)  ACROSS: sum{e in delta(S)} x_e >= 2
  (ii) SIDE  : sum{e in E(S)}     x_e <= |S| - 1
  ---------------------------------------------------------------------- */

#include "Decomp.h"
#include "DecompCut.h"
#include "UtilMacros.h"

#include <vector>
using namespace std;

/*---------------------------------------------------------------------------*/
class TSP_SubtourCut : public DecompCut {
public:
   enum storageType {VECTOR, BITSET, BOTH, NONE};
   enum cutType     {ACROSS, SIDE};

private:
   vector<int>       m_S;
   vector<bool>      m_inS;
   cutType           m_type;
   storageType       m_storage;
   int               m_nverts;

public:
   //these (pure virutal) methods are inherited from DecompCut
   void expandCutToRow(CoinPackedVector * row);
   //void expandCutToRow(DecompRow * row);
   void setBounds();

   //these (virutal) methods are inherited from DecompCut
   void print(ostream * os = &cout) const;
   bool isSame(const DecompCut * cut) const;

public:
   void init();
   void setCutType();
   void create_bitset();
   void create_vector();

public:
   TSP_SubtourCut(const vector<bool> & inS, 
		  const cutType        type = ACROSS){
      //assert(inS.count() < inS.size());
      //assert(inS.count() > 1);
      m_inS     = inS;
      m_storage = BITSET;
      m_nverts  = m_inS.size();
      m_type    = type;
      setBounds();
      init();
   };
  
   TSP_SubtourCut(const vector<bool> & inS, 
		  const vector<int>  & S,
		  const cutType        type){
      //assert(inS.count() < inS.size());
      //assert(inS.count() > 1);
      m_inS     = inS;
      m_S       = S;
      m_storage = BOTH;
      m_nverts  = m_inS.size();
      m_type    = type;
      setBounds();
      init();
   };

   TSP_SubtourCut(const vector<bool> & inS, 
		  const vector<int>  & S){
      //assert(inS.count() < inS.size());
      //assert(inS.count() > 1);
      m_inS     = inS;
      m_S       = S;
      m_storage = BOTH;
      m_nverts  = m_inS.size();
      setCutType();
      setBounds();
      init();
   };

   virtual ~TSP_SubtourCut() {}
};



   
  


#endif
