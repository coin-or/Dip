// $Id: VRP_GSECCut.hpp,v 1.2 2004/06/16 22:26:37 magh Exp $

/*-------------------------------------------------------------------------
  Author: Matthew Galati (magh@lehigh.edu)

  (c) Copyright 2004 Lehigh University. All Rights Reserved.

  This software is licensed under the Common Public License. Please see 
  accompanying file for terms.
  ---------------------------------------------------------------------------*/

#ifndef VRP_SUBTOUR_CUT
#define VRP_SUBTOUR_CUT
/*----------------------------------------------------------------------
  Type 1: VRP_GSEC_cut
 
  Generalized Subtour Elimination Constraints:

  (i) ACROSS:
  sum{e in delta(S)} x_e >=   2 * ceil( sum{i in S} d_i / C )

  (ii) SIDE:
  sum{e in E(S)    } x_e <= |S| - ceil( sum{i in S} d_i / C )

  (iii) SIDE_COMPL:
  sum{e in E(Shat)} + 0.5 x({0} : Shat) - 0.5 x({0} : S) <= |Shat| - k(S)

  if(|S| <= n/2) use (ii), else use (iii) for sparsest form

  Shat = V \ S (note, V does not include depot here)
  i.e., customers that are not in the set S
  ---------------------------------------------------------------------- */

#include "Decomp.h"
#include "DecompCut.h"

#include <vector>
using namespace std;

//TODO: switch to just use bool vector?
#include <boost/dynamic_bitset.hpp>
using namespace boost;

/*---------------------------------------------------------------------- */
class VRP_GSECCut : public DecompCut {
public:
   enum storageType {VECTOR, BITSET, BOTH, NONE};
   enum cutType     {ACROSS, SIDE, SIDE_COMPL};

private:
   /** Class id tag (for log / debugging). */
   const string m_classTag;


   vector<int>        m_S;
   dynamic_bitset<>   m_inS;
   cutType            m_type;
   storageType        m_storage;
   int                m_nverts;
   int                m_demandS;
   int                m_capacity;

public:
   //these (pure virutal) methods are inherited from DecompCut
   virtual void expandCutToRow(CoinPackedVector * row);
   virtual void setBounds();  
  
   //these (virutal) methods are inherited from DecompCut
   virtual void print(ostream * os = &cout) const;
   virtual bool isSame(const DecompCut * cut) const;
  
public:
   void setStorage();
   void create_bitset();
   void create_vector();
   void setCutType();
   void setDemand(const int * vertex_wt);
   const int getSize();
    
public:
   VRP_GSECCut(const dynamic_bitset<> & inS, 
	       const int              * vertex_wt,
	       const int                capacity,
	       const int                demandS = 0) :
      m_classTag("VRP-GSEC"),
      m_inS     (inS),
      m_storage (BITSET),
      m_nverts  (m_inS.size()),
      m_demandS (demandS),
      m_capacity(capacity)
   {
      setStorage();
      setCutType();
      setDemand(vertex_wt);
      setBounds();
   }

   virtual ~VRP_GSECCut() {}
};

#endif
