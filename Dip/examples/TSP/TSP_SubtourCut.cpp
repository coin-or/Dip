// $Id: TSP_SubtourCut.cpp,v 1.6 2004/03/03 01:00:19 magh Exp $

/*-------------------------------------------------------------------------
  Author: Matthew Galati (magh@lehigh.edu)

  (c) Copyright 2004 Lehigh University. All Rights Reserved.

  This software is licensed under the Common Public License. Please see 
  accompanying file for terms.
  ---------------------------------------------------------------------------*/

#include "TSP_SubtourCut.h"

#include "Decomp.h"
#include "UtilMacros.h"

/*-------------------------------------------------------------------------*/
void TSP_SubtourCut::init(){
   //setCutType();
   switch(m_storage){
   case VECTOR:
      create_bitset();
      break;
   case BITSET:
      create_vector();
      break;
   case BOTH:
      break;
   default:
      //throw exception
      assert(0);
      return;
   }
}

/*-------------------------------------------------------------------------*/
bool TSP_SubtourCut::isSame(const DecompCut * cut) const{
   const TSP_SubtourCut * sec_cut = dynamic_cast<const TSP_SubtourCut*>(cut);
   if(!sec_cut)
      return false;

   if(m_type != sec_cut->m_type)
      return false;
   switch(m_storage){
   case VECTOR:
      return m_S == sec_cut->m_S;
   case BITSET:
   case BOTH:
      return m_inS == sec_cut->m_inS;
   }
   return false;
}

/*-------------------------------------------------------------------------*/
void TSP_SubtourCut::setCutType(){
   //what does concorde do?
   //sec_type = getSize() <= (n_vertices - 1)/2 ? SIDE : ACROSS;
   m_type = ACROSS;
}

//another function that probably belongs in a utility class    
/*-------------------------------------------------------------------------*/
void TSP_SubtourCut::create_bitset(){
   //create bitset from vector
   for(int i = 0; i < m_nverts; i++)
      m_inS.push_back(false);//??
   for(vector<int>::iterator it = m_S.begin(); it != m_S.end(); it++)
      m_inS[*it] = true;
   m_storage = m_storage == VECTOR ? BOTH : BITSET;
}
  
//another function that probably belongs in a utility class    
/*-------------------------------------------------------------------------*/
void TSP_SubtourCut::create_vector(){
   //create vector from bistet
   //m_S.reserve(m_inS.count());//is this worth it? or is count costly?
   for(unsigned int u = 0; u < m_inS.size(); u++)
      if(m_inS[u]) 
	 m_S.push_back(u);
   m_storage = m_storage == BITSET ? BOTH : VECTOR;
}

/*-------------------------------------------------------------------------*/
void TSP_SubtourCut::setBounds(){
   switch(m_type){
   case ACROSS:
      setLowerBound(2.0);
      setUpperBound(DecompInf);
      break;
   case SIDE:
      setLowerBound(-DecompInf);
      setUpperBound(static_cast<int>(m_S.size()) - 1.0);
      break;
   default:
      assert(0);
   }   
}

/*-------------------------------------------------------------------------*/
void TSP_SubtourCut::expandCutToRow(CoinPackedVector * row){
   vector<int>    indices;
   vector<double> elements;
  
   switch(m_type){
   case ACROSS:
      {	
	 for(unsigned int i = 0; i < m_S.size(); i++){
	    for(int v = 0; v < m_nverts; v++){
	       if(m_inS[v] || m_S[i] == v) 
		  continue;
	       indices.push_back(UtilIndexU(m_S[i],v));	      
	    }
	 }
	 fill_n(back_inserter(elements), indices.size(), 1.0);
      }
      break;
    
   case SIDE:
      {
	 for(unsigned int i = 0; i < m_S.size(); i++)
	    for(unsigned int j = i + 1; j < m_S.size(); j++)
	       indices.push_back(UtilIndexU(m_S[i], m_S[j]));
	 fill_n(back_inserter(elements), indices.size(), 1.0);
      }
      break;
   default:
      cerr << "ERROR expandCutToRow sec_type" << endl;
      abort();
   }  
   row->setVector(static_cast<int>(indices.size()), 
		  &indices[0], &elements[0], false);
}

/*-------------------------------------------------------------------------*/
void TSP_SubtourCut::print(ostream * os) const{
   double lb, ub;
   switch(m_type){
   case ACROSS:
      (*os) << "ACROSS ";
      break;
   case SIDE:
      (*os) << "SIDE ";
      break;
   }

   switch(m_storage){
   case VECTOR:
      {
	 (*os) << "S: ";
	 for(vector<int>::const_iterator it = m_S.begin(); 
	     it != m_S.end(); it++)
	    (*os) << *it << " ";
      }
      break;
   case BITSET:
   case BOTH:
      {
	 (*os) << "S: ";
	 for(int i = 0; i < m_nverts; i++)
	    if(m_inS[i])
	       (*os) << i << " ";
      }
      break;
   default:
      cerr << "ERROR in print - BAD cut storage_type" << endl;
      abort();
   }
   lb = getLowerBound();
   ub = getUpperBound();
   if(lb > -DecompInf){
      (*os) << "\tm_lb\t" << lb;
   }
   else{
      (*os) << "\tm_lb\t-INF";
   }
   if(ub < DecompInf){
      (*os) << "\tm_ub\t" << ub;
   }
   else{
      (*os) << "\tm_ub\t INF";
   }
   (*os) << endl;
}
