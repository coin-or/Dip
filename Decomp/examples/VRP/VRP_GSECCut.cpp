// $Id: VRP_GSECCut.cpp,v 1.4 2004/07/28 20:49:10 magh Exp $

/*-------------------------------------------------------------------------
  Author: Matthew Galati (magh@lehigh.edu)

  (c) Copyright 2004 Lehigh University. All Rights Reserved.

  This software is licensed under the Common Public License. Please see 
  accompanying file for terms.
  ---------------------------------------------------------------------------*/

#include "VRP_GSECCut.h"
#include "UtilMacros.h"

//how to deal with logs

/*-------------------------------------------------------------------------*/
void VRP_GSECCut::expandCutToRow(CoinPackedVector * row){
   //UtilPrintFuncBegin(&cout, m_classTag,
   //	      "expandCutToRow()", 3, 2);

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
   row->setVector(indices.size(), &indices[0], &elements[0], false);

   //UtilPrintVector(indices);
   //UtilPrintVector(elements);

   //UtilPrintFuncEnd(&cout, m_classTag,
   //	    "expandCutToRow()", 3, 2);

}

/*-------------------------------------------------------------------------*/
void VRP_GSECCut::print(ostream * os) const{
   //DecompCut::print(os);
   switch(m_type){
   case ACROSS:
      (*os) << "ACROSS ";
      break;
   case SIDE:
      (*os) << "SIDE ";
      break;
   case SIDE_COMPL:
      (*os) << "SIDE_COMPL ";
      break;
   }

   switch(m_storage){
   case VECTOR:
      {
         vector<int>::const_iterator it;
	 (*os) << "S: ";
	 for(it = m_S.begin(); it != m_S.end(); it++)
	    (*os) << *it << " ";
      }
      break;
   case BITSET:
   case BOTH:
      {
         int i;
	 (*os) << "S: ";
	 for(i = 0; i < m_nverts; i++)
	    if(m_inS[i])
	       (*os) << i << " ";
      }
      break;
   case NONE:
   default:
      cerr << "ERROR in print - BAD cut storage_type" << endl;
      abort();
   }
   (*os) << endl;
}

/*-------------------------------------------------------------------------*/
void VRP_GSECCut::setBounds(){

   //(i)   ACROSS:     2 * ceil( sum{i in S} d_i / C )
   //(ii)  SIDE:       |S| - ceil( sum{i in S} d_i / C )
   //(iii) SIDE_COMPL: |Shat| - ceil( sum{i in S} d_i / C )
  
   //cout << "m_demandS: " << m_demandS << " m_cap: " << m_capacity << endl;
   int bin = 
      static_cast<int>(ceil( static_cast<double>(m_demandS) / m_capacity ));
   switch(m_type){
   case ACROSS:
      setLowerBound(2.0 * bin);
      setUpperBound(DecompInf);
      break;
   case SIDE:
      setLowerBound(-DecompInf);
      setUpperBound(getSize() - bin);
      break;
   case SIDE_COMPL:
      setLowerBound(-DecompInf);
      setUpperBound(((m_nverts - 1) - getSize()) - bin);
      break;
   default:
      cerr << "ERROR in setBounds - BAD cut type" << endl;
      abort();
   }

   /*
   cout << "in setBounds LB: ";
   if(getLowerBound() > -(DecompInf/2))
      cout << getLowerBound();
   else
      cout << "-inf";
   cout << "in setBounds UB: ";
   if(getUpperBound() < (DecompInf/2))
      cout << getUpperBound();
   else
      cout << "inf";
   */

}

/*-------------------------------------------------------------------------*/
bool VRP_GSECCut::isSame(const DecompCut * cut) const{
   const VRP_GSECCut * gsec_cut = dynamic_cast<const VRP_GSECCut*>(cut);
   if(!gsec_cut)
      return false;

   if(m_type != gsec_cut->m_type)
      return false;
   switch(m_storage){
   case VECTOR:
      return m_S == gsec_cut->m_S;
   case BITSET:
   case BOTH:
      return m_inS == gsec_cut->m_inS;
   case NONE:
      return false;
   }
   return false;
}

/*-------------------------------------------------------------------------*/
void VRP_GSECCut::setStorage(){
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
      m_storage = NONE;
      return;
   }
}

//another function that probably belongs in a utility class    
/*-------------------------------------------------------------------------*/
void VRP_GSECCut::create_bitset(){
   //create bitset from vector
   m_inS.resize(m_nverts);
   for(vector<int>::iterator it = m_S.begin(); it != m_S.end(); it++)
      m_inS.set(*it);
   m_storage = m_storage == VECTOR ? BOTH : BITSET;
}
  
//another function that probably belongs in a utility class    
/*-------------------------------------------------------------------------*/
void VRP_GSECCut::create_vector(){
   //create vector from bistet
   m_S.reserve(m_inS.count());//is this worth it? or is count costly?
   for(unsigned int u = 0; u < m_inS.size(); u++)
      if(m_inS[u]) 
	 m_S.push_back(u);
   m_storage = m_storage == BITSET ? BOTH : VECTOR;
}
/* SIDE_COMPL needs more thought... see Decomp's VRP_GSEC_cut */

/*-------------------------------------------------------------------------*/
/* TODO */
void VRP_GSECCut::setCutType(){  
   //m_type = ACROSS;
   m_type = SIDE;
}

/*-------------------------------------------------------------------------*/
void VRP_GSECCut::setDemand(const int * vertex_wt){
   if(m_demandS != 0) 
      return;
   switch(m_type){
   case VECTOR:
   case BOTH:
      for(unsigned int i = 0; i < m_S.size(); i++){
	 m_demandS += vertex_wt[m_S[i]]; 
         //printf("m_demandS=%d, i:%d m_S:%d vwt:%d\n",
	 //     m_demandS, i, m_S[i], vertex_wt[m_S[i]]);
      }
      break;
   case BITSET:
      for(unsigned int u = 0; u < m_inS.size(); u++){
	 if(m_inS[u]){
            m_demandS += vertex_wt[u];
            //printf("m_demandS=%d, u:%d vwt:%d\n",
	    //     m_demandS, u, vertex_wt[u]);                                           
         } 
      }
      break;
   }
}

/*-------------------------------------------------------------------------*/
const int VRP_GSECCut::getSize(){
   switch(m_type){
   case VECTOR:
   case BOTH:
      return m_S.size();
   case BITSET:
      return m_inS.count();
   }
   return 0;
}


