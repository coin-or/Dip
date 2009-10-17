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

//===========================================================================//
#include "MMKP_Instance.h"
#include "UtilMacrosDecomp.h"
//===========================================================================//

//===========================================================================//
void MMKP_Instance::readInstance(string & fileName){
   
   int      k, i, j, ij, numIJ, dummy;
   int      status = 0;
   ifstream is;
   
   //---
   //--- ftp://cermsem.univ-paris1.fr/pub/CERMSEM/hifi/MMKP/MMKP.html
   //---     NOTE: l[i] = l, for all i
   //---
   status = UtilOpenFile(is, fileName.c_str());   
   if(status)
      throw UtilException("Failed to read instance",
                          "readInstance", "MMKP_Instance");
   is >> m_nGroupRows
      >> m_nGroupCols
      >> m_nKnapRows;
   
   //---
   //--- allocate memory for capacity, value and weight
   //---
   numIJ      = m_nGroupCols * m_nGroupRows;   
   m_capacity = new double[m_nKnapRows];
   m_value    = new double[numIJ];
   m_weight   = new double*[m_nKnapRows];   
   if(!(m_capacity && m_value && m_weight))
      throw UtilExceptionMemory("readInstance", "MMKP_Instance");
	     
   for(k = 0; k < m_nKnapRows; k++){
      m_weight[k] = new double[numIJ];
      if(!m_weight[k])
         throw UtilExceptionMemory("readInstance", "MMKP_Instance");
      is >> m_capacity[k];
   }
   
   for(i = 0; i < m_nGroupRows; i++){
      is >> dummy;
      for(j = 0; j < m_nGroupCols; j++){
         ij = getIndexIJ(i,j);
         is >> m_value[ij];
         for(k = 0; k < m_nKnapRows; k++){
            is >> m_weight[k][ij];
         }
      }
   }
   
   is.close();
}

//===========================================================================//
void MMKP_Instance::readBestKnown(string & fileName,
                                  string & instanceName){

   ifstream is;
   string   instance;
   double   bestUpperBound;
   bool     isProvenOptimal;
   int      status  = 0;
   status = UtilOpenFile(is, fileName);
   if(status)
      throw UtilException("Failed to best-known file",
                          "readBestKnown", "MMKP_Instance");
   while(!is.eof()){
      is >> instance >> bestUpperBound >> isProvenOptimal;
      instance = UtilStrTrim(instance);
      if(instance == instanceName){
         if(isProvenOptimal)
            m_bestKnownLB = bestUpperBound;
         else
            m_bestKnownLB = -DecompInf;
         m_bestKnownUB     = bestUpperBound;
         m_isProvenOptimal = isProvenOptimal;
         break;
      }
   }   
}
