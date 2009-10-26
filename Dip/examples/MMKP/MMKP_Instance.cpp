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
void MMKP_Instance::readInstanceSimon(string & fileName){
   
   int      k, i, j, ij, numIJ, dummy;
   int      status = 0;
   string   dummyStr;
   char     dummyChr;
   char     dummyChrA[1000];
   ifstream is;
   

   /*
     # MMKP instances generated from definition file: tmp2.xml
     # serie_name    rp_hep_hec_strong
     # num_instances 1
     instance_id     1
     num_class       10
     num_dimension   5
     # class_id      1
     num_item        5
     1       2       2       3       3       3
     3       4       4       4       4       5
     5       6       6       6       6       6
     7       8       8       8       8       8
     9       10      10      10      10      10
     # class_id      2
     num_item        5
     ...
     c:      518     524     531     540     553
   */
   status = UtilOpenFile(is, fileName.c_str());   
   if(status)
      throw UtilException("Failed to read instance",
                          "readInstance", "MMKP_Instance");
   is.getline(dummyChrA, 1000);
   is.getline(dummyChrA, 1000);
   is.getline(dummyChrA, 1000);
   is.getline(dummyChrA, 1000);
   is >> dummyStr >> m_nGroupRows //num_class
      >> dummyStr >> m_nKnapRows; //num_dimension
   is.getline(dummyChrA, 1000);
   is.getline(dummyChrA, 1000);
   is >> dummyStr >> m_nGroupCols;
   printf("dummyStr = %s\n", dummyStr.c_str());
   printf("nGroupRows = %d\n", m_nGroupRows);
   printf("nKnapRows  = %d\n", m_nKnapRows);
   printf("nGroupCols = %d\n", m_nGroupCols);
   fflush(stdout);

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
   }
	       
   for(i = 0; i < m_nGroupRows; i++){
      for(j = 0; j < m_nGroupCols; j++){
         ij = getIndexIJ(i,j);
         is >> m_value[ij];
         printf("value[%d]: %g\n", ij, m_value[ij]);
         fflush(stdout);
         for(k = 0; k < m_nKnapRows; k++){
            is >> m_weight[k][ij];
            printf("weight[%d][%d]: %g\n", k, ij, m_weight[k][ij]);
            fflush(stdout);
         }
      }
      is >> dummyChr;
      if(dummyChr == '#'){
         is.getline(dummyChrA, 1000);
         is.getline(dummyChrA, 1000);
         printf("dummyChrA = %s\n", dummyChrA);
      }
      else{
         is >> dummyChr;
      }
   }
   

   printf("dummyChr = %c\n", dummyChr);
   for(k = 0; k < m_nKnapRows; k++){
      m_weight[k] = new double[numIJ];
      if(!m_weight[k])
         throw UtilExceptionMemory("readInstance", "MMKP_Instance");
      is >> m_capacity[k];
      printf("cap[k=%d]: %g\n", k, m_capacity[k]);
   }
      
   is.close();
}

//===========================================================================//
void MMKP_Instance::readInstance(string & fileName,
                                 string & dataFormat){
   
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
      if(dataFormat == "khan")
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
