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
#include "GAP_Instance.h"
#include "UtilMacrosDecomp.h"

//===========================================================================//
void GAP_Instance::readInstance(string & fileName){

   int      i, j, n_ij, indexIJ;
   ifstream is;

   //---
   //--- File format (.../Decomp/data/GAP)
   //---
   //---   agents = machines (m, i index)
   //---   jobs   = tasks    (n, j index)
   //---
   //--- number of machines (m), number of tasks (n)
   //---   for each machine i (i=1,...,m) in turn:
   //---     cost of allocating task j to machine i (j=1,...,n)
   //---   for each machine i (i=1,...,m) in turn:
   //---     resource consumed in allocating task j to machine i (j=1,...,n)
   //--- resource capacity of machine j (j=1,...,m)
   //---

   UtilOpenFile(is, fileName.c_str());
   
   is >> m_nMachines
      >> m_nTasks;
   
   //---
   //--- allocate memory for capacity, value and weight
   //---
   n_ij = m_nMachines * m_nTasks;
   
   m_capacity = new int[m_nMachines];
   m_profit   = new int[n_ij];
   m_weight   = new int[n_ij];
   if(!(m_capacity && m_profit && m_weight))
      throw UtilExceptionMemory("readInstance", "GAP_Instance");

             
   indexIJ = 0;
   for(i = 0; i < m_nMachines; i++){
      for(j = 0; j < m_nTasks; j++){
         is >> m_profit[indexIJ++];//TODO: bad name - since cost
      }
   }
   indexIJ = 0;
   for(i = 0; i < m_nMachines; i++){
      for(j = 0; j < m_nTasks; j++){
         is >> m_weight[indexIJ++];
      }
   }
   for(j = 0; j < m_nMachines; j++)
      is >> m_capacity[j];
   
   is.close();
}

//===========================================================================//
void GAP_Instance::readBestKnown(string & fileName,
				 string & instanceName){

   ifstream is;
   string   instance;
   double   bestUpperBound;
   bool     isProvenOptimal;
   int      status  = 0;
   status = UtilOpenFile(is, fileName);
   if(status)
      throw UtilException("Failed to best-known file",
                          "readBestKnown", "GAP_Instance");
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
