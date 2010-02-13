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
#include "MCF_Instance.h"
#include "UtilMacrosDecomp.h"
//===========================================================================//

//===========================================================================//
int MCF_Instance::readInstance(string & fileName,
                               bool     addDummyArcs){
   
   ifstream is;   
   int      status = UtilOpenFile(is, fileName.c_str());   
   if(status)
      throw UtilException("Failed to read instance",
                          "readInstance", "MCF_Instance");
   
   double sumweight        = 0;
   bool   size_read        = true;
   int    arcs_read        = 0;
   int    commodities_read = 0;;
   char   line[1000];
   char   name[1000];
   while(is.good()) {
      is.getline(line, 1000);
      if (is.gcount() >= 999) {
         cerr << "ERROR: Input file is incorrect. "
              << "A line more than 1000 characters is found." << endl;
         return 1;
      }
      switch (line[0]) {
      case 'p':
         if (sscanf(line, "p%s%i%i%i",
                    name, &m_numNodes, &m_numArcs, &m_numCommodities) != 4) {
            cerr << "ERROR: Input file is incorrect. (p line)" << endl;
            return 1;
         }
         m_problemName = name;
         m_arcs        = new arc[m_numArcs + 
                                 (addDummyArcs ? m_numCommodities : 0)];
         if(!m_arcs)
            throw UtilExceptionMemory("readInstance", "MCF_DecompApp");
         m_commodities = new commodity[m_numCommodities];
         if(!m_commodities)
            throw UtilExceptionMemory("readInstance", "MCF_DecompApp");
         break;
      case 'c':
         break;
      case 'd':
         if (sscanf(line, "d%i%i%i",
                    &m_commodities[commodities_read].source,
                    &m_commodities[commodities_read].sink,
                    &m_commodities[commodities_read].demand) != 3) {
            cerr << "ERROR: Input file is incorrect. (d line)" << endl;
            return 1;
         }
         ++commodities_read;
         break;
      case 'a':
         if (sscanf(line, "a%i%i%i%i%lf",
                    &m_arcs[arcs_read].tail,
                    &m_arcs[arcs_read].head,
                    &m_arcs[arcs_read].lb,
                    &m_arcs[arcs_read].ub,
                    &m_arcs[arcs_read].weight) != 5) {
            cerr << "Input file is incorrect. (a line)" << endl;
            return 1;
         }
         sumweight += fabs(m_arcs[arcs_read].weight);
         ++arcs_read;
         break;
      default:
         if (sscanf(line+1, "%s", name) <= 0) {
            cerr << "Input file is incorrect. (non-recognizable line)" << endl;
            return 1;
         }
         break;
      }
   }
   
   if (!size_read           || 
       arcs_read!=m_numArcs || 
       commodities_read!=m_numCommodities) {
      cerr << "Input file is incorrect."
           << " size_read=" << size_read 
           << " arcs_read=" << arcs_read
           << " commodities_read=" << commodities_read << endl;
      return 1;
   }
   
   if (addDummyArcs) {
      for (int i = 0; i < m_numCommodities; ++i) {
         m_arcs[m_numArcs].tail   = m_commodities[i].source;
         m_arcs[m_numArcs].head   = m_commodities[i].sink;
         m_arcs[m_numArcs].lb     = 0;
         m_arcs[m_numArcs].ub     = m_commodities[i].demand;
         m_arcs[m_numArcs].weight = sumweight+1;
         ++m_numArcs;
      }
   }
   is.close();
   return 0;
}

