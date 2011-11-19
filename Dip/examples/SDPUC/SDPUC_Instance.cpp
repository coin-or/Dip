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
#include "SDPUC_Instance.h"
#include "UtilMacrosDecomp.h"
//===========================================================================//

//===========================================================================//
int SDPUC_Instance::readInstance(string & fileName,
                               bool     addDummyArcs){
   
   ifstream is;   
   int      status = UtilOpenFile(is, fileName.c_str());   
   if(status)
      throw UtilException("Failed to read instance",
                          "readInstance", "MCF_Instance");
   
   double sumweight        = 0;
   bool   size_read        = true;
   int    arcs_read        = 0;
   int    nodes_read	   = 0;
   int	  ts_read   = 0;
   int	  nt = 0;
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
         if (sscanf(line, "p%s%i%i%i%i%i",
                    name, &m_numNodes, &m_numArcs, &m_numSwitchings, &m_numTimeseries, &m_numTimeperiods) != 6) {
            cerr << "ERROR: Input file is incorrect. (p line)" << endl;
            return 1;
         }
         m_problemName = name;
         m_arcs        = new arc[m_numArcs + (addDummyArcs ? 0 : 0)];
         if(!m_arcs)
            throw UtilExceptionMemory("readInstance", "MCF_DecompApp");
         m_nodes = new node[m_numNodes];
         if(!m_nodes)
            throw UtilExceptionMemory("readInstance", "MCF_DecompApp");
		 m_timeseries = new timeseries[m_numTimeseries];
         if(!m_timeseries)
            throw UtilExceptionMemory("readInstance", "MCF_DecompApp");

         break;
      case 'c':
         break;
	  case '#':
		 break;
      case 'n':
         if (sscanf(line, "n%i%lf%i",
                    &m_nodes[nodes_read].id,
                    &m_nodes[nodes_read].demand,
                    &m_nodes[nodes_read].tsdemand) != 3) {
            cerr << "ERROR: Input file is incorrect. (n line)" << endl;
            return 1;
         }
         ++nodes_read;
         break;
      case 'a':
         if (sscanf(line, "a%i%i%lf%lf%lf%lf%lf%lf%i%i%i%i",
                    &m_arcs[arcs_read].tail,
                    &m_arcs[arcs_read].head,
                    &m_arcs[arcs_read].lb,
                    &m_arcs[arcs_read].ub,
                    &m_arcs[arcs_read].weight,
					&m_arcs[arcs_read].mcost,
					&m_arcs[arcs_read].fcost1,
					&m_arcs[arcs_read].fcost2,
					&m_arcs[arcs_read].tscap,
					&m_arcs[arcs_read].tscost,
					&m_arcs[arcs_read].acline,
					&m_arcs[arcs_read].switchable
					) != 12) {
            cerr << "Input file is incorrect. (a line)" << endl;
            return 1;
         }
         sumweight += fabs(m_arcs[arcs_read].mcost);
         ++arcs_read;
         break;
	   case 't':
		  //cout << "ts_read=" << ts_read ;
		  //cout << " numTimeperiods=" << m_numTimeperiods << endl;
		  m_timeseries[ts_read].values = new double[m_numTimeperiods];
		/* if (sscanf(line, "t%i%lf%lf%lf%lf",
                    &m_timeseries[ts_read].id,
                    &m_timeseries[ts_read].values[0],
                    &m_timeseries[ts_read].values[1],
					&m_timeseries[ts_read].values[2],
					&m_timeseries[ts_read].values[3]
					) != 5) {
            cerr << "ERROR: Input file is incorrect. (t line) << " << line << endl;
            return 1;	
         }*/
		 
		 nt = 0;
		 char * pch;
		 //printf ("Splitting string \"%s\" into tokens:\n",line);
		 pch = strtok (line,"\t");  //stripping the initial 't'
		 //printf ("%s ",pch);
		 pch = strtok (NULL, "\t"); //timeseries id
		 m_timeseries[ts_read].id = atoi(pch);
		 //printf ("%s\n",pch);
		 while (pch != NULL && nt < m_numTimeperiods)
		 {
			
			pch = strtok (NULL, "\t");
			m_timeseries[ts_read].values[nt] = atof(pch);
			//printf ("%s\n",pch);
			nt++;
		 }



         ++ts_read;
         break;
      default:
         if (sscanf(line+1, "%s", name) <= 0) {
            cerr << "Input file is incorrect. (non-recognizable line)" << endl;
            return 1;
         }
         break;
      }
   }
   
   if (!size_read			    || 
       arcs_read  != m_numArcs  || 
       nodes_read != m_numNodes ||
	   ts_read	  != m_numTimeseries
	   ) {
      cerr << "Input file is incorrect."
           << " size_read=" << size_read 
           << " arcs_read=" << arcs_read
           << " nodes_read=" << nodes_read 
		   << " ts_read=" << ts_read << endl;
      return 1;
   }
   
   /*if (addDummyArcs) {
      for (int i = 0; i < m_numCommodities; ++i) {
         m_arcs[m_numArcs].tail   = m_commodities[i].source;
         m_arcs[m_numArcs].head   = m_commodities[i].sink;
         m_arcs[m_numArcs].lb     = 0;
         m_arcs[m_numArcs].ub     = m_commodities[i].demand;
         m_arcs[m_numArcs].weight = sumweight+1;
         ++m_numArcs;
      }
   }*/
   is.close();
   return 0;
}

