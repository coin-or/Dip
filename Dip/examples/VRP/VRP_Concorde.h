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

#ifndef VRP_CONCORDE_INCLUDED
#define VRP_CONCORDE_INCLUDED

//---
//--- Interface class to Concorde (TSP solver)
//---

// --------------------------------------------------------------------- //
#include "VRP_ConcordeI.h"

// --------------------------------------------------------------------- //
#define VRP_CONCORDE_DEBUG       0
#define VRP_CONCORDE_SOLVER_TINY 1
//#define VRP_CONCORDE_SOLVER_TINY 0

// --------------------------------------------------------------------- //
// --------------------------------------------------------------------- //
class ConcordeGraph {
 public:
   int    m_nEdges;
   int    m_nVerts;
   int  * m_edgeList;
   int  * m_edgeValue;

   //index on complete graph to index on complete graph minus d-to-d edges
   //int  * m_edgeMapCtoN;
   //index on complete graph minus d-to-d edges to index on complete graph
   //int  * m_edgeMapNtoC;

 public:
 ConcordeGraph() :
   m_nEdges        (0),
      m_nVerts     (0),
      m_edgeList   (0),
      m_edgeValue  (0)
      //m_edgeMapCtoN(0),
      //m_edgeMapNtoC(0)
{};
   
   void clear(){
      m_nEdges   = 0;
      m_nVerts   = 0;
      UTIL_DELARR(m_edgeList);
      UTIL_DELARR(m_edgeValue);
      //UTIL_DELARR(m_edgeMapCtoN);
      //UTIL_DELARR(m_edgeMapNtoC);
   }

   void setEdgeValue(const double * edgeValue){
      
   }

   void init(const int nVerts,
	     const int nEdges){
      m_nVerts = nVerts;
      m_nEdges = nEdges;
      if(nEdges > 0){
	 m_edgeList    = new int[2 * m_nEdges];
	 m_edgeValue   = new int[    m_nEdges];
         //m_edgeMapCtoN = new int[    m_nEdges];
         //m_edgeMapNtoC = new int[    m_nEdges];
         //UtilFillN(m_edgeMapCtoN, m_nEdges, -1);
         //UtilFillN(m_edgeMapNtoC, m_nEdges, -1);
	 CoinAssertHint(m_edgeList    && m_edgeValue,
                        //m_edgeMapCtoN && m_edgeMapNtoC, 
                        "Error: Out of Memory");
      }
   }
   
   ~ConcordeGraph() {
      clear();
   };
};

// --------------------------------------------------------------------- //
// --------------------------------------------------------------------- //
class VRP_Concorde {
 public:   
   const VRP_Instance * m_vrp;        //ptr to vrp instance
   ConcordeGraph        m_cgExpand;   //complete graph expanded
   int                * m_CCoptTour;  //store optimal tour out of concorde
   double               m_CCoptVal;   //store optimal obj  out of concorde
   //tmp storage of doubles (sz = num of original edges)
   double             * tmpDblNOrigE;   
   //tmp storage of doubles (sz = num of original edges)
   int                * tmpIntNOrigE;

 public:
 VRP_Concorde() :
   m_vrp          (0),
      m_cgExpand  ( ) ,
      m_CCoptTour (0),
      m_CCoptVal  ( ) ,
      tmpDblNOrigE( ) ,
      tmpIntNOrigE( ) {}
   
   ~VRP_Concorde() {
      UTIL_DELARR(m_CCoptTour);
      UTIL_DELARR(tmpDblNOrigE);
      UTIL_DELARR(tmpIntNOrigE);
   }
   
 public:

   /** @name Helper Functions */
   void init(const VRP_Instance * vrp) {
      m_vrp       = vrp;
      buildExpandedCompleteGraph();
      m_CCoptTour = new int[m_cgExpand.m_nVerts];
      CoinAssertHint(m_CCoptTour, "Error: Out of Memory");

      const UtilGraphLib & graphLib  = m_vrp->m_graphLib;
      tmpDblNOrigE = new double[graphLib.n_edges];
      tmpIntNOrigE = new int   [graphLib.n_edges];
      CoinAssertHint(tmpDblNOrigE && tmpIntNOrigE,
		     "Error: Out of Memory");
   }


   void buildExpandedCompleteGraph(){
      //---
      //--- To represent a VRP graph as a TSP instance we need to 
      //--- create (nRoutes-1) additional dummy depots and  their 
      //--- associated customer edges.
      //---
      //--- Original VRP graph has depot=0 and n_vertices-1 customers
      //---    depot => 0
      //---       C  => {1,..,n_vertices-1}
      //--- Now, we will shift the customers back to index:
      //---       C  => {0,..,n_vertices-2}
      //--- and make the depot indices 
      //---    depot => n_vertices-1, .... n_vertices - 1 + nRoutes
      //---
      //--- We do not want or need edges between dummy depots. But, it
      //--- actually makes the accounting easier to have them. So, we will
      //--- include them and treat this like a complete graph but make the
      //--- edge weights on depot-depot edges high (or just fix to 0).
      //---
      //--- NOTE: changed to remove depot-to-depot edges, now using
      //---   a mapping for the book-keeping
      //---
      int i, j, edgeIndex;
      const UtilGraphLib & graphLib  = m_vrp->m_graphLib;
      const int nRoutes              = m_vrp->m_numRoutes;
      int       nCustomers           = graphLib.n_vertices - 1;
      int       nVerts               = nCustomers + nRoutes;
      int       nEdges               = UtilNumEdgesU(nVerts);
      int     * edgeList             = NULL;
      //int     * edgeMapCtoN          = NULL;
      //int     * edgeMapNtoC          = NULL;
      
      //---
      //--- Assumption: a complete undirected graph, 
      //---   (i,j) = (j,i), i!=j (setup for i>j)      
      //--- Loop thru edges: i: 1 -> n-1, j: 0 -> i-1
      //--- Number of edges: m = (n^2 - n)/2 - (nR^2 - nR)/2
      //---
      m_cgExpand.init(nVerts, nEdges);
      edgeIndex   = 0;
      //edgeIndexC  = 0;
      edgeList    = m_cgExpand.m_edgeList;
      //edgeMapCtoN = m_cgExpand.m_edgeMapCtoN;
      //edgeMapNtoC = m_cgExpand.m_edgeMapNtoC;
      for(i = 1; i < nVerts; i++){
         for(j = 0; j < i; j++){
            //printf("(i=%d, j=%d) edgeIndexC:%d edgeIndex:%d", 
            //     i,j,edgeIndexC,edgeIndex);
            //if(i >= nCustomers && j >= nCustomers)
            // printf("  --> depot edge");
            //printf("\n");
            //edgeMapCtoN[edgeIndexC] = edgeIndex;
            //edgeMapNtoC[edgeIndex ] = edgeIndexC;
            //exclude depot-to-depot edges            
            //if(i < nCustomers || j < nCustomers){
               edgeList[2*edgeIndex  ]  = i;
               edgeList[2*edgeIndex+1]  = j;
               CoinAssert(edgeIndex < nEdges);                           
               edgeIndex++;
               //}
               //edgeIndexC++;
	 }
      }
      //TODO: or, open up less, except for m_edgeMapCtoN?
      //      m_cgExpand.m_nEdges = edgeIndex;
      //assert(m_cgExpand.m_nEdges == 
      //     (UtilNumEdgesU(nVerts) - UtilNumEdgesU(nRoutes)));
   }
   
   void setExpandedCost(const double * origGraphCost){
      //---
      //--- origGraphCost contains the costs based on original graph
      //---
      //--- Original VRP graph has depot=0 and n_vertices-1 customers
      //---    depot => 0
      //---       C  => {1,..,n_vertices-1}
      //---
      const UtilGraphLib & graphLib         = m_vrp->m_graphLib;
      const int            nOrigEdges       = graphLib.n_edges;
      const int            nOrigVerts       = graphLib.n_vertices;
      const int            nRoutes          = m_vrp->m_numRoutes;      
      const int            nEdges           = m_cgExpand.m_nEdges;
      const int            nCustomers       = nOrigVerts - 1;
      double             * origGraphCostDbl = tmpDblNOrigE;
      int                * origGraphCostInt = tmpIntNOrigE;
      int                * edgeValue        = m_cgExpand.m_edgeValue;
      int i;
      //---
      //--- offset the orig (reduced costs) costs so all positive
      //---
#if VRP_CONCORDE_DEBUG > 0  
      for(i = 0; i < nOrigEdges; i++){
      	 printf("origGraphCost[%d -> %s]: %10.5f\n", 
                i, UtilEdgeToStr(i).c_str(), origGraphCost[i]);
      }
#endif
      const double * minElement = min_element(origGraphCost, 
					      origGraphCost + nOrigEdges);
      const double   offset     = min(0.0, *minElement);
#if VRP_CONCORDE_DEBUG > 0
      printf("minElement = %10.5f offset = %10.5f\n", *minElement, -offset);
#endif
      memcpy(origGraphCostDbl, origGraphCost, nOrigEdges * sizeof(double));
      if(offset < 0.0){
	 UtilAddOffsetArr(nOrigEdges, -offset, origGraphCostDbl);
      }
      
      //---
      //--- scale to integers (do something very simple)
      //---    just multiply by 1.0e4 and then round
      //---    this should avoid any kind of overflow and
      //---    give accuracy at the 4th decimal place which should be
      //---    sufficient
      //--- NOTE: epsilon=1.0e-6, 1.0e-4 causes overflow issues
      //---
      //int scaleFactor = 
      UtilScaleDblToIntArr(nOrigEdges, 
                           origGraphCostDbl,
                           origGraphCostInt, 1.0e-3);
      //UtilScaleDblToIntArrCC(nOrigEdges, 
      //                                   origGraphCostDbl,
      //                                   origGraphCostInt);
#if VRP_CONCORDE_DEBUG > 0
      //printf("scaleFactor = %d\n", scaleFactor);      
      for(i = 0; i < nOrigEdges; i++){
         printf("origGraphCost[%d]: %10.5f -> %10d\n", 
                i, 
                origGraphCostDbl[i],
                origGraphCostInt[i]);
         //CoinAssert(origGraphCostInt[i] < bigCost);
      }
#endif

      //---
      //--- find max edge length for "bigM" on depot-to-depot edges
      //---   check for integer overflow
      //---
      int    bigCost    = 0;
      double bigCostDbl = 0.0;
      int    maxEdgeLen = *max_element(origGraphCostInt,
                                       origGraphCostInt+nOrigEdges);
      bigCostDbl = ((double)maxEdgeLen+1.0) * (double)m_cgExpand.m_nVerts;
      if((256*bigCostDbl) > (double)CCutil_MAXINT) {
         printf("WARNING: Large edge lengths!!!!\n");
         bigCost = CCutil_MAXINT / 256;         
      }
      else
         bigCost = (maxEdgeLen+1) * m_cgExpand.m_nVerts;
      
      while(maxEdgeLen > bigCost){
         printf("maxEdgeLen=%d > bigCost=%d !!! SCALE BY 10\n",
                maxEdgeLen, bigCost);
         //scale by 10
         for(i = 0; i < nOrigEdges; i++){
            origGraphCostInt[i] /= 10;//and int will trunc it
         }
         maxEdgeLen /= 10;
      }
      assert(maxEdgeLen < bigCost);
         
      //---
      //--- shift the customers back to index
      //---       C  => {0,..,n_vertices-2}
      //--- and make the depot indices 
      //---    depot => n_vertices-1, .... n_vertices - 1 + nRoutes
      //---			         
      //--- set all depot to depot costs to a big number (DON'T NEED)
      //---   and
      //--- duplicate the original depot to customer edges for dummies
      //---
      int origDepot = 0;      
      int origCost  = 0;
      int u, v, d, d1, d2, origIndex, newIndex, custIndex, depot;
      for(u = 1; u < nOrigVerts; u++){
	 for(v = 1; v < u; v++){
	    origIndex = UtilIndexU(u,v);
	    newIndex  = UtilIndexU(u-1, v-1);
#if VRP_CONCORDE_DEBUG > 0
	    printf("origIndex: [%d -> %s] to newIndex: [%d -> %s]\n",
                   origIndex, UtilEdgeToStr(origIndex).c_str(),
                   newIndex,  UtilEdgeToStr(newIndex).c_str());
#endif
            CoinAssert(origGraphCostInt[origIndex] >= 0);            
            CoinAssert(origGraphCostInt[origIndex] < bigCost);
	    edgeValue[newIndex] = origGraphCostInt[origIndex];
	 }
      }
      for(u = 1; u <= nCustomers; u++){
         //get the original cost of customer to depot
	 origIndex = UtilIndexU(u, origDepot);
         origCost  = origGraphCostInt[origIndex];
         //in the concorde graph, the customers are 0..c-1
	 custIndex = u - 1;
         //in the concorde graph, the depots    are c..c+r-1
	 for(d = 0; d < nRoutes; d++){
	    depot    = nCustomers + d; 
            //get the index in the concorde graph based on complete graph
	    //newIndex = edgeMapCtoN[UtilIndexU(custIndex, depot)];
            newIndex = UtilIndexU(custIndex, depot);
	    CoinAssert(newIndex < nEdges);
            //CoinAssert(origGraphCostInt[origIndex] < bigCost);
#if VRP_CONCORDE_DEBUG > 0
            printf("(u:%d,0) origIndex=%d origCost=%d d=%d cInd=%d nInd=%d\n",
                   u, origIndex, origCost, d, 
                   UtilIndexU(custIndex, depot), newIndex);
#endif
	    edgeValue[newIndex] = origCost;
	 }
      }
      for(d1 = 1; d1 < nRoutes; d1++){
	 for(d2 = 0; d2 < d1; d2++){
	    newIndex = UtilIndexU(nCustomers + d1,
				  nCustomers + d2);            
            edgeValue[newIndex] = bigCost;//bigCost (causing overflow)
	 }
      }

      //#if VRP_CONCORDE_DEBUG > 0
      //for(i = 0; i < nEdges; i++){
         //printf("newGraphCostInt[%d -> %s]: %d\n", 
         //     i, UtilEdgeToStr(i).c_str(), edgeValue[i]);
         //printf("newGraphCostInt[%d -> (%d,%d)]: %d\n", 
		//        i, 
      //     edgeList[2*i],
      //        edgeList[2*i+1], 
      //        edgeValue[i]);
      //}
      //#endif
   }

   int solveTSP(vector<int>    & vrpRouteInd,
		vector<double> & vrpRouteEls){


      //---
      //--- NOTE: using this requires linking with Cplex right now
      //---
      //--- int CCtsp_solve_dat (int           ncount, 
      //---                      CCdatagroup * indat, 
      //---                      int         * in_tour,
      //---                      int         * out_tour, 
      //---                      double      * in_val, 
      //---                      double      * optval, 
      //---		         int         * success,
      //---                      int         * foundtour, 
      //---			 char        * name, 
      //---			 double      * timebound, 
      //---			 int         * hit_timebound,
      //---			 int           silent, 
      //---			 CCrandstate * rstate);
      //---
      
      //---
      //--- Description:
      //--- SOLVES the TSP over the graph specfied in the edgelist.
      //---   - elist is an array giving the ends of the edges (in pairs)
      //---   - elen is an array giving the weights of the edges.
      //---   - in_tour gives a starting tour in node node node format (it can
      //---     be NULL)
      //---   - out_tour will return the optimal tour (it can be NULL, if it is
      //---     not NULL then it should point to an array of length at least
      //---     ncount.
      //---   - in_val can be used to specify an initial upperbound (it can be
      //---     NULL)
      //---   - optval will return the value of the optimal tour.
      //---   - success will be set to 1 if the run finished normally, 
      //---     and set to if the search was terminated early (by hitting 
      //---     some predefined limit)
      //---   - foundtour will be set to 1 if a tour has been found (if success
      //---     is 0, then it may not be the optimal tour)
      //---   - name specifes a char string that will be used to name various
      //---     files that are written during the branch and bound search 
      //---     (if it is NULL, then "noname" will be used - this will 
      //---     cause problems in a multithreaded program, so specify a 
      //---     distinct name in that case).
      //---   - silent will suppress most output if set to a nonzero value.
      //---
      int returnCode = 0;
      int silent     = 2;
      int success, foundtour, hit_timebound;
      
      CCrandstate rstate;
      int seed       = 1;
      CCutil_sprand(seed, &rstate);


      //do this once
      int    i, d1, d2;
      const  UtilGraphLib & graphLib = m_vrp->m_graphLib;
      int    nRoutes    = m_vrp->m_numRoutes;      
      int    nOrigVerts = graphLib.n_vertices;
      int    nCustomers = nOrigVerts - 1;
      double upbound    = DecompInf;
      double optval     = DecompInf;      

      //---
      //--- try tiny first, switch if needed... if fails
      //---   for now, full blown solver is overkill... ?
      //---
      int * lower      = new int[m_cgExpand.m_nEdges];
      int * upper      = new int[m_cgExpand.m_nEdges];
      int * xsol       = new int[m_cgExpand.m_nEdges];
      for(i = 0; i < m_cgExpand.m_nEdges; i++){
         lower[i] = 0;
         upper[i] = 1;
         xsol[i]  = 0;
      }
      for(d1 = 1; d1 < nRoutes; d1++){
	 for(d2 = 0; d2 < d1; d2++){
	    upper[UtilIndexU(nCustomers + d1, nCustomers + d2)] = 0;
         }
      }
      returnCode = CCtiny_bnc_msp (m_cgExpand.m_nVerts,
                                   m_cgExpand.m_nEdges,
                                   m_cgExpand.m_edgeList,
                                   m_cgExpand.m_edgeValue,
                                   -1,//depot
                                   lower, 
                                   upper, 
                                   &upbound, 
                                   CC_TINYTSP_MINIMIZE,
                                   &optval, 
                                   xsol,
                                   1,      //checkresult (make option)
                                   1000); //nodelimit

      assert(returnCode != CC_TINYTSP_ERROR);
      assert(returnCode != CC_TINYTSP_INFEASIBLE);
      bool useFullSolver = false;
      if(returnCode == CC_TINYTSP_SEARCHLIMITEXCEEDED){
	 printf("Search limit exceeded - let's go to full solver\n");
	 useFullSolver = true;
      }

      delete [] lower;
      delete [] upper;

      if(useFullSolver){
	 CCdatagroup dataGroup;
	 returnCode = CCutil_graph2dat_matrix(m_cgExpand.m_nVerts,      
					      m_cgExpand.m_nEdges,
					      m_cgExpand.m_edgeList,
					      m_cgExpand.m_edgeValue, 
					      -1, &dataGroup);      
	 char nameC[100];
	 string name = UtilStringRandom(10);
	 strcpy(nameC, name.c_str());
	 
#if VRP_CONCORDE_DEBUG > 0
	 printf("CONCORDE GRAPH:\n");
	 printf("   nVerts=%d nEdges=%d\n", 
		m_cgExpand.m_nVerts,
		m_cgExpand.m_nEdges);
	 for(i =  0; i < m_cgExpand.m_nEdges; i++){
	    printf("(%4d,%4d) val:%4d\n",
		   m_cgExpand.m_edgeList[2*i],
		   m_cgExpand.m_edgeList[2*i+1],
		   m_cgExpand.m_edgeValue[i]);
	 }
#endif
	 
	 returnCode = CCtsp_solve_dat(m_cgExpand.m_nVerts,
				      &dataGroup,
				      NULL, //in_tour
				      m_CCoptTour,
				      NULL, //in_val,
				      &m_CCoptVal,
				      &success,
				      &foundtour,
				      nameC,
				      NULL, //timebound
				      &hit_timebound,
				      silent,
				      &rstate);
	 CoinAssert((success == 1) && (foundtour == 1));
      }

#if VRP_CONCORDE_DEBUG > 0    
      printf("Optimal Route:\n");
#endif

      
      if(!useFullSolver){
	 vector<int> edgeListSol;
	 for(i = 0; i < m_cgExpand.m_nEdges; i++){
	    if(xsol[i]){
#if VRP_CONCORDE_DEBUG > 0
	       printf("(%d,%d) wt:%d\n", 
		      m_cgExpand.m_edgeList[2*i],
		      m_cgExpand.m_edgeList[2*i+1],
		      m_cgExpand.m_edgeValue[i]);
#endif
	       edgeListSol.push_back(m_cgExpand.m_edgeList[2*i]);
	       edgeListSol.push_back(m_cgExpand.m_edgeList[2*i+1]);               
	    }
	 }
	 createVrpRouteFromTspEdgeList(edgeListSol,
				       vrpRouteInd,
				       vrpRouteEls);
      }
      else{
	 createVrpRouteFromTspRoute(m_CCoptTour, 
				    m_cgExpand.m_nVerts,
				    vrpRouteInd,
				    vrpRouteEls);
      }
      return returnCode;
   }
   
   void createVrpRouteFromTspEdgeList(vector<int>    & tspEdgeList,
                                      vector<int>    & vrpRouteInd,
                                      vector<double> & vrpRouteEls){
      
      int i, u, v;
      int        nVerts               = m_vrp->m_graphLib.n_vertices;
      int        nCustomers           = m_vrp->m_graphLib.n_vertices - 1;
      int      * depotEdges           = tmpIntNOrigE;
      
      //---
      //--- keep track of depot edges separate for 1-routes (els=2.0)
      //---
      UtilFillN(depotEdges, nVerts, 0);
#if VRP_CONCORDE_DEBUG > 0    
      printf("Optimal Route:\n");
#endif
      for(i = 0; i < m_cgExpand.m_nVerts; i++){
         //for(i = 0; i < tspRouteLen-1; i++){
         u = tspEdgeList[2*i];
         v = tspEdgeList[2*i+1];
	 u = u >= nCustomers ? 0 : u+1;	 
	 v = v >= nCustomers ? 0 : v+1;
	 CoinAssert((u + v) > 0);
	 if(u == 0){
	    depotEdges[v]++; //(0,v) depot edge
	 }
	 else if(v == 0){
	    depotEdges[u]++; //(u,0) depot edge
	 }
	 else{
	    vrpRouteInd.push_back(UtilIndexU(u,v));
#if VRP_CONCORDE_DEBUG > 0    
            printf("(%d,%d) -> %d\n", u, v, UtilIndexU(u,v));
#endif	
	 } 
      }
      UtilFillN(vrpRouteEls, vrpRouteInd.size(), 1.0);

      for(v = 1; v < nVerts; v++){
	 if(depotEdges[v]){
	    CoinAssert(depotEdges[v] >= 0);
	    CoinAssert(depotEdges[v] <= 2);
	    vrpRouteInd.push_back(UtilIndexU(v, 0));
	    vrpRouteEls.push_back(static_cast<double>(depotEdges[v]));
#if VRP_CONCORDE_DEBUG > 0    
            printf("%d D(%d,%d) -> %d\n", 
                   depotEdges[v], v, 0, UtilIndexU(v,0));
#endif
	 }
      }
   }

   void createVrpRouteFromTspRoute(const int      * tspRoute,
				   const int        tspRouteLen,
				   vector<int>    & vrpRouteInd,
				   vector<double> & vrpRouteEls){
			       
      int i, u, v;
      int        nVerts               = m_vrp->m_graphLib.n_vertices;
      int        nCustomers           = m_vrp->m_graphLib.n_vertices - 1;
      int      * depotEdges           = tmpIntNOrigE;
      
      //---
      //--- keep track of depot edges separate for 1-routes (els=2.0)
      //---
      UtilFillN(depotEdges, nVerts, 0);
#if VRP_CONCORDE_DEBUG > 0    
      printf("Optimal Route Nodes:");
      for(i = 0; i < tspRouteLen; i++){
         printf("%d ", tspRoute[i]);
      }      
      printf("\nOptimal Route Edges:\n");
#endif
      for(i = 0; i < tspRouteLen-1; i++){
	 u = tspRoute[i]   >= nCustomers ? 0 : tspRoute[i]  +1;	 
	 v = tspRoute[i+1] >= nCustomers ? 0 : tspRoute[i+1]+1;
	 CoinAssert((u + v) > 0);
	 if(u == 0){
	    depotEdges[v]++; //(0,v) depot edge
	 }
	 else if(v == 0){
	    depotEdges[u]++; //(u,0) depot edge
	 }
	 else{
	    vrpRouteInd.push_back(UtilIndexU(u,v));
	 }	    
#if VRP_CONCORDE_DEBUG > 0    
	 printf("(%d,%d) -> %d\n", u, v, UtilIndexU(u,v));	 
#endif
      }
      u = tspRoute[tspRouteLen-1] >= nCustomers ? 0 
	 : tspRoute[tspRouteLen-1]+1;	 
      v = tspRoute[0] >= nCustomers ? 0 : tspRoute[0] + 1;
      CoinAssert((u + v) > 0);
      if(u == 0){
	 depotEdges[v]++; //(0,v) depot edge
      }
      else if(v == 0){
	 depotEdges[u]++; //(u,0) depot edge
      }
      else{
	 vrpRouteInd.push_back(UtilIndexU(u,v));
      }	
#if VRP_CONCORDE_DEBUG > 0        
      printf("(%d,%d) -> %d\n", u, v, UtilIndexU(u,v));
#endif
      UtilFillN(vrpRouteEls, vrpRouteInd.size(), 1.0);

      for(v = 1; v < nVerts; v++){
	 if(depotEdges[v]){
	    CoinAssert(depotEdges[v] >= 0);
	    CoinAssert(depotEdges[v] <= 2);
	    vrpRouteInd.push_back(UtilIndexU(v, 0));
	    vrpRouteEls.push_back(static_cast<double>(depotEdges[v]));
#if VRP_CONCORDE_DEBUG > 0    
            printf("%d D(%d,%d) -> %d\n", 
                   depotEdges[v], v, 0, UtilIndexU(v,0));
#endif
	 }
      }
   }

   void createTSPLIBFile(const string fileName){
      ofstream    os;
      int         i, j, edgeIndex;
      const int   nVerts    = m_cgExpand.m_nVerts;
      const int * edgeValue = m_cgExpand.m_edgeValue;

      printf("Open File %s\n", fileName.c_str());

      UtilOpenFile(os, fileName);
      os << "NAME: " << fileName << "\n"
	 << "TYPE: TSP\n"
	 << "DIMENSION: " << nVerts << "\n"
	 << "EDGE_WEIGHT_TYPE: EXPLICIT\n"
	 << "EDGE_WEIGHT_FORMAT: LOWER_DIAG_ROW\n"
	 << "EDGE_WEIGHT_SECTION\n";
      os << " 0\n";
      edgeIndex = 0;
      for(i = 1; i < nVerts; i++){
	 for(j = 0; j < i; j++){
	    os << edgeValue[edgeIndex] << " ";
	    edgeIndex++;
	 }
	 os << " 0\n";
      }
      os << "EOF\n";

      os.close();
   }


};

   
#endif
