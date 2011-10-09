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

#ifndef VRP_CONCORDEI_INCLUDED
#define VRP_CONCORDEI_INCLUDED

//---
//--- Interface class to Concorde (TSP solver)
//---     those prototypes from concorde lib that we use
//---
#define CC_TINYTSP_MINIMIZE            (1)
#define CC_TINYTSP_ERROR               -1
#define CC_TINYTSP_SEARCHLIMITEXCEEDED  1
#define CC_TINYTSP_INFEASIBLE           2
#define CCutil_MAXINT       (2147483647)

extern "C"{
   typedef struct CCrandstate {
      int a;
      int b;
      int arr[55];
   } CCrandstate;

   typedef struct CCdata_user {
      double  *x;
      double  *y;
   } CCdata_user;

   typedef struct CCdata_rhvector {
      int dist_00;
      int dist_01;
      int dist_02;
      int dist_12;
      int dist_22;
      double p;
      int rhlength;
      char *space;
      char **vectors;
   } CCdata_rhvector;

   typedef struct CCdatagroup {
      int    (*edgelen) (int i, int j, struct CCdatagroup *dat);
      double  *x;
      double  *y;
      double  *z;
      int    **adj;
      int     *adjspace;
      int    **len;
      int     *lenspace;
      int     *degree;
      int      norm;
      int      dsjrand_param;
      int      default_len;     /* for edges not in sparse graph   */
      int      sparse_ecount;   /* number of edges in sparse graph */
      double   gridsize;        /* for toroidal norm */
      double   dsjrand_factor;
      CCdata_rhvector rhdat;
      CCdata_user     userdat;
      int      ndepot;          /* used with the subdivision code   */
      int      orig_ncount;     /* just ncount-ndepot               */
      int     *depotcost;       /* cost from each node to the depot */
      int     *orig_names;      /* the nodes names from full problem */
   } CCdatagroup;
   
   void CCutil_sprand (int seed, CCrandstate *r);

   int CCtsp_solve_sparse (int           ncount, 
			   int           ecount, 
			   int         * elist, 
			   int         * elen,
			   int         * in_tour, 
			   int         * out_tour, 
			   double      * in_val, 
			   double      * optval,
			   int         * success, 
			   int         * foundtour, 
			   char        * name, 
			   double      * timebound,
			   int         * hit_timebound, 
			   int           silent, 
			   CCrandstate * rstate);

   int CCtsp_solve_dat (int           ncount, 
			CCdatagroup * indat, 
			int         * in_tour,
			int         * out_tour, 
			double      * in_val, 
			double      * optval, 
			int         * success,
			int         * foundtour, 
			char        * name, 
			double      * timebound, 
			int         * hit_timebound,
			int           silent, 
			CCrandstate * rstate);

   int CCutil_graph2dat_matrix (int           ncount,  
				int           ecount, 
				int         * elist, 
				int         * elen,
				int           defaultlen, 
				CCdatagroup * dat);

   int CCtiny_bnc_tsp (int           ncount, 
                       CCdatagroup * dat, 
                       double      * upbound,
                       double      * optval, 
                       int           nodelimit);

   int CCtiny_bnc_msp (int ncount, int ecount, int *elist, int *elen, int depot,
                       int *lower, int *upper, double *upperbound, int objsense,
                       double *optval, int *xsol, int checkresult, int searchlimit);
   


}
#endif
