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

#ifndef TSP_CONCORDEI_INCLUDED
#define TSP_CONCORDEI_INCLUDED

// --------------------------------------------------------------------- //
typedef struct CCtsp_skeleton {
   int   atomcount;
   int * atoms;
} CCtsp_skeleton;

typedef struct CCtsp_segment {
   int lo;
   int hi;
} CCtsp_segment;

typedef struct CCtsp_lpclique {
   int                    segcount;
   struct CCtsp_segment * nodes;
   int                    hashnext;
   int                    refcount;
} CCtsp_lpclique;

typedef struct CCtsp_lpcut_in {
   int                    cliquecount;
   int                    dominocount;
   int                    rhs;
   char                   sense;
   char                   branch;
   CCtsp_lpclique        *cliques;
   struct CCtsp_lpdomino        *dominos;
   CCtsp_skeleton         skel;
   struct CCtsp_lpcut_in *next;
   struct CCtsp_lpcut_in *prev;
} CCtsp_lpcut_in;

#define CC_FOREACH_NODE_IN_CLIQUE(i,c,tmp) \
    for(tmp=0;tmp<(c).segcount;tmp++) \
        for(i=(c).nodes[tmp].lo;i<=(c).nodes[tmp].hi;i++)




#if 0


typedef struct CCtsp_lpdomino {
   CCtsp_lpclique        sets[2];
   int                   hashnext;
   int                   refcount;
} CCtsp_lpdomino;

#endif


extern int CCtsp_exact_subtours (CCtsp_lpcut_in ** cuts, 
				 int             * cutcount, 
				 int               ncount,
				 int               ecount, 
				 int             * elist, 
				 double          * x);
extern void
CCtsp_free_lpcut_in (CCtsp_lpcut_in *c);

 
#endif
