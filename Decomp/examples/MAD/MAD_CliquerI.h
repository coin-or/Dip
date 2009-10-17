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

#ifndef MAD_CLIQUERI_INCLUDED
#define MAD_CLIQUERI_INCLUDED

// --------------------------------------------------------------------- //
#ifndef ELEMENTSIZE
typedef unsigned long int setelement;
# if (ULONG_MAX == 65535)
#  define ELEMENTSIZE 16
# elif (ULONG_MAX == 4294967295)
#  define ELEMENTSIZE 32
# else
#  define ELEMENTSIZE 64
# endif
#endif  /* !ELEMENTSIZE */
typedef setelement * set_t;

// --------------------------------------------------------------------- //
/*** Counting amount of 1 bits in a setelement ***/

/* Array for amount of 1 bits in a byte. */
static int set_bit_count[256] = {
   0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,
   1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
   1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
   1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
   3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
   1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
   3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
   3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
   3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
   4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8 };

/* The following macros assume that all higher bits are 0.
 * They may in some cases be useful also on with other ELEMENTSIZE's,
 * so we define them all.  */
#define SET_ELEMENT_BIT_COUNT_8(a)  (set_bit_count[(a)])
#define SET_ELEMENT_BIT_COUNT_16(a) (set_bit_count[(a)>>8] + \
				     set_bit_count[(a)&0xFF])
#define SET_ELEMENT_BIT_COUNT_32(a) (set_bit_count[(a)>>24] + \
				     set_bit_count[((a)>>16)&0xFF] + \
				     set_bit_count[((a)>>8)&0xFF] + \
				     set_bit_count[(a)&0xFF])
#define SET_ELEMENT_BIT_COUNT_64(a) (set_bit_count[(a)>>56] + \
				     set_bit_count[((a)>>48)&0xFF] + \
				     set_bit_count[((a)>>40)&0xFF] + \
				     set_bit_count[((a)>>32)&0xFF] + \
				     set_bit_count[((a)>>24)&0xFF] + \
				     set_bit_count[((a)>>16)&0xFF] + \
				     set_bit_count[((a)>>8)&0xFF] + \
				     set_bit_count[(a)&0xFF])
#if (ELEMENTSIZE==64)
# define SET_ELEMENT_BIT_COUNT(a) SET_ELEMENT_BIT_COUNT_64(a)
# define FULL_ELEMENT ((setelement)0xFFFFFFFFFFFFFFFF)
#elif (ELEMENTSIZE==32)
# define SET_ELEMENT_BIT_COUNT(a) SET_ELEMENT_BIT_COUNT_32(a)
# define FULL_ELEMENT ((setelement)0xFFFFFFFF)
#elif (ELEMENTSIZE==16)
# define SET_ELEMENT_BIT_COUNT(a) SET_ELEMENT_BIT_COUNT_16(a)
# define FULL_ELEMENT ((setelement)0xFFFF)
#else
# error "SET_ELEMENT_BIT_COUNT(a) not defined for current ELEMENTSIZE"
#endif

// --------------------------------------------------------------------- //
/*
 * Gives a value with bit x (counting from lsb up) set.
 *
 * Making this as a table might speed up things on some machines
 * (though on most modern machines it's faster to shift instead of
 * using memory).  Making it a macro makes it easy to change.
 */
#define SET_BIT_MASK(x) ((setelement)1<<(x))

/* Set handling macros */
#define SET_ELEMENT_CONTAINS(e,v)   ((e)&SET_BIT_MASK(v))
#define SET_ADD_ELEMENT(s,a) \
                       ((s)[(a)/ELEMENTSIZE] |= SET_BIT_MASK((a)%ELEMENTSIZE))
#define SET_DEL_ELEMENT(s,a) \
                       ((s)[(a)/ELEMENTSIZE] &= ~SET_BIT_MASK((a)%ELEMENTSIZE))
#define SET_CONTAINS_FAST(s,a) (SET_ELEMENT_CONTAINS((s)[(a)/ELEMENTSIZE], \
						      (a)%ELEMENTSIZE))
#define SET_CONTAINS(s,a) (((a)<SET_MAX_SIZE(s))?SET_CONTAINS_FAST(s,a):0)

/* Sets can hold values between 0,...,SET_MAX_SIZE(s)-1 */
#define SET_MAX_SIZE(s) ((s)[-1])

/* Sets consist of an array of SET_ARRAY_LENGTH(s) setelements */
#define SET_ARRAY_LENGTH(s) (((s)[-1]+ELEMENTSIZE-1)/ELEMENTSIZE)

/*
 * set_size()
 *
 * Returns the number of elements in set s.
 */
static inline int set_size(set_t s) {
   int count=0;
   setelement * c;
   for (c=s; c < s+SET_ARRAY_LENGTH(s); c++)
      count+=SET_ELEMENT_BIT_COUNT(*c);
   return count;
}


/*
 * set_new()
 *
 * Create a new set that can hold values in the range 0,...,size-1.
 */
static inline set_t set_new(int size) {
   int n;
   set_t s;
   n=(size/ELEMENTSIZE+1)+1;
   s=(setelement*)calloc(n,sizeof(setelement));
   s[0]=size;        
   return &(s[1]);
}

/*
 * set_free()
 *
 * Free the memory associated with set s.
 */
static inline void set_free(set_t s) {
   free(&(s[-1]));
}

/*
 * set_duplicate()
 *
 * Returns a newly allocated duplicate of set s.
 */
static inline set_t set_duplicate(set_t s) {
   set_t newT;

   newT=set_new(SET_MAX_SIZE(s));
   memcpy(newT,s,SET_ARRAY_LENGTH(s)*sizeof(setelement));
   return newT;
}

// --------------------------------------------------------------------- //
typedef struct _graph_t graph_t;
struct _graph_t {
   int     n;             /* Vertices numbered 0...n-1      */
   set_t * edges;         /* A list of n sets (the edges).  */
   int   * weights;       /* A list of n vertex weights.    */
};

// --------------------------------------------------------------------- //
#define GRAPH_ADD_EDGE(g,i,j) do {            \
	SET_ADD_ELEMENT((g)->edges[(i)],(j)); \
	SET_ADD_ELEMENT((g)->edges[(j)],(i)); \
} while (0)
#define GRAPH_DEL_EDGE(g,i,j) do {            \
	SET_DEL_ELEMENT((g)->edges[(i)],(j)); \
	SET_DEL_ELEMENT((g)->edges[(j)],(i)); \
} while (0)

#define GRAPH_IS_EDGE_FAST(g,i,j)  (SET_CONTAINS_FAST((g)->edges[(i)],(j)))

// --------------------------------------------------------------------- //
extern graph_t * graph_new  (int n);
extern void      graph_free (graph_t * g);
extern void      graph_print(graph_t * g);
extern int       graph_edge_count(graph_t *g);

static inline int graph_subgraph_weight(graph_t * g, set_t s) {
   unsigned int i,j;
   int count=0;
   setelement e;
   for (i=0; i<SET_ARRAY_LENGTH(s); i++) {
      if (s[i]) {
         e=s[i];
         for (j=0; j<ELEMENTSIZE; j++) {
            if (e&1)
               count+=g->weights[i*ELEMENTSIZE+j];
            e = e>>1;
         }
      }
   }
   return count;
}

// --------------------------------------------------------------------- //
typedef struct _clique_options clique_options;
struct _clique_options {
   int * (*reorder_function)(graph_t *, int);
   int * reorder_map;
   
   /* arguments:  level, n, max, user_time, system_time, opts */
   int (*time_function)(int, int, int, int, double, double, clique_options *);
   FILE * output;
   
   int   (*user_function)(set_t, graph_t *, clique_options *);
   void   * user_data;
   set_t  * clique_list;
   int      clique_list_length;
};
extern clique_options * clique_default_options;

// --------------------------------------------------------------------- //
/*
 * clique_find_all()
 *
 * Find all cliques with weight at least min_weight and at most max_weight.
 *
 *   g          - the graph
 *   min_weight - minimum weight of cliques to search for.  If min_weight==0,
 *                searches for maximum weight cliques.
 *   max_weight - maximum weight of cliques to search for.  If max_weight==0,
 *                no upper limit is used.  If min_weight==0, max_weight must
 *                also be 0.
 *   maximal    - require cliques to be maximal cliques
 *   opts       - time printing and clique storage options
 *
 * Returns the number of cliques found.  This can be less than the number
 * of cliques in the graph iff opts->time_function() or opts->user_function()
 * returns FALSE (request abort).
 *
 * The cliques found are stored in opts->clique_list[] and
 * opts->user_function() is called with them (if non-NULL).  The cliques
 * stored in opts->clique_list[] are newly allocated, and can be freed
 * by set_free().
 *
 * Note: Automatically uses clique_unweighted_find_all if all vertex
 *       weights are the same.
 */
int clique_find_all(graph_t        * g, 
                    int              min_weight, 
                    int              max_weight,
                    int              maximal, 
                    clique_options * opts);


/*
 * clique_find_single()
 *
 * Returns a clique with weight at least min_weight and at most max_weight.
 *
 *   g          - the graph
 *   min_weight - minimum weight of clique to search for.  If min_weight==0,
 *                searches for a maximum weight clique.
 *   max_weight - maximum weight of clique to search for.  If max_weight==0,
 *                no upper limit is used.  If min_weight==0, max_weight must
 *                also be 0.
 *   maximal    - require returned clique to be maximal
 *   opts       - time printing options
 *
 * Returns the set of vertices forming the clique, or NULL if a clique
 * of requested weight/maximality does not exist in the graph  (or if
 * opts->time_function() requests abort).
 *
 * The returned clique is newly allocated and can be freed by set_free().
 *
 * Note: Does NOT use opts->user_function() or opts->clique_list[].
 * Note: Automatically uses clique_unweighted_find_single if all vertex
 *       weights are the same.
 */
set_t clique_find_single(graph_t        * g,
                         int              min_weight,
                         int              max_weight,
			 int              maximal,
                         clique_options * opts);

/*
 * clique_unweighted_find_all()
 *
 * Find all cliques with size at least min_size and at most max_size.
 *
 *   g        - the graph
 *   min_size - minimum size of cliques to search for.  If min_size==0,
 *              searches for maximum cliques.
 *   max_size - maximum size of cliques to search for.  If max_size==0, no
 *              upper limit is used.  If min_size==0, this must also be 0.
 *   maximal  - require cliques to be maximal cliques
 *   opts     - time printing and clique storage options
 *
 * Returns the number of cliques found.  This can be less than the number
 * of cliques in the graph iff opts->time_function() or opts->user_function()
 * returns FALSE (request abort).
 *
 * The cliques found are stored in opts->clique_list[] and
 * opts->user_function() is called with them (if non-NULL).  The cliques
 * stored in opts->clique_list[] are newly allocated, and can be freed
 * by set_free().
 */
int clique_unweighted_find_all(graph_t        * g,
                               int              min_size,
                               int              max_size,
			       int              maximal,
                               clique_options * opts);

#endif
