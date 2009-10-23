
#ifndef MCKNAP_INCLUDED
#define MCKNAP_INCLUDED

/* ====================================================================== */
#define MCKNAP_RC_OK             0
#define MCKNAP_RC_INF            1
#define MCKNAP_RC_TRIVIAL_MAXSUM 2

/* ======================================================================
                                  definitions
   ====================================================================== */

#define TRACELEVEL  10                /* level of debug information */
#define START       1                /* first test to be run */
#define TESTS       100              /* last test to be run */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <stdarg.h>
/*#include <values.h>*/
#include <math.h>
#include <string.h>
/*#include <values.h>*/
#include <limits.h>
#include <malloc.h>
#define _INCLUDE_POSIX_SOURCE
#ifndef _MSC_VER
#include <sys/times.h>
#include <unistd.h>
#endif


/* ======================================================================
				   macros
   ====================================================================== */

#ifndef _MSC_VER
#define srand(x)    srand48(x)
#define random(x)   (lrand48() % (x))
#else
#define random(x)   (rand() % (x))
#endif

#define SYNC           5   /* when to switch to linear scan in binary scan */
#define MEDIMAX        15
#define MAXSTACK       100
#define MAXLIST        32
#define MAXVTYPE       ULONG_MAX

#define TRUE           1
#define FALSE          0

#define MAXIMIZE       1
#define MINIMIZE       0

#define DET(a1, a2, b1, b2)    ((a1) * (stype) (b2) - (a2) * (stype) (b1))
#define SWAPS(a,b)      { register itemset t; t=*(a); *(a)=*(b); *(b)=t; }
#define SWAPI(a,b)      { register itemrec t; t=*(a); *(a)=*(b); *(b)=t; }
#define SWAPO(a,b)      { register ordrec  t; t=*(a); *(a)=*(b); *(b)=t; }
#define SIZE(a)                          ((int) (((a)->lset+1)-(a)->fset))


/* ======================================================================
				 type declarations
   ====================================================================== */

typedef int           boolean; /* logical variable */
typedef int           ntype;   /* number of stages */


typedef int           itype;   /* item profits and weights */
/*typedef long          stype;*/   /* sum of pofit or weight */
typedef double          stype;   /* sum of pofit or weight */
typedef unsigned long vtype;   /* solution vector */


/* partial vector */
typedef struct {
  stype    psum;
  stype    wsum;
  vtype    vect;
} partvect;

/* item */
typedef struct {
   int i;/*MVG*/
   int j;/*MVG*/
   itype    psum;
   itype    wsum;
} itemrec;

/* set of partial vectors */
typedef struct {
  ntype    size;
  itemrec  *fset;
  itemrec  *lset;
  itemrec  *no;
  itemrec  f,l;
  boolean  used;
} itemset;

/* set of partial vectors */
typedef struct {
  ntype    size;
  partvect *fset;
  partvect *lset;
} partset;

/* set of itemsets */
typedef struct {
  itemset  *fset;
  itemset  *lset;
  ntype    size;
} isetset;

/* order record */
typedef struct {
  itype    dp;
  itype    dw;
  itemset  *ref;
} ordrec;

/* order interval */
typedef struct {
  ordrec   *f;
  ordrec   *l;
} ordintv;

/* order stack */
typedef struct {
  ordintv  intv[MAXSTACK];
  int      level;
  int      optim;
  ordrec   *first;
  ordrec   *last;
  ordrec   *i;
} ordstack;

/* solution record */
typedef struct {
  ntype    size;
  itemset  *set;
} solrec;

/* solution structure */
typedef struct {
  solrec   list[MAXLIST];
  ntype    size;
  stype    psum;
  stype    wsum;
  vtype    vect;
  vtype    vmax;
  ordrec   *a;
  ordrec   *b;
} solstruct;

typedef int (*funcptr) (const void *, const void *);


typedef struct { /* all problem information */
  ntype k;
  ntype n;
  int   type;
  itype range;

  stype capacity;	      /* capacity of knapsack */
  stype dantzig;              /* the dantzig upper bound */
  stype zstar;                /* optimal solution */
  stype summul;	  	      /* sum of multiplications */
  stype antmul;		      /* number of multiplications */
  stype maxmul;               /* max multiplied set */
  stype redusets;             /* sum of reduced sets */
  stype reduitems;            /* sum of items which are tested for reduce */
  stype redukill;             /* sum of tested items which were reduced */
  stype gap;                  /* current gap */
  stype partitions;           /* number of partitions */
  stype domikill;             /* number of dominated-kills */
  stype lpkill;               /* number of lp-kills */
  long  timepar;              /* time used for partitioning */
  long  timesort;             /* time used for sorting of gradients */
  long  time;                 /* time used for all solution */
  long  welldef;              /* is the found solution correct */
  long  checked;              /* optimal solution checked */
  long  iterates;             /* number of iterations to find optimal sol */
} allinfo;




extern int minmcknapSolve(int cap,
			  isetset * head,
			  itemrec * solRec,
			  stype   * minObj);
extern void visitems(itemset *d);
extern void inittrace(char *ext);



/* ======================================================================
				  global variables
   ====================================================================== */
/*MVG: TODO - not reentrant*/

/* #define DEBUG(x) x */
#define DEBUG(x)

/*solstruct solution;*/
/*solstruct optsol;*/

#endif
