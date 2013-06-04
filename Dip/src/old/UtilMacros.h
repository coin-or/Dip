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

#ifndef UTIL_MACROS_INCLUDED
#define UTIL_MACROS_INCLUDED

// =========================================================================
#include "DecompConfig.h"
#include "DecompPortable.h"
class DecompApp;
// =========================================================================
const string UtilSpaces   = " \t\r\n";
const double UtilEpsilon  = 1.0e-6;
// =========================================================================

// =========================================================================
// Memory Macros
// =========================================================================
#define UTIL_DELPTR(x) if(x) {delete    x; x = 0;}
#define UTIL_DELARR(x) if(x) {delete [] x; x = 0;}

// =========================================================================
// Debug Macros
// =========================================================================
#define UTIL_DEBUG(param, level, x) if(param > level) {x fflush(stdout);}
#define UTIL_DEBUG0(x)                                {x fflush(stdout);}

// ------------------------------------------------------------------------- //
template <class T> inline void
UtilPrintVector(const vector<T> & v,
                ostream*          os = &cout)
{
   typename vector<T>::const_iterator it;
   (*os) << "\n";

   for (it = v.begin(); it != v.end(); it++) {
      (*os) << *it << " ";
   }
}

// ------------------------------------------------------------------------- //
template <class T> inline void
UtilPrintList(const list<T> & v,
              ostream*        os = &cout)
{
   typename list<T>::const_iterator it;
   (*os) << "\n";

   for (it = v.begin(); it != v.end(); it++) {
      (*os) << *it << " ";
   }
}


// =========================================================================
// COIN Macros
// TODO: anything that depends on COIN should probably not be in util lib
// =========================================================================
CoinPackedVector* UtilPackedVectorFromDense(const int      len,
      const double* dense,
      const double   etol);
void UtilPackedVectorFromDense(const int          len,
                               const double*      dense,
                               const double       etol,
                               CoinPackedVector& v);

//TODO: now depends on DecompApp!? ugh... then belongs in DecompUtil,
// not Util...
void UtilPrintPackedVector(const CoinPackedVector& v,
                           ostream*                 os  = &cout,
                           DecompApp*               app = 0);

// =========================================================================
// Graph Macros
// =========================================================================

/* -------------------------------------------------------------------------
   --- Assumption: a complete undirected graph,
   ---   (i,j) = (j,i), i!=j (setup for i>j)

   --- Loop thru edges: i: 1 -> n-1, j: 0 -> i-1
   --- Number of edges: m = (n^2 - n)/2

   --- Get the edge index from (i,j):
   ---   INDEX_U(i,j) = i > j ? (i * (i-1) / 2) + j : (j * (j-1) / 2) + i

   --- Get (i,j) from the edge index:
   ---   index = (i * (i-1) / 2) + j
   ---   index = (j * (j-1) / 2) + i
   ---    ax^2 + bx + c = 0 -> x = (-b +- sqrt(b^2 - 4ac)) / 2a
   ---    j = index - (j * (j-1))/2
   ---    i = -1/2 + 1/2 sqrt(1 + 8 * index)

   --- Example: n = 5 (i>j)

   ---  0       1       2      3 = j
   ---  0
   ---  1  0 (1,0)
   ---  2  1 (2,0)  2 (2,1)
   ---  3  3 (3,0)  4 (3,1)  5 (3,2)
   ---  4  6 (4,0)  7 (4,1)  8 (4,2)  9 (4,3)

   --- For the directed version (see EWCP):

   --- Loop thru edges:
   ---  i: 1 -> n-1, j: 0 -> i-1 (i>j)
   ---  j: 1 -> n-1, i: 0 -> j-1 (i<j)

   --- Number of edges: m = (n^2 - n)/2

   ---  0       1        2        3       4 = j
   ---  0           10 (0,1) 11 (0,2) 13 (0,3) 16 (0,4)
   ---  1                    12 (1,2) 14 (1,3) 17 (1,4)
   ---  2                             15 (2,3) 18 (2,4)
   ---  3                                      19 (3,4)
   -------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/
inline int UtilNumEdgesU(const int n)
{
   return ((n * n) - n) / 2;
}

/*-----------------------------------------------------------------------*/
inline int UtilIndexU(const int i, const int j)
{
   return i > j ? (i * (i - 1) / 2) + j : (j * (j - 1) / 2) + i;
}

/*-----------------------------------------------------------------------*/
pair<int, int> UtilBothEndsU(const int index);

/*-----------------------------------------------------------------------*/
inline void UtilPrintEdge(const int   index,
                          ostream*    os = &cout)
{
   pair<int, int> uv = UtilBothEndsU(index);
   (*os) << "(" << setw(2) << uv.first << "," << setw(2) << uv.second << ") ";
}

// =========================================================================
// Fill-In Macros
// =========================================================================

#include "CoinHelperFunctions.hpp"

/*-----------------------------------------------------------------------*/
template <class T> inline void
UtilFillN(T* to, const int size, const T value)
{
   CoinFillN(to, size, value); //dep on COIN
}

/*-----------------------------------------------------------------------*/
template <class T> inline void
UtilFillN(vector<T> & v, const int size, const T value)
{
   std::fill_n(back_inserter(v), size, value);
}

/*-----------------------------------------------------------------------*/
inline void UtilIotaN(int*        first,
                      const int   size,
                      const int   init)
{
   int val = init + size;
   int ii;

   for (ii = size; ii-- != 0; ) {
      first[ii] = --val;
   }
}

/*-----------------------------------------------------------------------*/
//TODO: something faster?
inline void UtilIotaN(vector<int> & first,
                      const int     size,
                      const int     init)
{
   first.reserve(size);
   int i, val = init + size;

   for (i = init; i < val; i++) {
      first.push_back(i);
   }
}

// =========================================================================
// Random Numbers
// =========================================================================

/*-----------------------------------------------------------------------*/
inline double UtilURand(const double a, const double b)
{
   double rand01 = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
   return a + (rand01 * (b - a));
}

/*-----------------------------------------------------------------------*/
inline int UtilURand(const int a, const int b)
{
   double rand01 = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
   return a + static_cast<int>(rand01 * (b - a));
}

// =========================================================================
// Statistics
// =========================================================================
//?? make templates?
/*-----------------------------------------------------------------------*/
inline double UtilAve(const vector<double> & x)
{
   return std::accumulate(x.begin(), x.end(), 0.0) / x.size();
}

/*-----------------------------------------------------------------------*/
inline double UtilAve(const vector<int> & x)
{
   return std::accumulate(x.begin(), x.end(), 0.0) / x.size();
}

/*-----------------------------------------------------------------------*/
inline double UtilAve(const double* x,
                      const int      len)
{
   return std::accumulate(x, x + len, 0.0) / len;
}

// =========================================================================
// Other Macros
// =========================================================================

// ------------------------------------------------------------------------- //
inline double UtilFracPart(const double x)
{
   double floor_x     = floor(x);
   double floor_xplus = floor(x + 0.5);

   if (fabs(floor_xplus - x) < (UtilEpsilon * (fabs(floor_xplus) + 1.0))) {
      return 0.0;
   }

   return x - floor_x;
}

// ------------------------------------------------------------------------- //
int UtilScaleDblToIntArr(const int      arrLen,
                         const double* arrDbl,
                         int*           arrInt,
                         const double   oneDbl,
                         int*           oneInt,
                         const double   epstol = UtilEpsilon);

// ------------------------------------------------------------------------- //
int UtilScaleDblToIntArr(const int      arrLen,
                         const double* arrDbl,
                         int*           arrInt,
                         const double   epstol = UtilEpsilon);



// ------------------------------------------------------------------------- //
inline bool UtilIsZero(const double x,
                       const double etol = 1.0e-8)
{
   return fabs(x) < etol;
}

// ------------------------------------------------------------------------- //
inline string UtilIntToStr(const int i)
{
   stringstream ss;
   ss << i;
   return ss.str();
}

// ------------------------------------------------------------------------- //
template <class T>
void UtilDeleteVectorPtr(vector<T*> & vectorPtr,
                         typename vector<T*>::iterator first,
                         typename vector<T*>::iterator last)
{
   typename vector<T*>::iterator it;

   for (it = first; it != last; it++) {
      delete *it;
   }

   vectorPtr.erase(first, last);
}

// ------------------------------------------------------------------------- //
template <class T> void UtilDeleteVectorPtr(vector<T*> & vectorPtr)
{
   UtilDeleteVectorPtr(vectorPtr, vectorPtr.begin(), vectorPtr.end());
}

// ------------------------------------------------------------------------- //
template <class T> void UtilDeleteListPtr(list<T*> & listPtr,
      typename list<T*>::iterator first,
      typename list<T*>::iterator last)
{
   typename list<T*>::iterator it;

   for (it = first; it != last; it++) {
      delete *it;
   }

   listPtr.erase(first, last);
}

// ------------------------------------------------------------------------- //
template <class T> void UtilDeleteListPtr(list<T*> & listPtr)
{
   UtilDeleteListPtr(listPtr, listPtr.begin(), listPtr.end());
}

// ------------------------------------------------------------------------- //
inline bool UtilIsIntegral(const double x,
                           const double etol = 1.0e-10)
{
   return UtilIsZero(x - floor(x), etol) || UtilIsZero(ceil(x) - x, etol);
}

// ------------------------------------------------------------------------- //
template <class T> inline void UtilNegateArr(const int   arrLen,
      T*          arr)
{
   transform(arr, arr + arrLen, arr, negate<T>());
}

// ------------------------------------------------------------------------- //
template <class T>
struct AddOffset : public unary_function<T, T> {
   T m_n;
   T operator() (const T& k) {
      return k + m_n;
   }
   AddOffset(T n) : m_n (n) {};
};

// ------------------------------------------------------------------------- //
template <class T> inline void UtilAddOffsetArr(const int   arrLen,
      T           offset,
      T*          arr)
{
   transform(arr, arr + arrLen, arr, AddOffset<T>(offset));
}

// ------------------------------------------------------------------------- //
struct Perturb { //: public unary_function
   double m_randLB;
   double m_randUB;
   double operator() (const double& k) {
      return k + UtilURand(m_randLB, m_randUB);
   }
   Perturb(double randLB, double randUB) :
      m_randLB(randLB), m_randUB(randUB) {};
};

// ------------------------------------------------------------------------- //
inline void UtilPerturbCost(const int      seed,
                            const int      arrLen,
                            const double   randLB,
                            const double   randUB,
                            double*        arr)
{
   srand(seed);
   transform(arr, arr + arrLen, arr, Perturb(randLB, randUB));
}

// ------------------------------------------------------------------------- //
inline void UtilFlipRowLtoG(const int    len,
                            double*      els,
                            char&        sense,
                            double&      rhs)
{
   if (sense == 'L') {
      return;
   }

   if (sense == 'G') {
      //C++ row to negate? TODO
      for (int i = 0; i < len; i++) {
         els[i] = -els[i];
      }

      sense = 'L';
      rhs   = -rhs;
   }

   assert(0);
}

// ------------------------------------------------------------------------- //
inline void UtilBoundToSense(const double   lb,
                             const double   ub,
                             const double   inf,
                             char&          sense,
                             double&        rhs,
                             double&        range)
{
   range = 0.0;

   if (lb > -inf) {
      if (ub < inf) {
         rhs = ub;

         if (UtilIsZero(ub - lb)) {
            sense = 'E';
         } else {
            sense = 'R';
            range = ub - lb;
         }
      } else {
         sense = 'G';
         rhs = lb;
      }
   } else {
      if (ub < inf) {
         sense = 'L';
         rhs = ub;
      } else {
         sense = 'N';
         rhs = 0.0;
      }
   }
}

// ------------------------------------------------------------------------- //
inline void UtilSenseToBound(const char     sense,
                             const double   rhs,
                             const double   range,
                             const double   inf,
                             double&        lb,
                             double&        ub)
{
   switch (sense) {
   case 'E':
      lb = rhs;
      ub = rhs;
      break;
   case 'L':
      lb = -inf;
      ub = rhs;
      break;
   case 'G':
      lb = rhs;
      ub = inf;
      break;
   case 'R':
      lb = rhs - range;
      ub = rhs;
      break;
   case 'N':
      lb = -inf;
      ub = inf;
      break;
   }
}

// --------------------------------------------------------------------- //
template<class S, class T> class UtilIsGreaterThan {
public:
   bool operator()( const pair<S, T> & x,
                    const pair<S, T> & y) {
      return x.second > y.second;
   }
};

// --------------------------------------------------------------------- //
template<class S, class T> class UtilIsLessThan {
public:
   bool operator()( const pair<S, T> & x,
                    const pair<S, T> & y) {
      return x.second < y.second;
   }
};

#if 0
// ------------------------------------------------------------------------- //
class UtilIsLessThan {
public:
   bool operator()( const pair<int, double> & x,
                    const pair<int, double> & y) {
      return x.second < y.second;
   }
   bool operator()( const pair< pair<int, int>, double> & x,
                    const pair< pair<int, int>, double> & y) {
      return x.second < y.second;
   }
};

// ------------------------------------------------------------------------- //
class UtilIsGreaterThan {
public:
   bool operator()( const pair<int, double> & x,
                    const pair<int, double> & y) {
      return x.second > y.second;
   }
};
#endif

// ------------------------------------------------------------------------- //
inline string UtilDirSlash()
{
   string slash;
#if defined(_MSC_VER)
   slash = "\\";
#else
   slash = "/";
#endif
   return slash;
}

// ------------------------------------------------------------------------- //
inline void UtilOpenFile(ifstream&    fs,
                         const char* fileName) throw(CoinError)
{
   fs.open(fileName);

   if (!fs) {
      string errMessage = "Error: Filename = ";
      errMessage += fileName;
      errMessage += " failed to open.";
      CoinAssertHint(fs, errMessage.c_str());
   }
}

// =========================================================================
// String Macros
// =========================================================================

// ------------------------------------------------------------------------- //
//trims white space (as defined by UtilSpaces) in-place
inline string& UtilStrTrim(string&        s,
                           const string& t = UtilSpaces)
{
   if (s.size() == 0) {
      return s;
   }

   string::size_type pos = s.find_last_not_of(t);

   if (pos != string::npos) {
      s.erase(pos + 1);
      pos = s.find_first_not_of(' ');

      if (pos != string::npos) {
         s.erase(0, pos);
      }
   } else {
      s.erase(s.begin(), s.end());
   }

   return s;
}

// ------------------------------------------------------------------------- //
// returns a lower case version of the string in-place
inline string& UtilStrToLower(string& s)
{
   // Purify did not like this version:
   //   transform (s.begin(), s.end(), s.begin(), myToLower());
   if (s.size() == 0) {
      return s;
   }

   int i;

   for (i = 0; s[i] != '\0'; i++) {
      s[i] = static_cast<char>(tolower(s[i]));
   }

   return s;
}


// ------------------------------------------------------------------------- //
// returns an upper case version of the string in-place
inline string& UtilStrToUpper(string& s)
{
   // Purify did not like this version:
   //   transform (s.begin(), s.end(), s.begin(), myToUpper());
   int i;

   if (s.size() == 0) {
      return s;
   }

   for (i = 0; s[i] != '\0'; i++) {
      s[i] = static_cast<char>(toupper(s[i]));
   }

   return s;
}


#endif
