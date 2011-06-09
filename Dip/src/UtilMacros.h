//===========================================================================//
// This file is part of the DIP Solver Framework.                            //
//                                                                           //
// DIP is distributed under the Eclipse Public License as part of the        //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)       //
//                                                                           //
// Conceptual Design: Matthew Galati, SAS Institute Inc.                     //
//                    Ted Ralphs, Lehigh University                          //
//                                                                           //
// Copyright (C) 2002-2011, Lehigh University, Matthew Galati, Ted Ralphs    //
// All Rights Reserved.                                                      //
//===========================================================================//

#ifndef UTIL_MACROS_INCLUDED
#define UTIL_MACROS_INCLUDED

// =========================================================================
#include <cstdio>
#include <cassert>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <algorithm>
#include <functional>
#include <string>
#include <map>
#include <climits>
#include <cmath>
#include <cstring>
#include <ctime>
#include <memory>
using namespace std;

// =========================================================================
const string UtilSpaces   = " \t\r\n";
const double UtilEpsilon  = 1.0e-6;
const double UtilTooBig   = 1.0e20;
const double UtilSmallerThanTooBig = 1.0e19;
#ifndef INT_MAX
#define INT_MAX (static_cast<int>((~(static_cast<unsigned int>(0))) >> 1))
#endif

#ifndef round
#define round(x) floor(x+0.5)
#endif

// =========================================================================

// =========================================================================
// Util Error Codes
// =========================================================================
enum UtilStatus {
   UtilStatusOk = 0,
   UtilStatusFileIO 
};

// =========================================================================
// Memory Macros                                                            
// =========================================================================
#define UTIL_DELPTR(x) if(x) {delete    x; x = 0;}
#define UTIL_DELARR(x) if(x) {delete [] x; x = 0;}

// =========================================================================
// Debug Macros                                                             
// =========================================================================
#ifdef NDEBUG
//use with LogDebugLevel
#define UTIL_DEBUG(param, level, x)
//#define UTIL_DEBUG(param, level, x) if(param >= level) {x fflush(stdout);}
#else
#define UTIL_DEBUG(param, level, x) if(param >= level) {x fflush(stdout);}
#endif

//use with LogLevel
#define UTIL_MSG(param, level, x)   if(param >= level) {x fflush(stdout);}


// ------------------------------------------------------------------------- //
#ifndef NDEBUG
#define UtilAssert(expression,errorMsg,os) assert(expresssion)
#else
inline void UtilAssert(bool      expression,
		       string    errorMsg,
		       ostream * os){
   //---
   //--- this is a forced assertion (even when -NDEBUG)
   //---   
   if(!expression){
      (*os) << "ERROR:"  << errorMsg << endl;
      abort();
   }
}
#endif

// ------------------------------------------------------------------------- //
inline void UtilPrintParameter(ostream      * os,
			       const string & section,
			       const string & name,
			       const int      value){
   (*os) << left << setw(15) << section
	 << left << setw(25) << name
	 << setw(10) << value  << endl;
}

// ------------------------------------------------------------------------- //
inline void UtilPrintParameter(ostream      * os,
			       const string & section,
			       const string & name,
			       const double   value){
   (*os) << left << setw(15) << section
	 << left << setw(25) << name
	 << setw(10) << value  << endl;
}

// ------------------------------------------------------------------------- //
template <class T> inline void 
UtilPrintVector(const vector<T> & v,
                ostream         * os = &cout) {
   typename vector<T>::const_iterator it;
   for(it = v.begin(); it != v.end(); it++){
      (*os) << *it << " ";
   }
   (*os) << "\n";
}

// ------------------------------------------------------------------------- //
template <class T> inline void 
UtilPrintVector(const vector<T>      & v,
		const vector<string> & label,
                ostream              * os = &cout) {
   typename vector<T>::const_iterator it;
   for(it = v.begin(); it != v.end(); it++){
      (*os) << setw(5) << *it << " -> " 
	    << setw(25) << label[*it] << endl;
   }
}

// ------------------------------------------------------------------------- //
template <class T> inline void 
UtilPrintList(const list<T> & v,
              ostream       * os = &cout) {
   typename list<T>::const_iterator it;

   (*os) << "\n";
   for(it = v.begin(); it != v.end(); it++){
      (*os) << *it << " ";
   }
}

// =========================================================================
// Graph Macros                                                             
// =========================================================================

//TODO: should this be in a UtilGraph.h or something?
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

   ---          0       1       2      3 = j
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
inline int UtilNumEdgesU(const int n) { 
   return ((n * n) - n) / 2; 
}

/*-----------------------------------------------------------------------*/
inline int UtilIndexU(const int i, const int j){ 
   return i > j ? (i * (i-1) / 2) + j : (j * (j-1) / 2) + i;
}

/*-----------------------------------------------------------------------*/
pair<int,int> UtilBothEndsU(const int index);

/*-----------------------------------------------------------------------*/
inline void UtilPrintEdge(const int   index, 
                          ostream   * os = &cout){
   pair<int,int> uv = UtilBothEndsU(index);
   (*os) << "(" << setw(2) << uv.first << "," << setw(2) << uv.second << ") ";
}

/*-----------------------------------------------------------------------*/
inline string UtilEdgeToStr(const int index){
   stringstream  ss;
   pair<int,int> uv = UtilBothEndsU(index);
   ss << "(" << setw(2) << uv.first << "," << setw(2) << uv.second << ") ";
   return ss.str();
}

// =========================================================================
// Fill-In Macros                                                           
// =========================================================================

// =========================================================================
template <class T> inline void 
UtilFillN(T * to, const int size, const T value){
   int i;
   for(i = 0; i < size; i++)
      to[i] = value;
}

// =========================================================================
template <class T> inline void
UtilFillN(vector<T> & v, const int size, const T value){
   std::fill_n(back_inserter(v), size, value);
}

/*-----------------------------------------------------------------------*/
inline void UtilIotaN(int       * first, 
                      const int   size, 
                      const int   init){
   int val = init + size;
   int ii;
   for(ii = size; ii-- != 0; ){
      first[ii] = --val;
   }
}

/*-----------------------------------------------------------------------*/
inline void UtilIotaN(vector<int> & first, 
                      const int     size, 
                      const int     init){
   first.reserve(size);
   int i, val = init + size;
   for(i = init; i < val; i++){
      first.push_back(i);
   }
}

// =========================================================================
// Random Numbers                                                           
// =========================================================================

/*-----------------------------------------------------------------------*/
inline double UtilURand(const double a, const double b){
   double rand01 = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
   return a + (rand01 * (b - a));
} 

/*-----------------------------------------------------------------------*/
inline int UtilURand(const int a, const int b){
   double rand01 = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
   return a + static_cast<int>(rand01 * (b-a));
} 

/*-----------------------------------------------------------------------*/
inline double UtilNormRand(const double mean,
			   const double sigma){
   //http://mathworld.wolfram.com/Box-MullerTransformation.html
   double rand01a = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
   double rand01b = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
   const double pi = 3.14159265358979323846;
   double       z1 = sqrt(-2.0 * log(rand01a)) * cos(2.0 * pi * rand01b);
   return z1 * sigma + mean;
}

// =========================================================================
// Statistics                                                               
// =========================================================================

// ------------------------------------------------------------------------- //
inline double UtilAve(const vector<double> & x){
   return std::accumulate(x.begin(), x.end(), 0.0) / 
      static_cast<double>(x.size());
}

// ------------------------------------------------------------------------- //
inline double UtilAve(const vector<int> & x){
   return std::accumulate(x.begin(), x.end(), 0.0) / 
      static_cast<double>(x.size());
}

// ------------------------------------------------------------------------- //
inline double UtilAve(const double * x,
                      const int      len){
   return std::accumulate(x, x+len, 0.0) / static_cast<double>(len);
}

// =========================================================================
// String Macros                                                             
// =========================================================================

// ------------------------------------------------------------------------- //
inline void UtilStringTokenize(string const   & input,
			       string const   & delimiters,
			       vector<string> & tokens) {
   using namespace std;
   string::size_type last_pos = 0;
   string::size_type pos = 0;
   while(true){
      pos = input.find_first_of(delimiters, last_pos);
      if( pos == string::npos ){
	 tokens.push_back(input.substr(last_pos));
	 break;
      }
      else{
	 tokens.push_back(input.substr(last_pos, pos - last_pos));
	 last_pos = pos + 1;
      }
   }
}

// ------------------------------------------------------------------------- //
inline string UtilStringRandom(int iLength) {
   string strReturn;
   srand( (unsigned int)time(NULL) );
   for( int i = 0 ; i < iLength ; ++i ) {
      int iNumber;
      iNumber = rand()%122;
      if( 48 > iNumber )
         iNumber += 48;
      if( ( 57 < iNumber ) && ( 65 > iNumber ) )
         iNumber += 7;
      if( ( 90 < iNumber ) && ( 97 > iNumber ) )
         iNumber += 6;
      strReturn += (char)iNumber;
   }  
   srand(1);
   return strReturn;
} 

// ------------------------------------------------------------------------- //
//trims white space (as defined by UtilSpaces) in-place
inline string & UtilStrTrim(string       & s,
                            const string & t = UtilSpaces) {
   if(s.size() == 0)
      return s;
   string::size_type pos = s.find_last_not_of(t);
   if(pos != string::npos) {
      s.erase(pos + 1);
      pos = s.find_first_not_of(t);
      if(pos != string::npos) s.erase(0, pos);
   }
   else s.erase(s.begin(), s.end());
   return s;
}

// ------------------------------------------------------------------------- //
// returns a lower case version of the string in-place
inline string & UtilStrToLower(string & s) {
   // Purify did not like this version:
   //   transform (s.begin(), s.end(), s.begin(), myToLower());
   if(s.size() == 0)
      return s;
   int i;
   for (i = 0; s[i] != '\0'; i++)
      s[i] = static_cast<char>(tolower(s[i]));
   return s;
}


// ------------------------------------------------------------------------- //
// returns an upper case version of the string in-place
inline string & UtilStrToUpper(string & s) {
   // Purify did not like this version:
   //   transform (s.begin(), s.end(), s.begin(), myToUpper());   
   int i;
   if(s.size() == 0)
      return s;
   for (i = 0; s[i] != '\0'; i++)
      s[i] = static_cast<char>(toupper(s[i]));
   return s;
}

// =========================================================================
// Other Macros
// =========================================================================

// ------------------------------------------------------------------------- //
template <class T> inline int UtilGetSize(const vector<T> & vec){
  return static_cast<int>(vec.size());
}

// ------------------------------------------------------------------------- //
inline bool UtilIsInSet(const int   value,
                        const int * set,
                        const int   setSize){
   int  i;
   bool inSet = false;
   for(i = 0; i < setSize; i++){
      if(set[i] == value){
         inSet = true;
         break;
      }
   }
   return inSet;
}

// ------------------------------------------------------------------------- //
inline int UtilNumNonzeros(const double * x,
			   const int      len,
			   const double   etol = 1.0e-8){
   int i;
   int nzs = 0;
   for(i = 0; i < len; i++){
      if(fabs(x[i]) > etol)
	 nzs++;
   }
   return nzs;
}

// ------------------------------------------------------------------------- //
inline double UtilFracPart(const double x){
   double floor_x     = floor(x);
   double floor_xplus = floor(x + 0.5);
   if(fabs(floor_xplus - x) < (UtilEpsilon * (fabs(floor_xplus) + 1.0)))
      return 0.0;
   return x - floor_x;
}

// ------------------------------------------------------------------------- //
int UtilScaleDblToIntArr(const int      arrLen,
                         const double * arrDbl,
                         int          * arrInt,
                         const double   oneDbl,
                         int          * oneInt,
                         const double   epstol = UtilEpsilon);

// ------------------------------------------------------------------------- //
int UtilScaleDblToIntArr(const int      arrLen,
                         const double * arrDbl,
                         int          * arrInt,
                         const double   epstol = UtilEpsilon);

// ------------------------------------------------------------------------- //
inline bool UtilIsZero(const double x,
                       const double etol = 1.0e-8){
   return fabs(x) < etol;
}

// ------------------------------------------------------------------------- //
inline string UtilIntToStr(const int i){
   stringstream ss;
   ss << i;
   return ss.str();
}

// ------------------------------------------------------------------------- //
inline string UtilDblToStr(const double x, 
			   const int    precision = -1,
                           const double tooBig    = UtilSmallerThanTooBig){
   stringstream ss;
   if(fabs(x) > tooBig){
      if(x < 0)
	 ss << "-INF";
      else
	 ss << " INF";
   }
   else{
      if(precision >= 0){
	 ss << setiosflags(ios::fixed|ios::showpoint);   
	 ss << setprecision(precision);
      }
      ss << x;
   }
   return ss.str();
}

// ------------------------------------------------------------------------- //
inline void UtilPrintMemUsage(ostream  * os       = &cout,
			      int        logLevel = 0,
			      int        logLimit = 2){
// This doesn't build in gcc 4.5 (at least on MinGW)
#if 0
#if not defined(_MSC_VER)
   if(logLevel >= logLimit){
      struct mallinfo memInfo = mallinfo();
      double memUsage = static_cast<double>(memInfo.uordblks +
					    memInfo.hblkhd) / 1024.0;
      memUsage /= 1024.0;
      (*os) << "memUsage = " << UtilDblToStr(memUsage, 2) << " MB\n";
   }
#endif
#endif
}

// ------------------------------------------------------------------------- //
template <class T>
void UtilDeleteVectorPtr(vector<T*> & vectorPtr,
                         typename vector<T*>::iterator first,
                         typename vector<T*>::iterator last){
   
   typename vector<T*>::iterator it;
   for(it = first; it != last; it++)
      delete *it;
   vectorPtr.erase(first,last);
}

// ------------------------------------------------------------------------- //
template <class T> void UtilDeleteVectorPtr(vector<T*> & vectorPtr){
   UtilDeleteVectorPtr(vectorPtr, vectorPtr.begin(), vectorPtr.end());
}

// ------------------------------------------------------------------------- //
template <class T> void UtilDeleteListPtr(list<T*> & listPtr,
                                          typename list<T*>::iterator first,
                                          typename list<T*>::iterator last){
   
   typename list<T*>::iterator it;
   for(it = first; it != last; it++)
      delete *it;
   listPtr.erase(first,last);
}

// ------------------------------------------------------------------------- //
template <class T> void UtilDeleteListPtr(list<T*> & listPtr){
   UtilDeleteListPtr(listPtr, listPtr.begin(), listPtr.end());
}

// ------------------------------------------------------------------------- //
template <class S, class T>
void UtilDeleteMapPtr(map<S,T*> & mapPtr,
		      typename map<S,T*>::iterator first,
		      typename map<S,T*>::iterator last){  
   typename map<S,T*>::iterator it;
   for(it = first; it != last; it++)
      delete (*it).second;
   mapPtr.erase(first,last);
}

// ------------------------------------------------------------------------- //
template <class S, class T> void UtilDeleteMapPtr(map<S,T*> & mapPtr){
   UtilDeleteMapPtr(mapPtr, mapPtr.begin(), mapPtr.end());
}

// ------------------------------------------------------------------------- //
template <class S, class T>
void UtilDeleteMapVecPtr(map<S, vector<T*> > & mapPtr,
		      typename map<S, vector<T*> >::iterator first,
		      typename map<S, vector<T*> >::iterator last){  
   typename map<S, vector<T*> >::iterator it;
   for(it = first; it != last; it++)
      UtilDeleteVectorPtr((*it).second);
   mapPtr.erase(first,last);
}

// ------------------------------------------------------------------------- //
template <class S, class T> 
void UtilDeleteMapVecPtr(map<S,vector<T*> > & mapPtr){
   UtilDeleteMapVecPtr(mapPtr, mapPtr.begin(), mapPtr.end());
}

// ------------------------------------------------------------------------- //
inline bool UtilIsIntegral(const double x,
                           const double etol = 1.0e-10){
   return UtilIsZero(x - floor(x), etol) || UtilIsZero(ceil(x) - x, etol);
}

// ------------------------------------------------------------------------- //
inline bool UtilIsIntegral(const double * x,
                           const int      len,
                           const double   etol = 1.0e-10){
   int i;
   for(i = 0; i < len; i++){
      if(!UtilIsIntegral(x[i], etol))
         return false;
   }
   return true;
}

// ------------------------------------------------------------------------- //
template <class T> inline void UtilNegateArr(const int   arrLen,
                                             T         * arr){
   transform(arr, arr + arrLen, arr, negate<T>());
}

// ------------------------------------------------------------------------- //
template <class T> 
struct AddOffset : public unary_function<T,T>
{
   T m_n;   
   T operator() (const T & k) { return k + m_n; }
   AddOffset(T n) : m_n (n) {};
};

// ------------------------------------------------------------------------- //
template <class T> inline void UtilAddOffsetArr(const int   arrLen,
                                                T           offset,
                                                T         * arr){
   transform(arr, arr + arrLen, arr, AddOffset<T>(offset));
}

// ------------------------------------------------------------------------- //
struct Perturb //: public unary_function
{
    double m_randLB;
    double m_randUB;   
    double operator() (const double & k) {
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
			    double       * arr){
    srand(seed);
    transform(arr, arr + arrLen, arr, Perturb(randLB, randUB));
}

// ------------------------------------------------------------------------- //
inline void UtilFlipRowLtoG(const int    len,
                            double     * els,
                            char       & sense,
                            double     & rhs){
   if(sense == 'L')
      return;
  
   if(sense == 'G'){
      for(int i = 0; i < len; i++)
         els[i] = -els[i];
      sense = 'L';
      rhs   = -rhs;    
   }
   assert(0);
}

// ------------------------------------------------------------------------- //
inline void UtilBoundToSense(const double   lb, 
                             const double   ub,
                             const double   inf,
                             char         & sense, 
                             double       & rhs,
                             double       & range) {
  
   range = 0.0;
   if (lb > -inf) {
      if (ub < inf) {
         rhs = ub;
         if (UtilIsZero(ub-lb)) {
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
                             double       & lb, 
                             double       & ub) {

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
   bool operator()( const pair<S,T> & x, 
                    const pair<S,T> & y){
      return x.second > y.second;
   }
};

// --------------------------------------------------------------------- //
template<class S, class T> class UtilIsLessThan {
public:
   bool operator()( const pair<S,T> & x, 
                    const pair<S,T> & y){
      return x.second < y.second;
   }
};

// ------------------------------------------------------------------------- //
inline string UtilDirSlash(){
   string slash;
#if defined(_MSC_VER)
   slash = "\\"; 
#else
   slash = "/"; 
#endif
   return slash;
}

// ------------------------------------------------------------------------- //
inline int UtilOpenFile(ofstream   & os,
                        const char * fileName) {
   int status = UtilStatusOk;
   os.open(fileName);
   if(!os){
      string errMessage = "Error: Filename = ";
      errMessage += fileName;
      errMessage += " failed to open.";  
      cerr << errMessage.c_str() << endl; fflush(stdout);
      status = UtilStatusFileIO;
      assert(os);
   }
   return status;
}

// ------------------------------------------------------------------------- //
inline int UtilOpenFile(ifstream   & is,
                        const char * fileName) {
   int status = UtilStatusOk;
   is.open(fileName);
   if(!is){
      string errMessage = "Error: Filename = ";
      errMessage += fileName;
      errMessage += " failed to open.";
      cerr << errMessage.c_str() << endl; fflush(stdout);
      status = UtilStatusFileIO;
      assert(is);
   }
   return status;
}

// ------------------------------------------------------------------------- //
inline int UtilOpenFile(ofstream     & os,
                       const string & fileName) {
   return UtilOpenFile(os, fileName.c_str());
}

// ------------------------------------------------------------------------- //
inline int UtilOpenFile(ifstream     & is,
                        const string & fileName) {
   return UtilOpenFile(is, fileName.c_str());
}



#endif
