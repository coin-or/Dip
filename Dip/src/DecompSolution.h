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

#ifndef DECOMP_SOLUTION_INCLUDED
#define DECOMP_SOLUTION_INCLUDED

//TODO: make this public AlpsDecompSolution?
class DecompSolution {
protected:
   /** Length of solution (number of columns). */
   int      m_size;

   /** Solution values. */
   double * m_values;

   /** Quality of solution (bound wrt to objective). */
   double   m_quality;

public:
   /** @name Helper functions (public). */   

   /** Get length of solution. */
   inline const int getSize() const {return m_size;}
   
   /** Get solution values. */
   inline const double * getValues() const {return m_values;}
   
   /** Get quality of solution. */
   inline const double getQuality() const {return m_quality;}

public:
   /** Print solution. */
   virtual void print(int       precision = 2,
		      ostream & os        = cout) const {
      int i;
      os << setprecision(precision);
      os << setiosflags(ios::fixed|ios::showpoint)
         << setw(14);
      
      os << "-------------------------" << endl;
      os << "Quality = " << getQuality() << endl;
      for (i = 0; i < m_size; i++) {
         if (!UtilIsZero(m_values[i])){
	    os << setw(6) << i << " " << m_values[i] << endl;
         }
      }
      os << "-------------------------" << endl;
      os << resetiosflags(ios::fixed|ios::showpoint|ios::scientific); 
   }

   /** Print solution in MIPLIB2010 solution checker format. */
   virtual void print(const vector<string> & colNames,
		      int                    precision = 2,
		      ostream              & os        = cout) const {
      int i;
      os << setprecision(precision);
      os << setiosflags(ios::fixed|ios::showpoint);
      
      //os << "-------------------------" << endl;
      os << "=obj= " << getQuality() << endl;
      for (i = 0; i < m_size; i++) {
         if (!UtilIsZero(m_values[i])){
	    os << colNames[i] << " " << m_values[i] << endl;
         }
      }
      //os << "-------------------------" << endl;
      os << resetiosflags(ios::fixed|ios::showpoint|ios::scientific); 
   }
   
public:
   /** @name Copy Constructors */
   DecompSolution(const DecompSolution & source) :
      m_size(source.m_size),
      m_values(0),
      m_quality(source.m_quality)
      {
         m_values = new double[m_size];
         CoinAssertHint(m_values, "Error: Out of Memory"); 
         memcpy(m_values, source.m_values, m_size * sizeof(double));         
      }
   DecompSolution & operator=(const DecompSolution & rhs)
      {
         
         if(this != &rhs){
            m_size    = rhs.m_size;
            m_quality = rhs.m_quality;
            
            m_values = new double[m_size];
            CoinAssertHint(m_values, "Error: Out of Memory"); 
            memcpy(m_values, rhs.m_values, m_size * sizeof(double));         
         }
         return *this;
      }
   
 public:
   /** @name Constructor and Destructor */

   /** Default constructor. Takes size of solution. */
   DecompSolution() :
      m_size(0),
      m_values(0),
      m_quality(1e75)
   {
   }
   
   /** Constructor. */
   DecompSolution(const int      size,
                  const double * values,
                  const double   quality) :
      m_size(size),
      m_values(0),
      m_quality(quality)
   {
      CoinAssert(m_size > 0);
      m_values = new double[m_size];
      CoinAssertHint(m_values, "Error: Out of Memory"); 
      memcpy(m_values, values, m_size * sizeof(double));
   }

   DecompSolution(const int      size,
                  const double * values,
		  const double * cost) :
      m_size(size),
      m_values(0),
      m_quality(0.0)
   {
      CoinAssert(m_size > 0);
      m_values = new double[m_size];
      CoinAssertHint(m_values, "Error: Out of Memory"); 
      memcpy(m_values, values, m_size * sizeof(double));

      //---
      //--- calculate quality
      //---
      for(int i = 0; i < size; i++)
	 m_quality += cost[i] * values[i];
   }
   
   virtual ~DecompSolution() {  
      UTIL_DELARR(m_values);      
   };
};

#endif
