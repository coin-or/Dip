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

#ifndef VRP_CVRPSEP_INCLUDED
#define VRP_CVRPSEP_INCLUDED

//---
//--- Interface class to CVRPSEP (separation routines for CVRP)
//---   Note: all CVRSEP arrays start counting from 1 (not 0)
//---         for n customers, consider n+1 the depot node
//---

// --------------------------------------------------------------------- //
#include "cnstrmgr.h"
#include "capsep.h"

// --------------------------------------------------------------------- //
#include "VRP_GSECCut.h"
#include "VRP_Instance.h"

// --------------------------------------------------------------------- //
// --------------------------------------------------------------------- //
class CVRPsep_LPSol {
public:
   int      nEdges;
   double * EdgeX;
   int    * EdgeTail;
   int    * EdgeHead;

public:
   CVRPsep_LPSol() :
      nEdges  (0),
      EdgeX   (0),
      EdgeTail(0),      
      EdgeHead(0)   
   {};
   
   void clear(){
      UTIL_DELARR(EdgeX);
      UTIL_DELARR(EdgeTail);
      UTIL_DELARR(EdgeHead);
   }
   
   void init(const int n){
      nEdges = n;
      if(nEdges > 0){
	 clear(); //TODO: doing this each time is inefficient
	 EdgeX    = new double[nEdges+1];
	 EdgeTail = new int[nEdges+1];
	 EdgeHead = new int[nEdges+1];
	 CoinAssertHint(EdgeX && EdgeTail && EdgeHead, "Error: Out of Memory");
	 EdgeX[0]    = 0.0;
	 EdgeTail[0] = 0;
	 EdgeHead[0] = 0;
      }
   }
   
   ~CVRPsep_LPSol()
   {
      clear();
   };
};



// --------------------------------------------------------------------- //
// --------------------------------------------------------------------- //
class VRP_CVRPsep {
public:
   const VRP_Instance * m_vrp;     //ptr to vrp instance
   CVRPsep_LPSol        m_lpSol;   //storage of current lp solution
   CnstrMgrPointer      m_newCuts; //my new CVRPsep cuts
   CnstrMgrPointer      m_oldCuts; //my old CVRPsep cuts


public:
   VRP_CVRPsep(const int maxCuts = 500) :
      m_vrp        (0),
      m_lpSol      (),
      m_newCuts    (),
      m_oldCuts    ()
   {
      CMGR_CreateCMgr(&m_newCuts, maxCuts);
      CMGR_CreateCMgr(&m_oldCuts, 10*maxCuts);
   }
   ~VRP_CVRPsep() {
      CMGR_FreeMemCMgr(&m_newCuts);
      CMGR_FreeMemCMgr(&m_oldCuts);
   }
   
public:

   /** @name Helper Functions */
   void init(const VRP_Instance * vrp) {
      m_vrp = vrp;
      //TODO: could alloc cuts here and send in param for size
   }


   void buildLpSol(const double * x,
		   const int      nNzs,
		   const double   etol = 1.0e-8){
      int e;
      int ind = 1;
      int nEdges = m_vrp->m_graphLib.n_edges;
      int nVerts = m_vrp->m_graphLib.n_vertices;

      //---
      //--- free old memory, open new memory of size nszs
      //---    TODO: would be more efficient if we did not
      //---          do this every time - just create lpsol
      //---          vecs of max size (nEdges) and don't
      //---          bother clearing
      //---
      m_lpSol.init(nNzs);//allocates nNsz+1
      for(e = 0; e < nEdges; e++){
	 if(x[e] > etol){
	    pair<int,int> uv = UtilBothEndsU(e);
	    
	    //--- 
	    //--- CVRPSEP convention:
	    //---   1.) edgelist is 1 to NoOfEdges
            //????????? nodes start from - so does ours... 
	    //---   2.) depot is numbered nCustomers+1=nVerts
	    //---
	    m_lpSol.EdgeX[ind] = x[e];
	    m_lpSol.EdgeTail[ind]               
               = uv.second == 0 ? nVerts : uv.second;
	    m_lpSol.EdgeHead[ind]               
               = uv.first  == 0 ? nVerts : uv.first;
#if 0
            printf("XLP[%d,%d -> %d,%d] = %g\n",
		   uv.first, uv.second,
		   m_lpSol.EdgeTail[ind],
		   m_lpSol.EdgeHead[ind], x[e]);
	    UTIL_DEBUG(m_appParam.Log, 5,
		       (*m_osLog) << "XLP[" << uv.first << "," << uv.second;
		       (*m_osLog) << " -> " << m_lpSol.EdgeTail[ind] << ",";
		       (*m_osLog) << m_lpSol.EdgeHead[ind] << "]";
		       (*m_osLog) << " = " << x[e] << endl;
		       );
#endif
		assert(ind <= nNzs);
	    ind++;
	 }
      }
   }



   int sepCapacityCuts(const int maxCuts = 500){   

      const UtilGraphLib & graphLib = m_vrp->m_graphLib;
      const int            nVerts   = graphLib.n_vertices;
      const int          * demand   = graphLib.vertex_wt;
      const int            capacity = graphLib.capacity;
      
      //---
      //--- check to make sure m_newCuts is big enough for max size
      //---
      //CMGR_ExpandCMgr(m_newCuts, maxCuts);//??
      
      char   IntegerAndFeasible;
      double MaxViolation;

      m_newCuts->Size = 0;      //need to do?
      CAPSEP_SeparateCapCuts(nVerts - 1,
			     const_cast<int*>(demand),
			     capacity,
			     m_lpSol.nEdges,
			     m_lpSol.EdgeTail,
			     m_lpSol.EdgeHead, 
			     m_lpSol.EdgeX,
			     m_oldCuts,
			     maxCuts,
			     0.001, //EpsForIntegrality
			     &IntegerAndFeasible,
			     &MaxViolation,
			     m_newCuts);
#if 0
      printf("Found %d capacity cuts. MaxViolation=%g\n", 
             m_newCuts->Size, MaxViolation);
#endif
      return m_newCuts->Size;
   }


   void createVrpCuts(DecompCutList & newCuts){
      int i, j;
      const UtilGraphLib & graphLib = m_vrp->m_graphLib;
      const int            nVerts   = graphLib.n_vertices;
      const int          * demand   = graphLib.vertex_wt;
      const int            capacity = graphLib.capacity;
      
      
      for(i = 0; i < m_newCuts->Size; i++){
	 CnstrPointer CP = m_newCuts->CPL[i];
	 
	 switch(CP->CType){
	 case CMGR_CT_CAP: 
	    {
	       dynamic_bitset<> inS(nVerts);
	       for(j = 1; j <= CP->IntListSize; j++){
		  if(CP->IntList[j] == nVerts)//depot
		     inS.set(0);
		  else
		     inS.set(CP->IntList[j]);
	       }
	       VRP_GSECCut * cut = new VRP_GSECCut(inS, demand, capacity);
	       //cut->print();
	       //double lb = cut->getLowerBound();
	       double ub = cut->getUpperBound();
	       //printf("CUT lb=%g ub=%g\n", lb, ub);
	       //printf("CVRPSEP RHS=%g\n", CP->RHS);
	       assert(UtilIsZero(CP->RHS - ub));
	       newCuts.push_back(cut);	    
	    }
	    break;
	 default:
	    CoinAssert(0);	 
	 }
      }
      
      //---
      //--- move new cuts into old cut pool
      //---
      //---   NOTE: if the framework decides that some cut sould not
      //---         be entered (exceeds max number per pass or something)
      //---         then we don't really want this to go into old cuts yet
      //---
      for(i = 0; i < m_newCuts->Size; i++)
	 CMGR_MoveCnstr(m_newCuts, m_oldCuts, i, 0);
      //cout << "After Move m_newCuts size = " << m_newCuts->Size << endl;
      //cout << "After Move m_oldCuts size = " << m_oldCuts->Size << endl;
   }
};
   
#endif
