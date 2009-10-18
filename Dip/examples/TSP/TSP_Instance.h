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

#ifndef TSP_INSTANCE_INCLUDED
#define TSP_INSTANCE_INCLUDED

// --------------------------------------------------------------------- //
#include "UtilGraphLib.h"
#include "TSP_Concorde.h"
#include "TSP_Boost.h"

// --------------------------------------------------------------------- //
/*!
 * \class TSP_Instance
 * Storage of TSP instance data and utility methods.
 */
// --------------------------------------------------------------------- //

class TSP_Instance{
public:
   /** Data for an instance from TSPLIB. */
   UtilGraphLib    m_graphLib;

   //THINK: these next two are for algos not really input - 
   //       better as members of TSP_DecompApp?
   /** Interface class for Concorde methods. */
   TSP_Concorde    m_concorde;

   //** Interface class for Boost methods. */
   TSP_Boost       m_boost;

   /** The current support graph */
   //Graph           m_sg;

   /** The complete graph G\{m_vert} (needed by MODEL_ONETREE). */
   //Graph           m_cgV;       
   int             m_vert;

   //TODO: access methods
   
private:
   /** @name Copy Constructors */
   /** Disable the default copy constructors. */
   TSP_Instance(const TSP_Instance &);
   TSP_Instance & operator=(const TSP_Instance &);
   
public:
   /** @name Constructor and Destructor */
   TSP_Instance() :
      m_graphLib(),
      m_concorde(),
      m_boost   (),
      //m_sg      (),
      //m_cgV     (),
      m_vert    (0)
   {}
   ~TSP_Instance() {};
   
   
   
};

#endif
