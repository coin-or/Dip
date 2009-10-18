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

#ifndef VRP_INSTANCE_INCLUDED
#define VRP_INSTANCE_INCLUDED

// --------------------------------------------------------------------- //
#include "UtilGraphLib.h"

// --------------------------------------------------------------------- //
/*!
 * \class VRP_Instance
 * Storage of VRP instance data and utility methods.
 */
// --------------------------------------------------------------------- //

class VRP_Instance{
public:
   /** Data for an instance from VRPLIB. */
   UtilGraphLib    m_graphLib;
   int             m_numRoutes;
   
private:
   /** @name Copy Constructors */
   /** Disable the default copy constructors. */
   VRP_Instance(const VRP_Instance &);
   VRP_Instance & operator=(const VRP_Instance &);
   
public:
   /** @name Constructor and Destructor */
   VRP_Instance() :
      m_graphLib (),
      m_numRoutes(0)
   {}
   ~VRP_Instance() {};   
};

#endif
