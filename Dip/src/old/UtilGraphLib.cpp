//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)       //
//                                                                           //
// Copyright (C) 2002-2007, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//
/*--------------------------------------------------------------------------
  ---
  --- Author  : Matthew Galati
  --- Purpose : A class for storing data instances from TSPLIB and VRPLIB.

  --- Notes:
  --- read_data() - converted from Ted Ralph's Symphony/VRP C code
  --- assumes the use of a complete symmetric undirected graph,
  --- i.e., |E| = (|V| * |V| - |V|) / 2

  ---
-------------------------------------------------------------------------- */


#include "UtilGraphLib.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
using namespace std;

/* TODO - deal with logging, error handling, etc */
// ==========================================================================
void UtilGraphLib::read_data(const char* datafile)
{
   ifstream is(datafile); /* ??? */

   if (!is) {
      cerr << "UtilGraphLib::read_data failed to open " << datafile << endl;
      exit(1);
   }

   enum DIST {_EXPLICIT, _EUC_2D, _EUC_3D, _MAX_2D, _MAX_3D,
              _MAN_2D, _MAN_3D, _CEIL_2D, _GEO, _ATT
             };
   const int LENGTH      = 255;
   const int KEY_NUM     = 41;
   const int NCTYPE_NUM  = 3;
   const int WTYPE_NUM   = 10;
   const int WFORMAT_NUM = 10;
   const int DTYPE_NUM   = 3;
   const double MY_PI    = 3.141592;
   //This section lists the names of the possible fields in the data file
   static char keywords[KEY_NUM][22] = {
      "NAME", "NAME:", "TYPE", "TYPE:", "COMMENT", "COMMENT:",
      "DIMENSION", "DIMENSION:", "CAPACITY", "CAPACITY:",
      "EDGE_WEIGHT_TYPE", "EDGE_WEIGHT_TYPE:",
      "EDGE_WEIGHT_FORMAT", "EDGE_WEIGHT_FORMAT:",
      "DISPLAY_DATA_TYPE", "DISPLAY_DATA_TYPE:",
      "EDGE_WEIGHT_SECTION", "EDGE_WEIGHT_SECTION:",
      "DISPLAY_DATA_SECTION", "DISPLAY_DATA_SECTION:",
      "NODE_COORD_SECTION", "NODE_COORD_SECTION:",
      "NODE_COORD_TYPE", "NODE_COORD_TYPE:",
      "DEPOT_SECTION", "DEPOT_SECTION:",
      "CAPACITY_VOL", "CAPACITY_VOL:",
      "DEMAND_SECTION", "DEMAND_SECTION:",
      "TIME_WINDOW_SECTION", "TIME_WINDOW_SECTION:",
      "STANDTIME_SECTION", "STANDTIME_SECTION:",
      "PICKUP_SECTION", "PICKUP_SECTION:",
      "EOF", "EOF.", "", "", "NO_MORE_TYPE"
   };
   //This section lists the possible node coordinate data types
   static char nctypes[NCTYPE_NUM][14] = {"TWOD_COORDS", "THREED_COORDS",
                                          "NO_COORDS"
                                         };
   //This is a list of the possible data types for edge weights
   static char wtypes[WTYPE_NUM][9] = {
      "EXPLICIT", "EUC_2D", "EUC_3D",
      "MAX_2D", "MAX_3D", "MAN_2D", "MAN_3D", "CEIL_2D", "GEO", "ATT"
   };
   //This is a list of the possible formats that the edge weight matrix
   //could be given in if it is given explicitly
   static char wformats[WFORMAT_NUM][20] = {
      "UPPER_ROW", "LOWER_ROW", "UPPER_DIAG_ROW", "LOWER_DIAG_ROW",
      "UPPER_COL", "LOWER_COL", "UPPER_DIAG_COL", "LOWER_DIAG_COL",
      "FULL_MATRIX", "FUNCTION"
   };
   //This is a list of the various display data types
   static char dtypes[DTYPE_NUM][14] = {"COORD_DISPLAY", "TWOD_DISPLAY",
                                        "NO_DISPLAY"
                                       };
   FILE* f;

   //if ((err = fopen(&f, datafile, "r")) != 0){
   if ((f = fopen(datafile, "r")) == 0) {
      cerr << "UtilGraphLib::read_data ERROR I/O : reading " << datafile
           << ". Aborting." << endl;
      abort();
   }

   /*
     This loop reads in the next line of the data file and compares it
     to the list of possible keywords to determine what data will follow.
     It then reads the data into the appropriate field and iterates
   */
   char line1[LENGTH], line[LENGTH], key[30], tmp[80];
   int wtype = -1, wformat = -1, dtype = -1, nctype = -1;
   int depot, k, l, m, i, j, node, *coef2;
   double deg, min, x, y, fdummy;
   double coord_x, coord_y, coord_z;
   bool capacity_volume = false;

   while (0 != fgets( line1, LENGTH, f)) {
      strcpy(key, "");
      sscanf(line1, "%s", key); //read in next keyword
      int k;

      for (k = 0; k < KEY_NUM; k++)
         if (strcmp(keywords[k], key) == 0) {
            break;
         }

      if (k == KEY_NUM) {
         cerr << "UtilGraphLib::read_data ERROR I/O : unknown keyword "
              << key << ". Aborting." << endl;
         abort();
      }

      //This is a bit shift operation that divides k by 2 since in the list
      //of keywords, there are two possible formats for the keyword
      k >>= 1;

      if (strchr(line1, ':')) {
         strcpy(line, strchr(line1, ':') + 1);
      }

      switch (k) {
      case 0: { //NAME (set name)
         if (!sscanf(line, "%s", tmp)) {
            cerr << "UtilGraphLib::read_data ERROR I/O : reading NAME"
                 << tmp << ". Aborting." << endl;
            abort();
         }

         //cout << "PROBLEM NAME: \t\t" << tmp << endl;
         string tmp_str(tmp);
         this->name = tmp_str.substr(0, tmp_str.find_first_of("."));
         break;
      }
      case 3 : //DIMENSION (set n_vertices, n_edges)

         if (!sscanf(line, "%d", &k)) {
            cerr << "UtilGraphLib::read_data ERROR I/O : reading DIMENSION "
                 << k << ". Aborting." << endl;
            abort();
         }

         this->n_vertices = k;
         this->n_edges = ((n_vertices * n_vertices) - n_vertices) / 2;
         break;
      case 4 : //CAPACITY (set capacity)

         if (!sscanf(line, "%d", &k)) {
            cerr << "GrabLib:: read_data ERROR I/O : reading CAPACITY "
                 << k << ". Aborting." << endl;
            abort();
         }

         this->capacity = k;
         break;
      case 5 : //EDGE_WEIGHT_TYPE
         sscanf(line, "%s", tmp);

         for (wtype = 0; wtype < WTYPE_NUM; wtype++)
            if (strcmp(wtypes[wtype], tmp) == 0) {
               break;
            }

         if (wtype == WTYPE_NUM) {
            cerr << "GrabLib::read_data ERROR I/O : unknown weight type "
                 << tmp << ". Aborting." << endl;
            abort();
         }

         break;
      case 6 : //EDGE_WEIGHT_FORMAT
         sscanf(line, "%s", tmp);

         for (wformat = 0; wformat < WFORMAT_NUM; wformat++)
            if (strcmp(wformats[wformat], tmp) == 0) {
               break;
            }

         if (wformat == WFORMAT_NUM) {
            cerr << "UtilGraphLib::read_data ERROR I/O : unknown weight format "
                 << tmp << ". Aborting." << endl;
            abort();
         }

         break;
      case 7 : //DISPLAY_DATA_TYPE
         sscanf(line, "%s", tmp);

         for (dtype = 0; dtype < DTYPE_NUM; dtype++)
            if (strcmp(dtypes[dtype], tmp) == 0) {
               break;
            }

         if (dtype == DTYPE_NUM) {
            cerr << "UtilGraphLib::read_data ERROR I/O : unknown display type "
                 << tmp << ". Aborting." << endl;
            abort();
         }

         break;
      case 8: //EDGE_WEIGHT_SECTION (open memory and set edge_wt)

         if (wtype != _EXPLICIT) {
            break;
         }

         edge_wt = new int[n_edges];

         switch (wformat) {
         case 1 : //LOWER_ROW
         case 4 : //UPPER_COL
         case 3 : //LOWER_DIAG_ROW
         case 6 : //UPPER_DIAG_COL

            for (i = 0, coef2 = edge_wt; i < n_vertices; i++) {
               for (j = 0; j  < i; j++, coef2++) {
                  if (!fscanf(f, "%lf", &fdummy)) {
                     cerr << "UtilGraphLib::read_data ERROR I/O : not enough data "
                          << "-- DIMENSION or "
                          << "EDGE_WEIGHT_TYPE declared wrong. Aborting." << endl;
                     abort();
                  } else {
                     *coef2 = (int)fdummy;
                  }
               }

               if ((wformat == 3 || wformat == 6) && !fscanf(f, "%lf", &fdummy)) {
                  cerr << "UtilGraphLib::read_data ERROR I/O : not enough data "
                       << "-- DIMENSION or EDGE_WEIGHT_TYPE declared wrong. "
                       << "Aborting." << endl;
                  abort();
               }
            }

            if (fscanf(f, "%lf", &fdummy)) {
               cerr << "UtilGraphLib::read_data ERROR I/O : too much data "
                    << "-- DIMENSION or "
                    <<       "EDGE_WEIGHT_TYPE declared wrong. Aborting." << endl;
               abort();
            }

            break;
         case 0 : //UPPER_ROW
         case 5 : //LOWER_COL
         case 2 : //UPPER_DIAG_ROW
         case 7 : //LOWER_DIAG_COL

            for (i = 0, coef2 = edge_wt; i < n_vertices; i++) {
               if (wformat == 2 || wformat == 7)
                  if (!fscanf(f, "%lf", &fdummy)) {
                     cerr << "UtilGraphLib::read_data ERROR I/O : not enough data "
                          << "-- DIMENSION or EDGE_WEIGHT_TYPE declared wrong. "
                          << "Aborting." << endl;
                     abort();
                  }

               for (j = i + 1; j < n_vertices; j++) {
                  if (!fscanf(f, "%lf", &fdummy)) {
                     cerr << "UtilGraphLib::read_data ERROR I/O : not enough data "
                          << "-- DIMENSION or EDGE_WEIGHT_TYPE declared wrong. "
                          << "Aborting." << endl;
                     abort();
                  } else {
                     coef2[j*(j-1)/2+i] = (int)fdummy;
                  }
               }
            }

            if (fscanf(f, "%lf", &fdummy)) {
               cerr << "UtilGraphLib::read_data ERROR I/O : too much data "
                    << "-- DIMENSION or EDGE_WEIGHT_TYPE declared wrong. "
                    << "Aborting." << endl;
               abort();
            }

            break;
         case 8 : //FULL_MATRIX

            for (i = 0, coef2 = edge_wt; i < n_vertices; i++) {
               for (j = 0; j <= i; j++)
                  if (!fscanf(f, "%lf", &fdummy)) {
                     cerr << "UtilGraphLib::read_data ERROR I/O : not enough data "
                          << "-- DIMENSION or EDGE_WEIGHT_TYPE declared wrong. "
                          << "Aborting." << endl;
                     abort();
                  }

               for (j = i + 1; j < n_vertices; j++) {
                  if (!fscanf(f, "%lf", &fdummy)) {
                     cerr << "UtilGraphLib::read_data ERROR I/O : not enough data "
                          << "-- DIMENSION or EDGE_WEIGHT_TYPE declared wrong. "
                          << "Aborting." << endl;
                     abort();
                  }

                  coef2[j*(j-1)/2+i] = (int) fdummy;
               }
            }

            if (fscanf(f, "%lf", &fdummy)) {
               cerr << "UtilGraphLib::read_data ERROR I/O : too much data "
                    << "-- DIMENSION or EDGE_WEIGHT_TYPE declared wrong. "
                    << "Aborting." << endl;
               abort();
            }

            break;
         }

         break;
      case 9 : //DISPLAY_DATA_SECTION (open memory and set posx, posy)

         if (dtype != 1) {
            cerr << "UtilGraphLib::read_data ERROR I/O : DISPLAY_DATA_SECTION "
                 << "exists but not TWOD_DISPLAY. Aborting." << endl;
            abort();
         }

         posx = new int[n_vertices];
         posy = new int[n_vertices];

         for (i = 0; i < n_vertices; i++) {
            if ((k = fscanf(f, "%d%lf%lf", &node, &x, &y)) != 3) {
               cerr << "UtilGraphLib::read_data ERROR I/O : error reading "
                    << "DISPLAY_DATA" << endl;
               break;
            }

            posx[node - 1] = (int)(x + 0.5);
            posy[node - 1] = (int)(y + 0.5);
         }

         if (fscanf(f, "%lf", &fdummy)) {
            cerr << "UtilGraphLib::read_data ERROR I/O : too much display data"
                 << endl;
            break;
         }

         break;
      case 10 : //NODE_COORD_SECTION (open memory and set posx,
         //posy, coordx, coordy, coordz)
         if (nctype == -1) {
            nctype = 0;   //if not given: TWOD_COORDS
         }

         if (dtype == -1 && ((wtype == _EUC_2D) ||  //display type not defd
                             (wtype == _MAX_2D) ||  //yet && can disp
                             (wtype == _MAN_2D) ||
                             (wtype == _ATT) )) {
            dtype = 0;   //COORD_DISPLY
         }

         if (dtype == 0) {
            posx = new int[n_vertices];
            posy = new int[n_vertices];
         }

         coordx = new double[n_vertices];
         coordy = new double[n_vertices];

         if (nctype == 1) {
            coordz = new double[n_vertices];
         }

         for (i = 0; i < n_vertices; i++) {
            if (nctype == 0)          //TWOD_COORDS
               if (fscanf(f, "%d%lf%lf", &node, &coord_x, &coord_y) != 3) {
                  cerr << "UtilGraphLib::read_data ERROR I/O : error reading "
                       << "NODE_COORD. Aborting." << endl;
                  abort();
               }

            if (nctype == 1)          //THREED_COORDS
               if (fscanf(f, "%d%lf%lf%lf", &node, &coord_x, &coord_y, &coord_z) != 4) {
                  cerr << "UtilGraphLib::read_data ERROR I/O : error reading "
                       << "NODE_COORD. Aborting." << endl;
                  abort();
               }

            coordx[node-1] = coord_x;
            coordy[node-1] = coord_y;

            //since position is an integer and coord is
            //a double, round off here if dtype is EXPLICIT
            if (dtype == 0) {
               posx[node-1] = (int)coord_x;
               posy[node-1] = (int)coord_y;
            }

            if (nctype == 1) {
               coordz[node-1] = coord_z;
            }

            if (wtype == _GEO) { //GEO
               deg = (int)(coordx[node-1]);
               min = coordx[node-1] - deg;
               coordx[node-1] = MY_PI * (deg + 5.0 * min / 3.0 ) / 180.0;
               deg = (int)(coordy[node-1]);
               min = coordy[node-1] - deg;
               coordy[node-1] = MY_PI * (deg + 5.0 * min / 3.0 ) / 180.0;
            }
         }

         if (fscanf(f, "%d%lf%lf%lf", &node, &coord_x, &coord_y, &coord_z)) {
            cerr << "UtilGraphLib::read_data ERROR I/O: too much data in "
                 << "NODE_COORD. Aborting." << endl;
            abort();
         }

         break;
      case 11: //NODE_COORD_TYPE
         sscanf(line, "%s", tmp);

         for (nctype = 0; nctype < NCTYPE_NUM; nctype++)
            if (strcmp(nctypes[nctype], tmp) == 0) {
               break;
            }

         if (nctype == NCTYPE_NUM) {
            cerr << "UtilGraphLib::read_data ERROR I/O : unknown node_coord_type"
                 << tmp << ". Aborting." << endl;
            abort();
         }

         break;
      case 12: //DEPOT_SECTION
         fscanf(f, "%d", &k);

         if (k != 1) {
            cerr << "UtilGraphLib::read_data ERROR I/O : depot must be node 1."
                 << "Aborting." << endl;
            abort();
         }

         depot = k - 1;

         while (-1 != k) {
            fscanf(f, "%d", &k);
         }

         break;
      case 13: //CAPACITY_VOL
         sscanf(line, "%d", &k);
         capacity_volume = true;
         break;
      case 14: //DEMAND_SECTION
         vertex_wt = new int[n_vertices];

         for (i = 0; i < n_vertices; i++) {
            if (capacity_volume) {
               if (fscanf(f, "%d%d%d", &k, &l, &m) != 3) {
                  cerr << "UtilGraphLib::read_data ERROR I/O : error reading "
                       << "DEMAND_SECTION. Aborting." << endl;
                  abort();
               }
            } else if (fscanf(f, "%d%d", &k, &l) != 2) {
               cerr << "UtilGraphLib::read_data ERROR I/O : error reading "
                    << "DEMAND_SECTION. Aborting." << endl;
               abort();
            }

            vertex_wt[k-1]  = l;
         }

         if (fscanf(f, "%d%d", &k, &l)) {
            cerr << "UtilGraphLib::read_data ERROR I/O : too much data in "
                 << "DEMAND_SECTION. Aborting." << endl;
            abort();
         }

         break;
      case 18: //EOF
      default:
         break;
      }
   }

   if (f != stdin) {
      fclose(f);
   }

   //calculate all the distances explcitly and then use distance type EXPLICIT
   if (wtype != _EXPLICIT) {
      edge_wt = new int[n_edges];

      for (i = 1, k = 0; i < n_vertices; i++) {
         for (j = 0; j < i; j++) {
            edge_wt[k++] = compute_icost(wtype, i, j);
         }
      }
   }

   /*for(int i = 0; i < n_edges; i++){
     printf("\ngraphLib edge_wt[%d]: %d", i, edge_wt[i]);
     }*/
}

/**********************************************************************/
//This function computes the cost of the edge from va to vb
int UtilGraphLib::compute_icost(const int wtype, const int va, const int vb)
{
   double q1, q2, q3, dx, dy, dz;
   int cost = 0;
   const double RRR      = 6378.388;
   enum DIST {_EXPLICIT, _EUC_2D, _EUC_3D, _MAX_2D, _MAX_3D,
              _MAN_2D, _MAN_3D, _CEIL_2D, _GEO, _ATT
             };

   if (wtype == _GEO) {
      q1 = cos( coordy[va] - coordy[vb] );
      q2 = cos( coordx[va] - coordx[vb] );
      q3 = cos( coordx[va] + coordx[vb] );
      cost = (int) (RRR * acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
   } else {
      dx = coordx[va] - coordx[vb];
      dy = coordy[va] - coordy[vb];

      switch (wtype) {
      case _EUC_2D :
         cost = (int) floor( sqrt( dx * dx + dy * dy ) + 0.5);
         break;
      case _EUC_3D :
         dz = coordz[va] - coordz[vb];
         cost = (int) floor( sqrt( dx * dx + dy * dy + dz * dz) + 0.5);
         break;
      case _MAX_2D :
         cost = (int) fabs(dx);

         if (cost < fabs(dy)) {
            cost = (int) fabs(dy);
         }

         break;
      case _MAX_3D :
         dz = coordz[va] - coordz[vb];
         cost = (int) fabs(dx);

         if (cost < fabs(dy)) {
            cost = (int) fabs(dy);
         }

         if (cost < fabs(dz)) {
            cost = (int) fabs(dz);
         }

         break;
      case _MAN_2D :
         cost = (int) floor( dx + dy + 0.5 );
         break;
      case _MAN_3D :
         dz = coordz[va] - coordz[vb];
         cost = (int) floor( dx + dy + dz + 0.5 );
         break;
      case _CEIL_2D:
         cost = (int)ceil( sqrt( dx * dx + dy * dy ) + 0.5);
         break;
      case _ATT:
         cost = (int)( sqrt( (dx * dx + dy * dy ) / 10 ) + 1);
         break;
      }
   }

   return( cost );
}


