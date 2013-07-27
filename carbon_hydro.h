////////////////////////////////////////////////////////////////////////////////
// This file was originally written by Tianhong Yu.
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//! \file
//! \brief A header file for the thermodynamics evolution file.
////////////////////////////////////////////////////////////////////////////////

#ifndef CARBON_HYDRO
#define CARBON_HYDRO

//##############################################################################
// Includes.
//##############################################################################

#include <iostream>
#include <fstream>

#include "nnt/auxiliary.h"
#include "nnt/string_defs.h"
#include "nnt/math.h"

#include "user/network_utilities.h"

#include "my_network_utilities.h"

//##############################################################################
// Defines.
//##############################################################################

#define S_POWER                    "power"   // Temperature cooling power. 
#define S_OUTPUT_XML_FILE          "output xml file"

//##############################################################################
// Validation.  "no" = no validation, "yes" = validation.
//##############################################################################

#define VALIDATE       "no"

//##############################################################################
// Prototypes.
//##############################################################################

Libnucnet * get_nucnet( int, char ** );

void initialize_zone( nnt::Zone&, char ** );

void update_zone_properties( nnt::Zone& );

int set_zone( Libnucnet *, nnt::Zone&, char ** );

void get_trajectory_data( char * );

double compute_carbon_k1( nnt::Zone & );

void initialize_bin( nnt::Zone & );

void update_bin( nnt::Zone & );

#endif // CARBON_HYDRO
