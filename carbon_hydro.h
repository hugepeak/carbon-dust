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

#include "my_molecule_utilities.h"
#include "my_reaction_utilities.h"

//##############################################################################
// Defines.
//##############################################################################

#define S_POWER                    "power"   // Temperature cooling power. 
#define S_OUTPUT_XML_FILE          "output xml file"

#define S_RUN_DECADE               "run decade"
#define S_EXP_START                "exp start"
#define S_BASE                     "base"
#define S_EXP                      "exp"
#define S_PHOTON_END               "photon end"

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

#endif // CARBON_HYDRO
