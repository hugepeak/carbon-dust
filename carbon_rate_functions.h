//////////////////////////////////////////////////////////////////////////////
// This file was originally written by Bradley S. Meyer.
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

//##############################################################################
// Includes.
//##############################################################################

#include <iostream>
#include <Libnucnet.h>
#include <boost/lexical_cast.hpp>
#include <nnt/wrappers.h>
#include <nnt/string_defs.h>

//##############################################################################
// Defines.
//##############################################################################

#define S_ARRHENIUS_RATE                "arrhenius rate"
#define S_CARBON_PHOTOINVERSE_RATE      "carbon photoinverse rate"
#define S_CARBON_CONDENSATION_RATE      "carbon condensation rate"
#define S_COMPTON_ELECTRON_RATE         "compton electron rate"
#define S_YA                            "Ya"
#define S_INVERSE_REDUCED_MASS          "inverse reduced mass"
#define S_NUCLEON_NUMBER_PER_ATOM       "nucleon number per atom"  

//##############################################################################
// Validation.  "no" = no validation, "yes" = validation.
//##############################################################################

//#define VALIDATE       "no"

//##############################################################################
// Prototypes.
//##############################################################################

void register_carbon_rate_functions( Libnucnet__Reac * );

void update_carbon_rate_functions_data( nnt::Zone & );

double arrhenius_rate_function( 
  Libnucnet__Reaction *, double, void * 
);

double carbon_photoinverse_rate_function( 
  Libnucnet__Reaction *, double, void * 
);

double carbon_condensation_rate_function( 
  Libnucnet__Reaction *, double, void * 
);

double compton_electron_rate_function( 
  Libnucnet__Reaction *, double, void * 
);

int count_reactants( Libnucnet__Reaction__Element *, void * );

double carbon_compute_Ya( nnt::Zone& );

int carbon_compute_nucleon_number( Libnucnet__Species *, void * );

int carbon_compute_reduced_mass( Libnucnet__Reaction__Element *, void * );
