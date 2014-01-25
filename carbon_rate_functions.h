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
#define S_ARRHENIUS_INVERSE_RATE        "arrhenius inverse rate"
#define S_CARBON_CONDENSATION_RATE      "carbon condensation rate"
#define S_CARBON_CONDENSATION_INVERSE_RATE  \
                                        "carbon condensation inverse rate"
#define S_COMPTON_ELECTRON_RATE         "compton electron rate"
#define S_SINGLE_RATE                   "single rate"
#define S_YA                            "Ya"
#define S_INVERSE_REDUCED_MASS          "inverse reduced mass"
#define S_NUCLEON_NUMBER_PER_ATOM       "nucleon number per atom"  

#define S_CARBON_NUMBER                 "carbon number"

//##############################################################################
// Validation.  "no" = no validation, "yes" = validation.
//##############################################################################

//#define VALIDATE       "no"

//##############################################################################
// Prototypes.
//##############################################################################

void register_my_rate_functions( Libnucnet__Reac * );

void update_my_rate_functions_data( nnt::Zone & );

double arrhenius_rate_function( 
  Libnucnet__Reaction *, double, void * 
);

double arrhenius_inverse_rate_function( 
  Libnucnet__Reaction *, double, void * 
);

double carbon_condensation_rate_function( 
  Libnucnet__Reaction *, double, void * 
);

double carbon_condensation_inverse_rate_function( 
  Libnucnet__Reaction *, double, void * 
);

int carbon_number_callback(
  Libnucnet__Reaction__Element *, void *
);

double compton_electron_rate_function( 
  Libnucnet__Reaction *, double, void * 
);

double single_rate_function( 
  Libnucnet__Reaction *, double, void * 
);

int count_reactants( Libnucnet__Reaction__Element *, void * );

double carbon_compute_Ya( nnt::Zone& );

int carbon_compute_Ya_callback( Libnucnet__Species *, void * );

int compute_reduced_mass( Libnucnet__Reaction__Element *, void * );

int compute_bond_energy( Libnucnet__Reaction__Element *, void * );

