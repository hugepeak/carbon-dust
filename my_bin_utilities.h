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
//! \brief A header file for the bin evolution.
////////////////////////////////////////////////////////////////////////////////

#ifndef BIN_UTILITIES
#define BIN_UTILITIES

//##############################################################################
// Includes.
//##############################################################################

#include <iostream>
#include <fstream>

#include "nnt/auxiliary.h"
#include "nnt/string_defs.h"
#include "nnt/math.h"

#include "user/network_utilities.h"

#include "carbon_rate_functions.h"

#define S_BIN_EFFECTIVE_REVERSE "bin effective reverse"

//##############################################################################
// Prototypes.
//##############################################################################

void initialize_bin( nnt::Zone & );

void update_bin_properties( nnt::Zone & );

double compute_k1( nnt::Zone &, double );

void
bin_update_timestep(
  nnt::Zone &, double &, double, double, double
);

int print_bin_abundances( nnt::Zone & ); 

double compute_bin_sum( nnt::Zone & );

void
update_bin_rates(
  nnt::Zone &, WnMatrix *, gsl_vector *
);

void
evolve_bin( 
  nnt::Zone &, std::vector<double> &
);

double
check_bin_change(
  nnt::Zone &
);

std::vector<double>
get_bin_abundances(
  nnt::Zone &
);

int
compute_atom_numbers_in_bin(
  int, int
);

double
compute_effective_rate(
  nnt::Zone &, double, double
);

#endif // BIN_UTILITIES
