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
//! \brief A header file for the evolution routines.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes. 
//##############################################################################

#include "user/evolve.h"

#include "carbon_rate_functions.h"
#include "carbon_hydro.h"

int
evolve( nnt::Zone& );

double
compute_effective_rate(
  nnt::Zone &, double, double
);

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

std::pair<WnMatrix *, gsl_vector *>
get_evolution_matrix_and_vector( 
  nnt::Zone &
);
