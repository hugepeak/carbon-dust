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

/**
 * @brief A namespace for user-defined rate functions.
 */
namespace my_user
{

int
evolve( nnt::Zone& );

void
update_102_rates(
  WnMatrix *,
  gsl_vector *,
  nnt::Zone &
);

void
update_decade_rates(
  WnMatrix *,
  gsl_vector *,
  nnt::Zone &
);

size_t
get_gsl_vector_index_for_species(
  unsigned int,
  unsigned int
);

unsigned int
get_matrix_row_index_for_species(
  unsigned int,
  unsigned int
);

double
compute_effective_rate(
  double, double, double
);

} // namespace my_user
