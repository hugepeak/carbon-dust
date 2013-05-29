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
//! \brief Carbon molecule utility code. 
////////////////////////////////////////////////////////////////////////////////

#include "my_molecule_utilities.h"

namespace my_user
{

//##############################################################################
// add_default_molecules_to_nuc().
//##############################################################################

void
add_default_molecules_to_nuc( 
  Libnucnet__Nuc * p_nuc
)
{

  add_molecules_to_nuc( p_nuc, 1, 100 );

}

//##############################################################################
// add_molecules_to_nuc().
//##############################################################################

void
add_molecules_to_nuc(
  Libnucnet__Nuc * p_nuc,
  unsigned int i_lo,
  unsigned int i_hi
)
{

  gsl_vector * p_t9, * p_log10_partf;
  std::string s_state;

  p_t9 = gsl_vector_alloc( 1 );
  p_log10_partf = gsl_vector_alloc( 1 );

  gsl_vector_set( p_t9, 0, 1. );
  gsl_vector_set( p_log10_partf, 0, 0. );

  //============================================================================
  // Add Carbon and neutron. 
  //============================================================================

  s_state = "";

  for( unsigned int i = i_lo; i <= i_hi; i++ )
    add_molecule_to_nuc( p_nuc, i, i, s_state, p_t9, p_log10_partf );

  add_molecule_to_nuc( p_nuc, 0, 1, s_state, p_t9, p_log10_partf );

  gsl_vector_free( p_t9 );
  gsl_vector_free( p_log10_partf );

}

//##############################################################################
// add_molecule_to_nuc().
//##############################################################################

void
add_molecule_to_nuc(
  Libnucnet__Nuc * p_nuc,
  unsigned int i_Z,
  unsigned int i_A,
  std::string s_state,
  gsl_vector * p_t9,
  gsl_vector * p_log10_parf
)
{
 
  Libnucnet__Species * p_species;

  std::string s_source = "example";
  int i_state = 0;
  double d_mass_excess = 0.;
  double d_spin = 0.;

  i_state = 0;
     
  p_species = 
    Libnucnet__Species__new(
      i_Z,
      i_A,
      s_source.c_str(),
      i_state,
      s_state.c_str(),
      d_mass_excess,
      d_spin,
      p_t9,
      p_log10_parf
    );

  Libnucnet__Nuc__addSpecies( p_nuc, p_species );

}

}
