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

#include "carbon_molecule_utilities.h"

//##############################################################################
// add_default_molecules_to_nuc().
//##############################################################################

void
add_default_molecules_to_nuc( 
  Libnucnet__Nuc * p_nuc
)
{

  gsl_vector * p_t9, * p_log10_partf;
  unsigned int i_Z, i_A;
  std::string s_state;

  p_t9 = gsl_vector_alloc( 1 );
  p_log10_partf = gsl_vector_alloc( 1 );

  gsl_vector_set( p_t9, 0, 1. );
  gsl_vector_set( p_log10_partf, 0, 0. );

  //============================================================================
  // Add O, O2, CO.
  //============================================================================

  s_state = "";

  i_Z = 0;
  i_A = 1;

  add_molecule_to_nuc( p_nuc, i_Z, i_A, s_state, p_t9, p_log10_partf );

  i_Z = 0;
  i_A = 2;

  add_molecule_to_nuc( p_nuc, i_Z, i_A, s_state, p_t9, p_log10_partf );

  i_Z = 1;
  i_A = 2;

  add_molecule_to_nuc( p_nuc, i_Z, i_A, s_state, p_t9, p_log10_partf );

  //============================================================================
  // Add C and C+.
  //============================================================================

  i_Z = 1;
  i_A = 1;

  s_state = "g";

  add_molecule_to_nuc( p_nuc, i_Z, i_A, s_state, p_t9, p_log10_partf );

  s_state = "+";

  add_molecule_to_nuc( p_nuc, i_Z, i_A, s_state, p_t9, p_log10_partf );

  //============================================================================
  // Add He and He+.
  //============================================================================

  i_Z = 2;
  i_A = 4;

  s_state = "g";

  add_molecule_to_nuc( p_nuc, i_Z, i_A, s_state, p_t9, p_log10_partf );

  s_state = "+";

  add_molecule_to_nuc( p_nuc, i_Z, i_A, s_state, p_t9, p_log10_partf );

  //============================================================================
  // Add carbon chains and rings.
  //============================================================================

  add_carbon_molecules_to_nuc( p_nuc, 2, 8, 2, 8 );

  gsl_vector_free( p_t9 );
  gsl_vector_free( p_log10_partf );

}

//##############################################################################
// add_carbon_molecules_to_nuc().
//##############################################################################

void
add_carbon_molecules_to_nuc(
  Libnucnet__Nuc * p_nuc,
  unsigned int i_chain_lo,
  unsigned int i_chain_hi,
  unsigned int i_ring_lo,
  unsigned int i_ring_hi
)
{

  gsl_vector * p_t9, * p_log10_partf;
  std::string s_state;

  p_t9 = gsl_vector_alloc( 1 );
  p_log10_partf = gsl_vector_alloc( 1 );

  gsl_vector_set( p_t9, 0, 1. );
  gsl_vector_set( p_log10_partf, 0, 0. );

  //============================================================================
  // Add Carbon chains. 
  //============================================================================

  s_state = "c";

  for( unsigned int i = i_chain_lo; i <= i_chain_hi; i++ )
    add_molecule_to_nuc( p_nuc, i, i, s_state, p_t9, p_log10_partf );

  //============================================================================
  // Add Carbon rings. 
  //============================================================================

  s_state = "r";

  for( unsigned int i = i_ring_lo; i <= i_ring_hi; i++ ) 
    add_molecule_to_nuc( p_nuc, i, i, s_state, p_t9, p_log10_partf );

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

  if( s_state.length() )
    i_state = 1;
     
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

