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

  //============================================================================
  // Add carbon chains and rings.
  //============================================================================

  add_carbon_molecules_to_nuc( p_nuc, 2, 100, "c" );
  add_carbon_molecules_to_nuc( p_nuc, 2, 100, "r" );

  //============================================================================
  // Add O, O2, CO.
  //============================================================================

  add_oxygen_molecules_to_nuc( p_nuc );

  //============================================================================
  // Add C and C+.
  //============================================================================

  add_carbon_ion_molecules_to_nuc( p_nuc );

  //============================================================================
  // Add He and He+.
  //============================================================================

  add_helium_molecules_to_nuc( p_nuc );

}

//##############################################################################
// add_carbon_molecules_to_nuc().
//##############################################################################

void
add_carbon_molecules_to_nuc(
  Libnucnet__Nuc * p_nuc,
  unsigned int i_lo,
  unsigned int i_hi,
  std::string s_state 
)
{

  int i_state = 0;

  if( s_state.length() != 0 )
    i_state = 1;

  for( unsigned int i = i_lo; i <= i_hi; i++ )
    add_molecule_to_nuc( p_nuc, i, i, i_state, s_state );

}

//##############################################################################
// add_oxygen_molecules_to_nuc().
//##############################################################################

void
add_oxygen_molecules_to_nuc(
  Libnucnet__Nuc * p_nuc
)
{

  int i_state = 0;            // no oxygen ions for now.
  std::string s_state;
  unsigned i_Z, i_A;

  s_state = "";

  i_Z = 0;
  i_A = 1;

  add_molecule_to_nuc( p_nuc, i_Z, i_A, i_state, s_state );

  i_Z = 0;
  i_A = 2;

  add_molecule_to_nuc( p_nuc, i_Z, i_A, i_state, s_state );

  i_Z = 1;
  i_A = 2;

  add_molecule_to_nuc( p_nuc, i_Z, i_A, i_state, s_state );

}

//##############################################################################
// add_carbon_ion_molecules_to_nuc().
//##############################################################################

void
add_carbon_ion_molecules_to_nuc(
  Libnucnet__Nuc * p_nuc
)
{

  int i_state = 1;
  std::string s_state;

  unsigned i_Z = 1;
  unsigned i_A = 1;

  s_state = "g";

  add_molecule_to_nuc( p_nuc, i_Z, i_A, i_state, s_state );

  s_state = "+";

  add_molecule_to_nuc( p_nuc, i_Z, i_A, i_state, s_state );

}

//##############################################################################
// add_helium_molecules_to_nuc().
//##############################################################################

void
add_helium_molecules_to_nuc(
  Libnucnet__Nuc * p_nuc
)
{

  int i_state = 1;
  std::string s_state;

  unsigned i_Z = 2;
  unsigned i_A = 4;

  s_state = "g";

  add_molecule_to_nuc( p_nuc, i_Z, i_A, i_state, s_state );

  s_state = "+";

  add_molecule_to_nuc( p_nuc, i_Z, i_A, i_state, s_state );

}

//##############################################################################
// add_molecule_to_nuc().
//##############################################################################

void
add_molecule_to_nuc(
  Libnucnet__Nuc * p_nuc,
  unsigned int i_Z,
  unsigned int i_A,
  int i_state,
  std::string s_state
)
{
 
  Libnucnet__Species * p_species;
  std::string s_source = "example";
  double d_mass_excess = 0.;
  double d_spin = 0.;

  gsl_vector * p_t9, * p_log10_partf;

  p_t9 = gsl_vector_alloc( 1 );
  p_log10_partf = gsl_vector_alloc( 1 );

  gsl_vector_set( p_t9, 0, 1. );
  gsl_vector_set( p_log10_partf, 0, 0. );

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
      p_log10_partf
    );

  Libnucnet__Nuc__addSpecies( p_nuc, p_species );

  gsl_vector_free( p_t9 );
  gsl_vector_free( p_log10_partf );

}

//##############################################################################
// add_default_generic_molecules_to_nuc().
//##############################################################################

void
add_default_generic_molecules_to_nuc( 
  Libnucnet__Nuc * p_nuc
)
{

  add_generic_molecules_to_nuc( p_nuc, 1, (unsigned) 1e1 );

}

//##############################################################################
// add_generic_molecules_to_nuc().
//##############################################################################

void
add_generic_molecules_to_nuc(
  Libnucnet__Nuc * p_nuc,
  unsigned int i_lo,
  unsigned int i_hi
)
{

  std::string s_state = "";
  int i_state = 0;

  for( unsigned int i = i_lo; i <= i_hi; i++ )
    add_molecule_to_nuc( p_nuc, i, i, i_state, s_state );

}

