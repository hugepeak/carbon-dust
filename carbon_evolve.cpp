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
//! \brief Code for evolving a network.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes. 
//##############################################################################

#include "carbon_evolve.h"

/**
 * @brief A namespace for user-defined functions.
 */

namespace my_user
{

//##############################################################################
// evolve()
//##############################################################################

int
evolve( 
  nnt::Zone& zone
) {

  WnMatrix *p_matrix; 
  size_t i_iter;
  gsl_vector *p_y_old, *p_rhs, *p_sol, *p_work;
  double d_dt;
  std::pair<double,double> check;
  
  //==========================================================================
  // Evolve NSE + weak rates, if appropriate.
  //==========================================================================

  if( zone.hasProperty( nnt::s_USE_HI_T_EQUIL ) )
  {

    if( zone.getProperty( nnt::s_USE_HI_T_EQUIL ) == "yes" )
    {

      if(
        boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) >
        boost::lexical_cast<double>( zone.getProperty( nnt::s_HI_T9_EQUIL ) )
      )
      {
        user::evolve_nse_plus_weak_rates( zone );
        return 1;
      }

    }

  }

  //============================================================================
  // Get timestep. 
  //============================================================================

  d_dt = boost::lexical_cast<double>( zone.getProperty( nnt::s_DTIME ) );

  //============================================================================
  // Save the old abundances.
  //============================================================================

  p_y_old = Libnucnet__Zone__getAbundances( zone.getNucnetZone() );

  //============================================================================
  // Newton-Raphson Iterations.
  //============================================================================

  for( i_iter = 1; i_iter <= I_ITMAX; i_iter++ ) {

    //--------------------------------------------------------------------------
    // Get matrix and rhs vector.
    //--------------------------------------------------------------------------

    boost::tie( p_matrix, p_rhs ) = 
      user::get_evolution_matrix_and_vector( zone );

    //--------------------------------------------------------------------------
    // Update decade rates. 
    //--------------------------------------------------------------------------

    if( 
      zone.hasProperty( S_RUN_DECADE ) &&
      zone.getProperty( S_RUN_DECADE ) == "yes"
    ) 
      update_decade_rates( p_matrix, p_rhs, zone );

    //--------------------------------------------------------------------------
    // Add 1/dt to diagonal.
    //--------------------------------------------------------------------------

    WnMatrix__addValueToDiagonals(
      p_matrix,
      1.0 / boost::lexical_cast<double>( zone.getProperty( nnt::s_DTIME ) )
    );

    //--------------------------------------------------------------------------
    // Correct vector for iteration.
    //--------------------------------------------------------------------------

    p_work = Libnucnet__Zone__getAbundances( zone.getNucnetZone() );
    gsl_vector_sub( p_work, p_y_old );
    gsl_vector_scale( p_work, 1. / d_dt );
    gsl_vector_sub( p_rhs, p_work );

    gsl_vector_free( p_work );

    //--------------------------------------------------------------------------
    // If desired, change value of specific species to prescribed value.
    //--------------------------------------------------------------------------

    if( zone.hasProperty( nnt::s_SPECIFIC_SPECIES ) )
    {
      user::set_specific_species( zone, p_matrix, p_rhs, d_dt );
    }

    //--------------------------------------------------------------------------
    // Solve matrix equation.
    //--------------------------------------------------------------------------

    p_sol = user::solve_matrix_for_zone( zone, p_matrix, p_rhs );

    //--------------------------------------------------------------------------
    // Check solution.
    //--------------------------------------------------------------------------

    check = user::check_matrix_solution( zone, p_sol );

    //--------------------------------------------------------------------------
    // Update abundances.
    //--------------------------------------------------------------------------

    p_work = Libnucnet__Zone__getAbundances( zone.getNucnetZone() );

    gsl_vector_add( p_work, p_sol );

    Libnucnet__Zone__updateAbundances( zone.getNucnetZone(), p_work );

    gsl_vector_free( p_work );

    //--------------------------------------------------------------------------
    // Free matrix, p_rhs, and p_sol. Remember
    // Libnucnet__Zone__computeJacobianMatrix returns a new matrix and
    // Libnucnet__computeFlowVector and WnMatrix__solve return new gsl_vectors
    // each time they are called.
    //--------------------------------------------------------------------------

    WnMatrix__free( p_matrix ); 
    gsl_vector_free( p_rhs );
    gsl_vector_free( p_sol );

    //--------------------------------------------------------------------------
    // Exit iterations if converged.
    //--------------------------------------------------------------------------

    if( check.first < D_MIN ) break;

  }

  //==========================================================================
  // Update abundance changes.
  //==========================================================================

  p_work = Libnucnet__Zone__getAbundances( zone.getNucnetZone() );

  gsl_vector_sub( p_work, p_y_old );

  Libnucnet__Zone__updateAbundanceChanges( zone.getNucnetZone(), p_work );

  gsl_vector_free( p_work );
  
  //==========================================================================
  // Free allocated memory and return.
  //==========================================================================

  gsl_vector_free( p_y_old );

  return (int) i_iter;

}

//##############################################################################
// update_decade_rates().
//##############################################################################

void
update_decade_rates(
  WnMatrix * p_matrix,
  gsl_vector * p_rhs,
  nnt::Zone & zone
)
{

  if( 
    !zone.hasProperty( S_EXP_START ) ||
    !zone.hasProperty( S_BASE ) ||
    !zone.hasProperty( S_EXP )
  )
  {
    std::cerr << "not enough input in udpate decade rates" << std::endl;
    exit( EXIT_FAILURE );
  }

  unsigned int i_start =
    boost::lexical_cast<unsigned int>( zone.getProperty( S_EXP_START ) ); 
  unsigned int i_base = 
    boost::lexical_cast<unsigned int>( zone.getProperty( S_BASE ) ); 
  unsigned int i_exp =
    boost::lexical_cast<unsigned int>( zone.getProperty( S_EXP ) ); 

  double d_t9 = boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) );

  gsl_vector * p_abunds =
    Libnucnet__Zone__getAbundances(
      zone.getNucnetZone()
    );

  double d_na =
    carbon_compute_Ya( zone ) *
    boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ) *
    GSL_CONST_NUM_AVOGADRO;

  unsigned int i_size = 
    (unsigned int) WnMatrix__get_gsl_vector_size( p_rhs );

  //==========================================================================
  // Loop from start to start + exp.
  //==========================================================================

  for( unsigned int m = 0; m < i_exp; m++ )
  {
    
    double d_lo = 
      double( i_start * pow( i_base, m ) );
    double d_hi = 
      double( i_start * pow( i_base, m + 1 ) );

    double d_rate_m = 
      d_na *
      compute_effective_rate( 
        d_lo, d_hi, d_t9
      );

    double d_rate_1 =
      d_rate_m * ( d_hi - d_lo );

    size_t i_vector_m =
      get_gsl_vector_index_for_species( i_start, m );
    size_t i_vector_m1 =
      get_gsl_vector_index_for_species( i_start, m + 1 );

    unsigned int i_matrix_m =
      get_matrix_row_index_for_species( i_start, m );
    unsigned int i_matrix_m1 =
      get_matrix_row_index_for_species( i_start, m + 1 );
      
    //--------------------------------------------------------------------------
    // Matrix elements.
    //--------------------------------------------------------------------------
    
    WnMatrix__assignElement(
      p_matrix,
      i_size,
      i_size,
      d_rate_1 *
        gsl_vector_get( p_abunds, i_vector_m )
    );

    WnMatrix__assignElement(
      p_matrix,
      i_size,
      i_matrix_m,
      d_rate_1 *
        gsl_vector_get( p_abunds, i_size - 1 )
    );

    WnMatrix__assignElement(
      p_matrix,
      i_matrix_m,
      i_size,
      d_rate_m *
        gsl_vector_get( p_abunds, i_vector_m )
    );

    WnMatrix__assignElement(
      p_matrix,
      i_matrix_m,
      i_matrix_m,
      d_rate_m *
        gsl_vector_get( p_abunds, i_size - 1 )
    );

    WnMatrix__assignElement(
      p_matrix,
      i_matrix_m1,
      i_size, 
      -d_rate_m *
        gsl_vector_get( p_abunds, i_vector_m )
    );

    WnMatrix__assignElement(
      p_matrix,
      i_matrix_m1,
      i_matrix_m,
      -d_rate_m *
        gsl_vector_get( p_abunds, i_size - 1 )
    );

    //--------------------------------------------------------------------------
    // Rhs vector. 
    //--------------------------------------------------------------------------
    
    gsl_vector_set( 
      p_rhs,
      i_size - 1,
      gsl_vector_get( p_rhs, i_size - 1) -
        d_rate_1 *
        gsl_vector_get( p_abunds, i_size - 1 ) *
        gsl_vector_get( p_abunds, i_vector_m )
    );

    gsl_vector_set( 
      p_rhs,
      i_vector_m,
      gsl_vector_get( p_rhs, i_vector_m ) -
        d_rate_m *
        gsl_vector_get( p_abunds, i_size - 1 ) *
        gsl_vector_get( p_abunds, i_vector_m )
    );

    gsl_vector_set( 
      p_rhs,
      i_vector_m1,
      gsl_vector_get( p_rhs, i_vector_m1 ) +
        d_rate_m *
        gsl_vector_get( p_abunds, i_size - 1 ) *
        gsl_vector_get( p_abunds, i_vector_m )
    );

  }

  //==========================================================================
  // Free. 
  //==========================================================================
  
  gsl_vector_free( p_abunds );

}

//##############################################################################
// get_gsl_vector_index_for_species().
//##############################################################################

size_t
get_gsl_vector_index_for_species(
  unsigned int i_start,
  unsigned int m
)
{

  return
    size_t(
      i_start + m - 2
    );

}

//##############################################################################
// get_matrix_row_index_for_species().
//##############################################################################

unsigned int
get_matrix_row_index_for_species(
  unsigned int i_start,
  unsigned int m
)
{

  return 
    i_start + m - 1;

} 

//##############################################################################
// compute_effective_rate().
//##############################################################################

double
compute_effective_rate(
  double d_lo,
  double d_hi, 
  double d_t9
)
{

  double d_Rc = 1.7e-8;
  double d_A = 12.;

  return
    M_PI *
    gsl_pow_2( d_Rc ) *
    sqrt(
      2. *
      GSL_CONST_CGSM_BOLTZMANN *
      d_t9 * GSL_CONST_NUM_GIGA /
      (
        d_A / GSL_CONST_NUM_AVOGADRO
      )
    )
    /
    (
      3. *
      ( pow( d_hi, 1./3. ) - pow( d_lo, 1./3. ) )
    );

}

} // namespace my_user
