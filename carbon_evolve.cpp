////////////////////////////////////////////////////////////////////////////////
// This file was originally written by Bradley S. Meyer and Tianhong Yu.
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
  double check_bin = 0.;
  
  //============================================================================
  // Get timestep. 
  //============================================================================

  d_dt = zone.getProperty<double>( nnt::s_DTIME );

  //============================================================================
  // Save the old abundances.
  //============================================================================

  p_y_old = Libnucnet__Zone__getAbundances( zone.getNucnetZone() );

  //============================================================================
  // Save the old bin abundances.
  //============================================================================
    
  std::vector<double> v_bin_old;

  if( 
    zone.hasProperty( "run bin" ) &&
    zone.getProperty<std::string>( "run bin" ) == "yes"
  )
    v_bin_old = get_bin_abundances( zone );

  //============================================================================
  // Newton-Raphson Iterations.
  //============================================================================

  for( i_iter = 1; i_iter <= I_ITMAX; i_iter++ ) {

    //--------------------------------------------------------------------------
    // Get matrix and rhs vector.
    //--------------------------------------------------------------------------

    boost::tie( p_matrix, p_rhs ) = 
      get_evolution_matrix_and_vector( zone );

    //--------------------------------------------------------------------------
    // Add 1/dt to diagonal.
    //--------------------------------------------------------------------------

    WnMatrix__addValueToDiagonals(
      p_matrix,
      1.0 / zone.getProperty<double>( nnt::s_DTIME )
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
    // Update bin rates. 
    //--------------------------------------------------------------------------

    if( 
      zone.hasProperty( "run bin" ) &&
      zone.getProperty<std::string>( "run bin" ) == "yes"
    ) 
      update_bin_rates( zone, p_matrix, p_rhs );

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
    // Compute bin abundance changes.
    //--------------------------------------------------------------------------

    if( 
      zone.hasProperty( "run bin" ) &&
      zone.getProperty<std::string>( "run bin" ) == "yes"
    )
    {
 
      evolve_bin( zone, v_bin_old );

      check_bin = check_bin_change( zone );

    }

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

    if( check.first < D_MIN && check_bin < 1.e-8 ) break;

  }

  //==========================================================================
  // Update abundance changes.
  //==========================================================================

  p_work = Libnucnet__Zone__getAbundances( zone.getNucnetZone() );

  gsl_vector_sub( p_work, p_y_old );

  Libnucnet__Zone__updateAbundanceChanges( zone.getNucnetZone(), p_work );

  gsl_vector_free( p_work );
  
  //==========================================================================
  // Update bin abundance changes.
  //==========================================================================

  if( 
    zone.hasProperty( "run bin" ) &&
    zone.getProperty<std::string>( "run bin" ) == "yes"
  ) {
 
    std::vector<double> v_work = get_bin_abundances( zone );

    for( 
      size_t i = 1; 
      i <= zone.getProperty<size_t>( "bin number" );
      i++ 
    ) 
    {
  
      zone.updateProperty(
        "bin abundance change", 
        boost::lexical_cast<std::string>( i ),
        boost::lexical_cast<std::string>( v_work[i-1] - v_bin_old[i-1] )
      );

    }

  }

  //==========================================================================
  // Free allocated memory and return.
  //==========================================================================

  gsl_vector_free( p_y_old );

  return (int) i_iter;

}

//##############################################################################
// get_evolution_matrix_and_vector().
//##############################################################################

std::pair< WnMatrix *, gsl_vector *>
get_evolution_matrix_and_vector( nnt::Zone& zone )
{

  WnMatrix * p_matrix;
  gsl_vector * p_rhs;

  //--------------------------------------------------------------------------
  // Compute rates.
  //--------------------------------------------------------------------------

  Libnucnet__Zone__computeRates(
    zone.getNucnetZone(),
    zone.getProperty<double>( nnt::s_T9 ),
    zone.getProperty<double>( nnt::s_RHO )
  ); 

  //--------------------------------------------------------------------------
  // Get right hand side vector.
  //--------------------------------------------------------------------------

  p_rhs = Libnucnet__Zone__computeFlowVector( zone.getNucnetZone() );

  //--------------------------------------------------------------------------
  // Get Jacobian matrix.
  //--------------------------------------------------------------------------

  p_matrix = Libnucnet__Zone__computeJacobianMatrix( zone.getNucnetZone() );

  //--------------------------------------------------------------------------
  // Return.
  //--------------------------------------------------------------------------

  return std::make_pair( p_matrix, p_rhs );

}

