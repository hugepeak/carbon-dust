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

  //==========================================================================
  // Save the old bin abundances.
  //==========================================================================
    
  std::vector<double> v_bin_old;

  if( 
    zone.hasProperty( "run bin" ) && zone.getProperty( "run bin" ) == "yes"
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
    // Update bin rates. 
    //--------------------------------------------------------------------------

    if( 
      zone.hasProperty( "run bin" ) &&
      zone.getProperty( "run bin" ) == "yes"
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
    // Evolve bins. 
    //--------------------------------------------------------------------------

    if( 
      zone.hasProperty( "run bin" ) &&
      zone.getProperty( "run bin" ) == "yes"
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

    if( check.first < D_MIN && check_bin < D_MIN ) break;

  }

    if( i_iter > I_ITMAX )
      std::cout << "check " << check.first << " " << check_bin << std::endl;

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
    zone.getProperty( "run bin" ) == "yes"
  )
  {
 
    std::vector<double> v_work = get_bin_abundances( zone );

    for( 
      size_t i = 1; 
      i <= boost::lexical_cast<size_t>( zone.getProperty( "bin number" ) );
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
// update_bin_rates().
//##############################################################################

void
update_bin_rates(
  nnt::Zone & zone,
  WnMatrix * p_matrix,
  gsl_vector * p_rhs
)
{

  if( zone.getProperty( nnt::s_SOLVER ) != nnt::s_ARROW ) {
    std::cerr << "Only work for arrow solver for now" << std::endl;
    exit( EXIT_FAILURE );
  }

  double d_rate;
  int i_begin, i_end;
  int i_atom_number_diff;

  //==========================================================================
  // Get parameters.
  //==========================================================================
    
  int i_bin_start =
    boost::lexical_cast<int>( zone.getProperty( "bin start" ) ); 
  int i_bin_factor = 
    boost::lexical_cast<int>( zone.getProperty( "bin factor" ) ); 
  int i_bin_number =
    boost::lexical_cast<int>( zone.getProperty( "bin number" ) ); 

  if( i_bin_number == 0 )
    return;

  gsl_vector * p_abunds =
    Libnucnet__Zone__getAbundances(
      zone.getNucnetZone()
    );

  double d_na =
    carbon_compute_Ya( zone ) *
    boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ) *
    GSL_CONST_NUM_AVOGADRO;

  size_t i_size = WnMatrix__get_gsl_vector_size( p_rhs );

  //==========================================================================
  // bin_start + c -> b1 
  //==========================================================================
    
  i_begin = i_bin_start + 1;
  i_end  = i_bin_start * i_bin_factor;

  i_atom_number_diff =
    compute_atom_numbers_in_bin( 
      i_begin,
      i_end
    ) -
    i_bin_start;

  d_rate = 
    d_na *
    compute_effective_rate( 
      zone, i_bin_start
    );

  //--------------------------------------------------------------------------
  // Matrix elements.
  //--------------------------------------------------------------------------
  
  WnMatrix__assignElement(
    p_matrix,
    i_size - 1,
    i_size - 1,
    d_rate * gsl_vector_get( p_abunds, i_size - 1 )
  );

  WnMatrix__assignElement(
    p_matrix,
    i_size - 1,
    i_size,
    d_rate * gsl_vector_get( p_abunds, i_size - 2 )
  );

  WnMatrix__assignElement(
    p_matrix,
    i_size,
    i_size - 1,
    d_rate * gsl_vector_get( p_abunds, i_size - 1 ) * i_atom_number_diff
  );

  WnMatrix__assignElement(
    p_matrix,
    i_size,
    i_size,
    d_rate * gsl_vector_get( p_abunds, i_size - 2 ) * i_atom_number_diff
  );

  //--------------------------------------------------------------------------
  // Rhs vector. 
  //--------------------------------------------------------------------------
  
  gsl_vector_set( 
    p_rhs,
    i_size - 2,
    gsl_vector_get( p_rhs, i_size - 2 ) -
      d_rate *
      gsl_vector_get( p_abunds, i_size - 2 ) *
      gsl_vector_get( p_abunds, i_size - 1 )
  );

  gsl_vector_set( 
    p_rhs,
    i_size - 1,
    gsl_vector_get( p_rhs, i_size - 1 ) -
      d_rate *
      gsl_vector_get( p_abunds, i_size - 2 ) *
      gsl_vector_get( p_abunds, i_size - 1 ) * 
      i_atom_number_diff
  );

  //==========================================================================
  // Loop from b1 + c -> b2 to bN-1 + c -> bN.
  //==========================================================================
   
  for( int i = 1; i < i_bin_number; i++ )
  {

    i_begin = i_bin_start * (int)pow( i_bin_factor, i - 1 ) + 1;
    i_end = i_bin_start * (int)pow( i_bin_factor, i );

    d_rate = 
      d_na *
      compute_effective_rate( 
        zone, i_begin, i_end
      );

    i_atom_number_diff =
      compute_atom_numbers_in_bin(
        i_end + 1,
        i_end * i_bin_factor
      ) -
      compute_atom_numbers_in_bin(
        i_begin,
        i_end
      );

    WnMatrix__assignElement(
      p_matrix,
      i_size,
      i_size,
      d_rate * 
        boost::lexical_cast<double>(
          zone.getProperty( 
            "bin abundance", 
            boost::lexical_cast<std::string>( i )
          )
        ) *
        i_atom_number_diff
    );

    gsl_vector_set( 
      p_rhs,
      i_size - 1,
      gsl_vector_get( p_rhs, i_size - 1 ) -
        d_rate *
        boost::lexical_cast<double>(
          zone.getProperty( 
            "bin abundance", 
            boost::lexical_cast<std::string>( i )
          )
        ) *
        gsl_vector_get( p_abunds, i_size - 1 ) *
        i_atom_number_diff
    );

  }

  //==========================================================================
  // Free. 
  //==========================================================================
  
  gsl_vector_free( p_abunds );

}

//##############################################################################
// compute_effective_rate().
//##############################################################################

double
compute_effective_rate(
  nnt::Zone & zone,
  int i_begin
)
{

  if( zone.hasProperty( "k1" ) ) {

    return
      boost::lexical_cast<double>( zone.getProperty( "k1" ) ) *
      pow( (double) i_begin, 2./3. );

  } else {

    return
      compute_carbon_k1( zone ) *
      pow( (double) i_begin, 2./3. );
  }

}

double
compute_effective_rate(
  nnt::Zone & zone,
  int i_begin,
  int i_end 
)
{

  if( zone.hasProperty( "k1" ) ) {

/*
    return
      boost::lexical_cast<double>( zone.getProperty( "k1" ) ) /
      (
        3. *
        ( pow( i_end, 1./3. ) - pow( i_begin, 1./3. ) )
      );
*/
    return
      boost::lexical_cast<double>( zone.getProperty( "k1" ) ) *
      pow( i_end, 2./3. ) * i_end /
      (
        ( i_end + i_begin ) * ( i_end - i_begin + 1 ) / 2.
      );

  } else {

/*
    return
      compute_carbon_k1( zone ) /    
      (
        3. *
        ( pow( i_end, 1./3. ) - pow( i_begin, 1./3. ) )
      );
*/

    return
      compute_carbon_k1( zone ) *
      pow( i_end, 2./3. ) * i_end /
      (
        ( i_end + i_begin ) * ( i_end - i_begin + 1 ) / 2.
      );

  }

}

//##############################################################################
// evolve_bin().
//##############################################################################

void
evolve_bin(
  nnt::Zone & zone,
  std::vector<double> & v_old
)
{

  if( zone.getProperty( nnt::s_SOLVER ) != nnt::s_ARROW ) {
    std::cerr << "Only work for arrow solver for now" << std::endl;
    exit( EXIT_FAILURE );
  }

  int i_begin, i_end;

  int i_bin_start =
    boost::lexical_cast<int>( zone.getProperty( "bin start" ) ); 
  int i_bin_factor = 
    boost::lexical_cast<int>( zone.getProperty( "bin factor" ) ); 
  int i_bin_number =
    boost::lexical_cast<int>( zone.getProperty( "bin number" ) ); 

  if( i_bin_number == 0 )
    return;

  gsl_vector * p_abunds =
    Libnucnet__Zone__getAbundances( zone.getNucnetZone() );

  double d_na =
    carbon_compute_Ya( zone ) *
    boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ) *
    GSL_CONST_NUM_AVOGADRO;

  size_t i_size = WnMatrix__get_gsl_vector_size( p_abunds );

  //==========================================================================
  // "Allocate" new abundance vector.
  //==========================================================================
    
  std::vector<double> v_new( i_bin_number );

  //==========================================================================
  // Compute effective rates.
  //==========================================================================
    
  double d_rate_last =
    d_na *
    compute_effective_rate( 
      zone, i_bin_start
    );

  std::vector<double> v_rates;

  for( int i = 1; i < i_bin_number; i++ ) {

    i_begin = i_bin_start * (int)pow( i_bin_factor, i - 1 ) + 1;
    i_end = i_bin_start * (int)pow( i_bin_factor, i );

    v_rates.push_back( 
      d_na *
      compute_effective_rate( 
        zone, i_begin, i_end
      )
    );

  }

  if( i_bin_number == 1 ) {
 
    v_new[0] =
      d_rate_last *
      gsl_vector_get( p_abunds, i_size - 2 ) *
      gsl_vector_get( p_abunds, i_size - 1 ) *
      boost::lexical_cast<double>( zone.getProperty( nnt::s_DTIME ) )
      +
      v_old[0];

  } else {
    
    v_new[0] =
      (
        d_rate_last *
        gsl_vector_get( p_abunds, i_size - 2 ) *
        gsl_vector_get( p_abunds, i_size - 1 ) *
        boost::lexical_cast<double>( zone.getProperty( nnt::s_DTIME ) )
        +
        v_old[0]
      )
      /
      (
        1. +
        v_rates[0] *     
        gsl_vector_get( p_abunds, i_size - 1 ) *
        boost::lexical_cast<double>( zone.getProperty( nnt::s_DTIME ) )
      );
      
    for( int i = 1; i < i_bin_number - 1; i++ ) {
  
      v_new[i] = 
        (
          v_rates[i-1] *
          v_new[i-1] *
          gsl_vector_get( p_abunds, i_size - 1 ) *
          boost::lexical_cast<double>( zone.getProperty( nnt::s_DTIME ) )
          +
          v_old[i]
        )
        /
        (
          1. +
          v_rates[i] *     
          gsl_vector_get( p_abunds, i_size - 1 ) *
          boost::lexical_cast<double>( zone.getProperty( nnt::s_DTIME ) )
        );
    
    }
  
    v_new[i_bin_number - 1] =
      v_rates[i_bin_number - 2] *
      v_new[i_bin_number - 2] *
      gsl_vector_get( p_abunds, i_size - 1 ) *
      boost::lexical_cast<double>( zone.getProperty( nnt::s_DTIME ) )
      +
      v_old[i_bin_number - 1];

   } // i_bin_number > 1

  //==========================================================================
  // Update bin abundances and abundance changes.
  //==========================================================================
  
  std::vector<double> v_work = get_bin_abundances( zone );

  for( int i = 1; i <= i_bin_number; i++ ) 
  {

    zone.updateProperty(
      "bin abundance", 
      boost::lexical_cast<std::string>( i ),
      boost::lexical_cast<std::string>( v_new[i-1] )
    );

    zone.updateProperty(
      "bin abundance change", 
      boost::lexical_cast<std::string>( i ),
      boost::lexical_cast<std::string>( v_new[i-1] - v_work[i-1] )
    );

  }

  //==========================================================================
  // Free. 
  //==========================================================================
  
  gsl_vector_free( p_abunds );

}

//##############################################################################
// check_bin_change().
//##############################################################################

double
check_bin_change(
  nnt::Zone & zone
)
{

  double d_check = 0., d_checkT;

  for(
    size_t i = 1;
    i <= boost::lexical_cast<size_t>( zone.getProperty( "bin number" ) );
    i++
  )
  {

    if(
      boost::lexical_cast<double>( 
        zone.getProperty( 
          "bin abundance", 
          boost::lexical_cast<std::string>( i ) 
        )
      ) > 1.e-10
    )
    {

      d_checkT = 
        fabs(
          boost::lexical_cast<double>( 
            zone.getProperty( 
              "bin abundance change", 
              boost::lexical_cast<std::string>( i ) 
            )
          ) /   
          boost::lexical_cast<double>( 
            zone.getProperty( 
              "bin abundance", 
              boost::lexical_cast<std::string>( i ) 
            )
          ) 
        ); 

      if( d_checkT > d_check )
        d_check = d_checkT;

    }

  }

  return d_check;

}

//##############################################################################
// get_bin_abundances().
//##############################################################################

std::vector<double>
get_bin_abundances(
  nnt::Zone & zone
)
{

  std::vector<double> v_abunds;

  for( 
    size_t i = 1; 
    i <= boost::lexical_cast<size_t>( zone.getProperty( "bin number" ) ); 
    i++ 
  ) 
  {

    v_abunds.push_back( 
      boost::lexical_cast<double>(
        zone.getProperty( 
          "bin abundance", 
          boost::lexical_cast<std::string>( i ) 
        )
      )
    );

  }

  return v_abunds;

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
    boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ),
    boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) )
  ); 

  //--------------------------------------------------------------------------
  // Zero out small rates.
  //--------------------------------------------------------------------------

  if( zone.hasProperty( nnt::s_SMALL_RATES_THRESHOLD ) )
  {
    user::zero_out_small_rates(
      zone,
      boost::lexical_cast<double>(
        zone.getProperty( nnt::s_SMALL_RATES_THRESHOLD )
      )
    );
  }

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

//##############################################################################
// compute_atom_numbers_in_bin().
//##############################################################################

int
compute_atom_numbers_in_bin(
  int i_begin,
  int i_end
)
{

  return ( i_end + i_begin ) * ( i_end - i_begin + 1 ) / 2;

}

